package hospelhornbg_bioinformatics;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.GenomeBuild;

/*
 * UPDATES:
 * 
 * Initial Version (1.0.0) - December 19, 2017
 * 
 * 1.0.0 -> 1.1.0 | January 12, 2018
 * 	Javadoc
 * 	Updated to be compatible with other package 1.1.0 updates
 * 		SV converter - various getters and setters
 * 		BND pairing
 * 	Clean up VCF header parse storage - now can store ALT fields
 * 	Improved control of internal lists
 * 
 * 1.1.0 -> 1.1.1 | January 23, 2018
 * 	Minor change - added ability to add collection of variants.
 * 
 * 1.1.1 -> 1.1.2 | January 24, 2018
 * 	Minor change - SV bulk parser bug fixed (checks for mateID list size < 1 instead of null for non BND)
 * 
 * 1.1.2 -> 1.2.0 | February 20, 2018
 *  Update for GenomeBuild compatibility.
 * 
 * 1.2.0 -> 1.2.1 | February 22, 2018
 *  Added "makeSkeletonCopy" method that creates a new VariantPool instance with 
 *  the same header information, but no variants.
 *  Idea is that it can be used for filtering.
 *  
 * 1.2.1 -> 1.2.2 | February 26, 2018
 *  Updates genotype FORMAT fields of variants when added to be consistent with pool's.
 *  This is obviously problematic for variants in multiple pools, but whatever.
 *  
 * 1.2.2 -> 1.2.3 | March 5, 2018
 *  Added functions that both add definitions of INFO/FORMAT/ALT/FILTER fields and
 *  put them in active list.
 *  
 * 1.2.3 -> 1.2.4 | April 18, 2018
 *  BND variants with no partner are retained as single BNDs.
 *  These can be manually tossed later.
 *  
 * 1.2.4 -> 1.2.5 | April 19, 2018
 *  Now rejects duplicate genotype fields (quick solution to problem of double GT fields)
 *  
 * 1.2.5 -> 1.2.6 | June 22, 2018
 *  Added a method for counting the variants without having to make a list copy every time.
 *  
 * 1.2.6 -> 1.2.7 | July 12, 2018
 *  Should now be able to recognize BNDs with a CHR2 INFO field as TRA variants
 *  
 */

/*
 * Possible future improvements:
 * 
 */

/**
 * A collection of variants and variant pool metadata. Designed with the VCF format
 * in mind. Methods and implementation should make the kind of information access
 * needed for variant processing simpler.
 * @author Blythe Hospelhorn
 * @version 1.2.7
 * @since July 12, 2018
 */
public class VariantPool {
	
/* --- Constants --- */
	
	public static final int INFODEF_UNK = 0;
	public static final int INFODEF_INT = 1;
	public static final int INFODEF_FLOAT = 2;
	public static final int INFODEF_STRING = 3;
	public static final int INFODEF_FLAG = 4;
	
	public static final int INFODEF_NARGS_VARIABLE = -1;
	public static final int INFODEF_NARGS_PERALT = -2;
	public static final int INFODEF_NARGS_PERALLELE = -3;
	public static final int INFODEF_NARGS_PERGENOTYPE = -4;
	
	/* --- Instance Variables --- */
	
	private GenomeBuild genome;
	
	private List<Variant> varList;
	private List<String> sampleList;
	
	//private List<String> headerLines; //For VCF preservation
	
	private Map<String, InfoDefinition> infoFields;
	private Map<String, String> customAlts;
	private Map<String, String> filters;
	private Map<String, InfoDefinition> genotypeFields;
	
	//Eventually add custom genotype information structure?
	
	//Order lists
	private List<String> info_order;
	private List<String> filter_order;
	private List<String> geno_order;
	
	/* --- Construction/Parsing --- */
	
	/**
	 * Construct an empty variant pool knowing the approximate number of genotyped
	 * samples that this variant pool will have.
	 * @param samples Appoximate number of samples variant pool will refer to
	 */
	public VariantPool(int samples)
	{
		genome = null;
		infoFields = new HashMap<String, InfoDefinition>();
		varList = new LinkedList<Variant>();
		sampleList = new ArrayList<String>(samples);
		customAlts = new HashMap<String, String>();
		filters = new HashMap<String, String>();
		genotypeFields = new HashMap<String, InfoDefinition>();
		genotypeFields.put(Genotype.INFODEF_GT.getKey(), Genotype.INFODEF_GT);
		//headerLines = new LinkedList<String>();
		info_order = new LinkedList<String>();
		filter_order = new LinkedList<String>();
		geno_order = new LinkedList<String>();
		//geno_order.add(Genotype.INFODEF_GT.getKey());
	}
	
	/**
	 * Search the pool for variants that have SVTYPE tags and reparse them as StructuralVariants.
	 * This method also attempts to pair BND SVs using the SECONDARY, MATEID, and EVENT tags.
	 */
	public void castStructuralVariants()
	{
		//System.err.println(Thread.currentThread().getName() + " || VariantPool.castStructuralVariants || Method entered... ");
		List<Variant> nonSV = new LinkedList<Variant>();
		List<StructuralVariant> SV = new LinkedList<StructuralVariant>();
		
		for (Variant v : varList)
		{
			//Look for "SVTYPE" info field
			if (v.hasInfoField("SVTYPE"))
			{
				SV.add(new StructuralVariant(v));
			}
			else nonSV.add(v);
		}
		//System.err.println(Thread.currentThread().getName() + " || VariantPool.castStructuralVariants || SVs found: " + SV.size());
		//System.err.println(Thread.currentThread().getName() + " || VariantPool.castStructuralVariants || non-SVs found: " + nonSV.size());
		//Pair breakends
		List<StructuralVariant> bnd = new LinkedList<StructuralVariant>();
		List<StructuralVariant> otherSV = new LinkedList<StructuralVariant>();
		//StructuralVariant.resetMaxTotalEvidence();
		//StructuralVariant.resetMaxBEDEvidence();
		for (StructuralVariant sv : SV)
		{
			//System.out.println("VariantPool.castStructuralVariants || Structural variant found: " + sv.toString());
			if (sv.isSecondary())
			{
				//System.out.println("VariantPool.castStructuralVariants || Secondary variant found: " + sv.toString());
				try
				{
					BreakendPair bpair = new BreakendPair(sv, SV);
					bnd.add(bpair);
				}
				catch (IllegalArgumentException e)
				{
					otherSV.add(sv); //No partner found...
				}
			}
			else if (sv.getMateIDs().size() < 1) {
				otherSV.add(sv);
			}
			
		}
		//System.err.println(Thread.currentThread().getName() + " || VariantPool.castStructuralVariants || Standard SVs found: " + otherSV.size());
		//System.err.println(Thread.currentThread().getName() + " || VariantPool.castStructuralVariants || BND pairs found: " + bnd.size());
		
		//System.err.println(Thread.currentThread().getName() + " || VariantPool.castStructuralVariants || Casting TRA type SVs... ");
		List<Translocation> tra = new LinkedList<Translocation>();
		List<StructuralVariant> stdSV = new LinkedList<StructuralVariant>();
		for (StructuralVariant sv : otherSV)
		{
			if (sv.getType() == SVType.TRA || (sv.getType() == SVType.BND && (sv.getInfoEntry(Translocation.INFODEF_INFO_CHR2.getKey()) != null)))
			{
				Translocation tr = new Translocation(sv, genome);
				tra.add(tr);
			}
			else {
				sv.removeInfoField(Translocation.INFODEF_INFO_CHR2.getKey());
				stdSV.add(sv);
			}
		}
		
		//System.err.println(Thread.currentThread().getName() + " || VariantPool.castStructuralVariants || Non-TRA/BND SVs found: " + stdSV.size());
		//System.err.println(Thread.currentThread().getName() + " || VariantPool.castStructuralVariants || Translocations found: " + tra.size());
		
		varList.clear();
		varList.addAll(nonSV);
		varList.addAll(stdSV);
		varList.addAll(tra);
		varList.addAll(bnd);
		//System.err.println(Thread.currentThread().getName() + " || VariantPool.castStructuralVariants || Resorting variants... ");
		Collections.sort(varList);
		//System.err.println(Thread.currentThread().getName() + " || VariantPool.castStructuralVariants || Method returning... ");
	}
	
	/* --- Inner Classes --- */
	
	/**
	 * A simple class for holding information on the definition of a type of VCF style
	 * field, namely INFO and genotype FORMAT fields. This definition can be used
	 * to discern a field's type, the number of expected values, and can be easily
	 * formatted for a VCF header.
	 * @author Blythe Hospelhorn
	 * @version 1.0.0
	 * @since December 19, 2018
	 */
	public static class InfoDefinition
	{
		private String key;
		private int type;
		private String description;
		private int nEntries; //Negative means variable
		
		/**
		 * Create a new InfoDefinition and set its initial values.
		 * @param key The key string that will be used to access and recognize this field.
		 * @param type The type of value this field uses (see VariantPool constants).
		 * @param details A string containing the details about what this field is and how to use
		 * it. This details string is included in the VCF header record.
		 * @param entries The number of values that should be expected in this field. Non-finite
		 * or variable values can be represented with a negative number (see VariantPool constants).
		 */
		public InfoDefinition(String key, int type, String details, int entries)
		{
			this.key = key;
			this.type = type;
			this.description = details;
			nEntries = entries;
			if (this.type < 0 || this.type > 4) this.type = VariantPool.INFODEF_UNK;
		}
		
		/**
		 * Get the unique key string of the field defined by this object.
		 * @return Field key
		 */
		public String getKey()
		{
			return key;
		}
		
		/**
		 * Get the type of value this field contains.
		 * @return Value type as an int constant (see VariantPool constants)
		 */
		public int getType()
		{
			return type;
		}
		
		/**
		 * Get a String describing the type of value this field contains.
		 * @return An English string representation of the field type used in a VCF header.
		 */
		public String getTypeString()
		{
			switch (type)
			{
			case INFODEF_FLAG: return "Flag";
			case INFODEF_INT: return "Integer";
			case INFODEF_FLOAT: return "Float";
			case INFODEF_STRING: return "String";
			default: return "Unknown";	
			}
		}
		
		/**
		 * Get the description of this field as found in a VCF header record.
		 * @return Field description string
		 */
		public String getDescription()
		{
			return description;
		}
		
		/**
		 * Get the number of entries expected for a record of this field.
		 * @return A positive number representing a finite amount of that number,
		 * zero if the field functions as a flag, a negative number if the entry
		 * count is variable. See VariantPool constants to interpret negative values.
		 */
		public int getNumberEntries()
		{
			return nEntries;
		}
		
		/**
		 * Set the description to the specified string.
		 * @param s String to set description to.
		 */
		public void setDescription(String s)
		{
			description = s;
		}
		
		/**
		 * Return a string representation of the number of expected entries for this field. 
		 * Values at or above zero simply return that integer as a string. Variable values
		 * marked by a constant negative number are usually represented by a letter or other
		 * character.
		 * @return VCF style representation of expected entry count value.
		 */
		public String getNumberString()
		{
			if (nEntries >= 0) return Integer.toString(nEntries);
			switch (nEntries)
			{
			case INFODEF_NARGS_VARIABLE: return ".";
			case INFODEF_NARGS_PERALT: return "A";
			case INFODEF_NARGS_PERALLELE: return "R";
			case INFODEF_NARGS_PERGENOTYPE: return "G";
			}
			return "";
		}
		
		/**
		 * Parse a short VCF formatted string representing the expected number of entries
		 * for this field and return this value as an integer or a constant int representing
		 * a type of variable.
		 * @param number Expected number of entries as read from a VCF header line
		 * @return Integer (or reserved constant) representation of input value.
		 */
		public static int parseNumberString(String number)
		{
			try
			{
				int i = Integer.parseInt(number);
				return i;
			}
			catch(NumberFormatException e)
			{
				if (number.equals(".")) return INFODEF_NARGS_VARIABLE;
				if (number.equals("A")) return INFODEF_NARGS_PERALT;
				if (number.equals("R")) return INFODEF_NARGS_PERALLELE;
				if (number.equals("G")) return INFODEF_NARGS_PERGENOTYPE;
			}
			return -5;
		}
		
		/**
		 * Get a VCF formatted header record representation of this definition.
		 * @return A string, the piece that usually goes inside the angled brackets, containing
		 * the full VCF formatted definition.
		 */
		public String formatForVCFHeader()
		{
			String nString = getNumberString();
			return "ID=" + key + ",Number=" + nString + ",Type=" + getTypeString() + ",Description=\"" + description + "\"";
		}
		
		public boolean equals(Object o)
		{
			if (o == null) return false;
			if (this == o) return true;
			if (!(o instanceof InfoDefinition)) return false;
			InfoDefinition i = (InfoDefinition)o;
			return (i.getKey().equals(this.getKey()));
		}
		
		public int hashCode()
		{
			return this.getKey().hashCode();
		}
		
		public String toString()
		{
			return getKey() + ":" + getDescription();
		}
		
	}
	
	/* --- Getters --- */
	
	/**
	 * Get a list of the names of all the genotyped samples known by
	 * this variant pool.
	 * <br>The order of this list is set by and sets the order the sample
	 * columns appear in a VCF file.
	 * @return Ordered list of pool sample names.
	 */
	public List<String> getAllSamples()
	{
		List<String> all = new ArrayList<String>(sampleList.size());
		all.addAll(sampleList);
		return all;
	}
	
	/**
	 * Get a list of all the INFO keys "active" in this pool, in the 
	 * order they should appear, if present for a given variant, in a VCF record's
	 * INFO column.
	 * <br>Individual variants may have additional INFO keys not known by the pool as well, but
	 * those will not be found by this method.
	 * <br>The pool may contain definitions for other keys, but only keys included in this list
	 * are considered "active" and will be included in a VCF serialization of this pool.
	 * @return An ordered list of all of the pool's active INFO keys.
	 */
	public List<String> getOrderedInfoKeys()
	{
		List<String> list = new ArrayList<String>(info_order.size());
		list.addAll(info_order);
		return list;
	}
	
	/**
	 * Get a list of all the INFO definitions "active" in this pool, in the 
	 * order they should appear, if present for a given variant, in a VCF record's
	 * INFO column.
	 * <br>Individual variants may have additional INFO keys not known by the pool as well, but
	 * those will not be found by this method.
	 * <br>The pool may contain other definitions, but only keys included in this list
	 * are considered "active" and will be included in a VCF serialization of this pool.
	 * @return An ordered list of all definitions of the pool's active INFO keys.
	 */
	public List<InfoDefinition> getOrderedInfoDefs()
	{
		List<InfoDefinition> list = new ArrayList<InfoDefinition>(info_order.size());
		for (String key : info_order){
			InfoDefinition def = getInfoDef(key);
			if (def != null) list.add(def);
		}
		return list;
	}
	
	/**
	 * Get an unordered list of all unusual ALT field keys defined in a VCF header.
	 * @return List of ALT field keys recorded for this pool.
	 */
	public List<String> getAllAltKeys()
	{
		List<String> alts = new ArrayList<String>(customAlts.size());
		alts.addAll(customAlts.keySet());
		return alts;
	}
	
	/**
	 * Get an ordered list of all defined filter keys in use by this pool. The pool as
	 * well as individual variants may have additional filters defined, but only the
	 * filters with keys in this list will be included in a VCF generated from this pool.
	 * @return An ordered list of active filter keys for this pool.
	 */
	public List<String> getOrderedFilterKeys()
	{
		List<String> list = new ArrayList<String>(filter_order.size());
		list.addAll(filter_order);
		return list;
	}
	
	/**
	 * Get a list of all the FORMAT keys "active" in this pool, in the 
	 * order they should appear in a VCF header.
	 * <br>Individual variants and genotypes may have additional FORMAT keys not known by the pool as well, but
	 * those will not be found by this method.
	 * <br>The pool may contain definitions for other keys, but only keys included in this list
	 * are considered "active" and will be included in a VCF serialization of this pool.
	 * @return An ordered list of all of the pool's active FORMAT keys.
	 */
	public List<String> getOrderedGenotypeKeys()
	{
		List<String> list = new ArrayList<String>(geno_order.size());
		list.addAll(geno_order);
		return list;
	}
	
	/**
	 * Get the definition for an INFO field given a known key.
	 * @param key INFO field key
	 * @return InfoDefinition instance containing field definition information, null
	 * if the field referred to by the provided key has never been defined for this pool.
	 */
	public InfoDefinition getInfoDef(String key)
	{
		return infoFields.get(key);
	}
	
	/**
	 * Get the description for an ALT field given a known key.
	 * @param key ALT field key
	 * @return Description recorded for the ALT field with the given key, null
	 * if the field referred to by the provided key has never been defined for this pool.
	 */
	public String getAltDescription(String key)
	{
		return customAlts.get(key);
	}
	
	/**
	 * Get the description for a FILTER field given a known key.
	 * @param key FILTER field key
	 * @return Description recorded for the FILTER field with the given key, null
	 * if the field referred to by the provided key has never been defined for this pool.
	 */
	public String getFilterDescription(String key)
	{
		return filters.get(key);
	}
	
	/**
	 * Get the definition for a FORMAT field given a known key.
	 * @param key FORMAT field key
	 * @return InfoDefinition instance containing field definition information, null
	 * if the field referred to by the provided key has never been defined for this pool.
	 */
	public InfoDefinition getGenoDef(String key)
	{
		return genotypeFields.get(key);
	}
	
	/**
	 * Get an ordered list of all variants in this pool.
	 * @return List of all variants currently in pool.
	 */
	public List<Variant> getVariants()
	{
		List<Variant> all = new ArrayList<Variant>(varList.size());
		all.addAll(varList);
		return all;
	}
	
	/**
	 * Get the number of variants currently in the pool.
	 * @return Size of the internal variant list.
	 */
	public int countVariants()
	{
		return varList.size();
	}
	
	/**
	 * Scan all variants in the pool to form a set of all chromosomes referenced
	 * by any variant in the pool.
	 * @return Unordered set of unique chromosomes found in variants in the pool.
	 */
	public Set<Contig> getAllChromosomes()
	{
		Set<Contig> allChrom = new HashSet<Contig>();
		for (Variant v: varList)
		{
			allChrom.addAll(v.getAllChromosomes());
		}
		return allChrom;
	}
	
	/**
	 * Get the genome build recorded for this variant pool.
	 * <br>The genome build may be used for additional coordinate/contig
	 * information and can also be encoded into a VCF file.
	 * <br>VCF files with a recognized build name in their header should
	 * be parsed against that build automatically.
	 * @return GenomeBuild set for this Variant Pool
	 */
	public GenomeBuild getGenomeBuild()
	{
		return genome;
	}
	
	/* --- Setters --- */
	
	/**
	 * Add a variant to the pool. The variant will simply be tacked on to the end
	 * of the variant list, so if any order is desired, the internal list will have
	 * to be resorted.
	 * <br>Additionally, the pool does not screen for and automatically remove duplicates.
	 * @param v Variant to add. Must not be null.
	 */
	public void addVariant(Variant v)
	{
		if (v == null) return;
		varList.add(v);
		v.clearGenotypeFields();
		for (String k : geno_order) v.addGenotypeField(k);
	}
	
	/**
	 * Add a collection of variants to the pool. All variants will be tacked on
	 * to the end of the variant list, and duplicates will not be checked. If any order
	 * is desired, the internal list must be resorted with the sortVariants method.
	 * @param varSet Collection of variants to add to pool. Must not be null.
	 */
	public void addVariants(Collection<Variant> varSet)
	{
		if (varSet == null) return;
		//varList.addAll(varSet);
		for (Variant v : varSet) addVariant(v);
	}
	
	/**
	 * Sort the internal variant list. Variants are ordered as determined by the Variant
	 * class or subclass.
	 */
	public void sortVariants()
	{
		Collections.sort(varList);
	}
	
	/**
	 * Clear all variants from this pool.
	 */
	public void clearVariants()
	{
		varList.clear();
	}
	
	/**
	 * Add a sample to the end of the sample list.
	 * @param sampleName Name of the sample to add.
	 */
	public void addSample(String sampleName)
	{
		sampleList.add(sampleName);
	}
	
	/**
	 * Remove all sample references from this variant pool.
	 */
	public void clearSamples()
	{
		sampleList.clear();
	}
	
	/**
	 * Add an INFO field definition to the INFO dictionary. Note that this
	 * method DOES NOT add the field to the active key list; that must be done
	 * afterwards!
	 * <br>Keys that are null or empty are rejected. Definition cannot be null.
	 * @param key INFO field key.
	 * @param def Definition of INFO field.
	 */
	public void addInfoFieldDefinition(String key, InfoDefinition def)
	{
		//info_order.add(key);
		if (key == null || key.isEmpty()) return;
		if (def == null) return;
		infoFields.put(key, def);
	}
	
	/**
	 * If the INFO key has been defined for this variant, it will be added
	 * to the end of the list of active INFO keys. 
	 * This method does not check for or remove duplicates.
	 * @param key The INFO key to add to the list of fields to request of variants
	 * when this pool is serialized as a VCF.
	 */
	public void addInfoKeyToActiveList(String key)
	{
		if (key == null) return;
		if (key.isEmpty()) return;
		if (!infoFields.containsKey(key)) return;
		info_order.add(key);
	}
	
	/**
	 * Add an INFO field definition to the INFO dictionary AND add
	 * it to the active key list.
	 * @param key INFO field key.
	 * @param def Definition of INFO field.
	 */
	public void addInfoField(String key, InfoDefinition def)
	{
		addInfoFieldDefinition(key, def);
		addInfoKeyToActiveList(key);
	}
	
	/**
	 * Search for and remove all instances of the provided INFO key from
	 * the active INFO key list.
	 * @param key The INFO key to remove from the list of fields to request of variants
	 * when this pool is serialized as a VCF.
	 */
	public void removeInfoKeyFromActiveList(String key)
	{
		info_order.remove(key);
	}
	
	/**
	 * Clear the list of INFO keys of fields to request of variants 
	 * when this pool is serialized as a VCF.
	 */
	public void clearActiveInfoKeys()
	{
		info_order.clear();
	}
	
	/**
	 * Parse a VCF style string denoting the type of an annotation field (String, Integer, etc.)
	 * to an enumerated integer representing that type (see VariantPool constants).
	 * @param type String representing the field type
	 * @return Constant integer denoting field type
	 */
	public static int getInfoDefType(String type)
	{
		if (type.equals("Integer")) return INFODEF_INT;
		if (type.equals("Float")) return INFODEF_FLOAT;
		if (type.equals("String")) return INFODEF_STRING;
		if (type.equals("Flag")) return INFODEF_FLAG;
		return INFODEF_UNK;
	}
	
	/**
	 * Add a filter to this pool.
	 * @param key Name/key of FILTER to add.
	 * @param description Description of FILTER to add.
	 */
	public void addFilter(String key, String description)
	{
		filter_order.add(key);
		filters.put(key, description);
	}
	
	/**
	 * Clear the list of filters (FILTER fields) for this pool.
	 */
	public void clearFilters()
	{
		filter_order.clear();
	}

	/**
	 * Add a custom alternate allele definition (ALT field) to this pool.
	 * @param key ALT key of definition to add.
	 * @param description Description of field to add.
	 */
	public void addCustomAlt(String key, String description)
	{
		if (key == null | key.isEmpty()) return;
		if (description == null) return;
		customAlts.put(key, description);
	}
	
	/**
	 * Clear all custom ALT field definitions from this pool.
	 */
	public void clearCustomAlts()
	{
		customAlts.clear();
	}

	/**
	 * Add a genotype FORMAT field definition to the FORMAT dictionary. Note that this
	 * method DOES NOT add the field to the active key list; that must be done
	 * afterwards!
	 * <br>Keys that are null or empty are rejected. Definition cannot be null.
	 * @param key FORMAT field key.
	 * @param def Definition of FORMAT field.
	 */
	public void addFormatFieldDefinition(String key, InfoDefinition def)
	{
		if (key == null || key.isEmpty()) return;
		if (def == null) return;
		genotypeFields.put(key, def);
	}
	
	/**
	 * If the genotype FORMAT key has been defined for this variant, it will be added
	 * to the end of the list of active FORMAT keys. 
	 * This method does not check for or remove duplicates.
	 * @param key The FORMAT key to add to the list of fields to request of variants
	 * when this pool is serialized as a VCF.
	 */
	public void addFormatKeyToActiveList(String key)
	{
		if (key == null) return;
		if (key.isEmpty()) return;
		if (!genotypeFields.containsKey(key)) return;
		if (geno_order.contains(key)) return;
		geno_order.add(key);	
	}
	
	/**
	 * Add a genotype FORMAT field definition to the FORMAT dictionary AND add
	 * it to the active key list.
	 * @param key FORMAT field key.
	 * @param def Definition of FORMAT field.
	 */
	public void addFormatField(String key, InfoDefinition def)
	{
		addFormatFieldDefinition(key, def);
		addFormatKeyToActiveList(key);
	}
	
	/**
	 * Search for and remove all instances of the provided genotype FORMAT key from
	 * the active FORMAT key list.
	 * @param key The FORMAT key to remove from the list of fields to request of variants
	 * when this pool is serialized as a VCF.
	 */
	public void removeFormatKeyFromActiveList(String key)
	{
		geno_order.remove(key);
	}
	
	/**
	 * Clear the list of genotype FORMAT keys of fields to request of variants 
	 * when this pool is serialized as a VCF.
	 */
	public void clearActiveFormatKeys()
	{
		geno_order.clear();
	}
	
	/**
	 * Set a genome build (an object containing information on contig names
	 * and lengths for a genome reference build) for this pool.
	 * @param g Build to set for this pool.
	 */
	public void setGenomeBuild(GenomeBuild g)
	{
		this.genome = g;
	}
	
	/* --- Usage --- */
	
	/**
	 * Run this pool against a truth set pool and flag "confirmed" all variants in this pool
	 * that are equivalent to any variant in the truth set.
	 * @param truthset A variant pool containing a set of all confirmed variants to check
	 * this pool against.
	 * @param stringent_bothEnds Stringency setting: If on, both ends of the SV (not just one or the other) must be
	 * considered in range of each other.
	 * @param stringent_CI Stringency setting: If on, then the valid range around a breakpoint for considering
	 * two variants equivalent is the 95% confidence interval rather than the default 90% CI.
	 * @param stringent_overlap Stringency setting: If on, intervals at either end of the variant
	 * must encompass the actual endpoint calls themselves to be considered equivalent. Otherwise,
	 * a simple overlap of the two ranges is considered sufficient.
	 * @return The number of variants in this pool that have been flagged confirmed.
	 */
	public int confirm(VariantPool truthset, boolean stringent_bothEnds, boolean stringent_CI, boolean stringent_overlap)
	{
		final List<Variant> truth = truthset.getVariants();
		int confirmcount = 0;
		
		for (Variant t : truth)
		{
			for (Variant v : varList)
			{
				if (v.equivalent(t))
				{
					confirmcount++;
					v.confirm();
					break;
				}
				if (t instanceof StructuralVariant && v instanceof StructuralVariant)
				{
					StructuralVariant svt = (StructuralVariant)t;
					StructuralVariant svv = (StructuralVariant)v;
					if (svv.equivalent(svt, stringent_bothEnds, stringent_CI, stringent_overlap))
					{
						confirmcount++;
						v.confirm();
						break;
					}
				}
			}
		}
		
		return confirmcount;
	}
	
	/**
	 * Remove all "confirmed" flags from all variants in this pool.
	 */
	public void unconfirmAll()
	{
		for (Variant v : varList) v.unconfirm();
	}
	
	/**
	 * Generate a new VariantPool with the same metadata (such as INFO and ALT definitions),
	 * but with an empty variant set.
	 * @return A new VariantPool instance with copies of all the metadata in this instance,
	 * but without any variants.
	 */
	public VariantPool makeSkeletonCopy()
	{
		VariantPool copy = new VariantPool(this.sampleList.size());
		//Genome build
		copy.setGenomeBuild(this.getGenomeBuild());
		//Samples
		for (String s : sampleList) copy.addSample(s);
		//INFO
		for (String s : this.infoFields.keySet())copy.addInfoFieldDefinition(s, infoFields.get(s));
		copy.info_order.addAll(this.info_order);
		//FORMAT
		for (String s : this.genotypeFields.keySet())copy.addFormatFieldDefinition(s, genotypeFields.get(s));
		copy.geno_order.addAll(this.geno_order);
		//ALT
		for (String s : this.customAlts.keySet())copy.addCustomAlt(s, customAlts.get(s));
		//FILTER
		for (String s : this.filters.keySet())copy.addFilter(s, filters.get(s));
		copy.filter_order.addAll(this.filter_order);
		
		return copy;
	}
	
}
