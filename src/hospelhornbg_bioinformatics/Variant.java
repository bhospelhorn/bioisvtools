
package hospelhornbg_bioinformatics;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import hospelhornbg_bioinformatics.VariantPool.InfoDefinition;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.GeneFunc;
import hospelhornbg_genomeBuild.GeneSet;
import hospelhornbg_genomeBuild.GenomeBuild;
import hospelhornbg_segregation.CandidateFlags;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

/*
 * UPDATE NOTES
 * 
 * Initial version date: December 15, 2017 (1.0.0)
 * 
 * 1.0.0 -> 1.1.0 | January 10, 2018
 * 	Added some more setters for improved accessibility
 * 	Javadoc annotation
 *  Changed all internal arrays to ArrayLists (AltAlleles, Filters, genotype fields)
 *  
 * 1.1.0 -> 1.1.1 | January 24, 2018
 *  Added line generator for UCSC genome browser track BED line
 *  
 * 1.1.1 -> 1.1.2 | January 25, 2018
 * 	Turns out the quality field is a float?
 * 
 * 1.1.2 -> 1.1.3 | January 29, 2018
 * 	Minor bug fix with VCF serialization (null pointer to filters array if no filters... fix boolean funct)
 * 
 * 1.1.3 -> 1.1.4 | January 30, 2018
 * 	Minor bug fix with VCF serialization (genotype serialization)
 * 
 * 1.1.4 -> 1.1.5 | February 9, 2018
 * 	Minor bug fix with VCF serialization (Empty VarID)
 * 
 * 1.1.5 -> 1.2.0 | February 13, 2018
 * 	Updated to work with known genome builds. Chromosome field is now a Contig
 * 	object instead of a String.
 * 
 * 1.2.0 -> 1.2.1 | April 16, 2018
 * 	Added protected method to access list of genotyped samples
 * 
 * 1.2.1 -> 1.2.2 | April 19, 2018
 * 	Added boolean protector to VCF line parse constructor
 * 
 * 1.2.2 -> 1.2.3 | July 20, 2018
 * 	Added inRegion function
 * 
 * 1.2.3 -> 1.3.0 | August 10, 2018
 * 	Moved the GeneFunc instance variable from subclass StructuralVariant
 * 
 * 1.3.0 -> 1.3.1 | January 18, 2019
 * 	Added support vector field
 * 
 * 1.3.1 -> 1.3.2 | January 22, 2019
 * 	Added optional candidate flags map
 * 
 */

/*
 * Possible future improvements:
 * 	- Make an allele object so that special alt alleles can be better defined
 */

/**
 * An object representing a genetic variant. Instance variables, methods, and general interface
 * of this class were made with the VCF format in mind. 
 * <br>Because it uses INFO field style annotations, it can in theory take as many annotations as
 * memory allows - internally these are stored in a String-(type) map.
 * <br>Providing InfoDefinition instances allows for parsing of these INFO annotations into 
 * integers, floats, and flags. By default, they are all stored as either strings or flags.
 * @author Blythe Hospelhorn
 * @version 1.3.2
 * @since January 22, 2019
 *
 */
public class Variant implements Comparable<Variant>{
	
	/* --- Constants --- */
	
	/**
	 * Chromosome string is not of a standard recognizable type, or is unset.
	 */
	public static final int CHROMOSOME_UNKNOWN = 0;
	
	/**
	 * Chromosome string represents an autosome (positive integer value).
	 */
	public static final int CHROMOSOME_AUTOSOME = 1;
	
	/**
	 * Chromosome string represents a sex chromosome (X, Y, W, or Z).
	 */
	public static final int CHROMOSOME_SEXCHROMOSOME = 2;
	
	/**
	 * Chromosome string represents a mitochondrial chromosome (M, MT)
	 */
	public static final int CHROMOSOME_MITOCHONDRIAL = 3;
	
	/**
	 * Chromosome string represents a contig ID.
	 */
	public static final int CHROMOSOME_CONTIG = 4;
	
	/**
	 * Chromosome string represents a "compound" chromosome, such as that used to mark
	 * possible interchromosomal translocations.
	 */
	public static final int CHROMOSOME_COMPOUND = 5;

	/* --- Instance Variables --- */
	
	/**
	 * Switch set if variant is compared to a variant in a truth set and considered equivalent.
	 * If this value is true, the variant may be considered "confirmed."
	 */
	private boolean confirmed;
	
	/**
	 * The chromosome/contig the variant is mapped to.
	 */
	private Contig chromosome;
	/**
	 * The 1-based coordinate of the position of the first base in the variant along the chromosome.
	 */
	private int position;
	
	/**
	 * A field for giving the variant a name. This may be used to mark a variant found in a database.
	 */
	private String variantID;
	
	/**
	 * The reference allele for this variant (as a string).
	 */
	private String refAllele;
	/**
	 * A list of alternate alleles for this variant. Each alt should be a string.
	 */
	private ArrayList<String> altAlleles;
	
	/**
	 * A phred scaled quality score. More or less optional.
	 */
	private double quality;
	/**
	 * Whether the variant has passed all of the filters set for it.
	 */
	private boolean passedFilters;
	/**
	 * A list of the filters the variant was tried against and failed, represented as strings.
	 */
	private ArrayList<String> filters;

	/**
	 * A map of all string type INFO annotations for this variant. This map includes
	 * fields whose type is unknown and have been left as unparsed text.
	 */
	private Map<String, String[]> infoZ;
	/**
	 * A map of all (parsed) integer type INFO annotations for this variant.
	 */
	private Map<String, int[]> infoI;
	/**
	 * A map of all (parsed) float type INFO annotations for this variant.
	 */
	private Map<String, double[]> infoF;
	/**
	 * A set of all flag type INFO annotations for this variant.
	 */
	private Set<String> infoB;
	
	/**
	 * An ordered list of the genotype annotation field keys
	 * for the fields that should appear in each genotype column.
	 */
	private ArrayList<String> genoFields;
	/**
	 * A map of sample IDs (as strings) to genotypes. Sample IDs that are requested of this variant,
	 * but not in this map are assumed to be ref/ref.
	 */
	private Map<String, Genotype> genotypes;
	
	/**
	 * If noted, the coding effect the variant has.
	 */
	private GeneFunc function;
	
	/**
	 * A field for holding a support vector field. This field is not
	 * filled automatically from VCF parsing. It can be filled
	 * manually, though.
	 */
	private Set<String> support_set;
	
	/**
	 * Candidate flags map for optional use.
	 * This map is not instantiated normally (defaults to null), nor
	 * is is populated automatically by VCF parsing. 
	 * It's used by the VarTable parser, but not otherwise in order to 
	 * save memory.
	 */
	private Map<Integer, CandidateFlags> cflagsMap;
	
	/* --- Construction/Parsing --- */
	
	/**
	 * Construct a Variant by parsing information from a single line (variant record) of a VCF file.
	 * <br> This method REQUIRES the line be a line of VCF formatted text.
	 * <br> This method will throw an exception if fed a VCF header line. It expects a variant record.
	 * @param VCF_line Line, in plain ASCII text, from a VCF file representing the variant record to parse.
	 * @param samples An ordered list of the sample names (patient IDs) used in the VCF file so that the
	 * genotype columns can be properly read.
	 * @param genome An object detailing the coordinate system of the reference the calls in the
	 * supplied VCF were supposedly made against. If this value is null, new dummy contigs with the 
	 * chromosome names parsed from the VCF will be used instead.
	 * @throws UnsupportedFileTypeException If the line given cannot be parsed as a VCF variant.
	 */
	public Variant(String VCF_line, List<String> samples, GenomeBuild genome, boolean allowBuildModification) throws UnsupportedFileTypeException
	{
		//System.err.println(Thread.currentThread().getName() + " || Variant.<init> || Called! ");
		//System.err.println(Thread.currentThread().getName() + " || Variant.<init> || Sample list: ");
		//for (String s : samples) System.err.println("\t" + s);
		setDefaults();
		String[] fields = VCF_line.split("\t");
		if (fields == null) throw new FileBuffer.UnsupportedFileTypeException();
		if (fields.length < 8) throw new FileBuffer.UnsupportedFileTypeException();
		
		try
		{
			String chrom = fields[0];
			//Find chrom
			if (genome != null && !allowBuildModification) {
				//System.err.println(Thread.currentThread().getName() + " || Variant.<init> || Genome: " + genome.getBuildName());
				chromosome = genome.getContig(chrom);
				//Toss variant if can't find contig
				if (chromosome == null) {
					System.err.println(Thread.currentThread().getName() + " || Variant.<init> || Contig not recognized: " + chrom);
					throw new NullPointerException();
				}
			}
			else if (genome != null && allowBuildModification) 
			{
				chromosome = genome.getContig(chrom);
				if (chromosome == null)
				{
					chromosome = new Contig();
					chromosome.setUDPName(chrom);
					chromosome.setUCSCName(chrom);
					genome.addContig(chromosome);
				}
			}
			else
			{
				chromosome = new Contig();
				chromosome.setUDPName(chrom);
				chromosome.setUCSCName(chrom);
				//chromosome.setType(Contig.SORTCLASS_UNKNOWN);
			}
			position = Integer.parseInt(fields[1]);
			variantID = fields[2];
			if (variantID.equals(".")) variantID = ".";
			refAllele = fields[3]; if (refAllele.equals(".")) refAllele = ".";
			String[] stArr = fields[4].split(",");
			altAlleles = new ArrayList<String>(stArr.length);
			for (String s : stArr) altAlleles.add(s);
			if (!fields[5].equals(".")) quality = Double.parseDouble(fields[5]);
			if (!fields[6].equals("."))
			{
				if (fields[6].equals("PASS")) passedFilters = true;
				else
				{
					stArr = fields[6].split(";");
					filters = new ArrayList<String>(stArr.length);
					for (String s : stArr) filters.add(s);
				}
			}
			if (!fields[7].isEmpty() && !fields[7].equals("."))
			{
				String[] infofields = fields[7].split(";");
				for (String s : infofields)
				{
					if (s.indexOf("=") < 0)
					{
						//It's probably a flag
						//System.out.println("Variant<init> || Flag found: " + s);
						infoB.add(s);
					}
					else
					{
						String[] infopieces = s.split("=");
						if (infopieces == null) throw new UnsupportedFileTypeException();
						if (infopieces.length < 2) throw new UnsupportedFileTypeException();
						String[] infoVals = infopieces[1].split(",");
						infoZ.put(infopieces[0], infoVals);	
					}
				}
			}
			if (fields.length > 8 && samples != null)
			{
				//System.err.println(Thread.currentThread().getName() + " || Variant.<init> || More than 8 fields found!");
				int sampNum = fields.length - 9;
				//System.err.println(Thread.currentThread().getName() + " || Variant.<init> || Number of samples: " + sampNum);
				String format = fields[8];
				//System.err.println(Thread.currentThread().getName() + " || Variant.<init> || Format string: " + format);
				stArr = format.split(":");
				genoFields = new ArrayList<String>(stArr.length);
				for (String s : stArr) genoFields.add(s);
				//System.err.println(Thread.currentThread().getName() + " || Variant.<init> || format field count: " + genoFields.size());
				for (int i = 0; i < sampNum; i++)
				{
					//System.err.println(Thread.currentThread().getName() + " || Variant.<init> || genotyping sample " + i);
					String g = fields[9 + i];
					Genotype gt = new Genotype(format, g);
					genotypes.put(samples.get(i), gt);
				}
			}
			
		}
		catch (NumberFormatException | ArrayIndexOutOfBoundsException e)
		{
			e.printStackTrace();
			throw new UnsupportedFileTypeException();
		}
		
		String key = GeneSet.INFODEF_INFO_GFUNC.getKey();
		String func = getSingleStringInfoEntry(key);
		if (func != null)
		{
			removeInfoField(key);
			function = GeneFunc.getFunction(func);
		}
	
		//System.err.println(Thread.currentThread().getName() + " || Variant.<init> || Genotype map size: " + genotypes.size());
		
	}
 	
	/**
	 * PROTECTED constructor that constructs a new variant by REFERENCING the internal structures
	 * of an existing variant.
	 * This is used for converting an existing variant to an extended child type (such as a StructuralVariant).
	 * <br>Template variant should be deleted after this constructor is called, not used alongside its "copy."
	 * Otherwise, many changes made to one will be made to both!
	 * @param other The template variant to copy data references from.
	 */
	protected Variant(Variant other)
	{
		//REFERENCES, DOES NOT COPY. That's why it's protected - template is expected to be deleted!
		setDefaults();
		confirmed = other.confirmed;
		chromosome = other.chromosome;
		position = other.position;
		variantID = other.variantID;
		refAllele = other.refAllele;
		altAlleles = other.altAlleles;
		quality = other.quality;
		passedFilters = other.passedFilters;
		filters = other.filters;
		infoZ = other.infoZ;
		infoI = other.infoI;
		infoF = other.infoF;
		infoB = other.infoB;
		genoFields = other.genoFields;
		genotypes = other.genotypes;
		function = other.function;
	}
	
	/**
	 * Construct an empty variant.
	 * Leaving it empty may produce odd results.
	 * This constructor should be used primarily with parsers.
	 */
	public Variant()
	{
		setDefaults();
	}
	
	/**
	 * Set all fields to default "empty" values and instantiate collections.
	 * Arrays are left null.
	 */
	private void setDefaults()
	{
		confirmed = false;
		chromosome = null;
		position = -1;
		variantID = ".";
		refAllele = ".";
		altAlleles = new ArrayList<String>(4);
		genoFields = new ArrayList<String>(8); genoFields.add("GT");
		quality = -1.0;
		passedFilters = false;
		filters = new ArrayList<String>(4);
		infoZ = new HashMap<String, String[]>();
		infoI = new HashMap<String, int[]>();
		infoF = new HashMap<String, double[]>();
		infoB = new HashSet<String>();
		genotypes = new HashMap<String, Genotype>();
	}
	
	/**
	 * Request that a field value with a known key be treated as an integer
	 * instead of a string.
	 * This method searches for the key in the string map, and if found, removes
	 * it from the string map and attempts to parse the value as an integer.
	 * If it succeeds, it stores the key-value pair in the integer map instead.
	 * If it fails, it replaces the original key-(unparsed) value pair in the string map.
	 * Flags are not stored in the string map and will not be found by this method.
	 * @param fieldkey The string representation of the field key. In a VCF file, this is the
	 * string that comes before the equals sign in an INFO field annotation.
	 * @return True - If value was found and successfully parsed as an integer.
	 * <br>False - If value mapping to given key was not found or could not be interpreted as an integer.
	 */
	private boolean parseStringInfoField_toInt(String fieldkey)
	{
		String rawValue[] = infoZ.remove(fieldkey);
		if (rawValue == null) return false;
		try
		{
			int vals[] = new int[rawValue.length];
			for (int i = 0; i < rawValue.length; i++)
			{
				vals[i] = Integer.parseInt(rawValue[i]);
			}
			infoI.put(fieldkey, vals);
		}
		catch (NumberFormatException e)
		{
			//e.printStackTrace();
			infoZ.put(fieldkey, rawValue);
			return false;
		}
		
		return true;
	}
	
	/**
	 * Request that a field value with a known key be treated as an float (double)
	 * instead of a string.
	 * This method searches for the key in the string map, and if found, removes
	 * it from the string map and attempts to parse the value as a double.
	 * If it succeeds, it stores the key-value pair in the float map instead.
	 * If it fails, it replaces the original key-(unparsed) value pair in the string map.
	 * Flags are not stored in the string map and will not be found by this method.
	 * @param fieldkey The string representation of the field key. In a VCF file, this is the
	 * string that comes before the equals sign in an INFO field annotation.
	 * @return True - If value was found and successfully parsed as a double.
	 * <br>False - If value mapping to given key was not found or could not be interpreted as a float.
	 */
	private boolean parseStringInfoField_toFloat(String fieldkey)
	{
		String rawValue[] = infoZ.remove(fieldkey);
		if (rawValue == null) return false;
		try
		{
			double vals[] = new double[rawValue.length];
			for (int i = 0; i < rawValue.length; i++)
			{
				vals[i] = Double.parseDouble(rawValue[i]);
			}
			infoF.put(fieldkey, vals);
		}
		catch (NumberFormatException e)
		{
			//e.printStackTrace();
			infoZ.put(fieldkey, rawValue);
			return false;
		}
		
		return true;
	}
	
	/* --- Getters --- */
	
	/**
	 * Get whether this variant has been flagged confirmed by checking against a truth set.
	 * @return True - If flagged confirmed
	 * <br> False - Otherwise
	 */
	public boolean isConfirmed()
	{
		return confirmed;
	}
	
	/**
	 * Get the chromosome this variant is located on as a Contig object.
	 * If the chromosome is matched to a contig in a known build, additional information
	 * such as the length and alternative names can be obtained from this object as well.
	 * @return The chromosome this variant is mapped to.
	 */
	public Contig getChromosome()
	{
		return chromosome;
	}
	
	/**
	 * Get the chromosome this variant is located on as a String.
	 * @return Default String representation of chromosome.
	 */
	public String getChromosomeName()
	{
		return chromosome.getUDPName();
	}
	
	/**
	 * Get a collection of chromosomes that this variant touches. If this variant
	 * is not a translocation or interchromosomal fusion, this method will return a
	 * collection with only a single item: the same Contig that getChromosome() returns.
	 * @return A collection containing all chromosomes involved in this variant.
	 */
	public Collection<Contig> getAllChromosomes()
	{
		Set<Contig> cset = new HashSet<Contig>();
		cset.add(this.getChromosome());
		return cset;
	}
	
	/**
	 * Get the 1-based coordinate of the first base of this variant.
	 * @return An integer denoting the coordinate relative to the chromosome of the variant.
	 */
	public int getPosition()
	{
		return position;
	}
	
	/**
	 * Get the ID or name of the variant, if it has been set. The default value of this field,
	 * following VCF style, is "."
	 * This field should at no point be set to an empty string or null, so if this method returns such,
	 * it might be prudent to reset the ID to something more amiable to the VCF format.
	 * @return Variant ID as a string.
	 */
	public String getVarID()
	{
		return variantID;
	}
	
	/**
	 * Get the reference allele as a string.
	 * @return Reference allele string.
	 */
	public String getRefAllele()
	{
		return refAllele;
	}
	
	/**
	 * Get the number of alternate alleles recorded for this variant.
	 * @return Alternate allele count
	 */
	public int countAltAlleles()
	{
		if (altAlleles == null) return 0;
		return altAlleles.size();
	}
	
	/**
	 * Return an array listing all alternate alleles for this variant.
	 * This method returns an array representation of the internal list, not a reference.
	 * @return All alternate alleles.
	 * <br> null If the list is null or empty.
	 */
	public String[] getAllAltAlleles()
	{
		if (altAlleles == null) return null;
		if (altAlleles.isEmpty()) return null;
		String[] acopy = new String[altAlleles.size()];
		return altAlleles.toArray(acopy);
	}
	
	/**
	 * Get a particular alternate allele from its zero-based index, that is,
	 * its index in the array of just alt alleles. This is NOT the general
	 * allele index. The general allele index reserves zero for the ref allele.
	 * @param altIndex Index of alternate allele desired.
	 * @return Alternate allele string if index is valid.
	 * <br> null if index is invalid.
	 */
	public String getAltAllele(int altIndex)
	{
		if (altAlleles == null) return null;
		if (altIndex < 0) return null;
		if (altIndex >= altAlleles.size()) return null;
		return altAlleles.get(altIndex);
	}
	
	/**
	 * Get a particular allele from its zero-based index. This index
	 * should be the same as the one referred to in a standard genotype field.
	 * An index of zero will return the ref allele.
	 * @param index General index of allele desired.
	 * @return Allele string if index is valid.
	 * <br> null if index is invalid.
	 */
	public String getAllele(int index)
	{
		if (index < 0) return null;
		if (index == 0) return getRefAllele();
		index--;
		if (altAlleles == null) return null;
		if (index >= altAlleles.size()) return null;
		return altAlleles.get(index);
	}
	
	/**
	 * Get the phred quality score for this variant call.
	 * The default unset value is -1.
	 * @return The quality score, if set. -1 if it was never set.
	 */
	public double getQuality()
	{
		return quality;
	}

	/**
	 * Get whether this variant has any filters listed.
	 * If this variant was read from a VCF, if this method returns true,
	 * that probably means that the variant failed one or more filters.
	 * @return True - If this instance has a filter list stored.
	 * <br> False - If this instance does not have a filter list. Either it was
	 * tested against no filters, or it passed all filters.
	 */
	public boolean hasFilters()
	{
		if (passedFilters) return true;
		else return (filters != null);
	}
	
	/**
	 * Get whether this variant has passed all filters.
	 * @return True - If this variant has passed all filters. This can be
	 * marked internally or may be marked in the VCF FILTER field as "PASS".
	 * <br>False - Otherwise.
	 */
	public boolean passedAllFilters()
	{
		return passedFilters;
	}
	
	/**
	 * Get a list (array) of strings denoting all filters marked for this
	 * variant, usually a list of the filters this variant has failed.
	 * @return An array (a copy) of the failed filters list.
	 */
	public String[] getFiltersFailed()
	{
		if (filters == null) return null;
		if (filters.isEmpty()) return null;
		String[] acopy = new String[filters.size()];
		return filters.toArray(acopy);
	}
	
	/**
	 * Get the value(s) corresponding to an INFO field annotation with the provided key.
	 * This method returns a string representation of the values, regardless of their supposed
	 * type. This method searches all fields, not just string fields.
	 * @param key INFO field key string.
	 * @return String array containing found values. If there is only one value recorded, an array
	 * will still be returned, just with a length of one.
	 * <br> This method returns null if the value mapped to the key is null or if the key
	 * could not be found.
	 */
	public String[] getInfoEntry(String key)
	{
		//return this.getInfoEntryDirect(key);
		String[] vals = getInfoEntryDirect(key);
		if (vals != null) return vals;
		if(fieldMap == null) populateFieldMap();
		Field f = fieldMap.get(key);
		if (f == null) return null;
		return f.getAll();
	}
	
	/**
	 * Get the first value corresponding to an INFO field annotation with the provided key.
	 * This method returns a string representation of the values, regardless of their supposed
	 * type. This method searches all fields, not just string fields.
	 * <br>If this field only has one value, then this method should return that value.
	 * @param key INFO field key string.
	 * @return String corresponding to the first value in this field.
	 * <br> This method returns null if the value mapped to the key is null or if the key
	 * could not be found.
	 */
	public String getSingleInfoEntry(String key)
	{
		//return getSingleInfoEntryDirect(key);
		String val = getSingleInfoEntryDirect(key);
		if (val != null) return val;
		if(fieldMap == null) populateFieldMap();
		Field f = fieldMap.get(key);
		if (f == null) return null;
		return f.getFirst();
	}
	
	/**
	 * Does the same thing as getInfoEntry, but always retrieves from the Variant ancestor maps,
	 * ignoring any downstream overrides.
	 * Turns out the overrides were trying to look for info in the maps and getting pointed back to themselves.
	 * Stack overflow! Oops!
	 * @param key INFO field key string.
	 * @return String array containing found values. If there is only one value recorded, an array
	 * will still be returned, just with a length of one.
	 * <br> This method returns null if the value mapped to the key is null or if the key
	 * could not be found.
	 */
	protected final String[] getInfoEntryDirect(String key)
	{
		if (infoZ.containsKey(key)) return infoZ.get(key);
		if (infoB.contains(key))
		{
			String[] flag = {""};
			return flag;
		}
		if (infoI.containsKey(key))
		{
			int[] val = infoI.get(key);
			String[] vals = new String[val.length];
			for (int i = 0; i < val.length; i++)
			{
				vals[i] = Integer.toString(val[i]);
			}
			return vals;
		}
		if (infoF.containsKey(key))
		{
			double[] val = infoF.get(key);
			String[] vals = new String[val.length];
			for (int i = 0; i < val.length; i++)
			{
				vals[i] = Double.toString(val[i]);
			}
			return vals;
		}
		return null;
	}
	
	/**
	 * Does the same thing as getSingleInfoEntry, but always retrieves from the Variant ancestor maps,
	 * ignoring any downstream overrides.
	 * Turns out the overrides were trying to look for info in the maps and getting pointed back to themselves.
	 * Stack overflow! Oops!
	 * @param key INFO field key string.
	 * @return String corresponding to the first value in this field.
	 * <br> This method returns null if the value mapped to the key is null or if the key
	 * could not be found.
	 */
	protected final String getSingleInfoEntryDirect(String key)
	{
		String[] allval = getInfoEntryDirect(key);
		if (allval == null) return null;
		if (allval.length < 1) return null;
		return allval[0];
	}
	
	/**
	 * Get the value(s) corresponding to a String type INFO field annotation with the provided key.
	 * This method only finds values whose field is denoted as a String field or are unparsed as
	 * something else.
	 * @param fieldKey INFO field key string. Must correlate to a String type or unparsed field.
	 * @return String array containing found values. If there is only one value recorded, an array
	 * will still be returned, just with a length of one.
	 * <br> This method returns null if the value mapped to the key is null or if the key does
	 * not map to any string field value.
	 */
	public final String[] getStringInfoEntry(String fieldKey)
	{
		return infoZ.get(fieldKey);
	}
	
	/**
	 * Get the first value corresponding to a String type INFO field annotation with the provided key.
	 * This method only finds values whose field is denoted as a String field or are unparsed as
	 * something else.
	 * <br>If this field only has one value, then this method should return that value.
	 * @param key INFO field key string. Must correlate to a String type or unparsed field.
	 * @return String corresponding to the first value in this field.
	 * <br> This method returns null if the value mapped to the key is null or if the key does
	 * not map to any string field value.
	 */
	public String getSingleStringInfoEntry(String key)
	{
		String[] allval = getStringInfoEntry(key);
		if (allval == null) return null;
		if (allval.length < 1) return null;
		return allval[0];
	}
	
	/**
	 * Get the value(s) corresponding to an Integer type INFO field annotation with the provided key.
	 * This method only searches fields that have been explicitly denoted as Integer fields.
	 * @param fieldKey INFO field key string. Must correlate to an Integer type field.
	 * @return int array containing found values. If there is only one value recorded, an array
	 * will still be returned, just with a length of one.
	 * <br> This method returns null if the value mapped to the key is null or if the key does
	 * not map to any int field value.
	 */
	public final int[] getIntInfoEntry(String fieldKey)
	{
		//Check string fields if not already parse into ints and floats
		int[] vals = infoI.get(fieldKey);
		if (vals != null) return vals;
		
		if (parseStringInfoField_toInt(fieldKey))
		{
			return infoI.get(fieldKey);
		}
		
		return null;
	}
	
	/**
	 * Get the first value corresponding to an Integer type INFO field annotation with the provided key.
	 * This method only searches fields that have been explicitly denoted as Integer fields.
	 * <br>If this field only has one value, then this method should return that value.
	 * @param key INFO field key string. Must correlate to an Integer type field.
	 * @return The first (int) value in this field.
	 * <br> This method returns -1 if the key is null or the key was not mapped to any Integer
	 * field value.
	 */
	public int getSingleIntInfoEntry(String key)
	{
		int[] allval = getIntInfoEntry(key);
		if (allval == null) return -1;
		if (allval.length < 1) return -1;
		return allval[0];
	}
	
	/**
	 * Get the value(s) corresponding to a Float type INFO field annotation with the provided key.
	 * This method only searches fields that have been explicitly denoted as Float fields.
	 * @param fieldKey INFO field key string. Must correlate to an Float type field.
	 * @return double array containing found values. If there is only one value recorded, an array
	 * will still be returned, just with a length of one.
	 * <br> This method returns null if the value mapped to the key is null or if the key does
	 * not map to any Float field value.
	 */
	public final double[] getFloatInfoEntry(String fieldKey)
	{
		double[] vals = infoF.get(fieldKey);
		if (vals != null) return vals;
		
		if (parseStringInfoField_toFloat(fieldKey))
		{
			return infoF.get(fieldKey);
		}
		
		return null;
	}

	/**
	 * Get the first value corresponding to a Float type INFO field annotation with the provided key.
	 * This method only searches fields that have been explicitly denoted as Float fields.
	 * <br>If this field only has one value, then this method should return that value.
	 * @param key INFO field key string. Must correlate to an Float type field.
	 * @return The first (double) value in this field.
	 * <br> This method returns NaN if the key is null or the key was not mapped to any Float
	 * field value.
	 */
	public double getSingleFloatInfoEntry(String key)
	{
		double[] allval = getFloatInfoEntry(key);
		if (allval == null) return Double.NaN;
		if (allval.length < 1) return Double.NaN;
		return allval[0];
	}
	
	/**
	 * Check for an INFO field flag. This method checks whether this variant has been flagged
	 * using the provided key.
	 * @param fieldKey The INFO field key corresponding to the flag. In a VCF file, the presence
	 * of this key in the INFO field flags it as true, whereas the absence flags it as false.
	 * @return True - If flag was found.
	 * <br>False - If flag was not found.
	 */
	public boolean getInfoFlag(String fieldKey)
	{
		//for (String k : infoB) if (k.equals(fieldKey)) return true;
		//return false;
		return infoB.contains(fieldKey);
	}
	
	/**
	 * Get an ordered array of the genotype keys requested for all genotypes of this variant.
	 * @return The list of genotype fields in the order they should appear in a VCF file
	 * as an array of strings. The strings are the genotype field keys.
	 */
	public String[] getOrderedGenotypeFields()
	{
		if (genoFields == null) return null;
		if (genoFields.isEmpty()) return null;
		String[] acopy = new String[genoFields.size()];
		return genoFields.toArray(acopy);
	}
	
	/**
	 * Retrieve the variant Genotype for the sample with the provided name/ID string.
	 * @param sampleName Name or ID of sample genotype requested.
	 * @return Genotype object containing genotype information if found for requested sample
	 * and variant. 
	 * <br> null If there is no genotype recorded for that sample.
	 */
	public Genotype getSampleGenotype(String sampleName)
	{
		return genotypes.get(sampleName);
	}
	
	/**
	 * Generate a set of all of the INFO field keys currently being used to annotate this
	 * variant. This set includes keys corresponding to all types of INFO fields.
	 * @return A set of all INFO keys currently annotating this variant.
	 */
	public Set<String> getAllInfoKeys()
	{
		Set<String> allKeys = new HashSet<String>();
		allKeys.addAll(infoB);
		allKeys.addAll(infoZ.keySet());
		allKeys.addAll(infoI.keySet());
		allKeys.addAll(infoF.keySet());
		return allKeys;
	}

	/**
	 * Check whether this variant has an annotation with the provided INFO field key.
	 * @param fieldKey INFO field key to check.
	 * @return True - If annotation was found. This method returns true as long as the key
	 * is mapped to a value, even if this value is default or null.
	 * <br> False - If the provided key was not found mapped to any value.
	 */
	public boolean hasInfoField(String fieldKey)
	{
		if (infoB.contains(fieldKey)) return true;
		if (infoZ.containsKey(fieldKey)) return true;
		if (infoI.containsKey(fieldKey)) return true;
		if (infoF.containsKey(fieldKey)) return true;
		return false;
	}
	
	/**
	 * Get the difference in length between the reference allele and the alternate
	 * allele closest in length to the reference allele. This method should always
	 * return a positive number, regardless of whether the reference or the chosen
	 * alternate is longer.
	 * @return The smallest length difference between the ref allele and each alt allele.
	 * <br> 0 if there are no alt alleles.
	 */
	public int getSmallestAbsoluteLength()
	{
		int min = Integer.MAX_VALUE;
		int nAlt = this.countAltAlleles();
		if (nAlt < 1) return 0;
		for (int i = 0; i < nAlt; i++)
		{
			int len = getAltAllele(i).length() - getRefAllele().length();
			int abslen = Math.abs(len) + 1;
			if (abslen < min) min = abslen;
		}
		return min;
	}
	
	/**
	 * Get the difference in length between the reference allele and the alternate
	 * allele with the most different length. This method should always
	 * return a positive number, regardless of whether the reference or the chosen
	 * alternate is longer.
	 * @return The largest length difference between the ref allele and each alt allele.
	 * <br> 0 if there are no alt alleles.
	 */
	public int getLargestAbsoluteLength()
	{
		int max = 0;
		int nAlt = this.countAltAlleles();
		if (nAlt < 1) return 0;
		for (int i = 0; i < nAlt; i++)
		{
			int len = getAltAllele(i).length() - getRefAllele().length();
			int abslen = Math.abs(len) + 1;
			if (abslen > max) max = abslen;
		}
		return max;
	}
	
	/**
	 * Get a set containing the name (key) of every sample there is a genotype for
	 * in this variant record.
	 * @return Set of strings representing the names of all samples genotyped for this variant.
	 */
	protected Set<String> getAllGenotypedSamples()
	{
		return genotypes.keySet();
	}
	
	/**
	 * Get the GFUNC field value directly from the INFO map. If it isn't in the INFO map,
	 * this function will return null, even if the function is noted in the instance variable.
	 * @return Gene function string (from bioisvtools, following refSeq guidelines?) or null
	 * if field is not stored as an INFO entry.
	 */
	public String getGeneFuncINFO()
	{
		return getSingleInfoEntryDirect(GeneSet.INFODEF_INFO_GFUNC.getKey());
	}
	
	/**
	 * Get the location effect of the variant as a GeneFunc enum. First, the structural variant's
	 * instance variable is checked, then if nothing is there, the INFO map is checked.
	 * @return A GeneFunc enum if the variant location effect is noted, null otherwise.
	 */
	public GeneFunc getGeneFunction()
	{
		if (function != null) return function;
		String f = getGeneFuncINFO();
		if (f == null) return null;
		removeInfoField(GeneSet.INFODEF_INFO_GFUNC.getKey());
		function = GeneFunc.getFunction(f);
		return function;
	}
	
	/**
	 * Check whether the provided string matches a field
	 * marked as "supporting" this variant.
	 * @param s String to check support
	 * @return True if this variant has a support record denoted
	 * by the provided string. False otherwise.
	 */
	public boolean supportMarked(String s)
	{
		if(support_set == null) return false;
		return support_set.contains(s);
	}
	
	/**
	 * Get the candidate flags map (used by the VarTable parser as an
	 * intermediate before Candidates are created) for this Variant.
	 * @return Map of int transcriptUIDs to CandidateFlags objects, if 
	 * in use. Null otherwise.
	 */
	public Map<Integer, CandidateFlags> getCandidateFlagsMap()
	{
		return this.cflagsMap;
	}
	
	/* --- Setters --- */
	
	/**
	 * Set the variant chromosome.
	 * <br>If not using a GenomeBuild, then simply create a new Contig object
	 * and set the UDPname of the Contig to the desired chromosome name.
	 * @param chrom New chromosome.
	 * <br>If this argument is null, the method will return without changing anything.
	 */
 	public void setChromosome(Contig chrom)
	{
		if (chrom == null) return;
		chromosome = chrom;
	}
	
	/**
	 * Set the chromosome position of this variant.
	 * <br> One-based coordinate system.
	 * <br> Unset value is -1.
	 * <br> Must provide a positive number. Negative numbers (except -1) 
	 * will cause the method to return without doing anything.
	 * @param pos New position. Must be a positive number or -1.
	 */
	public void setPosition(int pos)
	{
		if (pos < -1) return;
		position = pos;
	}

	/**
	 * Set the variant name/ID to the given string.
	 * @param ID String to set variant name/ID to.
	 */
	public void setVariantName(String ID)
	{
		this.variantID = ID;
	}
	
	/**
	 * Set the reference allele to the given string.
	 * <br>Null strings are invalid.
	 * <br>Default value is "."
	 * @param ref String to set reference allele to.
	 */
	public void setRefAllele(String ref)
	{
		if (ref == null) return;
		if (ref.isEmpty()) return;
		refAllele = ref;
	}
	
	/**
	 * Add an alternate allele with the given string.
	 * <br>Note: Adding too many new alleles may cause memory reallocation due to the fact
	 * that the list is backed by an array.
	 * @param alt String representation of allele to add.
	 */
	public void addAltAllele(String alt)
	{
		if (altAlleles == null)
		{
			altAlleles = new ArrayList<String>(4);
		}
		altAlleles.add(alt);
	}
	
	/**
	 * Remove an alt allele from the variant using its index.
	 * @param index Index of allele to remove.
	 * @throws IndexOutOfBoundsException If index is invalid.
	 */
	public void removeAltAllele(int index)
	{
		altAlleles.remove(index);
	}
	
	/**
	 * Delete all alt alleles for this variant.
	 */
	public void clearAltAlleles()
	{
		altAlleles.clear();
	}
	
	/**
	 * Directly set the phred scaled quality score.
	 * -1 is used to represent an unset value.
	 * @param qual The quality score to set.
	 */
	public void setQuality(double qual)
	{
		quality = qual;
	}
	
	/**
	 * Set the filters passed flag. If this flag is set, it is assumed
	 * that this variant has passed all filters, regardless of whether
	 * the failed filter list is populated.
	 * @param b Filters passed flag
	 */
	public void setFilterPass(boolean b)
	{
		passedFilters = b;
	}
	
	/**
	 * Add the name of a failed filter to the end of the list of filters failed.
	 * If the filters passed flag is set, the list of failed filters is summarily ignored.
	 * @param filterName Name of the filter failed by this variant.
	 */
	public void addFailedFilter(String filterName)
	{
		if (filterName == null) return;
		if (filterName.isEmpty()) return;
		filters.add(filterName);
	}
	
	/**
	 * Search for a filter by the provided name and remove it if found.
	 * @param filterName Name of filter to remove.
	 */
	public void removeFailedFilter(String filterName)
	{
		filters.remove(filterName);
	}
	
	/**
	 * Clear the list of failed filters.
	 */
	public void clearFilterList()
	{
		filters.clear();
	}
	
	/**
	 * Reparse all INFO fields known to be numerical from strings by providing a map of keys of interest
	 * and their INFO field definitions.
	 * @param infoDefs A map of keys to reclassify and the INFO field definitions corresponding to them. The
	 * definitions tell the parser what type to attempt to read the field values as.
	 */
	public void defineInfoFields(Map<String, InfoDefinition> infoDefs)
	{
		Set<String> infoStrings = infoDefs.keySet();
		for (String key : infoStrings)
		{
			if (!infoZ.containsKey(key)) continue;
			InfoDefinition def = infoDefs.get(key);
			if (def.getType() == VariantPool.INFODEF_INT) parseStringInfoField_toInt(key);
			else if (def.getType() == VariantPool.INFODEF_FLOAT) parseStringInfoField_toFloat(key);
		}
	}

	/**
	 * Add an INFO flag with the provided key.
	 * @param key INFO field key for the desired flag to add.
	 */
	public void addInfoFlag(String key)
	{
		if (key == null || key.isEmpty()) return;
		infoB.add(key);
	}
	
	/**
	 * Remove an INFO flag with the given key, if it is present.
	 * <br>IMPORTANT: This method only removes boolean flags. It will not
	 * check for INFO field annotations of any other type!
	 * @param key The INFO field key of the flag to remove.
	 */
	public void removeInfoFlag(String key)
	{
		if (key == null || key.isEmpty()) return;
		infoB.remove(key);
	}

	/**
	 * Remove an INFO field annotation with the given key.
	 * <br>This method can remove fields of any type.
	 * @param key The INFO field key of the annotation to remove.
	 */
	public void removeInfoField(String key)
	{
		if (key == null || key.isEmpty()) return;
		infoB.remove(key);
		infoZ.remove(key);
		infoI.remove(key);
		infoF.remove(key);
	}
	
	/**
	 * Add an INFO annotation as a String type field.
	 * Annotation will be added as a String type field regardless of what it
	 * is intended to be. InfoDefinition is required to determine its intended type.
	 * @param key INFO field key.
	 * @param value Single value to set for the field with the provided key.
	 */
	public void addInfoField(String key, String value)
	{
		if (key == null || key.isEmpty()) return;
		String[] vals = new String[1];
		vals[0] = value;
		infoZ.put(key, vals);
	}	
	
	/**
	 * Add an INFO annotation as a String type field.
	 * Annotation will be added as a String type field regardless of what it
	 * is intended to be. InfoDefinition is required to determine its intended type.
	 * @param key INFO field key.
	 * @param values Array of values to set for the field with the provided key.
	 */
	public void addInfoField(String key, String[] values)
	{
		if (key == null || key.isEmpty()) return;
		infoZ.put(key, values);
	}
	
	/**
	 * Add a single value INFO annotation and set key and type according
	 * to the provided InfoDefinition.
	 * @param value Unparsed string representation of value to use.
	 * @param def Definition of INFO field. Must include key and type.
	 */
	public void addInfoField(String value, InfoDefinition def)
	{
		if (def == null) return;
		String[] vals = new String[1];
		vals[0] = value;
		addInfoField(vals, def);
	}
	
	/**
	 * Add a multi value INFO annotation and set key and type according
	 * to the provided InfoDefinition.
	 * If parsing one or more values fails, annotation will be added
	 * as a string type.
	 * @param values Unparsed string representations of values to use.
	 * @param def Definition of INFO field. Must include key and type.
	 */
	public void addInfoField(String[] values, InfoDefinition def)
	{
		if (def == null || values == null || values.length < 1) return;
		String key = def.getKey();
		int type = def.getType();
		if (type == VariantPool.INFODEF_INT)
		{
			int[] iVals = new int[values.length];
			try
			{
				for (int i = 0; i < values.length; i++)
				{
					iVals[i] = Integer.parseInt(values[i]);
				}
				infoI.put(key, iVals);
				return;
			}
			catch (NumberFormatException e)
			{
				infoZ.put(key, values);
				return;
			}
		}
		else if (type == VariantPool.INFODEF_FLOAT)
		{
			double[] fVals = new double[values.length];
			try
			{
				for (int i = 0; i < values.length; i++)
				{
					fVals[i] = Double.parseDouble(values[i]);
				}
				infoF.put(key, fVals);
				return;
			}
			catch (NumberFormatException e)
			{
				infoZ.put(key, values);
				return;
			}
		}
		infoZ.put(key, values);
	}
	
	/**
	 * Add an INFO annotation as an Integer type field.
	 * @param key INFO field key.
	 * @param value Single value to set for the field with the provided key.
	 */
	public void addInfoField(String key, int value)
	{
		if (key == null || key.isEmpty()) return;
		int[] vals = new int[1];
		vals[0] = value;
		infoI.put(key, vals);
	}
	
	/**
	 * Add an INFO annotation as an Integer type field.
	 * @param key INFO field key.
	 * @param values Array of values to set for the field with the provided key.
	 */
	public void addInfoField(String key, int[] values)
	{
		if (key == null || key.isEmpty()) return;
		infoI.put(key, values);
	}
	
	/**
	 * Add an INFO annotation as a Float type field.
	 * @param key INFO field key.
	 * @param value Single value to set for the field with the provided key.
	 */
	public void addInfoField(String key, double value)
	{
		if (key == null || key.isEmpty()) return;
		double[] vals = new double[1];
		vals[0] = value;
		infoF.put(key, vals);
	}
	
	/**
	 * Add an INFO annotation as a Float type field.
	 * @param key INFO field key.
	 * @param values Array of values to set for the field with the provided key.
	 */
	public void addInfoField(String key, double[] values)
	{
		if (key == null || key.isEmpty()) return;
		infoF.put(key, values);
	}
	
	/**
	 * Add information about a sample genotype.
	 * If a sample with the given name already has a recorded variant for
	 * this genotype, this method will overwrite that record with the new one.
	 * @param sampleName Name of the sample to make genotype to. If this argument is
	 * null or empty, the method will return without adding anything.
	 * @param gt Genotype object to record for sample.
	 */
	public void addGenotype(String sampleName, Genotype gt)
	{
		if (sampleName == null) return;
		if (sampleName.isEmpty()) return;
		genotypes.put(sampleName, gt);
	}
	
	/**
	 * Add a genotype field to the end of the list of genotype fields
	 * to include should this variant be written out to a file.
	 * @param fieldKey Key of the genotype field to include.
	 */
	public void addGenotypeField(String fieldKey)
	{
		if (fieldKey == null) return;
		if (fieldKey.isEmpty()) return;
		genoFields.add(fieldKey);
	}
	
	/**
	 * Search the genotype field list and remove the field from the list
	 * if it is present. Note that removing it from the ordered list will make it
	 * lose its place.
	 * @param fieldKey Key of the genotype field to remove.
	 */
	public void removeGenotypeField(String fieldKey)
	{
		genoFields.remove(fieldKey);
	}
	
	/**
	 * Clear the list of genotype fields.
	 */
	public void clearGenotypeFields()
	{
		genoFields.clear();
	}

	/**
	 * Flag this variant as confirmed.
	 */
	public void confirm()
	{
		this.confirmed = true;
	}
	
	/**
	 * Reset variant confirmed flag to false.
	 */
	public void unconfirm()
	{
		this.confirmed = false;
	}
	
	/**
	 * Set the location effect of the variant as a GeneFunc enum. First, the variant's
	 * instance variable is set, then an INFO field is checked for and removed.
	 */
	public void setGeneFunction(GeneFunc eff)
	{
		function = eff;
		removeInfoField(GeneSet.INFODEF_INFO_GFUNC.getKey());
	}
	
	/**
	 * Add a "support" mark correlated to the provided string.
	 * @param supportString String to mark as supporting evidence for
	 * this variant.
	 */
	public void addSupportMark(String supportString)
	{
		if (this.support_set == null) this.support_set = new HashSet<String>();
		this.support_set.add(supportString);
	}
	
	/**
	 * Remove the "support" mark correlated to the provided string.
	 * @param supportString String to remove as supporting evidence for
	 * this variant.
	 */
	public void removeSupportMark(String supportString)
	{
		if (this.support_set == null) return;
		this.support_set.remove(supportString);
	}
	
	/**
	 * Note a transcriptUID, CandidateFlags pairing.
	 * @param transcriptUID Integer UID of gene transcript associated with the
	 * candidate the flags apply to.
	 * @param cf CandidateFlags object containing the flags.
	 */
	public void addCandidateFlagsNotation(int transcriptUID, CandidateFlags cf)
	{
		if (cflagsMap == null) cflagsMap = new HashMap<Integer, CandidateFlags>();
		cflagsMap.put(transcriptUID, cf);
	}
	
	/**
	 * Delete the CandidateFlags map from this Variant to free up memory when
	 * it is no longer in use.
	 */
	public void clearCandidateFlags()
	{
		cflagsMap = null;
	}
	
	/* --- Comparing --- */
	
	/**
	 * Get whether this variant is on/ associated with the given chromosome
	 * described by the provided string.
	 * @param chrom Chromosome to check against.
	 * @return True - If this variant is on/ associated with chrom.
	 * <br>False - If it is not.
	 */
	public boolean isOnChromosome(Contig chrom)
	{
		return chrom.equals(this.getChromosome());
	}
	
	/**
	 * Get whether the chromosome string provided is a "compound" - that
	 * is to say, whether it is a standin representing more than one chromosome.
	 * @param c Chromosome string to check.
	 * @return True - If chromosome string represents multiple chromosomes.
	 * False - If chromosome string represents only one chromosome.
	 */
	private static boolean chrIsCompound(String c)
	{
		return (c.indexOf(':') >= 0);
	}
	
	/**
	 * Get the sorting weight of an alphabetical chromosome string
	 * (such as a sex chromosome or mitochondrial chromosome).
	 * @param c Chromosome string to check.
	 * @return Relative weight of chromosome string used for sorting.
	 */
	private static int chrAlphaStringWeight(String c)
	{
		//Higher number means lower priority
			/*
			 * Compound = -1
			 * X = 1
			 * Y = 2
			 * W = 3
			 * Z = 4
			 * M/MT = 5
			 * Other = 6
			 */
		if (chrIsCompound(c)) return -1;
		if (c.equals("X")) return 1;
		if (c.equals("Y")) return 2;
		if (c.equals("W")) return 3;
		if (c.equals("Z")) return 4;
		if (c.equals("M")) return 5;
		if (c.equals("MT")) return 6;
		return 7;
	}
	
	/**
	 * Compare two chromosome strings to determine which should come
	 * first in sorting.
	 * @param c1 First chromosome to check
	 * @param c2 Second chromosome to check
	 * @return A positive number if c1 is greater (comes after c2 in order).
	 * <br>A negative number if c2 is greater.
	 * <br> Zero if they are equal.
	 */
	public static int compareChromosomeStrings(String c1, String c2)
	{
		//Return positive # if c1 is greater, negative # if c2 is greater
		if (c1.equals(c2)) return 0;
		boolean c1c = chrIsCompound(c1);
		boolean c2c = chrIsCompound(c2);
		if (c1c && !c2c) return 1;
		if (!c1c && c2c) return -1;
		if (c1c && c2c)
		{
			String[] c1p = c1.split(":");
			String[] c2p = c2.split(":");
			int compare1 = compareChromosomeStrings(c1p[0], c2p[0]);
			if (compare1 != 0) return compare1;
			int compare2 = compareChromosomeStrings(c1p[1], c2p[1]);
			return compare2;
		}
		try
		{
			int cn1 = Integer.parseInt(c1);
			try
			{
				int cn2 = Integer.parseInt(c2);
				//Both are integers
				return cn1 - cn2; //Negative if cn2 is bigger, positive if cn1 is bigger
			}
			catch (NumberFormatException e2)
			{
				//1 is an integer, but 2 is not
				return -1; //Integers < alpha (integers come first)
			}
			
		}
		catch (NumberFormatException e1)
		{
			try
			{
				Integer.parseInt(c2);
				//c2 is integer, but c1 is not
				return 1;
			}
			catch (NumberFormatException e2)
			{
				//Neither are integers
				int w1 = chrAlphaStringWeight(c1);
				int w2 = chrAlphaStringWeight(c2);
				if (w1 > 6 && w2 > 6) return c1.compareTo(c2);
				else return w1 - w2;
			}
		}
	}

	/**
	 * Get an enumeration (int) representing the type of the chromosome
	 * represented by the given string.
	 * @param c Chromosome to check.
	 * @return An int representing the chromosome type.
	 * @see Variant.CHROMOSOME_UNKNOWN
	 * <br>Variant.CHROMOSOME_AUTOSOME
	 * <br>Variant.CHROMOSOME_SEXCHROMOSOME
	 * <br>Variant.CHROMOSOME_MITOCHONDRIAL
	 * <br>Variant.CHROMOSOME_CONTIG
	 * <br>Variant.CHROMOSOME_COMPOUND
	 */
	public static int getChromosomeType(String c)
	{
		try
		{
			Integer.parseInt(c);
			return Variant.CHROMOSOME_AUTOSOME;
		}
		catch(NumberFormatException e)
		{
			int w = chrAlphaStringWeight(c);
			switch (w)
			{
			case -1: return Variant.CHROMOSOME_COMPOUND;
			case 1: return Variant.CHROMOSOME_SEXCHROMOSOME;
			case 2: return Variant.CHROMOSOME_SEXCHROMOSOME;
			case 3: return Variant.CHROMOSOME_SEXCHROMOSOME;
			case 4: return Variant.CHROMOSOME_SEXCHROMOSOME;
			case 5: return Variant.CHROMOSOME_MITOCHONDRIAL;
			case 6: return Variant.CHROMOSOME_CONTIG;
			}
		}
		return Variant.CHROMOSOME_UNKNOWN;
	}

	/**
	 * Get whether the chromosome associated with this variant is a single,
	 * known chromosome.
	 * @return True - If the chromosome associated with this variant is an autosome,
	 * sex chromosome, or mitochondrial chromosome.
	 * <br>False - If it is a contig, compound, or unknown.
	 */
	public boolean isChromosome()
	{
		int cType = this.getChromosome().getType();
		switch (cType)
		{
		case Contig.SORTCLASS_AUTOSOME: return true;
		case Contig.SORTCLASS_SEXCHROM: return true;
		case Contig.SORTCLASS_CONTIG: return false;
		case Contig.SORTCLASS_MITO: return true;
		case Contig.SORTCLASS_COMPOUND: return false;
		case Contig.SORTCLASS_UNKNOWN: return false;
		}
		return false;
	}
	
	/**
	 * Get whether the chromosome associated with this variant is a contig.
	 * @return True - If the chromosome coordinate refers to a contig.
	 * <br>False - If the chromosome coordinate refers to a known chromosome,
	 * a compound, or an unknown.
	 */
	public boolean isContig()
	{
		int cType = this.getChromosome().getType();
		switch (cType)
		{
		case Contig.SORTCLASS_AUTOSOME: return false;
		case Contig.SORTCLASS_SEXCHROM: return false;
		case Contig.SORTCLASS_CONTIG: return true;
		case Contig.SORTCLASS_MITO: return false;
		case Contig.SORTCLASS_COMPOUND: return false;
		case Contig.SORTCLASS_UNKNOWN: return false;
		}
		return false;
	}
	
	public boolean equals(Object o)
	{
		if (o == this) return true;
		if (o == null) return false;
		if (!(o instanceof Variant)) return false;
		Variant v = (Variant)o;
		if (!this.getChromosome().equals(v.getChromosome())) return false;
		if (this.getPosition() != v.getPosition()) return false;
		if (!this.getRefAllele().equals(v.getRefAllele())) return false;
		return true;
	}
	
	/**
	 * Compares another variant to this variant and determines whether they might be 
	 * considered to represent the same event, even if the data backing them are not
	 * perfectly identical, and thus would not qualify as "equal."
	 * @param v Variant to compare to this one.
	 * @return True - If the two variants might be considered equivalent.
	 * <br> False - Otherwise.
	 */
	public boolean equivalent(Variant v)
	{
		return this.equals(v);
	}
		
	public int compareTo(Variant o)
	{
		if (o == null) return 1;
		if (o == this) return 0;
		int chrComp = this.getChromosome().compareTo(o.getChromosome());
		if (chrComp != 0) return chrComp;
		if (this.getPosition() != o.getPosition()) return this.getPosition() - o.getPosition();
		if (!this.getVarID().equals(o.getVarID())) return this.getVarID().compareTo(o.getVarID());
		if (!this.getRefAllele().equals(o.getRefAllele())) return this.getRefAllele().length() - o.getRefAllele().length();
		if (this.countAltAlleles() != o.countAltAlleles()) return this.countAltAlleles() - o.countAltAlleles();
		return 0;
	}
	
	/**
	 * Check to see if the variant falls in a specified genomic region.
	 * @param c Contig of region to check.
	 * @param start Start coordinate of region.
	 * @param end End coordinate of region.
	 * @param anyend Whether this variant should be considered in the region if
	 * only one if its endpoints (assuming it has multiple endpoints) is in or not.
	 * <br>True = This variant is considered inside the region if any of its endpoints
	 * are.
	 * @return True if this variant falls inside the region as defined by the provided
	 * arguments. False if not.
	 */
	public boolean inRegion(Contig c, int start, int end, boolean anyend)
	{
		if (!this.getChromosome().equals(c)) return false;
		if (this.getPosition() < start) return false;
		if (this.getPosition() > end) return false;
		return true;
	}
	
	/* --- Parsed Fields as INFO Access --- */
	
	private HashMap<String, Field> fieldMap;
	
	private interface Field
	{
		public String get();
		public String[] getAll();
		public String getFirst();
		public void set(String value);
		public boolean hasFlag();
	}
	
	private void populateFieldMap()
	{
		fieldMap = new HashMap<String, Field>();
		
		fieldMap.put(GeneSet.INFODEF_INFO_GFUNC.getKey(), new Field(){
			public String get()
			{
				GeneFunc f = getGeneFunction();
				if (f == null) return null;
				return f.toString();
			}
			public String[] getAll()
			{
				GeneFunc f = getGeneFunction();
				if (f == null) return null;
				String[] all = new String[1];
				all[0] = f.toString();
				return all;
			}
			public String getFirst()
			{
				return get();
			}			
			public void set(String value)
			{
				function = GeneFunc.getFunction(value);
				//If function is non-null, then INFO text is ignored anyway
				if (value == null) removeInfoField(GeneSet.INFODEF_INFO_GFUNC.getKey());
			}
			public boolean hasFlag()
			{
				return false;
			}
		});
	}
	
	/* --- View --- */
	
	public String toString()
	{
		String s = "";
		s += getChromosome();
		s += "::" + getPosition();
		s += "|" + getRefAllele();
		s += "->";
		if (countAltAlleles() > 0)
		{
			String[] alts = getAllAltAlleles();
			for (String a : alts)
			{
				s += a + ",";
			}
			s = s.substring(0, s.length() - 1);
		}
		else s += ".";
		return s;
	}

	/**
	 * Get a string representing the genotype of a sample for viewing. This string
	 * is not in VCF format, and only shows the allele calls.
	 * @param sampleName Name of the sample to retrieve the genotype for.
	 * @return A string representing the allele calls for the requested sample, if
	 * the sample has a genotype recorded.
	 * <br> An empty string if sample record could not be found or sample string was empty.
	 */
	public String getSampleGenotypeString(String sampleName)
	{
		//System.err.println(Thread.currentThread().getName() + " || Variant.getSampleGenotypeString || Called!");
		if (sampleName == null || sampleName.isEmpty()) return "";
		//System.err.println(Thread.currentThread().getName() + " || Variant.getSampleGenotypeString || Passed check 1...");
		//System.err.println(Thread.currentThread().getName() + " || Variant.getSampleGenotypeString || Genotypes for this variant: " + this.genotypes.size());
		Genotype g = getSampleGenotype(sampleName);
		if (g == null) return "";
		//System.err.println(Thread.currentThread().getName() + " || Variant.getSampleGenotypeString || Passed check 2...");
		int[] salleles = g.getAlleles();
		if (salleles.length <= 0) return "";
		//System.err.println(Thread.currentThread().getName() + " || Variant.getSampleGenotypeString || Passed check 3...");
		String s = "";
		for (int i = 0; i < salleles.length; i++)
		{
			int a = salleles[i];
			if (a == 0)
			{
				s += this.getRefAllele();
			}
			else if (a < 0)
			{
				s += ".";
			}
			else if (a > this.countAltAlleles())
			{
				s += "ERR";
			}
			else
			{
				s += this.getAltAllele(a - 1);
			}
			if (i < salleles.length - 1) s += " / ";
		}
		//System.err.println(Thread.currentThread().getName() + " || Variant.getSampleGenotypeString || Returning " + s);
		return s;
	}
	
	/* --- Serialization --- */
	
	/**
	 * Output the data describing this variant as a line in VCF format.
	 * @param orderedSamples List of samples to included genotype columns for in the order
	 * to be included.
	 * @param orderedInfoFields List of INFO fields to include in the order they are to be
	 * included. If a key is not found for this variant, the field will simply be skipped.
	 * @return A VCF format variant record as a String.
	 */
	public String toVCFLine(List<String> orderedSamples, List<InfoDefinition> orderedInfoFields)
	{
		String s = "";
		//CHROM
		s += chromosome.getUDPName() + "\t";
		//POS
		s += getPosition() + "\t";
		//ID
		if (!getVarID().isEmpty()) s += getVarID() + "\t";
		else s += ".\t";
		//REF
		s += getRefAllele() + "\t";
		//ALT
		if (countAltAlleles() > 0)
		{
			String[] alts = getAllAltAlleles();
			for (String a : alts)
			{
				s += a + ",";
			}
			s = s.substring(0, s.length() - 1);
		}
		else s += ".";
		s += "\t";
		//QUAL
		if (getQuality() >= 0) s += Double.toString(getQuality());
		else s += ".";
		s += "\t";
		//FILTER
		if (!hasFilters()) s += ".";
		else if (passedAllFilters()) s += "PASS";
		else
		{
			String[] filters = getFiltersFailed();
			if (filters != null)
			{
				for (String f : filters)
				{
					s += f + ";";
				}
				s = s.substring(0, s.length() - 1);
			}
			else s += ".";
		}
		s += "\t";
		
		//INFO
		if (orderedInfoFields != null)
		{
			int nFields = orderedInfoFields.size();
			int c = 0;
			for (InfoDefinition def : orderedInfoFields)
			{
				//Is i a flag?
				if (def.getType() == VariantPool.INFODEF_FLAG)
				{
					boolean flag = this.getInfoFlag(def.getKey());
					if (flag){
						s += def.getKey();
						if (c < nFields - 1) s += ";";
					}
				}
				else
				{
					String[] vals = this.getInfoEntry(def.getKey());
					if (vals != null && vals.length > 0)
					{
						s += def.getKey() + "=";
						for (int v = 0; v < vals.length; v++)
						{
							s += vals[v];
							if (v < vals.length - 1) s += ",";
						}
						if (c < nFields - 1) s += ";";
					}
				}
				c++;
			}
		}
		else s += ".";
		s += "\t";
		
		//Format field & Genotypes
		//System.err.println("Variant.toVCFLine || Now serializing genotype information... ");
		if (orderedSamples != null)
		{
			//System.err.println("Variant.toVCFLine || Sample list is not null.");
			if (!orderedSamples.isEmpty() && genoFields != null)
			{
				//System.err.println("Variant.toVCFLine || Sample list is not empty, and field list is not null.");
				String formatString = "";
				for (int i = 0; i < genoFields.size(); i++)
				{
					formatString += genoFields.get(i);
					if (i < genoFields.size() - 1) formatString += ":";
				}
				s += formatString + "\t";
				//System.err.println("Variant.toVCFLine || Format string: " + formatString);
				
				for (int i = 0; i < orderedSamples.size(); i++)
				{
					String sample = orderedSamples.get(i);
					String genotype = genotypes.get(sample).toVCFField(genoFields);
					//System.err.println("Variant.toVCFLine || Genotype string for " + sample + ": " + genotype);
					if (i < orderedSamples.size() - 1) s += genotype + "\t";
					else s += genotype;
				}
			}
		}
	
		//System.err.println("Variant.toVCFLine || Returning...");
		//System.err.println("\t" + s);
		return s;
	}
	
	/**
	 * Output the data describing this variant as a line formatted for a BED to input as a track
	 * in the UCSC genome browser.
	 * @param sample Name of the sample this track will be generated for. This is required in case genotype
	 * information needs to be accessed.
	 * @return A single line as a string (with no newline character) representing the record in BED12 format.
	 */
	public String toViewerBEDLine(String sample)
	{
		//Single block
		Genotype g = this.getSampleGenotype(sample);
		if (g == null) return null;
		double altpercent = g.getPercentAlt();
		if (altpercent <= 0.0) return null;
		String s = "";
		String cname = getChromosome().getUCSCName();
		if (cname == null) return null;
		s += cname + "\t"; //chrom
		s += (getPosition() - 1) + "\t";
		s += getPosition() + "\t";
		s += getVarID() + "\t";
		/*
		int score = getQuality();
		score *= 200;
		if (score > 1000) score = 1000;
		if (score >= 0) s += Integer.toString(score);
		else s += "1000"; //score
		*/
		//Instead, I'm going to use the score for zygosity...
		int score = (int)Math.round(altpercent * 10.0);
		s += score + "\t"; //score
		s += ".\t"; //strand
		s += (getPosition() - 1) + "\t";
		s += getPosition() + "\t";
		s += UCSCGVBED.COLOR_SNV + "\t";
		s += "1\t";
		s += "1\t";
		s += "0";
		
		return s;
	}
	
}
