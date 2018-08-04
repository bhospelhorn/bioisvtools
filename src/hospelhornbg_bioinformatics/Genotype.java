package hospelhornbg_bioinformatics;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import hospelhornbg_bioinformatics.VariantPool.InfoDefinition;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

/*
 * UPDATE NOTES
 * 
 * Initial version date: December 14, 2017 (1.0.0)
 * 
 * 1.0.0 -> 1.1.0 | January 11, 2018
 * 	Javadoc
 * 	Added some common CNV/SV fields
 * 	Added constant InfoDefinitions for common fields
 * 	Reconciled parsed instance variable fields and unparsed mapped fields
 * 	Added getters and setters
 * 	Changed internal implementation of filter list to an ArrayList from an array
 * 	Added array size control (needed for setter capability)
 * 
 * 1.1.0 -> 1.1.1 | January 25, 2018
 * 	Added a few methods and constants for easier querying on zygosity
 * 
 * 1.1.1 -> 1.1.2 | January 29, 2018
 * 	Bug fix - null pointer in set max copy number/alleles if corresponding array doesn't exist (then don't need to copy, silly!)
 * 
 * 1.1.2 -> 1.1.3 | February 9, 2018
 * 	Bug fix - missed that "|" is a regex special character, so genotypes with phasing were not parsed.
 * 	Changed default GT field to ./. instead of .
 * 
 * 1.1.3 -> 1.1.4 | March 5, 2018
 * 	Added an allele setting that takes a raw VCF GT string.
 * 
 * 1.1.4 -> 1.1.5 | April 24, 2018
 * 	Small modification to allow for it to take genotype fields where there are fewer fields than requested
 * 
 * 1.1.5 -> 1.1.6 | July 10, 2018
 * 	The alt percent calculation now ignores the "." allele (-1). It doesn't count them towards the total at all.
 * 	Also added a genotypeUnknown function to determine if all alleles are "."
 * 
 * 1.1.6 -> 1.1.7 | July 26, 2018
 * 	Added fun little "isHomozygous" function for the segregation framework.
 * 
 * 1.1.7 -> 1.1.8 | August 1, 2018
 * 	Added "hasAllele" function
 * 
 */

/*
 * Possible future improvements:
 * 	- Doesn't appear to have special capabilities for flags... might want to add that, if needed
 */

/**
 * Class to contain information on a genotype - usually for a given 
 * sample and variant.
 * @author Blythe Hospelhorn (blythe.hospelhorn@nih.gov)
 * @version 1.1.8
 * @since August 1, 2018
 */
public class Genotype {
	
	/* --- Constants --- */
	
	//From the VCF format specification version 4.2
	
	/**
	 * ID = GT
	 * <br>Number = 1
	 * <br>Type = String
	 * <br>Description = "Genotype"
	 */
	public static final InfoDefinition INFODEF_GT = new InfoDefinition("GT", VariantPool.INFODEF_STRING, "Genotype", 1);
	
	/**
	 * ID = GQ
	 * <br>Number = 1
	 * <br>Type = Integer
	 * <br>Description = "Genotype Quality"
	 */
	public static final InfoDefinition INFODEF_GQ = new InfoDefinition("GQ", VariantPool.INFODEF_INT, "Genotype Quality", 1);
	
	/**
	 * ID = DP
	 * <br>Number = 1
	 * <br>Type = Integer
	 * <br>Description = "Read Depth"
	 */
	public static final InfoDefinition INFODEF_DP = new InfoDefinition("DP", VariantPool.INFODEF_INT, "Read Depth", 1);
	
	/**
	 * ID = HQ
	 * <br>Number = 2
	 * <br>Type = Integer
	 * <br>Description = "Haplotype Quality"
	 */
	public static final InfoDefinition INFODEF_HQ = new InfoDefinition("HQ", VariantPool.INFODEF_INT, "Haplotype Quality", 2);
	
	/**
	 * ID = FT
	 * <br>Number = .
	 * <br>Type = String
	 * <br>Description = "Filters"
	 */
	public static final InfoDefinition INFODEF_FT = new InfoDefinition("FT", VariantPool.INFODEF_STRING, "Filters", VariantPool.INFODEF_NARGS_VARIABLE);
	
	/**
	 * ID = GL
	 * <br>Number = G
	 * <br>Type = Float
	 * <br>Description = "Genotype Likelihood"
	 */
	public static final InfoDefinition INFODEF_GL = new InfoDefinition("GL", VariantPool.INFODEF_FLOAT, "Genotype Likelihood", VariantPool.INFODEF_NARGS_PERGENOTYPE);
	
	/**
	 * ID = GLE
	 * <br>Number = 1
	 * <br>Type = String
	 * <br>Description = "Genotype Likelihoods of Heterogeneous Ploidy"
	 */
	public static final InfoDefinition INFODEF_GLE = new InfoDefinition("GLE", VariantPool.INFODEF_STRING, "Genotype Likelihoods of Heterogeneous Ploidy", 1);
	
	/**
	 * ID = PL
	 * <br>Number = G
	 * <br>Type = Integer
	 * <br>Description = "Phred-Scaled Genotype Likelihood"
	 */
	public static final InfoDefinition INFODEF_PL = new InfoDefinition("PL", VariantPool.INFODEF_INT, "Phred-Scaled Genotype Likelihood", VariantPool.INFODEF_NARGS_PERGENOTYPE);
	
	/**
	 * ID = GP
	 * <br>Number = G
	 * <br>Type = Float
	 * <br>Description = "Phred-Scaled Genotype Posterior Probabilities"
	 */
	public static final InfoDefinition INFODEF_GP = new InfoDefinition("GP", VariantPool.INFODEF_FLOAT, "Phred-Scaled Genotype Posterior Probabilities", VariantPool.INFODEF_NARGS_PERGENOTYPE);
	
	/**
	 * ID = PS
	 * <br>Number = 1
	 * <br>Type = Integer
	 * <br>Description = "Phase Set"
	 */
	public static final InfoDefinition INFODEF_PS = new InfoDefinition("PS", VariantPool.INFODEF_INT, "Phase Set", 1);
	
	/**
	 * ID = PQ
	 * <br>Number = 1
	 * <br>Type = Integer
	 * <br>Description = "Phasing Quality"
	 */
	public static final InfoDefinition INFODEF_PQ = new InfoDefinition("PQ", VariantPool.INFODEF_INT, "Phasing Quality", 1);
	
	/**
	 * ID = EC
	 * <br>Number = A
	 * <br>Type = Integer
	 * <br>Description = "Expected Alternate Allele Counts"
	 */
	public static final InfoDefinition INFODEF_EC = new InfoDefinition("EC", VariantPool.INFODEF_INT, "Expected Alternate Allele Counts", VariantPool.INFODEF_NARGS_PERALT);
	
	/**
	 * ID = MQ
	 * <br>Number = 1
	 * <br>Type = Integer
	 * <br>Description = "RMS Mapping Quality"
	 */
	public static final InfoDefinition INFODEF_MQ = new InfoDefinition("MQ", VariantPool.INFODEF_INT, "RMS Mapping Quality", 1);
	
	/**
	 * ID = CN
	 * <br>Number = 1
	 * <br>Type = Integer
	 * <br>Description = "Copy number genotype for imprecise events"
	 */
	public static final InfoDefinition INFODEF_CN = new InfoDefinition("CN", VariantPool.INFODEF_INT, "Copy number genotype for imprecise events", 1);
	
	/**
	 * ID = CNQ
	 * <br>Number = 1
	 * <br>Type = Float
	 * <br>Description = "Copy number genotype quality for imprecise events"
	 */
	public static final InfoDefinition INFODEF_CNQ = new InfoDefinition("CNQ", VariantPool.INFODEF_FLOAT, "Copy number genotype quality for imprecise events", 1);
	
	/**
	 * ID = CNL
	 * <br>Number = .
	 * <br>Type = Float
	 * <br>Description = "Copy number genotype likelihood for imprecise events"
	 */
	public static final InfoDefinition INFODEF_CNL = new InfoDefinition("CNL", VariantPool.INFODEF_FLOAT, "Copy number genotype likelihood for imprecise events", VariantPool.INFODEF_NARGS_VARIABLE);
	
	/**
	 * ID = NQ
	 * <br>Number = 1
	 * <br>Type = Integer
	 * <br>Description = "Phred style probability score that the variant is novel"
	 */
	public static final InfoDefinition INFODEF_NQ = new InfoDefinition("NQ", VariantPool.INFODEF_INT, "Phred style probability score that the variant is novel", 1);
	
	/**
	 * ID = HAP
	 * <br>Number = 1
	 * <br>Type = Integer
	 * <br>Description = "Unique haplotype identifier"
	 */
	public static final InfoDefinition INFODEF_HAP = new InfoDefinition("HAP", VariantPool.INFODEF_INT, "Unique haplotype identifier", 1);
	
	/**
	 * ID = AHAP
	 * <br>Number = 1
	 * <br>Type = Integer
	 * <br>Description = "Unique identifier of ancestral haplotype"
	 */
	public static final InfoDefinition INFODEF_AHAP = new InfoDefinition("AHAP", VariantPool.INFODEF_INT, "Unique identifier of ancestral haplotype", 1);
	
	
	public static final int ZYGOSITY_HOMOREF = 0;
	public static final int ZYGOSITY_HOMOALT = 1;
	public static final int ZYGOSITY_HETERORA = 2;
	public static final int ZYGOSITY_HETEROAA = 3;
	public static final int ZYGOSITY_UNKNOWN = 4;
	public static final int ZYGOSITY_HOMOREF_CNV = 5;
	public static final int ZYGOSITY_HOMOALT_CNV = 6;
	public static final int ZYGOSITY_HETERORA_CNV = 7;
	public static final int ZYGOSITY_HETEROAA_CNV = 8;
	public static final int ZYGOSITY_HETERORAA_CNV = 9;
	
	/* --- Instance Variables --- */
	
	//private boolean isEmpty;
	
	private int numberAlleles; //Used to size GL(PL, GP) and EC arrays
	private int maxCopyNumber; //Used to size CNL array
	
	private int[] alleles; //GT
	private boolean genotypePhased; //GT
	
	private int readDepth; //DP
	
	private boolean passedFilters; //FT
	private ArrayList<String> filters; //FT
	
	private double[] scaledGenotypeLikelihoods; //GL, PL, GP
	private Map<String, Double> heteroPloidyLikelihoods; //GLE
	
	private int conditionalQuality; // GQ
	private int[] haploQualities; //HQ - phred scaled
	
	private int phaseSet; //PS
	private int phaseQuality; //PQ - phred scaled
	
	private int[] expectedAltCounts; // EC
	private int mappingQuality; //MQ
	
		/*CNV/SV fields...*/
	private int copynumber; //CN
	private double cnQuality; //CNQ
	private double[] cnLikelihood; //CNL
	private int novelProbability; //NQ
	private int haplotypeIdent; //HAP
	private int ancestorHTID; //AHAP
	
	private Map<String, String> unparsedFieldMap; //For unrecognized values
	private Map<String, Field> fieldMap;
	private static Map<String, InfoDefinition> defMap;
	
	/* --- Construction --- */
	
	/**
	 * Construct an empty genotype with no custom fields and all
	 * common fields set to default values.
	 */
	public Genotype()
	{
		setDefaults();
		populateMaps();
	}
	
	/**
	 * Construct a genotype from a string in a VCF sample genotype field and an array of 
	 * strings representing the order and keys for each sub-field in the string.
	 * @param VCF_keys Ordered array of key strings describing fields present
	 * @param VCF_field Raw VCF formatted genotype string
	 * @throws UnsupportedFileTypeException If there is a parsing error.
	 */
	public Genotype(String[] VCF_keys, String VCF_field) throws UnsupportedFileTypeException
	{
		setDefaults();
		populateMaps();
		if (VCF_field == null) return;
		if (VCF_field.isEmpty()) return;
		if (VCF_field.equals(".")) return;
		
		String[] fields = VCF_field.split(":");
		//System.err.println("Genotype.<init> || DEBUG: VCF_keys.length = " + VCF_keys.length);
		//for (String s : VCF_keys) System.err.println("Genotype.<init> || \t" + s);
		//System.err.println("Genotype.<init> || DEBUG: fields.length = " + fields.length);
		//for (String s : fields) System.err.println("Genotype.<init> || \t" + s);
		//if (fields.length != VCF_keys.length) throw new UnsupportedFileTypeException();
		for (int i = 0; i < VCF_keys.length; i++)
		{
			if (fields.length> i)
			{
				unparsedFieldMap.put(VCF_keys[i], fields[i]);
			}
		}
		
		tryParseCommonFields();
		tryParseCommonCNVFields();
		//isEmpty = false;
	}
	
	/**
	 * Construct a genotype from a string in a VCF sample genotype field and a 
	 * VCF format field string.
	 * @param VCF_format Raw VCF format field string
	 * @param VCF_field Raw VCF formatted genotype string
	 * @throws UnsupportedFileTypeException If there is a parsing error.
	 */
	public Genotype(String VCF_format, String VCF_field) throws UnsupportedFileTypeException
	{
		this(VCF_format.split(":"), VCF_field);
	}
	
	private void setDefaults()
	{
		//isEmpty = true;
		numberAlleles = 0;
		maxCopyNumber = 0;
		alleles = new int[2];
		alleles[0] = -1;
		alleles[1] = -1;
		genotypePhased = false;
		readDepth = -1;
		passedFilters = false;
		filters = new ArrayList<String>(4);
		scaledGenotypeLikelihoods = null;
		heteroPloidyLikelihoods = new HashMap<String, Double>();
		conditionalQuality = -1;
		haploQualities = null;
		phaseSet = -1;
		phaseQuality = -1;
		expectedAltCounts = null;
		mappingQuality = -1;
		copynumber = -1;
		cnQuality = -1.0;
		cnLikelihood = null;
		unparsedFieldMap = new HashMap<String, String>();
		novelProbability = -1;
		haplotypeIdent = -1;
		ancestorHTID = -1;
	}

	private void populateMaps()
	{
		fieldMap = new HashMap<String, Field>();
		
		fieldMap.put(Genotype.INFODEF_AHAP.getKey(), new Field(){

			public String get() 
			{
				if (ancestorHTID < 0) return null;
				return Integer.toString(ancestorHTID);
			}

			public void set(String value) {
				try {
					parseAHAP(value);
				} catch (UnsupportedFileTypeException e) {
					//e.printStackTrace();
					return;
				}	
			}
			
		});
		
		fieldMap.put(Genotype.INFODEF_CN.getKey(), new Field(){

			public String get() 
			{
				if (copynumber < 0) return null;
				return Integer.toString(copynumber);
			}

			public void set(String value) {
				try {
					parseCN(value);
				} catch (UnsupportedFileTypeException e) {
					//e.printStackTrace();
					return;
				}	
			}
			
		});
		
		fieldMap.put(Genotype.INFODEF_CNL.getKey(), new Field(){

			public String get() 
			{
				if (cnLikelihood == null) return null;
				String s = "";
				for (int i = 0; i < cnLikelihood.length; i++)
				{
					s += Double.toString(cnLikelihood[i]);
					if (i < cnLikelihood.length - 1) s += ",";
				}
				return s;
			}

			public void set(String value) {
				try {
					parseCNL(value);
				} catch (UnsupportedFileTypeException e) {
					//e.printStackTrace();
					return;
				}	
			}
			
		});
		
		fieldMap.put(Genotype.INFODEF_CNQ.getKey(), new Field(){

			public String get() 
			{
				if (cnQuality < 0.0) return null;
				return Double.toString(cnQuality);
			}

			public void set(String value) {
				try {
					parseCNQ(value);
				} catch (UnsupportedFileTypeException e) {
					//e.printStackTrace();
					return;
				}	
			}
			
		});
		
		fieldMap.put(Genotype.INFODEF_DP.getKey(), new Field(){

			public String get() 
			{
				if (readDepth < 0) return null;
				return Integer.toString(readDepth);
			}

			public void set(String value) {
				try {
					parseDP(value);
				} catch (UnsupportedFileTypeException e) {
					//e.printStackTrace();
					return;
				}	
			}
			
		});
		
		fieldMap.put(Genotype.INFODEF_EC.getKey(), new Field(){

			public String get() 
			{
				if (expectedAltCounts == null) return null;
				String s = "";
				for (int i = 0; i < expectedAltCounts.length; i++)
				{
					s += Integer.toString(expectedAltCounts[i]);
					if (i < expectedAltCounts.length - 1) s += ",";
				}
				return s;
			}

			public void set(String value) {
				try {
					parseEC(value);
				} catch (UnsupportedFileTypeException e) {
					//e.printStackTrace();
					return;
				}	
			}
			
		});
		
		fieldMap.put(Genotype.INFODEF_FT.getKey(), new Field(){

			public String get() 
			{
				if (passedFilters) return "PASS";
				if (filters == null) return null;
				if (filters.isEmpty()) return null;
				String s = "";
				int sz = filters.size();
				for (int i = 0; i < sz; i++)
				{
					s += filters.get(i);
					if (i < sz - 1) s += ";";
				}
				return s;
			}

			public void set(String value) {
				parseFT(value);
			}
			
		});
		
		fieldMap.put(Genotype.INFODEF_GL.getKey(), new Field(){

			public String get() 
			{
				if (scaledGenotypeLikelihoods == null) return null;
				String s = "";
				for (int i = 0; i < scaledGenotypeLikelihoods.length; i++)
				{
					s += Double.toString(scaledGenotypeLikelihoods[i]);
					if (i < scaledGenotypeLikelihoods.length - 1) s += ",";
				}
				return s;
			}

			public void set(String value) {
				try {
					parseGL(value);
				} catch (UnsupportedFileTypeException e) {
					//e.printStackTrace();
					return;
				}	
			}
			
		});
		
		fieldMap.put(Genotype.INFODEF_GLE.getKey(), new Field(){

			public String get() 
			{
				if (heteroPloidyLikelihoods == null) return null;
				int size = heteroPloidyLikelihoods.size();
				List<String> gleKeys = new ArrayList<String>(size);
				gleKeys.addAll(heteroPloidyLikelihoods.keySet());
				Collections.sort(gleKeys);
				String s = "";
				for (int i = 0; i < size; i++)
				{
					String key = gleKeys.get(i);
					double val = heteroPloidyLikelihoods.get(key);
					s += key + ":" + val;
					if (i < size - 1) s += ",";
				}
				return s;
			}

			public void set(String value) {
				try {
					parseGLE(value);
				} catch (UnsupportedFileTypeException e) {
					//e.printStackTrace();
					return;
				}	
			}
			
		});
		
		fieldMap.put(Genotype.INFODEF_GP.getKey(), new Field(){

			public String get() 
			{
				if (scaledGenotypeLikelihoods == null) return null;
				String s = "";
				for (int i = 0; i < scaledGenotypeLikelihoods.length; i++)
				{
					double scaledFloat = scaledGenotypeLikelihoods[i];
					double phredVal = scaledFloat * (-10.0);
					s += Double.toString(phredVal);
					if (i < scaledGenotypeLikelihoods.length - 1) s += ",";
				}
				return s;
			}

			public void set(String value) {
				try {
					parseGP(value);
				} catch (UnsupportedFileTypeException e) {
					//e.printStackTrace();
					return;
				}	
			}
			
		});
		
		fieldMap.put(Genotype.INFODEF_GQ.getKey(), new Field(){

			public String get() 
			{
				if (conditionalQuality < 0) return null;
				return Integer.toString(conditionalQuality);
			}

			public void set(String value) {
				try {
					parseGQ(value);
				} catch (UnsupportedFileTypeException e) {
					//e.printStackTrace();
					return;
				}	
			}
			
		});
		
		fieldMap.put(Genotype.INFODEF_GT.getKey(), new Field(){

			public String get() 
			{
				if (alleles == null) return null;
				String s = "";
				for (int i = 0; i < alleles.length; i++)
				{
					if (alleles[i] >= 0) s += alleles[i];
					else s += ".";
					if (i < alleles.length - 1)
					{
						if (genotypePhased) s += "|";
						else s += "/";
					}
				}
				return s;
			}

			public void set(String value) {
				try {
					parseGT(value);
				} catch (UnsupportedFileTypeException e) {
					//e.printStackTrace();
					return;
				}	
			}
			
		});
		
		fieldMap.put(Genotype.INFODEF_HAP.getKey(), new Field(){

			public String get() 
			{
				if (haplotypeIdent == -1) return null;
				return Integer.toString(haplotypeIdent);
			}

			public void set(String value) {
				try {
					parseHAP(value);
				} catch (UnsupportedFileTypeException e) {
					//e.printStackTrace();
					return;
				}	
			}
			
		});
	
		fieldMap.put(Genotype.INFODEF_HQ.getKey(), new Field(){

			public String get() 
			{
				if (haploQualities == null) return null;
				if (haploQualities.length != 2)
				{
					haploQualities = null;
					return null;
				}
				String s = "";
				for (int i = 0; i < haploQualities.length; i++)
				{
					s += haploQualities[i];
					if (i < haploQualities.length - 1) s += ",";
				}
				return s;
			}

			public void set(String value) {
				try {
					parseHQ(value);
				} catch (UnsupportedFileTypeException e) {
					//e.printStackTrace();
					return;
				}	
			}
			
		});
		
		fieldMap.put(Genotype.INFODEF_MQ.getKey(), new Field(){

			public String get() 
			{
				if (mappingQuality == -1) return null;
				return Integer.toString(mappingQuality);
			}

			public void set(String value) {
				try {
					parseMQ(value);
				} catch (UnsupportedFileTypeException e) {
					//e.printStackTrace();
					return;
				}	
			}
			
		});
	
		fieldMap.put(Genotype.INFODEF_NQ.getKey(), new Field(){

			public String get() 
			{
				if (novelProbability == -1) return null;
				return Integer.toString(novelProbability);
			}

			public void set(String value) {
				try {
					parseNQ(value);
				} catch (UnsupportedFileTypeException e) {
					//e.printStackTrace();
					return;
				}	
			}
			
		});
	
		fieldMap.put(Genotype.INFODEF_PL.getKey(), new Field(){

			public String get() 
			{
				if (scaledGenotypeLikelihoods == null) return null;
				String s = "";
				for (int i = 0; i < scaledGenotypeLikelihoods.length; i++)
				{
					double scaledFloat = scaledGenotypeLikelihoods[i];
					int phredVal = Utilities.floatScaleToPhred(scaledFloat);
					s += Integer.toString(phredVal);
					if (i < scaledGenotypeLikelihoods.length - 1) s += ",";
				}
				return s;
			}

			public void set(String value) {
				try {
					parsePL(value);
				} catch (UnsupportedFileTypeException e) {
					//e.printStackTrace();
					return;
				}	
			}
			
		});
		
		fieldMap.put(Genotype.INFODEF_PQ.getKey(), new Field(){

			public String get() 
			{
				if (phaseQuality == -1) return null;
				return Integer.toString(phaseQuality);
			}

			public void set(String value) {
				try {
					parsePQ(value);
				} catch (UnsupportedFileTypeException e) {
					//e.printStackTrace();
					return;
				}	
			}
			
		});
	
		fieldMap.put(Genotype.INFODEF_PS.getKey(), new Field(){

			public String get() 
			{
				if (phaseSet == -1) return null;
				return Integer.toString(phaseSet);
			}

			public void set(String value) {
				try {
					parsePS(value);
				} catch (UnsupportedFileTypeException e) {
					//e.printStackTrace();
					return;
				}	
			}
			
		});
	
		if (defMap == null) populateDefMap();
	}
	
	private static void populateDefMap()
	{
		defMap = new HashMap<String, InfoDefinition>();
		defMap.put(Genotype.INFODEF_AHAP.getKey(), INFODEF_AHAP);
		defMap.put(Genotype.INFODEF_CN.getKey(), INFODEF_CN);
		defMap.put(Genotype.INFODEF_CNL.getKey(), INFODEF_CNL);
		defMap.put(Genotype.INFODEF_CNQ.getKey(), INFODEF_CNQ);
		defMap.put(Genotype.INFODEF_DP.getKey(), INFODEF_DP);
		defMap.put(Genotype.INFODEF_EC.getKey(), INFODEF_EC);
		defMap.put(Genotype.INFODEF_FT.getKey(), INFODEF_FT);
		defMap.put(Genotype.INFODEF_GL.getKey(), INFODEF_GL);
		defMap.put(Genotype.INFODEF_GLE.getKey(), INFODEF_GLE);
		defMap.put(Genotype.INFODEF_GP.getKey(), INFODEF_GP);
		defMap.put(Genotype.INFODEF_GQ.getKey(), INFODEF_GQ);
		defMap.put(Genotype.INFODEF_GT.getKey(), INFODEF_GT);
		defMap.put(Genotype.INFODEF_HAP.getKey(), INFODEF_HAP);
		defMap.put(Genotype.INFODEF_HQ.getKey(), INFODEF_HQ);
		defMap.put(Genotype.INFODEF_MQ.getKey(), INFODEF_MQ);
		defMap.put(Genotype.INFODEF_NQ.getKey(), INFODEF_NQ);
		defMap.put(Genotype.INFODEF_PL.getKey(), INFODEF_PL);
		defMap.put(Genotype.INFODEF_PQ.getKey(), INFODEF_PQ);
		defMap.put(Genotype.INFODEF_PS.getKey(), INFODEF_PS);
	}
	
	private void tryParseCommonFields() throws UnsupportedFileTypeException
	{
		parseGT(unparsedFieldMap.remove(Genotype.INFODEF_GT.getKey()));
		parseDP(unparsedFieldMap.remove(Genotype.INFODEF_DP.getKey()));
		parseFT(unparsedFieldMap.remove(Genotype.INFODEF_FT.getKey()));
		parseGL(unparsedFieldMap.remove(Genotype.INFODEF_GL.getKey()));
		parseGLE(unparsedFieldMap.remove(Genotype.INFODEF_GLE.getKey()));
		parsePL(unparsedFieldMap.remove(Genotype.INFODEF_PL.getKey()));
		parseGP(unparsedFieldMap.remove(Genotype.INFODEF_GP.getKey()));
		parseGQ(unparsedFieldMap.remove(Genotype.INFODEF_GQ.getKey()));
		parseHQ(unparsedFieldMap.remove(Genotype.INFODEF_HQ.getKey()));
		parsePS(unparsedFieldMap.remove(Genotype.INFODEF_PS.getKey()));
		parsePQ(unparsedFieldMap.remove(Genotype.INFODEF_PQ.getKey()));
		parseEC(unparsedFieldMap.remove(Genotype.INFODEF_EC.getKey()));
		parseMQ(unparsedFieldMap.remove(Genotype.INFODEF_MQ.getKey()));
	}
	
	private void tryParseCommonCNVFields() throws UnsupportedFileTypeException
	{
		parseCN(unparsedFieldMap.remove(Genotype.INFODEF_CN.getKey()));
		parseCNQ(unparsedFieldMap.remove(Genotype.INFODEF_CNQ.getKey()));
		parseCNL(unparsedFieldMap.remove(Genotype.INFODEF_CNL.getKey()));
		parseNQ(unparsedFieldMap.remove(Genotype.INFODEF_NQ.getKey()));
		parseHAP(unparsedFieldMap.remove(Genotype.INFODEF_HAP.getKey()));
		parseAHAP(unparsedFieldMap.remove(Genotype.INFODEF_AHAP.getKey()));
	}
	
	/* --- Field Parsers --- */
	
	private void parseGT(String val) throws UnsupportedFileTypeException
	{
		//String gt = fieldMap.get("GT");
		if (val != null)
		{
			if (val.equals(".")) return;
			String sep = "|";
			int lineind = val.indexOf(sep);
			if (lineind >= 0) {
				genotypePhased = true;
				//System.err.println("Genotype.parseGT || Genotype phased. Sep char at " + lineind);
				sep = "\\|";
			}
			else sep = "/";
			String[] aStr = val.split(sep);
			//System.err.println("Genotype.parseGT || sep = " + sep);
			//System.err.println("Genotype.parseGT || alleles found: " + aStr.length);
			//System.err.println("Genotype.parseGT || aStr: ");
			alleles = new int[aStr.length];
			for (int i = 0; i < aStr.length; i++)
			{
				//System.err.println("\t" + aStr[i]);
				try
				{
					int n = Integer.parseInt(aStr[i]);
					alleles[i] = n;	
				}
				catch (NumberFormatException e)
				{
					if (aStr[i].equals(".")) alleles[i] = -1;
					else
					{
						System.err.println("Genotype.parseGT || Full input: " + val);
						e.printStackTrace();
						throw new UnsupportedFileTypeException();
					}
				}
			}
		}
	}
	
	private void parseDP(String val) throws UnsupportedFileTypeException
	{
		//String dp = fieldMap.get("DP");
		if (val != null)
		{
			if (val.equals(".")) return;
			try
			{
				readDepth = Integer.parseInt(val);
			}
			catch (NumberFormatException e)
			{
				e.printStackTrace();
				throw new UnsupportedFileTypeException();
			}
		}
	}
	
	private void parseFT(String val)
	{
		//String ft = fieldMap.get("FT");
		if (val != null)
		{
			if (val.equals(".")) return;
			if (val.equals("PASS"))
			{
				this.passedFilters = true;
				return;
			}
			String[] flist = val.split(";");
			filters = new ArrayList<String>(flist.length + 2);
			for (String s : flist) filters.add(s);
		}
	}
	
	private void parseGL(String val) throws UnsupportedFileTypeException
	{
		//String gl = fieldMap.get("GL");
		if (val != null)
		{
			if (val.equals(".")) return;
			String[] vals = val.split(",");
			scaledGenotypeLikelihoods = new double[vals.length];
			for (int i = 0; i < vals.length; i++)
			{
				try
				{
					scaledGenotypeLikelihoods[i] = Double.parseDouble(vals[i]);	
				}
				catch (NumberFormatException e)
				{
					e.printStackTrace();
					scaledGenotypeLikelihoods = null;
					//numberAlleles = 0;
					throw new UnsupportedFileTypeException();
				}
			}
			numberAlleles = comboToNumber(scaledGenotypeLikelihoods.length);
		}
	}
	
	private void parseGLE(String val) throws UnsupportedFileTypeException
	{
		//String gle = fieldMap.get("GLE");
		if (val != null)
		{
			if (val.equals(".")) return;
			if (val.indexOf("=") >= 0) val = val.substring(val.indexOf("=") + 1);
			String[] vals = val.split(",");
			for (int i = 0; i < vals.length; i++)
			{
				try
				{
					String[] subs = vals[i].split(":");
					String gt = subs[0];
					double d = Double.parseDouble(subs[1]);
					heteroPloidyLikelihoods.put(gt, d);
				}
				catch (NumberFormatException | NullPointerException | ArrayIndexOutOfBoundsException e)
				{
					e.printStackTrace();
					throw new UnsupportedFileTypeException();
				}
			}
		}
	}
	
	private void parsePL(String val) throws UnsupportedFileTypeException
	{
		//String pl = fieldMap.get("PL");
		if (val != null)
		{
			if (val.equals(".")) return;
			String[] vals = val.split(",");
			scaledGenotypeLikelihoods = new double[vals.length];
			for (int i = 0; i < vals.length; i++)
			{
				try
				{
					int phredVal = Integer.parseInt(vals[i]);
					double scaledFloat = Utilities.phredToFloatScale(phredVal);
					scaledGenotypeLikelihoods[i] = scaledFloat;
				}
				catch (NumberFormatException e)
				{
					e.printStackTrace();
					scaledGenotypeLikelihoods = null;
					//numberAlleles = 0;
					throw new UnsupportedFileTypeException();
				}
			}
			numberAlleles = comboToNumber(scaledGenotypeLikelihoods.length);
		}
	}
	
	private void parseGP(String val) throws UnsupportedFileTypeException
	{
		//String gp = fieldMap.get("GP");
		if (val != null)
		{
			if (val.equals(".")) return;
			String[] vals = val.split(",");
			scaledGenotypeLikelihoods = new double[vals.length];
			for (int i = 0; i < vals.length; i++)
			{
				try
				{
					double phredVal = Double.parseDouble(vals[i]);
					double scaledFloat = phredVal / (-10.0);
					scaledGenotypeLikelihoods[i] = scaledFloat;
				}
				catch (NumberFormatException e)
				{
					e.printStackTrace();
					scaledGenotypeLikelihoods = null;
					//numberAlleles = 0;
					throw new UnsupportedFileTypeException();
				}
			}
			numberAlleles = comboToNumber(scaledGenotypeLikelihoods.length);
		}
	}
	
	private void parseGQ(String val) throws UnsupportedFileTypeException
	{
		//String gq = fieldMap.get("GQ");
		if (val != null)
		{
			if (val.equals(".")) return;
			try
			{
				conditionalQuality = Integer.parseInt(val);	
			}
			catch (NumberFormatException e)
			{
				e.printStackTrace();
				throw new UnsupportedFileTypeException();
			}
		}
	}
	
	private void parseHQ(String val) throws UnsupportedFileTypeException
	{
		//String hq = fieldMap.get("HQ");
		if (val != null)
		{
			if (val.equals(".")) return;
			String[] vals = val.split(",");
			haploQualities = new int[vals.length];
			for (int i = 0; i < vals.length; i++)
			{
				try
				{
					haploQualities[i] = Integer.parseInt(vals[i]);
				}
				catch (NumberFormatException e)
				{
					e.printStackTrace();
					throw new UnsupportedFileTypeException();
				}
			}
		}
	}
	
	private void parsePS(String val) throws UnsupportedFileTypeException
	{
		//String ps = fieldMap.get("PS");
		if (val != null)
		{
			if (val.equals(".")) return;
			try
			{
				phaseSet = Integer.parseInt(val);	
			}
			catch (NumberFormatException e)
			{
				e.printStackTrace();
				throw new UnsupportedFileTypeException();
			}
		}
	}
	
	private void parsePQ(String val) throws UnsupportedFileTypeException
	{
		//String pq = fieldMap.get("PQ");
		if (val != null)
		{
			if (val.equals(".")) return;
			try
			{
				phaseQuality = Integer.parseInt(val);	
			}
			catch (NumberFormatException e)
			{
				e.printStackTrace();
				throw new UnsupportedFileTypeException();
			}
		}
	}
	
	private void parseEC(String val) throws UnsupportedFileTypeException
	{
		//String ec = fieldMap.get("EC");
		if (val != null)
		{
			if (val.equals(".")) return;
			String[] vals = val.split(",");
			expectedAltCounts = new int[vals.length];
			for (int i = 0; i < vals.length; i++)
			{
				try
				{
					expectedAltCounts[i] = Integer.parseInt(vals[i]);
				}
				catch (NumberFormatException e)
				{
					e.printStackTrace();
					expectedAltCounts = null;
					throw new UnsupportedFileTypeException();
				}
			}
			numberAlleles = expectedAltCounts.length + 1;
		}
	}
	
	private void parseMQ(String val) throws UnsupportedFileTypeException
	{
		//String mq = fieldMap.get("MQ");
		if (val != null)
		{
			if (val.equals(".")) return;
			try
			{
				mappingQuality = Integer.parseInt(val);	
			}
			catch (NumberFormatException e)
			{
				e.printStackTrace();
				throw new UnsupportedFileTypeException();
			}
		}
	}
	
	private void parseCN(String val) throws UnsupportedFileTypeException
	{
		//String cn = fieldMap.get(Genotype.INFODEF_CN.getKey());
		if (val != null)
		{
			if (val.equals(".")) return;
			try
			{
				copynumber = Integer.parseInt(val);
			}
			catch (NumberFormatException e)
			{
				e.printStackTrace();
				throw new UnsupportedFileTypeException();
			}
		}
	}
	
	private void parseCNQ(String val) throws UnsupportedFileTypeException
	{
		//String val = fieldMap.get(Genotype.INFODEF_CNQ.getKey());
		if (val != null)
		{
			if (val.equals(".")) return;
			try
			{
				cnQuality = Double.parseDouble(val);
			}
			catch (NumberFormatException e)
			{
				e.printStackTrace();
				throw new UnsupportedFileTypeException();
			}
		}
	}
	
	private void parseCNL(String val) throws UnsupportedFileTypeException
	{
		//String val = fieldMap.get(Genotype.INFODEF_CNL.getKey());
		if (val != null)
		{
			if (val.equals(".")) return;
			String[] vals = val.split(",");
			cnLikelihood = new double[vals.length];
			for (int i = 0; i < vals.length; i++)
			{
				try
				{
					cnLikelihood[i] = Double.parseDouble(vals[i]);
				}
				catch (NumberFormatException e)
				{
					e.printStackTrace();
					cnLikelihood = null;
					maxCopyNumber = 0;
					throw new UnsupportedFileTypeException();
				}
			}
			maxCopyNumber = cnLikelihood.length - 1;
		}
	}
	
	private void parseNQ(String val) throws UnsupportedFileTypeException
	{
		//String val = fieldMap.get(Genotype.INFODEF_NQ.getKey());
		if (val != null)
		{
			if (val.equals(".")) return;
			try
			{
				novelProbability = Integer.parseInt(val);
			}
			catch (NumberFormatException e)
			{
				e.printStackTrace();
				throw new UnsupportedFileTypeException();
			}
		}
	}
	
	private void parseHAP(String val) throws UnsupportedFileTypeException
	{
		//String val = fieldMap.get(Genotype.INFODEF_HAP.getKey());
		if (val != null)
		{
			if (val.equals(".")) return;
			try
			{
				haplotypeIdent = Integer.parseInt(val);
			}
			catch (NumberFormatException e)
			{
				e.printStackTrace();
				throw new UnsupportedFileTypeException();
			}
		}
	}
	
	private void parseAHAP(String val) throws UnsupportedFileTypeException
	{
		//String val = fieldMap.get(Genotype.INFODEF_AHAP.getKey());
		if (val != null)
		{
			if (val.equals(".")) return;
			try
			{
				ancestorHTID = Integer.parseInt(val);
			}
			catch (NumberFormatException e)
			{
				e.printStackTrace();
				throw new UnsupportedFileTypeException();
			}
		}
	}
	
	/* --- Field Management --- */
	
	private interface Field
	{
		public String get();
		public void set(String value);
	}
	
	/* --- Getters --- */
	
	/**
	 * Get the list of allele calls (by general allele index) for this
	 * genotype.
	 * <br>An allele index of zero indicates the reference allele. Anything greater
	 * than zero indicates an alt allele.
	 * <br>To get the alt allele index, subtract one from the general allele index.
	 * @return An int array listing the general allele indices of the called alleles, if
	 * they have been recorded.
	 * <br> This method returns null if no allele calls have been made.
	 */
	public int[] getAlleles()
	{
		if (alleles == null) return null;
		int[] acopy = new int[alleles.length];
		for (int i = 0; i < alleles.length; i++) acopy[i] = alleles[i];
		return acopy;
	}
	
	/**
	 * Get whether the genotype call is phased or not.
	 * @return True - If the genotype is phased.
	 * <br>False - If the genotype is not phased.
	 */
	public boolean isPhased()
	{
		return genotypePhased;
	}
	
	/**
	 * Get the read depth recorded at the position this genotype is called for.
	 * @return Read depth at sample position.
	 */
	public int getReadDepth()
	{
		return readDepth;
	}
	
	/**
	 * Get whether this genotype has any filter records, passed or failed.
	 * @return True - If this genotype is either marked has having passed all filters,
	 * or if there are records of failed filters.
	 * <br>False - If no trace of filtering was found for this genotype.
	 */
	public boolean hasFilters()
	{
		if (passedFilters) return true;
		else return (filters != null);
	}
	
	/**
	 * Get whether this genotype has been marked as passing all filters it
	 * was put through.
	 * @return True - If genotype has passed all filters.
	 * <br> False - Otherwise.
	 */
	public boolean passedFilters()
	{
		return passedFilters;
	}
	
	/**
	 * Get a list of filters (names/IDs of filters) that this genotype
	 * has failed.
	 * @return A string array containing the IDs of all filters failed by this
	 * genotype, or null if no such filters exist.
	 */
	public String[] failedFilters()
	{
		if (filters == null) return null;
		if (filters.isEmpty()) return null;
		int sz = filters.size();
		String[] acopy = new String[sz];
		for (int i = 0; i < sz; i++) acopy[i] = filters.get(i);
		return acopy;
	}
	
	private static int getAlleleComboIndex(int a1, int a2)
	{
		return (a1 * (a1 + 1))/2 + a2;
	}
	
	/**
	 * Given two general allele indices, get the likelihood of the particular
	 * diploid genotype generated from those two alleles.
	 * @param a1 Allele index 1.
	 * @param a2 Allele index 2.
	 * @return The likelihood as a scaled double value of the provided genotype, or 0.0
	 * if the allele indices are invalid or if there are no likelihoods recorded.
	 */
 	public double getGenotypeLikelihood(int a1, int a2)
	{
		if (scaledGenotypeLikelihoods == null) return 0.0;
		if (a1 >= numberAlleles) return 0.0;
		if (a2 >= numberAlleles) return 0.0;
		if (a1 > a2)
		{
			int temp = a1;
			a1 = a2;
			a2 = temp;
		}
		int i = getAlleleComboIndex(a1, a2);
		if (i < 0 || i >= scaledGenotypeLikelihoods.length) return 0.0;
		return scaledGenotypeLikelihoods[i];
	}
	
	/**
	 * Get the likelihood of a genotype of heterogeneous ploidy. This field
	 * is useful when the copy number is heterogeneous or anything other than two.
	 * <br>The genotype must be specified with a string following the VCF genotype standard.
	 * <br>This requires knowledge of the general allele index/indices of interest.
	 * General allele indices must be separated using a forward slash (/) if unphased, or a line
	 * (|) if phased. The string must match perfectly to the string recorded for the genotype.
	 * @param genotype A string representing the genotype of interest.
	 * @return The likelihood of the given genotype as a scaled double value. If the genotype
	 * provided is not on record, this method will return 0.0.
	 */
	public double getHeterogeneousPloidyGenotypeLikelihood(String genotype)
	{
		Double d = heteroPloidyLikelihoods.get(genotype);
		if (d == null) return 0.0;
		return d;
	}
	
	/**
	 * Phred-scaled genotype quality score. Logarithmically scaled probability
	 * that the genotype call is incorrect.
	 * @return Genotype quality score
	 */
	public int getConditionalGenotypeQuality()
	{
		return conditionalQuality;
	}
	
	/**
	 * Get the quality score of the first haplotype, if recorded.
	 * @return Phred-scaled haplotype quality score, or -1 if not found.
	 */
	public int getHaplotypeQuality_1()
	{
		if (haploQualities == null) return -1;
		if (haploQualities.length < 1) return -1;
		return haploQualities[0];
	}
	
	/**
	 * Get the quality score of the second haplotype, if recorded.
	 * @return Phred-scaled haplotype quality score, or -1 if not found.
	 */
	public int getHaplotypeQuality_2()
	{
		if (haploQualities == null) return -1;
		if (haploQualities.length < 2) return -1;
		return haploQualities[1];
	}
	
	/**
	 * ID of the set of phased genotypes to which this genotype, if phased, belongs.
	 * This value should never be negative, nor should it exceed the maximum value that
	 * can be represented by a signed 32-bit integer.
	 * This field is ignored by most parsers if the genotype is unphased.
	 * @return The phase set ID number.
	 */
	public int getPhaseSet()
	{
		return phaseSet;
	}
	
	/**
	 * Get the phasing quality - the probability that the alleles are ordered incorrectly
	 * relative to other members of the phase set.
	 * <br>This appears to be reserved for VCF 4.2 standard usage, but has no meaning as of now.
	 * @return Phred-scaled phasing quality score.
	 */
	public int getPhaseQuality()
	{
		return phaseQuality;
	}
	
	/**
	 * Get the expected allele count, if recorded, for a given alternate allele.
	 * @param altIndex Index of the alternate allele (zero-based alt allele index, or
	 * general allele index - 1).
	 * @return Expected allele count, if found. This method returns -1 if no count could be
	 * found for the provided allele, or the index is invalid.
	 */
	public int getExpectedAltCount(int altIndex)
	{
		if (expectedAltCounts == null) return -1;
		if (altIndex + 1 >= this.numberAlleles) return -1;
		if (expectedAltCounts.length <= altIndex) return -1;
		if (altIndex < 0) return -1;
		return expectedAltCounts[altIndex];
	}
	
	/**
	 * Get the recorded mapping quality at the position of the variant this genotype
	 * is associated with.
	 * @return Phred-scaled mapping quality, or -1 if unknown.
	 */
	public int getMappingQuality()
	{
		return mappingQuality;
	}
	
	/**
	 * Get the gene copy number for a copy number variant.
	 * @return Copy number, if on record. -1 if variant is not a CNV or CN is unknown.
	 */
	public int getCopyNumber()
	{
		return copynumber;
	}
	
	/**
	 * Get quality score for copy number call if genotype is for a copy number variant.
	 * @return Phred-scaled copy number call quality score. -1 if genotype is not for
	 * a CNV or copy number quality is unknown. 
	 */
	public double getCopyNumberQuality()
	{
		return this.getCopyNumberQuality();
	}
	
	/**
	 * Get the log scaled likelihood that the genotype has the given copy number.
	 * @param copyNumber Possible copy number of interest.
	 * @return Copy number likelihood, if known. 0.0 if unknown.
	 */
	public double getCopyNumberLikelihood(int copyNumber)
	{
		if (copyNumber < 0) return 0.0;
		if (copyNumber > maxCopyNumber) return 0.0;
		if (cnLikelihood == null) return 0.0;
		if (cnLikelihood.length <= copyNumber) return 0.0;
		return cnLikelihood[copyNumber];
	}

	/**
	 * Get the probability that the variant is novel. This field is generally
	 * defined for CNVs/ SVs.
	 * @return Phred-scaled probability as an integer, if known. If unknown, -1.
	 */
	public int getNovelProbability()
	{
		return novelProbability;
	}

	/**
	 * Get the ID number for the haplotype, if known and set.
	 * @return Haplotype ID number as an integer. -1 if unknown.
	 */
	public int getHaplotypeID()
	{
		return this.haplotypeIdent;
	}
	
	/**
	 * Get the ID number for the ancestral haplotype, if known and set.
	 * @return Ancestral haplotype ID number, if on record. -1 if unknown.
	 */
	public int getAncestralHaplotypeID()
	{
		return this.ancestorHTID;
	}
	
	/**
	 * Get a VCF formatted string representation of any field recorded for this genotype,
	 * including parsed fields, from the key.
	 * <br>This function is not recommended for regular access to pre-parsed common fields as it will 
	 * generate a string version of the native value.
	 * It's intent for common fields is use in VCF formatted serialization.
	 * @param fieldKey Genotype field key
	 * @return The field value as a string, if present. Null if no value for the field is
	 * on record.
	 */
	public String getField(String fieldKey)
	{
		//First check unparsed fields.
		String val = unparsedFieldMap.get(fieldKey);
		if (val != null) return val;
		
		//Check parsed fields
		Field f = fieldMap.get(fieldKey);
		if (f == null) return null;
		return f.get();
	}
	
	/* --- Setters --- */
	
	/**
	 * Set any genotype field to the given string. Parsing will be attempted if the key
	 * provided is for a common field accessible through another setter. If the given 
	 * String cannot be parsed correctly, this method will ultimately change nothing.
	 * <br>Use of this method for setting is recommended for fields that are not currently
	 * defined as instance variables.
	 * @param fieldKey Key of genotype field to set value of.
	 * @param value String representation of value to set.
	 */
	public void setField(String fieldKey, String value)
	{
		if (fieldKey == null || fieldKey.isEmpty()) return;
		//See if parsed field.
		Field f = this.fieldMap.get(fieldKey);
		if (f == null) this.unparsedFieldMap.put(fieldKey, value);
		else f.set(value);
	}

	/**
	 * Set the called alleles for this genotype using a VCF style GT field string.
	 * @param gtString GT format string to parse into genotype alleles.
	 */
	public boolean setAlleles(String gtString)
	{
		try 
		{
			this.parseGT(gtString);
			return true;
		} 
		catch (UnsupportedFileTypeException e) 
		{
			return false;
		}
	}
	
	/**
	 * Set the called alleles for this genotype to the specified values.
	 * Array should represent a list of the general allele indices of the alleles
	 * to set. A general allele index of zero represents the reference allele - anything
	 * higher is an alternate allele.
	 * <br>This method rejects any negative values, but it does not check for an index
	 * upper bound.
	 * @param newAlleles Array of general allele indicies to set.
	 */
	public void setAlleles(final int[] newAlleles)
	{
		if (newAlleles == null) return;
		if (newAlleles.length < 1) return;
		int[] acopy = new int[newAlleles.length];
		for (int i = 0; i < newAlleles.length; i++)
		{
			int a = newAlleles[i];
			if (a < 0) return;
			acopy[i] = a;
		}
		alleles = acopy;
	}

	/**
	 * Set the genotype phasing flag.
	 * @param phased Whether this genotype is phased.
	 */
	public void setPhased(boolean phased)
	{
		this.genotypePhased = phased;
	}
	
	/**
	 * Set the read depth at the position of the variant associated with this genotype.
	 * @param value Integer representing the read depth. -1 if unknown.
	 */
	public void setReadDepth(int value)
	{
		readDepth = value;
	}
	
	/**
	 * Manually set the filters passed flag. If this flag is set, it is assumed that this
	 * genotype has passed all filters, regardless of the contents of its filter list.
	 * @param b Desired flag setting
	 */
	public void setFilterPassFlag(boolean b)
	{
		this.passedFilters = b;
	}
	
	/**
	 * Add a filter by name or ID to the genotype's list of filters. This list is usually
	 * used to list filters that it has failed.
	 * @param filterName Name or ID of filter to add.
	 */
	public void addFilter(String filterName)
	{
		if (filters == null) filters = new ArrayList<String>(4);
		filters.add(filterName);
	}
	
	/**
	 * Remove a filter from the genotype's filter list with a given name.
	 * @param filterName Name of filter to remove.
	 */
	public void removeFilter(String filterName)
	{
		if (filters == null) filters = new ArrayList<String>(4);
		filters.remove(filterName);
	}
	
	/**
	 * Clear this genotype's list of filters. This list is usually used to store
	 * the names of filters the genotype or variant has failed.
	 */
	public void clearFilters()
	{
		if (filters == null) filters = new ArrayList<String>(4);
		filters.clear();
	}
	
	/**
	 * Set the total number of possible alleles (reference and alternate) recorded for the variant
	 * this genotype is associated with.
	 * <br>This controls internal array sizes so that probabilities of various allele combinations
	 * can be set properly.
	 * <br>If a negative number is provided to this function, it will return without changing anything.
	 * A value of zero indicates that the field is unset.
	 * <br>IMPORTANT: If the number is set lower than its current value, data may be lost.
	 * @param n Number of alleles to set value to. Must be greater than or equal to 0.
	 */
	public void setNumberAlleles(int n)
	{
		if (n < 0) return;
		numberAlleles = n;
		if (numberAlleles == 0)
		{
			scaledGenotypeLikelihoods = null;
			return;
		}
		if (numberAlleles > 0)
		{
			int targetSize = numberToCombo(n);
			double[] temp = new double[targetSize];
			if (scaledGenotypeLikelihoods != null)
			{
				for (int i = 0; i < targetSize; i++)
				{
					if (i < scaledGenotypeLikelihoods.length) temp[i] = this.scaledGenotypeLikelihoods[i];
					else temp[i] = 0.0;
				}	
			}
			scaledGenotypeLikelihoods = temp;
			int[] temp2 = new int[numberAlleles];
			if (expectedAltCounts != null)
			{
				for (int i = 0; i < targetSize; i++)
				{
					if (i < expectedAltCounts.length) temp2[i] = this.expectedAltCounts[i];
					else temp2[i] = 0;
				}	
			}
			expectedAltCounts = temp2;
		}
	}
	
	/**
	 * Set the likelihood, as a log 10 scaled probability, that the combination of the
	 * two given alleles (general allele indices) is the true genotype.
	 * <br>Allele indices must be valid. A general index of zero indicates the reference allele.
	 * The index cannot be less than zero, and cannot exceed the set number of alleles.
	 * @param a1 - General index of first allele
	 * @param a2 - General index of second allele
	 * @param scaledValue - Likelihood value to set
	 */
	public void setScaledGenotypeLikelihood(int a1, int a2, double scaledValue)
	{
		if (a1 < 0 || a2 < 0) return;
		if (a1 >= numberAlleles || a2 >= numberAlleles) return;
		if (a1 > a2)
		{
			int temp = a1;
			a1 = a2;
			a2 = temp;
		}
		int index = getAlleleComboIndex(a1, a2);
		if (index >= scaledGenotypeLikelihoods.length || index < 0) return;
		
		scaledGenotypeLikelihoods[index] = scaledValue;
	}
	
	/**
	 * Set the probability that the combination of the
	 * two given alleles (general allele indices) is the true genotype.
	 * The stored likelihood will be a log 10 scaled version of this probability. Because
	 * floating-point values can only be so precise, transferring information via probabilities
	 * may cause data loss.
	 * <br>Allele indices must be valid. A general index of zero indicates the reference allele.
	 * The index cannot be less than zero, and cannot exceed the set number of alleles.
	 * @param a1 - General index of first allele
	 * @param a2 - General index of second allele
	 * @param p - Probability to set
	 */
	public void setGenotypeProbability(int a1, int a2, double p)
	{
		if (a1 < 0 || a2 < 0) return;
		if (a1 >= numberAlleles || a2 >= numberAlleles) return;
		if (a1 > a2)
		{
			int temp = a1;
			a1 = a2;
			a2 = temp;
		}
		int index = getAlleleComboIndex(a1, a2);
		if (index >= scaledGenotypeLikelihoods.length || index < 0) return;
		
		scaledGenotypeLikelihoods[index] = Utilities.probToFloatScale(p);
	}
	
	/**
	 * Set the likelihood, as a phred scaled probability, that the combination of the
	 * two given alleles (general allele indices) is the true genotype.
	 * <br>Allele indices must be valid. A general index of zero indicates the reference allele.
	 * The index cannot be less than zero, and cannot exceed the set number of alleles.
	 * @param a1 - General index of first allele
	 * @param a2 - General index of second allele
	 * @param phred - Rounded phred scaled probability value
	 */
	public void setGenotypeLikelihoodPhred(int a1, int a2, int phred)
	{
		if (a1 < 0 || a2 < 0) return;
		if (a1 >= numberAlleles || a2 >= numberAlleles) return;
		if (a1 > a2)
		{
			int temp = a1;
			a1 = a2;
			a2 = temp;
		}
		int index = getAlleleComboIndex(a1, a2);
		if (index >= scaledGenotypeLikelihoods.length || index < 0) return;
		
		scaledGenotypeLikelihoods[index] = Utilities.phredToFloatScale(phred);
	}

	private String generateHeteroGenoMapKey(List<Integer> alleles)
	{
		if (alleles == null) return null;
		if (alleles.isEmpty()) return null;
		Collections.sort(alleles);
		String key = "";
		int sz = alleles.size();
		for (int n = 0; n < sz; n++)
		{
			Integer i = alleles.get(n);
			if (i == null) return null;
			if (i < 0) return null;
			if (i >= numberAlleles) return null;
			key += i.toString();
			if (i < sz - 1)
			{
				if (isPhased()) key += "|";
				else key += "/";
			}
		}
		return key;
	}
	
	/**
	 * Set a genotype likelihood for a non-diploid genotype. Provide the alleles as a list
	 * of general allele indices. Indices must be valid - if an index exceeds the number of
	 * possible alleles noted for this genotype or its variant, this method will return
	 * without changing anything.
	 * @param alleles List of allele indices making up the genotype to set value for.
	 * @param scaledValue Log 10 scaled probability of the genotype being the true genotype.
	 */
	public void setHeterogeneousGenotypeLikelihood(List<Integer> alleles, double scaledValue)
	{
		String key = generateHeteroGenoMapKey(alleles);
		if (key == null) return;
		heteroPloidyLikelihoods.put(key, scaledValue);
	}
	
	/**
	 * Set a genotype likelihood for a non-diploid genotype. Provide the alleles as a list
	 * of general allele indices. Indices must be valid - if an index exceeds the number of
	 * possible alleles noted for this genotype or its variant, this method will return
	 * without changing anything.
	 * @param alleles List of allele indices making up the genotype to set value for.
	 * @param probability Probability of the genotype being the true genotype.
	 */
	public void setHeterogeneousGenotypeProbability(List<Integer> alleles, double probability)
	{
		String key = generateHeteroGenoMapKey(alleles);
		if (key == null) return;
		double scaledValue = Utilities.probToFloatScale(probability);
		heteroPloidyLikelihoods.put(key, scaledValue);
	}
	
	/**
	 * Set a genotype likelihood for a non-diploid genotype. Provide the alleles as a list
	 * of general allele indices. Indices must be valid - if an index exceeds the number of
	 * possible alleles noted for this genotype or its variant, this method will return
	 * without changing anything.
	 * @param alleles List of allele indices making up the genotype to set value for.
	 * @param phred Phred scaled probability of the genotype being the true genotype.
	 */
	public void setHeterogeneousGenotypeLikelihood(List<Integer> alleles, int phred)
	{
		String key = generateHeteroGenoMapKey(alleles);
		if (key == null) return;
		double scaledValue = Utilities.phredToFloatScale(phred);
		heteroPloidyLikelihoods.put(key, scaledValue);
	}
	
	/**
	 * Set the conditional genotype quality, or the phred-scaled probability that the genotype
	 * call is incorrect.
	 * <br>A value of -1 will be treated as if the variable is unset.
	 * @param phred Phred-scaled quality score
	 */
	public void setGenotypeQuality(int phred)
	{
		this.conditionalQuality = phred;
	}
	
	/**
	 * Set the conditional genotype quality, or the phred-scaled probability that the genotype
	 * call is incorrect.
	 * <br>A phred value of -1 will be treated as if the variable is unset.
	 * @param probability Raw probability that the call is wrong.
	 */
	public void setGenotypeQuality(double probability)
	{
		int phred = Utilities.probToPhred(probability);
		this.conditionalQuality = phred;
	}

	/**
	 * Set the phred-scaled haplotype quality for one or the other haplotype.
	 * @param phred Phred score: -1 will be treated as unset.
	 * @param haplo Which haplotype to set the score for - False = 1, True = 2
	 */
	public void setHaplotypeQuality(int phred, boolean haplo)
	{
		if (this.haploQualities == null) haploQualities = new int[2];
		if (haplo) haploQualities[1] = phred;
		else haploQualities[0] = phred;
	}
	
	/**
	 * Set the ID of the phase set, if this genotype is phased. The ID can be any
	 * positive integer. A value of -1 is considered unset.
	 * This can be set even if the genotype is not phased, but will be ignored.
	 * @param setID ID of genotype phase set
	 */
	public void setPhaseSetID(int setID)
	{
		if (setID < -1) return;
		this.phaseSet = setID;
	}
	
	/**
	 * Set the phred-scaled quality of the genotype phasing. As of VCF standard 4.2, this
	 * field isn't actually fully defined, but it can be set nonetheless.
	 * <br>A value of -1 is treated as unset.
	 * @param phred Phred-scaled phasing quality
	 */
	public void setPhaseQuality(int phred)
	{
		this.phaseQuality = phred;
	}

	/**
	 * Set the expected count for an alternate allele given its alt allele index (general
	 * allele index minus one - its index relative to a list of only alt alleles).
	 * <br>A count value of -1 is treated as unset for that allele.
	 * <br>Negative counts less than -1 will be rejected. The alt allele index must be valid,
	 * and not exceed the total number of alleles known by this genotype object.
	 * @param altAlleleIndex Alternate allele index of allele to set count for.
	 * @param count Expected allele count.
	 */
	public void setExpectedAlternateAlleleCount(int altAlleleIndex, int count)
	{
		if (count < -1) return;
		if (altAlleleIndex < 0) return;
		if (altAlleleIndex + 1 >= this.numberAlleles) return;
		expectedAltCounts[altAlleleIndex] = count;
	}
	
	/**
	 * Set the mapping quality at the position of the variant associated with this genotype.
	 * The value provided should be a phred-scaled integer.
	 * <br>A value of -1 is treated as unset.
	 * @param phred Positional mapping quality
	 */
	public void setMappingQuality(int phred)
	{
		this.mappingQuality = phred;
	}
	
	/**
	 * Set the maximum copy number tracked for this genotype/variant. Copy numbers set
	 * cannot exceed this value. The purpose of this value is for probability recording
	 * and memory management.
	 * A value of 0 is treated as unset, and a value of 4 is default.
	 * Reducing this value from its current state may cause loss of information regarding
	 * likelihoods of possible copy numbers.
	 * @param CNMAX New maximum copy number to track.
	 */
	public void setMaximumCopyNumber(int CNMAX)
	{
		if (CNMAX < 0) return;
		this.maxCopyNumber = CNMAX;
		int len = maxCopyNumber + 1;
		double[] temp = new double[len];
		if (cnLikelihood != null)
		{
			for (int i = 0; i < len; i++)
			{
				if (i >= cnLikelihood.length) temp[i] = 0.0;
				else temp[i] = cnLikelihood[i];
			}	
		}
		cnLikelihood = temp;
	}
	
	/**
	 * Set the copy number of the genotype for a CNV/SV. This field can be set regardless of
	 * whether the associated variant is a CNV or SV. 
	 * A value of -1 is treated as unset. Values lower than -1 will cause the function to
	 * return without changing anything.
	 * Values above CNMAX will increase CNMAX and cause reallocation of memory accordingly.
	 * @param CN Copy number to set.
	 */
	public void setCopyNumber(int CN)
	{
		if (CN < -1) return;
		if (CN > maxCopyNumber) setMaximumCopyNumber(CN);
		copynumber = CN;
	}
	
	/**
	 * Set the call quality of the copy number call using an unrounded phred-scaled
	 * double precision float.
	 * <br>A value of -1.0 is treated as unset.
	 * @param phred Phred-scaled probability that the copy number call is incorrect.
	 */
	public void setCopyNumberCallQuality(double phred)
	{
		this.cnQuality = phred;
	}
	
	/**
	 * Set the likelihood, as a log 10 scaled probability, that any given CN call
	 * is the true genotype.
	 * @param CN Copy number call of interest.
	 * @param scaledValue Scaled probability to set.
	 */
	public void setCopyNumberCallLikelihood(int CN, double scaledValue)
	{
		if (CN < 0) return;
		if (CN > maxCopyNumber) setMaximumCopyNumber(CN);
		cnLikelihood[CN] = scaledValue;
	}
	
	/**
	 * Set the likelihood that any given CN call is the true genotype.
	 * @param CN Copy number call of interest.
	 * @param phredScore Phred scaled probability.
	 */
	public void setCopyNumberCallLikelihood(int CN, int phredScore)
	{
		if (CN < 0) return;
		if (CN > maxCopyNumber) setMaximumCopyNumber(CN);
		double scaledValue = Utilities.phredToFloatScale(phredScore);
		cnLikelihood[CN] = scaledValue;
	}
	
	/**
	 * Set the likelihood that any given CN call is the true genotype.
	 * @param CN Copy number call of interest
	 * @param p Probability
	 */
	public void setCopyNumberCallProbability(int CN, double p)
	{
		if (CN < 0) return;
		if (CN > maxCopyNumber) setMaximumCopyNumber(CN);
		double scaledValue = Utilities.probToFloatScale(p);
		cnLikelihood[CN] = scaledValue;
	}
	
	/**
	 * Set the phred-scaled probability that the variant is novel.
	 * <br>A value of -1 is treated as unset.
	 * @param p Phred-scaled probability score
	 */
	public void setNovelProbability(int phred)
	{
		this.novelProbability = phred;
	}
	
	/**
	 * Set the unique haplotype ID.
	 * <br>A value of -1 is treated as unset.
	 * @param ID ID number to set
	 */
	public void setHaplotypeID(int ID)
	{
		this.haplotypeIdent = ID;
	}
	
	/**
	 * Set the unique ancestral haplotype ID.
	 * <br>A value of -1 is treated as unset.
	 * @param ID ID number to set
	 */
	public void setAncestralHaplotypeID(int ID)
	{
		this.ancestorHTID = ID;
	}
	
	/* --- Zygosity Shortcuts --- */
	
	/**
	 * Evaluate the allele calls for this genotype and determine the relative
	 * zygosity. See class constants for possible values.
	 * @return Integer value representing a zygosity type. See class constants for valid values.
	 */
	public int getZygosity()
	{
		if (this.alleles == null || alleles.length < 1) return ZYGOSITY_UNKNOWN;
		if (alleles.length < 2)
		{
			int a = alleles[0];
			if (a == 0) return Genotype.ZYGOSITY_HOMOREF_CNV;
			else if (a > 0) return Genotype.ZYGOSITY_HOMOALT_CNV;
			else return ZYGOSITY_UNKNOWN;
		}
		else if (alleles.length > 2)
		{
			boolean hasr = false;
			Set<Integer> alts = new HashSet<Integer>();
			for (int i : alleles)
			{
				if (i == 0) hasr = true;
				if (i > 0) alts.add(i);
			}
			if (alts.isEmpty())
			{
				if (hasr) return Genotype.ZYGOSITY_HOMOREF_CNV;
				else return Genotype.ZYGOSITY_UNKNOWN;
			}
			else
			{
				if (hasr)
				{
					if(alts.size() > 1) return Genotype.ZYGOSITY_HETERORAA_CNV;
					else return Genotype.ZYGOSITY_HETERORA_CNV;
				}
				else
				{
					if (alts.size() > 1) return Genotype.ZYGOSITY_HETEROAA_CNV;
					else return Genotype.ZYGOSITY_HOMOALT_CNV;
				}
			}
		}
		else
		{
			int a1 = alleles[0];
			int a2 = alleles[1];
			if (a1 != a2)
			{
				if (a1 == 0 || a2 == 0) return Genotype.ZYGOSITY_HETERORA;
				else return Genotype.ZYGOSITY_HETEROAA;
			}
			else
			{
				if (a1 == 0) return Genotype.ZYGOSITY_HOMOREF;
				else return Genotype.ZYGOSITY_HOMOALT;
			}
		}
		
		//return ZYGOSITY_UNKNOWN;
	}
	
	/**
	 * Get the number of alleles that have been called for this genotype.
	 * This should resemble the copy number, if that has been noted, but it may not.
	 * To be sure what the copy number is, in the case of a CNV call, it might be better to 
	 * call getCopyNumber().
	 * @return The number of alleles called for this genotype.
	 */
	public int getAlleleCallCount()
	{
		//Should produce something similar to copy number if alleles have been called!
		if (this.alleles == null) return 0;
		return this.alleles.length;
	}
	
	/**
	 * Get an approximate percentage of alleles called that are alternate alleles.
	 * @return Decimal percentage of alleles that are alternate alleles.
	 */
	public double getPercentAlt()
	{
		if (alleles == null) return 0;
		int t = alleles.length;
		if (alleles.length < 1) return 0;
		
		int acount = 0;
		for (int a : alleles) {
			if (a > 0) acount++;
			if (a < 0) t--;
		}
		
		if (t <= 0) return 0.0;
		double percent = (double)acount / (double)t;
		percent *= 100.0;
		
		return percent;
	}
	
	/**
	 * Get whether the genotype is unknown (usually denoted "./." in a VCF.
	 * @return True - If there are no alleles stored or all alleles are "undefined". False - If there is at least
	 * one valid allele.
	 */
	public boolean isGenotypeUnknown()
	{
		if (this.alleles == null) return true;
		if (this.alleles.length < 1) return true;
		for (int i : alleles)
		{
			if (i >= 0) return false;
		}
		return true;
	}
	
	/**
	 * Get whether the genotype is a homozygous genotype or not.
	 * @return True if...
	 * <br>1. The genotype is homozygous. OR...
	 * <br>2. The genotype is unknown. OR...
	 * <br>3. The genotype is CNV1.
	 * <br>Otherwise, false.
	 */
	public boolean isHomozygous()
	{
		if (this.alleles == null || alleles.length < 1) return true; //I mean... it's not het?
		if (isGenotypeUnknown()) return true; //Ehhhh...
		int a = alleles[0];
		int j = 1;
		while (a < 0 && j < alleles.length)
		{
			a = alleles[j];
			j++;
		}
		if (a < 0) return true; //Should be caught by geno unknown, BUT just in case...
		for (int i : alleles)
		{
			if (i < 0) continue; //Ignore unknown alleles
			if (i != a) return false;
		}
		return true;
	}
	
	/**
	 * Check whether this genotype has a particular allele (given an integer index >= 0)
	 * @param allele Allele to check
	 * @return True if this genotype has that allele, false if it does not or the provided
	 * allele is negative.
	 */
	public boolean hasAllele(int allele)
	{
		if (allele < 0) return false;
		if (alleles == null) return false;
		if (alleles.length < 1) return false;
		for (int i : alleles)
		{
			if (i == allele) return true;
		}
		return false;
	}
	
	/* --- Definition Mapping --- */
	
	/**
	 * Get the remainder of a definition for a standard genotype field, knowing its key.
	 * @param key The key to the field of interest
	 * @return The InfoDefinition for that field, if known. If field is not standard and
	 * not recorded in this class, this method returns null.
	 */
	public static InfoDefinition getStandardDefinition(String key)
	{
		return defMap.get(key);
	}
	
	/* --- Array Sizing --- */
	
	private static int numberToCombo(int number)
	{
		int tot = 0;
		while (number > 0)
		{
			tot += number;
			number--;
		}
		return tot;
	}
	
	private static int comboToNumber(int combo)
	{
		int n = 1;
		while (combo > 0)
		{
			combo -= n;
			n++;
		}
		return (n - 1);
	}
	
	/* --- Serialization --- */
	
	/**
	 * Get a VCF formatted String representation of this genotype with
	 * the requested fields in the order provided by the argument collection.
	 * Any requested field for which there is no value will simply be recorded as "."
	 * @param fields Ordered keys of fields to include
	 * @return VCF formatted genotype field
	 */
	public String toVCFField(String[] fields)
	{
		if (fields == null) return "";
		if (fields.length < 1) return "";
		String s = "";
		for (int i = 0; i < fields.length; i++)
		{
			String val = getField(fields[i]);
			if (val == null) s += ".";
			else s += val;
			if (i < fields.length - 1) s += ":";
		}
		return s;
	}
	
	/**
	 * Get a VCF formatted String representation of this genotype with
	 * the requested fields in the order provided by the argument collection.
	 * Any requested field for which there is no value will simply be recorded as "."
	 * @param fields Ordered keys of fields to include
	 * @return VCF formatted genotype field
	 */
	public String toVCFField(List<String> fields)
	{
		if (fields == null) return "";
		if (fields.isEmpty()) return "";
		String s = "";
		for (int i = 0; i < fields.size(); i++)
		{
			String val = getField(fields.get(i));
			if (val == null) s += ".";
			else s += val;
			if (i < fields.size() - 1) s += ":";
		}
		return s;
	}
	
}
