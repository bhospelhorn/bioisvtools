package hospelhornbg_bioinformatics;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import hospelhornbg_bioinformatics.VariantPool.InfoDefinition;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

/*
 * UPDATE NOTES
 * 
 * Initial version date: December 20, 2017 (1.0.0)
 * 
 * 1.0.0 -> 1.1.0 | January 12, 2018
 * 	Made defo constructor public instead of protected.
 * 	Added and updated setters
 * 	Removed/ protected setters for LUMPY set fields - these should only be editable by LUMPY, of course!
 * 	Removed INFO field annotations for fields that are parsed as instance variables
 * 	Javadoc annotation
 * 	Cleaned out Lumpy score calculation - this has been delegated to another class.
 * 	Made available standard InfoDefs for VCF writing
 * 	Used convoluted field getter/setter system like with Genotype to override string based INFO field access for parsed
 * 		instance variable fields
 * 	Added BED and BEDPE serializers
 * 	Removed BEDPE parsing (it's just too vague)
 * 
 * 1.1.0 -> 1.1.1 | January 24, 2018
 * 	Added null checks for parsing SV INFO fields - wasn't handling properly when field wasn't present.
 * 	Added genome browser BED12 serializer
 * 
 * 1.1.1 -> 1.1.2 | February 6, 2018
 * 	Tweaked field getters for CIEND values to return null if there is no END set.
 * 
 * 1.1.2 -> 1.1.3 | February 9, 2018
 * 	If you attempt to set the prob curve to null, it won't evaluate its contents.
 * 
 * 1.1.3 -> 1.2.0 | February 20, 2018
 * 	Update for Contig/GenomeBuild addition.
 * 
 * 1.2.0 -> 1.2.1 | February 21, 2018
 * 	Added static function that automatically adds standard SV defs to a VariantPool
 * 	Altered BED trackmaker line generator to not throw out homozygous ref automatically.
 * 
 * 1.2.1 -> 1.2.2 | February 26, 2018
 * 	Definition maps can't be accessed if they haven't been instantiated/populated!
 * 
 * 1.2.2 -> 1.3.0 | March 5, 2018
 * 	Update for SURVIVOR compatibility - added TRA type.
 * 
 * 1.3.0 -> 1.3.1 | March 26, 2018
 * 	The get minimum and maximum lengths (for multiallelic variants) used string length to calculate size. This does not apply for SVs - 
 * should use SVlen instead!
 * 
 * 1.3.1 -> 1.4.0 | April 9, 2018
 * 	Added the gene function field for GeneFunc enum (slower parsing, but faster filtering!)
 * 	This is the location effect (eg. exonic, splicing...)
 * 
 * 1.4.0 -> 1.4.1 | April 16, 2018
 * 	Added split to/merge from BND pairs for VCF processing.
 * 
 * 1.4.1 -> 1.4.2 | April 19, 2018
 * 	Compatibility with Variant version 1.2.2
 * 	Added some null pointer checks on the merge/split methods
 * 
 * 1.4.2 -> 1.4.3 | April 23, 2018
 * 	Better updating of SVLEN
 * 
 * 1.4.3 -> 1.4.4 | May 2, 2018
 * 	Added a function for setting gene func
 * 
 * 1.4.4 -> 1.4.5 | July 12, 2018
 * 	Change in Genotype class was forcing the BED viewer function to toss unknown genotypes.
 * 	No more!
 * 
 * 1.4.5 -> 1.4.6 | July 20, 2018
 * 	Added inRegion method
 * 
 * 1.4.6 -> 1.5.0 | August 10, 2018
 * 	Moved the GeneFunc variable to the superclass
 * 
 * 1.5.0 -> 1.5.1 | August 17, 2018
 * 	Added getEndChromosome function
 */

/*
 * Possible future improvements:
 * 	- Mate list lists references to Variant objects instead of ID strings?
 * 	- There is a lot of code repetition when it comes to similar fields. Slim that down.
 */

/**
 * Container object extending the standard Variant to include information and methods
 * for easier processing of Structural Variants.
 * @author Blythe Hospelhorn
 * @version 1.5.1
 * @since August 17, 2018
 *
 */
public class StructuralVariant extends Variant implements Comparable<Variant>{
	
	/* --- Constants --- */
	
	//Definition maps
	private static Map<String, String> altDefMap;
	private static Map<String, InfoDefinition> infoDefMap;
	
	//From the VCF format specification version 4.2
	
	/**
	 * Alt allele type representing a deletion.
	 */
	public static final String INFODEF_ALT_DEL = "DEL";
	
	/**
	 * Alt allele type representing a de novo insertion.
	 */
	public static final String INFODEF_ALT_INS = "INS";
	
	/**
	 * Alt allele type representing a duplication.
	 */
	public static final String INFODEF_ALT_DUP = "DUP";
	
	/**
	 * Alt allele type representing an inversion.
	 */
	public static final String INFODEF_ALT_INV = "INV";
	
	/**
	 * Alt allele type representing a general copy number variation.
	 */
	public static final String INFODEF_ALT_CNV = "CNV";
	
	/**
	 * Alt allele type representing a tandem duplication.
	 */
	public static final String INFODEF_ALT_TANDEM = "DUP:TANDEM";
	
	/**
	 * Alt allele type representing a mobile element deletion.
	 */
	public static final String INFODEF_ALT_DELME = "DEL:ME";
	
	/**
	 * Alt allele type representing a mobile element insertion.
	 */
	public static final String INFODEF_ALT_INSME = "INS:ME";
	
	/**
	 * ID = SVTYPE
	 * <br>Number = 1
	 * <br>Type = String
	 * <br>Description = "Type of structural variant"
	 */
	public static final InfoDefinition INFODEF_INFO_SVTYPE = new InfoDefinition("SVTYPE", VariantPool.INFODEF_STRING, "Type of structural variant", 1);
	
	/**
	 * ID = SVLEN
	 * <br>Number = A
	 * <br>Type = Integer
	 * <br>Description = "Difference in length between REF and ALT alleles"
	 */
	public static final InfoDefinition INFODEF_INFO_SVLEN = new InfoDefinition("SVLEN", VariantPool.INFODEF_INT, "Difference in length between REF and ALT alleles", VariantPool.INFODEF_NARGS_PERALT);
	
	/**
	 * ID = END
	 * <br>Number = 1
	 * <br>Type = Integer
	 * <br>Description = "End position of the variant described in this record"
	 */
	public static final InfoDefinition INFODEF_INFO_END = new InfoDefinition("END", VariantPool.INFODEF_INT, "End position of the variant described in this record", 1);
	
	/**
	 * ID = PRECISE
	 * <br>Number = 0
	 * <br>Type = Flag
	 * <br>Description = "Precise structural variation"
	 */
	public static final InfoDefinition INFODEF_INFO_PRECISE = new InfoDefinition("PRECISE", VariantPool.INFODEF_FLAG, "Precise structural variation", 0);
	
	/**
	 * ID = IMPRECISE
	 * <br>Number = 0
	 * <br>Type = Flag
	 * <br>Description = "Imprecise structural variation"
	 */
	public static final InfoDefinition INFODEF_INFO_IMPRECISE = new InfoDefinition("IMPRECISE", VariantPool.INFODEF_FLAG, "Imprecise structural variation", 0);
	
	/**
	 * ID = NOVEL
	 * <br>Number = 0
	 * <br>Type = Flag
	 * <br>Description = "Indicates a novel structural variation"
	 */
	public static final InfoDefinition INFODEF_INFO_NOVEL = new InfoDefinition("NOVEL", VariantPool.INFODEF_FLAG, "Indicates a novel structural variation", 0);
	
	/**
	 * ID = CIPOS
	 * <br>Number = 2
	 * <br>Type = Integer
	 * <br>Description = "Confidence interval around POS for imprecise variants"
	 */
	public static final InfoDefinition INFODEF_INFO_CIPOS = new InfoDefinition("CIPOS", VariantPool.INFODEF_INT, "Confidence interval around POS for imprecise variants", 2);
	
	/**
	 * ID = CIEND
	 * <br>Number = 2
	 * <br>Type = Integer
	 * <br>Description = "Confidence interval around END for imprecise variants"
	 */
	public static final InfoDefinition INFODEF_INFO_CIEND = new InfoDefinition("CIEND", VariantPool.INFODEF_INT, "Confidence interval around END for imprecise variants", 2);
	
	/**
	 * ID = CIPOS95
	 * <br>Number = 2
	 * <br>Type = Integer
	 * <br>Description = "Confidence interval (95%) around POS for imprecise variants"
	 */
	public static final InfoDefinition INFODEF_INFO_CIPOS95 = new InfoDefinition("CIPOS95", VariantPool.INFODEF_INT, "Confidence interval (95%) around POS for imprecise variants", 2);
	
	/**
	 * ID = CIEND95
	 * <br>Number = 2
	 * <br>Type = Integer
	 * <br>Description = "Confidence interval (95%) around END for imprecise variants"
	 */
	public static final InfoDefinition INFODEF_INFO_CIEND95 = new InfoDefinition("CIEND95", VariantPool.INFODEF_INT, "Confidence interval (95%) around END for imprecise variants", 2);
	
	/**
	 * ID = MATEID
	 * <br>Number = .
	 * <br>Type = String
	 * <br>Description = "ID of mate breakends"
	 */
	public static final InfoDefinition INFODEF_INFO_MATEID = new InfoDefinition("MATEID", VariantPool.INFODEF_STRING, "ID of mate breakends", VariantPool.INFODEF_NARGS_VARIABLE);
	
	/**
	 * ID = EVENT
	 * <br>Number = 1
	 * <br>Type = String
	 * <br>Description = "ID of event associated to breakend"
	 */
	public static final InfoDefinition INFODEF_INFO_EVENT = new InfoDefinition("EVENT", VariantPool.INFODEF_STRING, "ID of event associated to breakend", 1);
	
	/**
	 * ID = MAPQ
	 * <br>Number = 1
	 * <br>Type = Integer
	 * <br>Description = "Median mapping quality of paired-ends"
	 */
	public static final InfoDefinition INFODEF_INFO_MAPQ = new InfoDefinition("MAPQ", VariantPool.INFODEF_INT, "Median mapping quality of paired-ends", 1);
	
	/**
	 * ID = RE
	 * <br>Number = 1
	 * <br>Type = Integer
	 * <br>Description = "Read support"
	 */
	public static final InfoDefinition INFODEF_INFO_RE = new InfoDefinition("RE", VariantPool.INFODEF_INT, "Read support", 1);
	
	/**
	 * ID = SVMETHOD
	 * <br>Number = 1
	 * <br>Type = String
	 * <br>Description = "Type of approach used to detect SV"
	 */
	public static final InfoDefinition INFODEF_INFO_SVMETHOD = new InfoDefinition("SVMETHOD", VariantPool.INFODEF_STRING, "Type of approach used to detect SV", 1);
	
	/**
	 * ID = TRUESVTYPE
	 * <br>Number = 1
	 * <br>Type = String
	 * <br>Description = "Type of structural variant before BND split"
	 */
	public static final InfoDefinition INFODEF_INFO_TRUESVTYPE = new InfoDefinition("TRUESVTYPE", VariantPool.INFODEF_STRING, "Type of structural variant before BND split", 1);
	
	
	//LUMPY additional fields
	
	/**
	 * ID = STRANDS
	 * <br>Number = .
	 * <br>Type = String
	 * <br>Description = "Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)"
	 */
	public static final InfoDefinition INFODEF_INFO_STRANDS = new InfoDefinition("STRANDS", VariantPool.INFODEF_STRING, "Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)", VariantPool.INFODEF_NARGS_VARIABLE);
	
	/**
	 * ID = SECONDARY
	 * <br>Number = 0
	 * <br>Type = Flag
	 * <br>Description = "Secondary breakend in a multi-line variant"
	 */
	public static final InfoDefinition INFODEF_INFO_SECONDARY = new InfoDefinition("SECONDARY", VariantPool.INFODEF_FLAG, "Secondary breakend in a multi-line variant", 0);
	
	/**
	 * ID = SU
	 * <br>Number = .
	 * <br>Type = Integer
	 * <br>Description = "Number of pieces of evidence supporting the variant across all samples"
	 */
	public static final InfoDefinition INFODEF_INFO_SU = new InfoDefinition("SU", VariantPool.INFODEF_INT, "Number of pieces of evidence supporting the variant across all samples", VariantPool.INFODEF_NARGS_VARIABLE);
	
	/**
	 * ID = PE
	 * <br>Number = .
	 * <br>Type = Integer
	 * <br>Description = "Number of paired-end reads supporting the variant across all samples"
	 */
	public static final InfoDefinition INFODEF_INFO_PE = new InfoDefinition("PE", VariantPool.INFODEF_INT, "Number of paired-end reads supporting the variant across all samples", VariantPool.INFODEF_NARGS_VARIABLE);
	
	/**
	 * ID = SR
	 * <br>Number = .
	 * <br>Type = Integer
	 * <br>Description = "Number of split reads supporting the variant across all samples"
	 */
	public static final InfoDefinition INFODEF_INFO_SR = new InfoDefinition("SR", VariantPool.INFODEF_INT, "Number of split reads supporting the variant across all samples", VariantPool.INFODEF_NARGS_VARIABLE);
	
	/**
	 * ID = BD
	 * <br>Number = .
	 * <br>Type = Integer
	 * <br>Description = "Amount of BED evidence supporting the variant across all samples"
	 */
	public static final InfoDefinition INFODEF_INFO_BD = new InfoDefinition("BD", VariantPool.INFODEF_INT, "Amount of BED evidence supporting the variant across all samples", VariantPool.INFODEF_NARGS_VARIABLE);
	
	/**
	 * ID = EV
	 * <br>Number = .
	 * <br>Type = String
	 * <br>Description = "Type of LUMPY evidence contributing to the variant call"
	 */
	public static final InfoDefinition INFODEF_INFO_EV = new InfoDefinition("EV", VariantPool.INFODEF_STRING, "Type of LUMPY evidence contributing to the variant call", VariantPool.INFODEF_NARGS_VARIABLE);
	
	/**
	 * ID = PRPOS
	 * <br>Number = .
	 * <br>Type = Float
	 * <br>Description = "LUMPY probability curve of the POS breakend"
	 */
	public static final InfoDefinition INFODEF_INFO_PRPOS = new InfoDefinition("PRPOS", VariantPool.INFODEF_FLOAT, "LUMPY probability curve of the POS breakend", VariantPool.INFODEF_NARGS_VARIABLE);
	
	/**
	 * ID = PREND
	 * <br>Number = .
	 * <br>Type = Float
	 * <br>Description = "LUMPY probability curve of the END breakend"
	 */
	public static final InfoDefinition INFODEF_INFO_PREND = new InfoDefinition("PREND", VariantPool.INFODEF_FLOAT, "LUMPY probability curve of the END breakend", VariantPool.INFODEF_NARGS_VARIABLE);
	
	// SURVIVOR additional fields
	
	/**
	 * ID = SUPP
	 * <br>Number = 1
	 * <br>Type = Integer
	 * <br>Description = "Number of samples supporting call"
	 */
	public static final InfoDefinition INFODEF_INFO_SUPP = new InfoDefinition("SUPP", VariantPool.INFODEF_INT, "Number of samples supporting call", 1);
	
	/**
	 * ID = SUPP_VEC
	 * <br>Number = 1
	 * <br>Type = String
	 * <br>Description = "Boolean string representing which samples provided evidence of call"
	 */
	public static final InfoDefinition INFODEF_INFO_SUPP_VEC = new InfoDefinition("SUPP_VEC", VariantPool.INFODEF_STRING, "Boolean string representing which samples provided evidence of call", 1);
	
	/**
	 * ID = AVGLEN
	 * <br>Number = 1
	 * <br>Type = Float
	 * <br>Description = "Average length"
	 */
	public static final InfoDefinition INFODEF_INFO_AVGLEN = new InfoDefinition("AVGLEN", VariantPool.INFODEF_FLOAT, "Average length", 1);
	
	/* --- Instance Variables --- */
	
	private int end;
	
	private SVType type;
	private int[] SVlen;
	
	private boolean imprecise;
	
		//All of these are offset from position and end, not true values
	private int CIPosLow;
	private int CIPosHi;
	private int CIEndLow;
	private int CIEndHi;
	private int CI95PosLow;
	private int CI95PosHi;
	private int CI95EndLow;
	private int CI95EndHi;
	
	private ArrayList<String> mateid;
	private String event;
	
		//LUMPY Stuff
	
	//private String strands;
	private boolean secondary;
	
	private int[] evidence_SU;
	private int[] evidence_PE;
	private int[] evidence_SR;
	private int[] evidence_BD;
	
	private double[] prPos;
	private double[] prEnd;
	
	private double highestPRPOS;
	private double highestPREND;
	
	
	/* --- Construction --- */
	
	/**
	 * Construct an empty structural variant based upon an empty variant.
	 * All values will be set to their defaults.
	 */
	public StructuralVariant()
	{
		super();
		setDefaults();
		populateFieldMap();
	}
	
	/**
	 * Construct a structural variant from a VCF record and an ordered list of samples.
	 * @param VCF_line The ASCII VCF record line.
	 * @param samples List of samples in the order they appear in the VCF.
	 * @throws UnsupportedFileTypeException If an INFO field needed to determine the natural of the
	 * structural variant could not be parsed properly.
	 */
	public StructuralVariant(String VCF_line, List<String> samples, GenomeBuild genome, boolean allowBuildModification) throws UnsupportedFileTypeException
	{
		super(VCF_line, samples, genome, allowBuildModification);
		setDefaults();
		parseStandardSVFields();
		parseLUMPYSVFields();
		populateFieldMap();
	}
	
	/**
	 * Create a structural variant by merging a BreakendPair that was split by
	 * this same class.
	 * @param mergePair BreakendPair variant to merge into a single SV record.
	 */
	public StructuralVariant(BreakendPair mergePair)
	{
		StructuralVariant sv1 = mergePair.getBNDVariant(false);
		StructuralVariant sv2 = mergePair.getBNDVariant(true);
		
		//Check for split tag
		String origType = sv1.getSingleInfoEntry(StructuralVariant.INFODEF_INFO_TRUESVTYPE.getKey());
		if (origType == null) throw new IllegalArgumentException();
		SVType t = SVType.getType(origType);
		if (t == null || t == SVType.OTHER) throw new IllegalArgumentException();
		
		//Positions & Ends
		this.setChromosome(sv1.getChromosome());
		this.setPosition(sv1.getPosition());
		this.setEndPosition(sv2.getPosition());
		
		//Other variant instance properties
		this.setVariantName(sv1.getEventID());
		this.setRefAllele(sv1.getRefAllele());
		this.addAltAllele(sv1.getAltAllele(0));
		this.setFilterPass(sv1.passedAllFilters());
		this.setQuality(sv1.getQuality());
		String[] filters = this.getFiltersFailed();
		if(filters != null)
		{
			for (String s : filters) this.addFailedFilter(s);
		}
		
		//SV Basics
		this.setType(t);
		this.setSVLength(sv1.getSVLength());
		this.setImprecise(sv1.isImprecise());
		
		//CIs
		this.setCIDiff(sv1.getCIDiff(false, false, false), false, false, false);
		this.setCIDiff(sv1.getCIDiff(false, false, true), false, false, true);
		this.setCIDiff(sv1.getCIDiff(false, true, false), false, true, false);
		this.setCIDiff(sv1.getCIDiff(false, true, true), false, true, true);
		this.setCIDiff(sv2.getCIDiff(false, false, false), true, false, false);
		this.setCIDiff(sv2.getCIDiff(false, false, true), true, false, true);
		this.setCIDiff(sv2.getCIDiff(false, true, false), true, true, false);
		this.setCIDiff(sv2.getCIDiff(false, true, true), true, true, true);
		
		//MateID and Event ID
			//Not copied
		
		//Other sv instance properties
		this.setSecondary(false);
		this.evidence_SU = sv1.evidence_SU;
		this.evidence_SR = sv1.evidence_SR;
		this.evidence_PE = sv1.evidence_PE;
		this.evidence_BD = sv1.evidence_BD;
		this.prPos = sv1.prPos;
		this.highestPRPOS = sv1.highestPRPOS;
		this.prEnd = sv2.prPos;
		this.highestPREND = sv2.highestPREND;
		this.setGeneFunction(sv1.getGeneFunction());
		//this.function = sv1.function;
		
		//Info fields
		Set<String> keys = getAllInfoKeys();
		if (keys != null)
		{
			for(String k : keys) this.addInfoField(k, sv1.getInfoEntry(k));
		}
		
		//Genotypes
		String[] gfields = getOrderedGenotypeFields();
		if(gfields != null)
		{
			for (String g : gfields) this.addGenotypeField(g);	
		}
		Set<String> samps = this.getAllGenotypedSamples();
		if(samps != null)
		{
			for (String s : samps) this.addGenotype(s, sv1.getSampleGenotype(s));	
		}
		
	}
	
	/**
	 * PROTECTED constructor for constructing a Structural Variant from an existing Variant
	 * by copying the REFERENCES to all internal structures.
	 * @param source Variant to copy
	 */
	protected StructuralVariant(Variant source)
	{
		super(source);
		setDefaults();
		parseStandardSVFields();
		parseLUMPYSVFields();
		populateFieldMap();
	}
	
	/**
	 * PROTECTED constructor for constructing a new Structural Variant from an existing StructuralVariant
	 * by copying the REFERENCES to all internal structures.
	 * @param source Variant to copy
	 */
	protected StructuralVariant(StructuralVariant source)
	{
		super(source);
		end = source.end;
		type = source.type;
		SVlen = source.SVlen;
		imprecise = source.imprecise;
		CIPosLow = source.CIPosLow;
		CIPosHi = source.CIPosHi;
		CIEndLow = source.CIEndLow;
		CIEndHi = source.CIEndHi;
		CI95PosLow = source.CI95PosLow;
		CI95PosHi = source.CI95PosHi;
		CI95EndLow = source.CI95EndLow;
		CI95EndHi = source.CI95EndHi;
		mateid = source.mateid;
		event = source.event;
		secondary = source.secondary;
		evidence_SU = source.evidence_SU;
		evidence_PE = source.evidence_PE;
		evidence_SR = source.evidence_SR;
		evidence_BD = source.evidence_BD;
		prPos = source.prPos;
		prEnd = source.prEnd;
		highestPRPOS = source.highestPRPOS;
		highestPREND = source.highestPREND;
		
		populateFieldMap();
	}
	
	/**
	 * Construct a structural variant from a BED line, previously knowing the variant type.
	 * This parser only recognizes the standard BED fields.
	 * @param BED_line Record from a BED file as an ASCII line
	 * @param varType Type of Structural Variant that the record represents
	 * @throws UnsupportedFileTypeException If there is an error parsing the line like a
	 * BED record.
	 */
	public StructuralVariant(String BED_line, SVType varType, GenomeBuild genome) throws UnsupportedFileTypeException
	{
		setDefaults();
		parseBEDLine(BED_line, varType, genome);
		populateFieldMap();
	}
	
	private void setDefaults()
	{
		end = -1;
		type = SVType.OTHER;
		SVlen = new int[1];
		SVlen[0] = 0;
		//strands = ".";
		imprecise = false;
		CIPosLow = 0;
		CIPosHi = 0;
		CIEndLow = 0;
		CIEndHi = 0;
		CI95PosLow = 0;
		CI95PosHi = 0;
		CI95EndLow = 0;
		CI95EndHi = 0;
		
		mateid = new ArrayList<String>(4);
		event = "";
		secondary = false;
		
		evidence_SU = null;
		evidence_PE = null;
		evidence_SR = null;
		evidence_BD = null;
		
		prPos = null;
		prEnd = null;
		
		highestPRPOS = -1;
		highestPREND = -1;
		
		//function = null;
	}
	
	/* --- Parsing --- */
	
	private void parseBEDLine(String BED_line, SVType varType, GenomeBuild genome) throws UnsupportedFileTypeException
	{
		String[] fields = BED_line.split("\t");
		String chrom = fields[0];
		String stPos = fields[1];
		String edPos = fields[2];
		//String name = "";
		//String strand = "";
		//if (fields.length >= 4) name = fields[3];
		//if (fields.length >= 6) strand = fields[5];
		super.setChromosome(genome.getContig(chrom));
		try
		{
			int st = Integer.parseInt(stPos) + 1;
			int ed = Integer.parseInt(edPos);
			super.setPosition(st);
			end = ed;
		}
		catch (NumberFormatException e)
		{
			throw new FileBuffer.UnsupportedFileTypeException();
		}
		//strands = strand;
		type = varType;
		SVlen = new int[1];
		SVlen[0] = end - super.getPosition();
		
		super.addInfoField("SVTYPE", type.toString());
		super.addInfoField("SVLEN", Integer.toString(SVlen[0]));
		super.addInfoField("END", Integer.toString(end));
		//if (!strands.isEmpty()) super.addInfoField("STRANDS", strands);
	}

	private void parseStandardSVFields()
	{
		
		String key = INFODEF_INFO_SVTYPE.getKey();
		type = SVType.getType(super.getSingleStringInfoEntry(key));
		if (type == null) type = SVType.OTHER;
		super.removeInfoField(key);
		
		key = INFODEF_INFO_SVLEN.getKey();
		int aAll = super.countAltAlleles();
		if (aAll < 1) aAll = 1;
		final int[] svlen = super.getIntInfoEntry(key);
		if (svlen == null || svlen.length < 1)
		{
			SVlen = new int[aAll];
			SVlen[0] = -1;
		}
		else
		{
			SVlen = new int[svlen.length];
			for (int i = 0; i < svlen.length; i++) SVlen[i] = svlen[i];
		}
		super.removeInfoField(key);
		
		key = INFODEF_INFO_END.getKey();
		this.setEndPosition(super.getSingleIntInfoEntry(key));
		//end = super.getSingleIntInfoEntry(key);
		super.removeInfoField(key);
		
		key = INFODEF_INFO_IMPRECISE.getKey();
		imprecise = super.getInfoFlag(key);
		super.removeInfoFlag(key);
		
		int[] ivals;
		
		key = INFODEF_INFO_CIPOS.getKey();
		ivals = super.getIntInfoEntry(key);
		if (ivals != null && ivals.length == 2) {
			CIPosLow = ivals[0];
			CIPosHi = ivals[1];
		}
		super.removeInfoField(key);
		
		key = INFODEF_INFO_CIEND.getKey();
		ivals = super.getIntInfoEntry(key);
		if (ivals != null && ivals.length == 2) {
			CIEndLow = ivals[0];
			CIEndHi = ivals[1];
		}
		super.removeInfoField(key);
		
		key = INFODEF_INFO_CIPOS95.getKey();
		ivals = super.getIntInfoEntry(key);
		if (ivals != null && ivals.length == 2) {
			CI95PosLow = ivals[0];
			CI95PosHi = ivals[1];
		}
		super.removeInfoField(key);
		
		key = INFODEF_INFO_CIEND95.getKey();
		ivals = super.getIntInfoEntry(key);
		if (ivals != null && ivals.length == 2) {
			CI95EndLow = ivals[0];
			CI95EndHi = ivals[1];
		}
		super.removeInfoField(key);
		
		key = INFODEF_INFO_MATEID.getKey();
		String[] arr = super.getStringInfoEntry(key);
		if (arr != null) for (String s : arr) mateid.add(s);
		super.removeInfoField(key);
		
		key = INFODEF_INFO_EVENT.getKey();
		event = super.getSingleStringInfoEntry(key);
		super.removeInfoField(key);
		
		/*key = GeneSet.INFODEF_INFO_GFUNC.getKey();
		String func = super.getSingleStringInfoEntry(key);
		if (func != null)
		{
			super.removeInfoField(key);
			function = GeneFunc.getFunction(func);
		}*/
		
	}
	
	private void parseLUMPYSVFields()
	{
		String key = INFODEF_INFO_SU.getKey();
		evidence_SU = super.getIntInfoEntry(key);
		super.removeInfoField(key);
		
		key = INFODEF_INFO_PE.getKey();
		evidence_PE = super.getIntInfoEntry(key);
		super.removeInfoField(key);
		
		key = INFODEF_INFO_SR.getKey();
		evidence_SR = super.getIntInfoEntry(key);
		super.removeInfoField(key);
		
		key = INFODEF_INFO_BD.getKey();
		evidence_BD = super.getIntInfoEntry(key);
		super.removeInfoField(key);
		
		key = INFODEF_INFO_SECONDARY.getKey();
		secondary = super.getInfoFlag(key);
		super.removeInfoFlag(key);
		
		key = INFODEF_INFO_PRPOS.getKey();
		prPos = super.getFloatInfoEntry(key);
		if (prPos != null)
		{
			for (int i = 0; i < prPos.length; i++)
			{
				if (prPos[i] > highestPRPOS) highestPRPOS = prPos[i];
			}	
		}
		super.removeInfoField(key);
		
		key = INFODEF_INFO_PREND.getKey();
		prEnd = super.getFloatInfoEntry(key);
		if (prEnd != null)
		{
			for (int i = 0; i < prEnd.length; i++)
			{
				if (prEnd[i] > highestPREND) highestPREND = prEnd[i];
			}	
		}
		super.removeInfoField(key);
	
	}
	
	/* --- Getters --- */
	
	/**
	 * Get the contig the "end" breakpoint lies on. For standard SVs, this is the same
	 * contig as the main one. For translocations and BND pairs, this may be a different chromosome.
	 * @return Contig of end breakpoint.
	 */
	public Contig getEndChromosome()
	{
		return this.getChromosome();
	}
	
	/**
	 * Get the end position of the structural variant.
	 * @return End position of structural variant in one-based coordinates.
	 */
	public int getEndPosition()
	{
		return end;
	}
	
	/**
	 * Get the type of the structural variant. Default is SVType.OTHER...
	 * @return Structural variant type as an SVType enum.
	 */
	public SVType getType()
	{
		return type;
	}
	
	/**
	 * Get the difference in length between the reference allele and the primary
	 * alternate allele.
	 * @return Int representing the difference in length between the ref and first alt allele.
	 * If there are no alt alleles, this will simply return 0.
	 */
	public int getSVLength()
	{
		if (SVlen == null) return 0;
		return SVlen[0];
	}
	
	/**
	 * Get the difference in length between the reference allele and a given alternate
	 * allele. Index is the relative alternate allele index (general allele index - 1).
	 * If index is invalid, this method simply returns 0.
	 * @param altAlleleIndex Alt allele index of SV alt allele to examine length of.
	 * @return Length difference between ref and specified alt allele, or 0 if undetermined.
	 */
	public int getSVLength(int altAlleleIndex)
	{
		if (altAlleleIndex < 0) return 0;
		if (SVlen == null) return 0;
		if (SVlen.length < 1) return 0;
		if (altAlleleIndex >= SVlen.length) return 0;
		return SVlen[altAlleleIndex];
	}

	/**
	 * Get the absolute value of the SVLEN of the primary alternate allele relative to the 
	 * reference allele.
	 * @return Absolute value of (getSVLength())
	 */
	public int getAbsoluteSVLength()
	{
		return Math.abs(this.getSVLength());
	}
	
	/**
	 * Get whether this structural variant has been flagged as imprecise.
	 * If imprecise, then the breakpoints have not been determined at basepair
	 * resolution and the ends are likely accompanied by confidence intervals.
	 * @return True - If this variant is imprecise.
	 * <br>False - If this variant is precise.
	 */
	public boolean isImprecise()
	{
		return imprecise;
	}
	
	/**
	 * Get the one-based chromosome coordinate position of one end of a recorded
	 * confidence interval. Which end of which CI to retrieve is set by the bool
	 * arguments.
	 * @param variantEnd End of the SV to examine CI of. True = END, False = POS
	 * @param narrow Whether to use the tighter (95%) CI. True = CI95, False = CI90
	 * @param top Which end of the CI to get the position of. True = higher, False = lower
	 * @return The CI position value requested as a position coordinate.
	 */
	public int getCIPosition(boolean variantEnd, boolean narrow, boolean top)
	{
		if (variantEnd)
		{
			if (narrow)
			{
				if (top) return end + CI95EndHi;
				else return end + CI95EndLow;
			}
			else
			{
				if (top) return end + CIEndHi;
				else return end + CIEndLow;
			}
		}
		else
		{
			if (narrow)
			{
				if (top) return getPosition() + CI95PosHi;
				else return getPosition() + CI95PosLow;
			}
			else
			{
				if (top) return getPosition() + CIPosHi;
				else return getPosition() + CIPosLow;
			}
		}
	}

	/**
	 * Get the difference between the SV marked end and the end of one of the confidence intervals.
	 * Which end of which CI to retrieve is set by the bool arguments.
	 * @param variantEnd End of the SV to examine CI of. True = END, False = POS
	 * @param narrow Whether to use the tighter (95%) CI. True = CI95, False = CI90
	 * @param top Which end of the CI to get the position of. True = higher, False = lower
	 * @return The CI position value as a difference between the requested CI boundary and
	 * the SV end (POS or END).
	 */
	public int getCIDiff(boolean variantEnd, boolean narrow, boolean top)
	{
		if (variantEnd)
		{
			if (narrow)
			{
				if (top) return CI95EndHi;
				else return CI95EndLow;
			}
			else
			{
				if (top) return CIEndHi;
				else return CIEndLow;
			}
		}
		else
		{
			if (narrow)
			{
				if (top) return CI95PosHi;
				else return CI95PosLow;
			}
			else
			{
				if (top) return CIPosHi;
				else return CIPosLow;
			}
		}
	}

	/**
	 * Get the IDs of all variants that have been "paired" to this one.
	 * @return A List of all mate IDs recorded for this variant, or null
	 * if there are none.
	 */
	public List<String> getMateIDs()
	{
		if (mateid == null) return null;
		List<String> idlist = new ArrayList<String>(mateid.size());
		for (String s : mateid) idlist.add(s);
		return idlist;
	}

	/**
	 * Get the break event ID for this variant. The ID denotes the break event that this
	 * variant is a part of, especially if this variant is a BND type and involves joining
	 * two distant or interchromosomal positions.
	 * @return String representing the break event ID for this variant. This string may be
	 * empty or null if the event ID is not set.
	 */
	public String getEventID()
	{
		return event;
	}
	
	/**
	 * Get whether this SV has been flagged as secondary - that is, the second 
	 * in a breakend pair. If this SV is secondary, it should probably be of type BND, 
	 * and it should probably have a partner somewhere. It may not have all the data
	 * required by itself (such as END) to use as a regular SV.
	 * @return True - If this variant is secondary
	 * <br>False - Otherwise
	 */
	public boolean isSecondary()
	{
		return secondary;
	}
	
	/**
	 * Get the total evidence count for the primary sample evidence counts have
	 * been recorded for in this variant.
	 * <br>A value of -1 indicates that this field is not set for any sample.
	 * @return LUMPY total evidence count for primary sample
	 */
	public int getTotalEvidenceCount()
	{
		if (evidence_SU == null) return -1;
		return evidence_SU[0];
	}

	/**
	 * Get the total evidence count for the sample with the specified sample index.
	 * The sample index is the relative ordered position the counts for this sample
	 * were recorded. At the moment, this value must be known by the user or code
	 * calling it; it is not parsed from a VCF.
	 * <br>A value of -1 indicates that this field is not set for the specified sample.
	 * @param sampleIndex
	 * @return LUMPY total evidence count for sample (sampleIndex) if known.
	 */
	public int getTotalEvidenceCount(int sampleIndex)
	{
		if (evidence_SU == null) return -1;
		if (sampleIndex < 0) return -1;
		if (sampleIndex > evidence_SU.length) return -1;
		return evidence_SU[sampleIndex];
	}
	
	/**
	 * Get the paired-end evidence count for the primary sample evidence counts have
	 * been recorded for in this variant.
	 * <br>A value of -1 indicates that this field is not set for any sample.
	 * @return LUMPY paired-end evidence count for primary sample
	 */
	public int getPEEvidenceCount()
	{
		if (evidence_PE == null) return -1;
		return evidence_PE[0];
	}
	
	/**
	 * Get the paired-end evidence count for the sample with the specified sample index.
	 * The sample index is the relative ordered position the counts for this sample
	 * were recorded. At the moment, this value must be known by the user or code
	 * calling it; it is not parsed from a VCF.
	 * <br>A value of -1 indicates that this field is not set for the specified sample.
	 * @param sampleIndex
	 * @return LUMPY paired-end evidence count for sample (sampleIndex) if known.
	 */
	public int getPEEvidenceCount(int sampleIndex)
	{
		if (evidence_PE == null) return -1;
		if (sampleIndex < 0) return -1;
		if (sampleIndex > evidence_PE.length) return -1;
		return evidence_PE[sampleIndex];
	}
	
	/**
	 * Get the split read evidence count for the primary sample evidence counts have
	 * been recorded for in this variant.
	 * <br>A value of -1 indicates that this field is not set for any sample.
	 * @return LUMPY split read evidence count for primary sample
	 */
	public int getSREvidenceCount()
	{
		if (evidence_SR == null) return -1;
		return evidence_SR[0];
	}
	
	/**
	 * Get the split read evidence count for the sample with the specified sample index.
	 * The sample index is the relative ordered position the counts for this sample
	 * were recorded. At the moment, this value must be known by the user or code
	 * calling it; it is not parsed from a VCF.
	 * <br>A value of -1 indicates that this field is not set for the specified sample.
	 * @param sampleIndex
	 * @return LUMPY split read evidence count for sample (sampleIndex) if known.
	 */
	public int getSREvidenceCount(int sampleIndex)
	{
		if (evidence_SR == null) return -1;
		if (sampleIndex < 0) return -1;
		if (sampleIndex > evidence_SR.length) return -1;
		return evidence_SR[sampleIndex];
	}
	
	/**
	 * Get the BED evidence count for the primary sample evidence counts have
	 * been recorded for in this variant.
	 * <br>A value of -1 indicates that this field is not set for any sample.
	 * @return LUMPY BED evidence count for primary sample
	 */
	public int getBDEvidenceCount()
	{
		if (evidence_BD == null) return -1;
		return evidence_BD[0];
	}
	
	/**
	 * Get the BED evidence count for the sample with the specified sample index.
	 * The sample index is the relative ordered position the counts for this sample
	 * were recorded. At the moment, this value must be known by the user or code
	 * calling it; it is not parsed from a VCF.
	 * <br>A value of -1 indicates that this field is not set for the specified sample.
	 * @param sampleIndex
	 * @return LUMPY BED evidence count for sample (sampleIndex) if known.
	 */
	public int getBDEvidenceCount(int sampleIndex)
	{
		if (evidence_BD == null) return -1;
		if (sampleIndex < 0) return -1;
		if (sampleIndex > evidence_BD.length) return -1;
		return evidence_BD[sampleIndex];
	}
	
	/**
	 * Get a copy of the LUMPY probability curve at one of the variant ends.
	 * If there is no curve recorded at that end, this method will return null.
	 * @param variantEnd End of SV to get the breakpoint probability curve of. True = END, False = POS
	 * @return Probability curve at breakpoint at SV end, if present. 
	 * <br>null if there is no curve recorded.
	 */
	public double[] getProbabilityCurve(boolean variantEnd)
	{
		double[] acopy;
		if (!variantEnd)
		{
			if (prPos == null) return null;
			acopy = new double[prPos.length];
			for (int i = 0; i < prPos.length; i++) acopy[i] = prPos[i];
		}
		else
		{
			if (prEnd == null) return null;
			acopy = new double[prEnd.length];
			for (int i = 0; i < prEnd.length; i++) acopy[i] = prEnd[i];
		}
		return acopy;
	}
	
	/**
	 * Get the highest value in the POS probability curve, if there is a POS breakpoint
	 * probability curve generated by the LUMPY SV caller.
	 * @return Max point in the POS probability curve, if present. A value of -1.0 indicates
	 * that this value is unset.
	 */
	public double getPrPosPeak()
	{
		return highestPRPOS;
	}
	
	/**
	 * Get the highest value in the END probability curve, if there is a END breakpoint
	 * probability curve generated by the LUMPY SV caller.
	 * @return Max point in the END probability curve, if present. A value of -1.0 indicates
	 * that this value is unset.
	 */
	public double getPrEndPeak()
	{
		return highestPREND;
	}
	
	/**
	 * Get the higher value of the peaks of the two end breakpoint probability curves.
	 * If there are no probability curves, this method returns -1.0.
	 * @return Value of probability curve peak with the greater value.
	 */
	public double getHigherProbabilityPeak()
	{
		if (highestPRPOS > highestPREND) return highestPRPOS;
		else return highestPREND;
	}
	
	/**
	 * Get the lower value of the peaks of the two end breakpoint probability curves.
	 * If there are no probability curves, this method returns -1.0.
	 * @return Value of probability curve peak with the lower value.
	 */
	public double getLowerProbabilityPeak()
	{
		if (highestPRPOS < highestPREND) return highestPRPOS;
		else return highestPREND;
	}
	
	public int getSmallestAbsoluteLength()
	{
		int min = Integer.MAX_VALUE;
		int nAlt = this.SVlen.length;
		if (nAlt < 1) return 0;
		for (int i = 0; i < nAlt; i++)
		{
			int len = SVlen[i];
			int abslen = Math.abs(len);
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
		int nAlt = this.SVlen.length;
		if (nAlt < 1) return 0;
		for (int i = 0; i < nAlt; i++)
		{
			int len = SVlen[i];
			int abslen = Math.abs(len);
			if (abslen > max) max = abslen;
		}
		return max;
	}
	
	/* --- Setters --- */
	
	/**
	 * Set the END position, in one-based coordinates, for this variant.
	 * This method also updates SVLEN.
	 * @param pos One-based chromosome coordinate of the END position.
	 */
	public void setEndPosition(int pos)
	{
		end = pos;
		if (SVlen == null) SVlen = new int[1];
		SVlen[0] = end - this.getPosition();
		if (type == SVType.DEL) SVlen[0] *= -1;
	}

	/**
	 * Set the structural variant type for this variant (deletion, duplication, etc.) using
	 * an SVType enumeration.
	 * @param t Type to set
	 */
	public void setType(SVType t)
	{
		type = t;
	}
	
	/**
	 * Manually set the SV length relative to the primary alternate allele. In most
	 * cases, this will suffice as SVs should really be recorded one at a time.
	 * The value entered via this method will take precedence over any values automatically
	 * calculated from the POS and END points. Use with caution.
	 * @param len SVLEN (length) value to set for the primary alternate allele of this variant.
	 * This is also used as the default SVLEN value.
	 */
	public void setSVLength(int len)
	{
		if (SVlen == null){
			int aAll = super.countAltAlleles();
			if (aAll < 1) aAll = 1;
			SVlen = new int[aAll];
		}
		SVlen[0] = len;
	}
	
	/**
	 * Manually set the SV length relative to the specified alternate allele.
	 * The value entered via this method will take precedence over any values automatically
	 * calculated from the POS and END points. Use with caution.
	 * If the alt allele index is invalid, this method will return without changing anything.
	 * @param altAllele The alternate allele index (general allele index - 1) of the alternate
	 * allele to set the SVLEN relative to.
	 * @param len SVLEN (length) value to set for the selected alternate allele of this variant.
	 */
	public void setSVLength(int altAllele, int len)
	{
		if (altAllele < 0) return;
		int aAll = super.countAltAlleles();
		if (altAllele >= aAll) return;
		if (SVlen == null){
			if (aAll < 1) aAll = 1;
			SVlen = new int[aAll];
		}
		SVlen[altAllele] = len;
	}
	
	/**
	 * Set the variant imprecise flag. When this flag is set, the variant is marked as
	 * imprecise, and it is assumed that the breakends are not known at bp resoultion.
	 * @param b Boolean value to set for the flag.
	 */
	public void setImprecise(boolean b)
	{
		imprecise = b;
	}
	
	/**
	 * Set a confidence interval value by specifying which end of which confidence interval
	 * to set and the difference between it and the called position.
	 * @param offset Value to set - the difference between the CI end specified and the called
	 * position.
	 * @param variantEnd End of the SV to set CI value for. True = END, False = POS
	 * @param narrow Whether to set the tighter (95%) CI. True = CI95, False = CI90
	 * @param top Which end of the CI to set the position of. True = higher, False = lower
	 */
	public void setCIDiff(int offset, boolean variantEnd, boolean narrow, boolean top)
	{
		if (variantEnd)
		{
			if (narrow)
			{
				if (top) CI95EndHi = offset;
				else CI95EndLow = offset;
			}
			else
			{
				if (top) CIEndHi = offset;
				else CIEndLow = offset;
			}
		}
		else
		{
			if (narrow)
			{
				if (top) CI95PosHi = offset;
				else CI95PosLow = offset;
			}
			else
			{
				if (top) CIPosHi = offset;
				else CIPosLow = offset;
			}
		}
	}
	
	/**
	 * Add a "mate" variant to this variant. This field is usually for complex events where
	 * multiple BND records are involved, but in can in theory be used for other things as
	 * well to mark mutiple variants as being related or part of the same event.
	 * @param mateID ID string of variant record to add to mate list.
	 */
	public void addMate(String mateID)
	{
		if (mateid == null) mateid = new ArrayList<String>(4);
		mateid.add(mateID);
	}
	
	/**
	 * Search this variant's list of mate variants and if there is one by the given ID, remove
	 * it as a mate.
	 * @param mateID Variant ID of mate to remove mate list.
	 */
	public void removeMate(String mateID)
	{
		if (mateid == null) mateid = new ArrayList<String>(4);
		mateid.remove(mateID);
	}
	
	/**
	 * Clear all variant IDs from this variant's mate list.
	 */
	public void clearMates()
	{
		if (mateid == null) mateid = new ArrayList<String>(4);
		mateid.clear();
	}
	
	/**
	 * Set the ID of the SV event - that is, the single break or biological event that
	 * encompasses the SV group of this variant and its mates.
	 * @param ID Event ID. An empty or null string is interpreted as unset.
	 */
	public void setEventID(String ID)
	{
		event = ID;
	}
	
	/**
	 * Set the secondary variant flag. When this flag is set and the variant has at least
	 * one mate, it is assumed to be the "secondary" variant of a breakend pair or group
	 * for a single event.
	 * <br>This flag is useful, but more or less optional.
	 * @param b Secondary flag to set.
	 */
	public void setSecondary(boolean b)
	{
		secondary = b;
	}
	
	private void reevaluateProbPeak(boolean variantEnd)
	{
		double[] curve;
		if (variantEnd) curve = prEnd;
		else curve = prPos;
		
		double max = -1.0;
		
		for (int i = 0; i < curve.length; i++)
		{
			if (curve[i] > max) max = curve[i];
		}
		
		if (variantEnd) highestPREND = max;
		else highestPRPOS = max;
	}
	
	/**
	 * Manually replace a LUMPY probability curve. This method is primarily intended for use
	 * by the BreakendPair class.
	 * @param curve Array to set as the curve. This argument will be referenced, not copied.
	 * @param variantEnd End of the SV to set curve for. True = END, False = POS
	 */
	protected void setProbabilityCurve(double[] curve, boolean variantEnd)
	{
		//if (curve == null) return;
		if (variantEnd) prEnd = curve;
		else prPos = curve;
		if (curve != null) reevaluateProbPeak(variantEnd);
		else
		{
			if (variantEnd) this.highestPREND = 0.0;
			else this.highestPRPOS = 0.0;
		}
	}
	
	public void addAltAllele(String alt)
	{
		super.addAltAllele(alt);
		int aAll = super.countAltAlleles();
		if (SVlen.length < aAll)
		{
			int[] copy = new int[aAll + 2];
			for (int i = 0; i < SVlen.length; i++) copy[i] = SVlen[i];
			SVlen = copy;
		}
	}
	
	public void removeAltAllele(int index)
	{
		super.removeAltAllele(index);
		for (int i = index; i < SVlen.length - 1; i++)
		{
			SVlen[i] = SVlen[i+1];
		}
		SVlen[SVlen.length - 1] = 0;
	}
	
	public void clearAltAlleles()
	{
		super.clearAltAlleles();
		for (int i = 0; i < SVlen.length; i++) SVlen[i] = 0;
	}
	
	/* --- Comparing --- */
	
	public boolean equals(Object o)
	{
		if (o == this) return true;
		if (o == null) return false;
		if (!(o instanceof StructuralVariant)) return false;
		StructuralVariant v = (StructuralVariant)o;
		if (!this.getChromosome().equals(v.getChromosome())) return false;
		if (this.getPosition() != v.getPosition()) return false;
		if (this.getType() != v.getType()) return false;
		if (this.getSVLength() != v.getSVLength()) return false;
		return true;
	}
	
	private boolean pointInRange(int point, boolean variantEnd, boolean stringent_CI)
	{
		int rangeBot = getCIPosition(variantEnd, stringent_CI, false);
		int rangeTop = getCIPosition(variantEnd, stringent_CI, true);
		return (point >= rangeBot && point <= rangeTop);
	}
	
	private boolean rangesOverlap(int t1, int b1, int t2, int b2)
	{
		int bcount = 0;
		if (t1 <= t2 && t1 >= b2) bcount++;
		if (b1 <= t2 && b1 >= b2) bcount++;
		if (t2 <= t1 && t2 >= b1) bcount++;
		if (b2 <= t1 && b2 >= b1) bcount++;
		return (bcount == 2);
	}
	
	private boolean rangesOverlap(StructuralVariant v, boolean stringent_CI, boolean stringent_overlap, boolean end)
	{
		if(stringent_overlap)
		{
			if (!this.pointInRange(v.getPosition(), end, stringent_CI)) return false;
			if (!v.pointInRange(this.getPosition(), end, stringent_CI)) return false;
		}
		else
		{
			int t1 = this.getCIPosition(end, stringent_CI, true);
			int b1 = this.getCIPosition(end, stringent_CI, false);
			int t2 = v.getCIPosition(end, stringent_CI, true);
			int b2 = v.getCIPosition(end, stringent_CI, false);
			if (!rangesOverlap(t1, b1, t2, b2)) return false;
		}
		return true;
	}
	
	/**
	 * Compares another variant to this variant and determines whether they might be 
	 * considered to represent the same event, even if the data backing them are not
	 * perfectly identical, and thus would not qualify as "equal."
	 * @param v Variant to compare to this one.
	 * @param stringent_bothEnds Stringency setting: If on, both ends of the SV (not just one or the other) must be
	 * considered in range of each other.
	 * @param stringent_CI Stringency setting: If on, then the valid range around a breakpoint for considering
	 * two variants equivalent is the 95% confidence interval rather than the default 90% CI.
	 * @param stringent_overlap Stringency setting: If on, intervals at either end of the variant
	 * must encompass the actual endpoint calls themselves to be considered equivalent. Otherwise,
	 * a simple overlap of the two ranges is considered sufficient.
	 * @return True - If the two variants might be considered equivalent.
	 * <br> False - Otherwise.
	 */
	public boolean equivalent(Variant v, boolean stringent_bothEnds, boolean stringent_CI, boolean stringent_overlap)
	{
		//Both ends - both ends must be in range of each other
		// CI - if on, use the CI95
		// Overlap - if off, breakpoint ranges need only overlap. If on, breakpoints themselves must be in range of other.
		if (this.equals(v)) return true;
		if (!this.getChromosome().equals(v.getChromosome())) return false;
		if (!(v instanceof StructuralVariant)) return false;
		StructuralVariant sv = (StructuralVariant)v;
		if (this.getType() != sv.getType()) return false;
		if (stringent_bothEnds)
		{
			return (rangesOverlap(sv, stringent_CI, stringent_overlap, false) && rangesOverlap(sv, stringent_CI, stringent_overlap, true));
		}
		else
		{
			return (rangesOverlap(sv, stringent_CI, stringent_overlap, false) || rangesOverlap(sv, stringent_CI, stringent_overlap, true));
		}
	}
		
	public int compareTo(Variant o)
	{
		if (o == null) return 1;
		if (o == this) return 0;
		int chrComp = 0;
		try
		{
			chrComp = this.getChromosome().compareTo(o.getChromosome());
		}
		catch(NullPointerException e)
		{
			System.err.println("StructuralVariant.compareTo || Issue comparing contigs...");
			System.err.println("StructuralVariant.compareTo || This variant:");
			System.err.println(this.toVCFLine(null, null));
			System.err.println("StructuralVariant.compareTo || Other variant:");
			System.err.println(o.toVCFLine(null, null));
			throw new NullPointerException();
		}
		//int chrComp = this.getChromosome().compareTo(o.getChromosome());
		if (chrComp != 0) return chrComp;
		if (this.getPosition() != o.getPosition()) return this.getPosition() - o.getPosition();
		if (!(o instanceof StructuralVariant)) return 1;
		return 0;
	}
	
	public boolean inRegion(Contig c, int start, int end, boolean anyend)
	{
		if (!this.getChromosome().equals(c)) return false;
		boolean s = this.getPosition() >= start && this.getPosition() <= end;
		boolean e = this.getEndPosition() >= start && this.getEndPosition() <= end;
		if (anyend) return s | e;
		return s && e;
	}
	
	/* --- Definitions --- */
	
	private static void populateDefMaps()
	{
		altDefMap = new HashMap<String, String>();
		altDefMap.put(INFODEF_ALT_DEL, "Deletion");
		altDefMap.put(INFODEF_ALT_DELME, "Deletion of mobile element");
		altDefMap.put(INFODEF_ALT_DUP, "Duplication");
		altDefMap.put(INFODEF_ALT_INS, "Insertion of novel sequence");
		altDefMap.put(INFODEF_ALT_INSME, "Insertion of mobile element");
		altDefMap.put(INFODEF_ALT_INV, "Inversion");
		altDefMap.put(INFODEF_ALT_TANDEM, "Tandem duplication");
		altDefMap.put(INFODEF_ALT_CNV, "Copy number variable region");
		
		infoDefMap = new HashMap<String, InfoDefinition>();
		infoDefMap.put(INFODEF_INFO_SVTYPE.getKey(), INFODEF_INFO_SVTYPE);
		infoDefMap.put(INFODEF_INFO_SVLEN.getKey(), INFODEF_INFO_SVLEN);
		infoDefMap.put(INFODEF_INFO_END.getKey(), INFODEF_INFO_END);
		infoDefMap.put(INFODEF_INFO_IMPRECISE.getKey(), INFODEF_INFO_IMPRECISE);
		infoDefMap.put(INFODEF_INFO_NOVEL.getKey(), INFODEF_INFO_NOVEL);
		infoDefMap.put(INFODEF_INFO_MATEID.getKey(), INFODEF_INFO_MATEID);
		infoDefMap.put(INFODEF_INFO_EVENT.getKey(), INFODEF_INFO_EVENT);
		infoDefMap.put(INFODEF_INFO_CIPOS.getKey(), INFODEF_INFO_CIPOS);
		infoDefMap.put(INFODEF_INFO_CIEND.getKey(), INFODEF_INFO_CIEND);
		infoDefMap.put(INFODEF_INFO_CIPOS95.getKey(), INFODEF_INFO_CIPOS95);
		infoDefMap.put(INFODEF_INFO_CIEND95.getKey(), INFODEF_INFO_CIEND95);
		
		infoDefMap.put(INFODEF_INFO_SECONDARY.getKey(), INFODEF_INFO_SECONDARY);
		infoDefMap.put(INFODEF_INFO_STRANDS.getKey(), INFODEF_INFO_STRANDS);
		infoDefMap.put(INFODEF_INFO_SU.getKey(), INFODEF_INFO_SU);
		infoDefMap.put(INFODEF_INFO_PE.getKey(), INFODEF_INFO_PE);
		infoDefMap.put(INFODEF_INFO_SR.getKey(), INFODEF_INFO_SR);
		infoDefMap.put(INFODEF_INFO_BD.getKey(), INFODEF_INFO_BD);
		infoDefMap.put(INFODEF_INFO_EV.getKey(), INFODEF_INFO_EV);
		infoDefMap.put(INFODEF_INFO_PRPOS.getKey(), INFODEF_INFO_PRPOS);
		infoDefMap.put(INFODEF_INFO_PREND.getKey(), INFODEF_INFO_PREND);
	}

	/**
	 * Get the InfoDefinition for standard SV and LUMPY INFO fields, knowing the key.
	 * @param key Key to field of interest.
	 * @return InfoDefinition if the key is to a known SV field. null if the key is unknown
	 * or not the key to one of the statically defined SV fields.
	 */
	public static InfoDefinition getStandardDefinition(String key)
	{
		if (infoDefMap == null) populateDefMaps();
		return infoDefMap.get(key);
	}
	
	/**
	 * Get the Alt allele description for a VCF header line describing a standard SV alt allele.
	 * @param key Alt allele key
	 * @return String containing the description of the alt allele field, as is required in
	 * the VCF header and by extension, likely many VCF parsers. This function returns null
	 * if the key is unknown.
	 */
	public static String getAltDefDescription(String key)
	{
		if (altDefMap == null) populateDefMaps();
		return altDefMap.get(key);
	}
	
	/**
	 * Add standard INFO and ALT definitions used commonly in VCF SV representation
	 * to a variant pool.
	 * @param pool Pool to add definitions to.
	 * @param CI95 Whether to include CI95 definitions.
	 */
	public static void addStandardDefs(VariantPool pool, boolean CI95)
	{
		if (pool == null) return;
		if (altDefMap == null) populateDefMaps();
		
		InfoDefinition def = INFODEF_INFO_SVTYPE;
		pool.addInfoFieldDefinition(def.getKey(), def);
		pool.addInfoKeyToActiveList(def.getKey());
		
		def = INFODEF_INFO_SVLEN;
		pool.addInfoFieldDefinition(def.getKey(), def);
		pool.addInfoKeyToActiveList(def.getKey());
		
		def = INFODEF_INFO_END;
		pool.addInfoFieldDefinition(def.getKey(), def);
		pool.addInfoKeyToActiveList(def.getKey());
		
		def = INFODEF_INFO_IMPRECISE;
		pool.addInfoFieldDefinition(def.getKey(), def);
		pool.addInfoKeyToActiveList(def.getKey());
		
		def = INFODEF_INFO_CIPOS;
		pool.addInfoFieldDefinition(def.getKey(), def);
		pool.addInfoKeyToActiveList(def.getKey());
		
		def = INFODEF_INFO_CIEND;
		pool.addInfoFieldDefinition(def.getKey(), def);
		pool.addInfoKeyToActiveList(def.getKey());
		
		if (CI95)
		{
			def = INFODEF_INFO_CIPOS95;
			pool.addInfoFieldDefinition(def.getKey(), def);
			pool.addInfoKeyToActiveList(def.getKey());
			
			def = INFODEF_INFO_CIEND95;
			pool.addInfoFieldDefinition(def.getKey(), def);
			pool.addInfoKeyToActiveList(def.getKey());
		}
		
		def = INFODEF_INFO_MATEID;
		pool.addInfoFieldDefinition(def.getKey(), def);
		pool.addInfoKeyToActiveList(def.getKey());
		def = INFODEF_INFO_EVENT;
		pool.addInfoFieldDefinition(def.getKey(), def);
		pool.addInfoKeyToActiveList(def.getKey());
		def = INFODEF_INFO_SECONDARY;
		pool.addInfoFieldDefinition(def.getKey(), def);
		pool.addInfoKeyToActiveList(def.getKey());
		
		String alt = INFODEF_ALT_DEL;
		pool.addCustomAlt(alt, altDefMap.get(alt));
		alt = INFODEF_ALT_DUP;
		pool.addCustomAlt(alt, altDefMap.get(alt));
		alt = INFODEF_ALT_INS;
		pool.addCustomAlt(alt, altDefMap.get(alt));
		alt = INFODEF_ALT_INV;
		pool.addCustomAlt(alt, altDefMap.get(alt));
		alt = INFODEF_ALT_CNV;
		pool.addCustomAlt(alt, altDefMap.get(alt));
		alt = INFODEF_ALT_TANDEM;
		pool.addCustomAlt(alt, altDefMap.get(alt));
		alt = INFODEF_ALT_DELME;
		pool.addCustomAlt(alt, altDefMap.get(alt));
		alt = INFODEF_ALT_INSME;
		pool.addCustomAlt(alt, altDefMap.get(alt));
		
		
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
		
		fieldMap.put(INFODEF_INFO_SVTYPE.getKey(), new Field(){
			public String get()
			{
				if (type == null) return null;
				return type.toString();
			}
			public String[] getAll()
			{
				if (type == null) return null;
				String[] all = new String[1];
				all[0] = get();
				return all;
			}
			public String getFirst()
			{
				return get();
			}
			public void set(String value)
			{
				SVType t = SVType.getType(value);
				if (t == null) t = SVType.OTHER;
				type = t;
			}
			public boolean hasFlag()
			{
				return false;
			}
		});
		
		fieldMap.put(INFODEF_INFO_SVLEN.getKey(), new Field(){
			public String get()
			{
				String s = "";
				int aAll = countAltAlleles();
				if (aAll < 1) return null;
				for (int i = 0; i < aAll; i++)
				{
					s += SVlen[i];
					if (i < aAll - 1) s += ",";
				}
				return s;
			}
			public String[] getAll()
			{
				int aAll = countAltAlleles();
				if (aAll < 1) return null;
				String[] all = new String[aAll];
				for (int i = 0; i < aAll; i++)
				{
					all [i]= Integer.toString(SVlen[i]);
				}
				return all;
			}
			public String getFirst()
			{
				if (SVlen == null) return null;
				return Integer.toString(SVlen[0]);
			}			
			public void set(String value)
			{
				String[] fields = value.split(",");
				if (fields == null) return;
				int aAll = countAltAlleles();
				int sz = aAll;
				if (fields.length > aAll) sz = fields.length;
				int[] arr = new int[sz];
				for (int i = 0; i < fields.length; i++)
				{
					try{arr[i] = Integer.parseInt(fields[i]);}
					catch(NumberFormatException e) {return;}
				}
				SVlen = arr;
			}
			public boolean hasFlag()
			{
				return false;
			}
		});

		fieldMap.put(INFODEF_INFO_END.getKey(), new Field(){
			public String get()
			{
				if (end < 0) return null;
				return Integer.toString(end);
			}
			public String[] getAll()
			{
				if (end < 0) return null;
				String[] all = new String[1];
				all[0] = get();
				return all;
			}
			public String getFirst()
			{
				return get();
			}			
			public void set(String value)
			{
				try{int i = Integer.parseInt(value); end = i;}
				catch(NumberFormatException e){return;}
			}
			public boolean hasFlag()
			{
				return false;
			}
		});
		
		fieldMap.put(INFODEF_INFO_IMPRECISE.getKey(), new Field(){
			public String get()
			{
				if (!imprecise) return null;
				return "";
			}
			public String[] getAll()
			{
				if (!imprecise) return null;
				String[] all = new String[1];
				all[0] = get();
				return all;
			}
			public String getFirst()
			{
				return get();
			}			
			public void set(String value)
			{
				imprecise = (value != null);
			}
			public boolean hasFlag()
			{
				return imprecise;
			}
		});
		
		fieldMap.put(INFODEF_INFO_CIPOS.getKey(), new Field(){
			public String get()
			{
				int bot = getCIDiff(false, false, false);
				int top = getCIDiff(false, false, true);
				return bot + "," + top;
			}
			public String[] getAll()
			{
				String[] all = new String[2];
				all[0] = Integer.toString(getCIDiff(false, false, false));
				all[1] = Integer.toString(getCIDiff(false, false, true));
				return all;
			}
			public String getFirst()
			{
				return Integer.toString(getCIDiff(false, false, false));
			}				
			public void set(String value)
			{
				String[] fields = value.split(",");
				if (fields == null) return;
				if (fields.length != 2) return;
				try{
					int i = Integer.parseInt(fields[0]);
					int j = Integer.parseInt(fields[1]);
					setCIDiff(i, false, false, false);
					setCIDiff(j, false, false, true);
				}
				catch(NumberFormatException e) {return;}
			}
			public boolean hasFlag()
			{
				return false;
			}
		});
		
		fieldMap.put(INFODEF_INFO_CIEND.getKey(), new Field(){
			public String get()
			{
				if (end < 0) return null;
				int bot = getCIDiff(true, false, false);
				int top = getCIDiff(true, false, true);
				return bot + "," + top;
			}
			public String[] getAll()
			{
				if (end < 0) return null;
				String[] all = new String[2];
				all[0] = Integer.toString(getCIDiff(true, false, false));
				all[1] = Integer.toString(getCIDiff(true, false, true));
				return all;
			}
			public String getFirst()
			{
				if (end < 0) return null;
				return Integer.toString(getCIDiff(true, false, false));
			}					
			public void set(String value)
			{
				String[] fields = value.split(",");
				if (fields == null) return;
				if (fields.length != 2) return;
				try{
					int i = Integer.parseInt(fields[0]);
					int j = Integer.parseInt(fields[1]);
					setCIDiff(i, true, false, false);
					setCIDiff(j, true, false, true);
				}
				catch(NumberFormatException e) {return;}
			}
			public boolean hasFlag()
			{
				return false;
			}
		});
		
		fieldMap.put(INFODEF_INFO_CIPOS95.getKey(), new Field(){
			public String get()
			{
				int bot = getCIDiff(false, true, false);
				int top = getCIDiff(false, true, true);
				return bot + "," + top;
			}
			public void set(String value)
			{
				String[] fields = value.split(",");
				if (fields == null) return;
				if (fields.length != 2) return;
				try{
					int i = Integer.parseInt(fields[0]);
					int j = Integer.parseInt(fields[1]);
					setCIDiff(i, false, true, false);
					setCIDiff(j, false, true, true);
				}
				catch(NumberFormatException e) {return;}
			}
			public String[] getAll()
			{
				String[] all = new String[2];
				all[0] = Integer.toString(getCIDiff(false, true, false));
				all[1] = Integer.toString(getCIDiff(false, true, true));
				return all;
			}
			public String getFirst()
			{
				return Integer.toString(getCIDiff(false, true, false));
			}		
			public boolean hasFlag()
			{
				return false;
			}
		});
		
		fieldMap.put(INFODEF_INFO_CIEND95.getKey(), new Field(){
			public String get()
			{
				if (end < 0) return null;
				int bot = getCIDiff(true, true, false);
				int top = getCIDiff(true, true, true);
				return bot + "," + top;
			}
			public String[] getAll()
			{
				if (end < 0) return null;
				String[] all = new String[2];
				all[0] = Integer.toString(getCIDiff(true, true, false));
				all[1] = Integer.toString(getCIDiff(true, true, true));
				return all;
			}
			public String getFirst()
			{
				if (end < 0) return null;
				return Integer.toString(getCIDiff(true, true, false));
			}					
			public void set(String value)
			{
				String[] fields = value.split(",");
				if (fields == null) return;
				if (fields.length != 2) return;
				try{
					int i = Integer.parseInt(fields[0]);
					int j = Integer.parseInt(fields[1]);
					setCIDiff(i, true, true, false);
					setCIDiff(j, true, true, true);
				}
				catch(NumberFormatException e) {return;}
			}
			public boolean hasFlag()
			{
				return false;
			}
		});

		fieldMap.put(INFODEF_INFO_MATEID.getKey(), new Field(){
			public String get()
			{
				if (mateid == null) return null;
				if (mateid.isEmpty()) return null;
				String s = "";
				int sz = mateid.size();
				for (int i = 0; i < sz; i++)
				{
					s += mateid.get(i);
					if (i < sz - 1) s += ",";
				}
				return s;
			}
			public String[] getAll()
			{
				if (mateid == null) return null;
				if (mateid.isEmpty()) return null;
				String[] all = new String[mateid.size()];
				return mateid.toArray(all);
			}
			public String getFirst()
			{
				if (mateid == null) return null;
				if (mateid.isEmpty()) return null;
				return mateid.get(0);
			}		
			public void set(String value)
			{
				String[] fields = value.split(",");
				if (fields == null) return;
				if (mateid == null) mateid = new ArrayList<String>(fields.length + 2);
				for (int i = 0; i < fields.length; i++)
				{
					mateid.add(fields[i]);
				}
			}
			public boolean hasFlag()
			{
				return false;
			}
		});

		fieldMap.put(INFODEF_INFO_EVENT.getKey(), new Field(){
			public String get()
			{
				if (event == null) return null;
				if (event.isEmpty()) return null;
				return event;
			}
			public String[] getAll()
			{
				if (event == null) return null;
				if (event.isEmpty()) return null;
				String[] all = new String[1];
				all[0] = get();
				return all;
			}
			public String getFirst()
			{
				return get();
			}
			public void set(String value)
			{
				event = value;
			}
			public boolean hasFlag()
			{
				return false;
			}
		});
		
		fieldMap.put(INFODEF_INFO_SECONDARY.getKey(), new Field(){
			public String get()
			{
				if (!secondary) return null;
				return "";
			}
			public String[] getAll()
			{
				if (!secondary) return null;
				String[] all = new String[1];
				all[0] = get();
				return all;
			}
			public String getFirst()
			{
				return get();
			}			
			public void set(String value)
			{
				secondary = (value != null);
			}
			public boolean hasFlag()
			{
				return secondary;
			}
		});

		fieldMap.put(INFODEF_INFO_SU.getKey(), new Field(){
			public String get()
			{
				if (evidence_SU == null) return null;
				String s = "";
				int len = evidence_SU.length;
				for (int i = 0; i < len; i++)
				{
					s += evidence_SU[i];
					if (i < len - 1) s += ",";
				}
				return s;
			}
			public String[] getAll()
			{
				if (evidence_SU == null) return null;
				int len = evidence_SU.length;
				String[] all = new String[len];
				for (int i = 0; i < len; i++)
				{
					all [i]= Integer.toString(evidence_SU[i]);
				}
				return all;
			}
			public String getFirst()
			{
				if (evidence_SU == null) return null;
				return Integer.toString(evidence_SU[0]);
			}			
			public void set(String value)
			{
				String[] fields = value.split(",");
				if (fields == null) return;
				int sz = fields.length;
				if (sz < 1) return;
				int[] arr = new int[sz];
				for (int i = 0; i < fields.length; i++)
				{
					try{arr[i] = Integer.parseInt(fields[i]);}
					catch(NumberFormatException e) {return;}
				}
				evidence_SU = arr;
			}
			public boolean hasFlag()
			{
				return false;
			}
		});

		fieldMap.put(INFODEF_INFO_PE.getKey(), new Field(){
			public String get()
			{
				if (evidence_PE == null) return null;
				String s = "";
				int len = evidence_PE.length;
				for (int i = 0; i < len; i++)
				{
					s += evidence_PE[i];
					if (i < len - 1) s += ",";
				}
				return s;
			}
			public String[] getAll()
			{
				if (evidence_PE == null) return null;
				int len = evidence_PE.length;
				String[] all = new String[len];
				for (int i = 0; i < len; i++)
				{
					all [i]= Integer.toString(evidence_PE[i]);
				}
				return all;
			}
			public String getFirst()
			{
				if (evidence_PE == null) return null;
				return Integer.toString(evidence_PE[0]);
			}			
			public void set(String value)
			{
				String[] fields = value.split(",");
				if (fields == null) return;
				int sz = fields.length;
				if (sz < 1) return;
				int[] arr = new int[sz];
				for (int i = 0; i < fields.length; i++)
				{
					try{arr[i] = Integer.parseInt(fields[i]);}
					catch(NumberFormatException e) {return;}
				}
				evidence_PE = arr;
			}
			public boolean hasFlag()
			{
				return false;
			}
		});
		
		fieldMap.put(INFODEF_INFO_SR.getKey(), new Field(){
			public String get()
			{
				if (evidence_SR == null) return null;
				String s = "";
				int len = evidence_SR.length;
				for (int i = 0; i < len; i++)
				{
					s += evidence_SR[i];
					if (i < len - 1) s += ",";
				}
				return s;
			}
			public String[] getAll()
			{
				if (evidence_SR == null) return null;
				int len = evidence_SR.length;
				String[] all = new String[len];
				for (int i = 0; i < len; i++)
				{
					all [i]= Integer.toString(evidence_SR[i]);
				}
				return all;
			}
			public String getFirst()
			{
				if (evidence_SR == null) return null;
				return Integer.toString(evidence_SR[0]);
			}			
			public void set(String value)
			{
				String[] fields = value.split(",");
				if (fields == null) return;
				int sz = fields.length;
				if (sz < 1) return;
				int[] arr = new int[sz];
				for (int i = 0; i < fields.length; i++)
				{
					try{arr[i] = Integer.parseInt(fields[i]);}
					catch(NumberFormatException e) {return;}
				}
				evidence_SR = arr;
			}
			public boolean hasFlag()
			{
				return false;
			}
		});
		
		fieldMap.put(INFODEF_INFO_BD.getKey(), new Field(){
			public String get()
			{
				if (evidence_BD == null) return null;
				String s = "";
				int len = evidence_BD.length;
				for (int i = 0; i < len; i++)
				{
					s += evidence_BD[i];
					if (i < len - 1) s += ",";
				}
				return s;
			}
			public String[] getAll()
			{
				if (evidence_BD == null) return null;
				int len = evidence_BD.length;
				String[] all = new String[len];
				for (int i = 0; i < len; i++)
				{
					all [i]= Integer.toString(evidence_BD[i]);
				}
				return all;
			}
			public String getFirst()
			{
				if (evidence_BD == null) return null;
				return Integer.toString(evidence_BD[0]);
			}			
			public void set(String value)
			{
				String[] fields = value.split(",");
				if (fields == null) return;
				int sz = fields.length;
				if (sz < 1) return;
				int[] arr = new int[sz];
				for (int i = 0; i < fields.length; i++)
				{
					try{arr[i] = Integer.parseInt(fields[i]);}
					catch(NumberFormatException e) {return;}
				}
				evidence_BD = arr;
			}
			public boolean hasFlag()
			{
				return false;
			}
		});

		fieldMap.put(INFODEF_INFO_PRPOS.getKey(), new Field(){
			public String get()
			{
				if (prPos == null) return null;
				String s = "";
				int len = prPos.length;
				for (int i = 0; i < len; i++)
				{
					s += prPos[i];
					if (i < len - 1) s += ",";
				}
				return s;
			}
			public String[] getAll()
			{
				if (prPos == null) return null;
				int len = prPos.length;
				String[] all = new String[len];
				for (int i = 0; i < len; i++)
				{
					all [i]= Double.toString(prPos[i]);
				}
				return all;
			}
			public String getFirst()
			{
				if (prPos == null) return null;
				return Double.toString(prPos[0]);
			}			
			public void set(String value)
			{
				String[] fields = value.split(",");
				if (fields == null) return;
				int sz = fields.length;
				if (sz < 1) return;
				double[] arr = new double[sz];
				for (int i = 0; i < fields.length; i++)
				{
					try{arr[i] = Double.parseDouble(fields[i]);}
					catch(NumberFormatException e) {return;}
				}
				prPos = arr;
			}
			public boolean hasFlag()
			{
				return false;
			}
		});
		
		fieldMap.put(INFODEF_INFO_PREND.getKey(), new Field(){
			public String get()
			{
				if (prEnd == null) return null;
				String s = "";
				int len = prEnd.length;
				for (int i = 0; i < len; i++)
				{
					s += prEnd[i];
					if (i < len - 1) s += ",";
				}
				return s;
			}
			public String[] getAll()
			{
				if (prEnd == null) return null;
				int len = prEnd.length;
				String[] all = new String[len];
				for (int i = 0; i < len; i++)
				{
					all [i]= Double.toString(prEnd[i]);
				}
				return all;
			}
			public String getFirst()
			{
				if (prEnd == null) return null;
				return Double.toString(prEnd[0]);
			}			
			public void set(String value)
			{
				String[] fields = value.split(",");
				if (fields == null) return;
				int sz = fields.length;
				if (sz < 1) return;
				double[] arr = new double[sz];
				for (int i = 0; i < fields.length; i++)
				{
					try{arr[i] = Double.parseDouble(fields[i]);}
					catch(NumberFormatException e) {return;}
				}
				prEnd = arr;
			}
			public boolean hasFlag()
			{
				return false;
			}
		});

	}

	public String[] getInfoEntry(String key)
	{
		//First check standard, then check SV specific
		String[] vals = super.getInfoEntry(key);
		if (vals != null) return vals;
		Field f = fieldMap.get(key);
		if (f == null) return null;
		return f.getAll();
	}
	
	public String getSingleInfoEntry(String key)
	{
		String val = super.getSingleInfoEntry(key);
		if (val != null) return val;
		Field f = fieldMap.get(key);
		if (f == null) return null;
		return f.getFirst();
	}
	
	public boolean getInfoFlag(String fieldKey)
	{
		boolean flag = super.getInfoFlag(fieldKey);
		if (flag) return flag;
		Field f = fieldMap.get(fieldKey);
		if (f == null) return false;
		return f.hasFlag();
	}
	
	/* --- Serialization --- */
	
	/**
	 * Get a representation of the structural variant as a BED formatted record.
	 * <br>IMPORTANT: BED records are minimalistic and this serializer only includes
	 * the chromosome, start position, end position, and variant ID string as currently set
	 * in this instance. Information may be lost if you don't know what you are doing!
	 * @return BED formatted variant record representing this structural variant.
	 */
	public String toBEDLine()
	{
		//chrom	start	end	name	score	strand	(other things)
		return this.getChromosome() + "\t" 
			 + this.getPosition() + "\t" 
			 + this.getEndPosition() + "\t" 
			 + this.getVarID();
	}
	
	/**
	 * Get a representation of the structural variant as a BEDPE formatted record.
	 * <br>IMPORTANT: BED and BEDPE records are minimalistic. The BEDPE format is designed
	 * to represent two regions at a time - in this case, it is used to mark possible ranges of
	 * the breakends of imprecise variants.
	 * Information may be lost if you don't know what you are doing!
	 * @return BEDPE formatted variant record representing this structural variant.
	 */
	public String toBEDPELine()
	{
		//chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2	(custom)
		String s = "";
		//chrom1
		s += getChromosome() + "\t";
		//start1
		s += this.getCIPosition(false, false, false) + "\t";
		//end1
		s += this.getCIPosition(false, false, true) + "\t";
		//chrom2
		s += getChromosome() + "\t";
		//start2
		s += this.getCIPosition(true, false, false) + "\t";
		//end2
		s += this.getCIPosition(true, false, true) + "\t";
		//name
		s += this.getVarID();
		return s;
	}
	
	public String toViewerBEDLine(String sample)
	{
		Genotype g = this.getSampleGenotype(sample);
		if (g == null) return null;
		double altpercent = g.getPercentAlt(); //TODO: Make sure we don't throw out non-genotyped!
		if (altpercent <= 0.0) {
			if (!g.isGenotypeUnknown()) return null;
			else altpercent = 100.0;
		}
		
		String s = "";
		//s += "chr" + getChromosome() + "\t"; //chrom
		String cname = this.getChromosome().getUCSCName();
		if (cname == null) return null;
		s += cname + "\t"; //chrom
		int st = getCIPosition(false, true, false) - 1;
		int ed = getCIPosition(true, true, true);
		int clen = (int)this.getChromosome().getLength();
		if (st < 0) st = 0;
		if (st > clen) return null; //Then the genome build read with this variant is NOT the one it is mapped to!
		if (ed > clen) ed = clen; //Probably incorrectly mapped, but in case it's a CI thing, we'll let it slide.
		s += st + "\t"; //start
		s += ed + "\t"; //end
		s += getVarID() + "\t"; //name

		int score = (int)Math.round(altpercent * 10.0);
		s += score + "\t"; //score
		
		s += ".\t"; //strand
		s += (getPosition() - 1) + "\t"; //thickstart
		s += getEndPosition() + "\t"; //thickend
		
		if (getType() == SVType.CNV)
		{
			//Attempt to find a genotype for our sample.
			int cn = g.getCopyNumber();
			if (cn < 0) s += UCSCGVBED.COLOR_CNVUNK + "\t";
			else
			{
				if (cn < 2) s += UCSCGVBED.COLOR_DEL + "\t";
				else if (cn > 2) s += UCSCGVBED.COLOR_DUP + "\t";
				else s += UCSCGVBED.COLOR_SNV + "\t";
			}
		}
		else s += UCSCGVBED.getSVColor(getType()) + "\t"; //color
		
		//TODO: If a block size is zero, then don't make that block!
		int sz1 = getCIPosition(false, true, true) - getCIPosition(false, true, false);
		int sz2 = ed - getCIPosition(true, true, false);
		int st2 = getCIPosition(true, true, false);
		if (sz1 <= 0)
		{
			sz1 = 1;
		}
		if (sz2 <= 0)
		{
			st2--;
			sz2 = 1;
		}
		
		s += "2\t"; //blockcount
		s += sz1 + "," + sz2 + "\t"; //blocksizes
		s += "0," + (st2 - st); //block starts
		
		
		return s;
	}
	
	/* --- Split & Merge --- */
	
	/**
	 * Split this structural variant into a Breakend Pair (a pair of BND
	 * structural variants). This may be needed for some VCF processing for
	 * VCF parsers that do not read the END coordinate.
	 * @return BreakendPair object containing two BND variants representing each
	 * end of this variant.
	 */
	public BreakendPair split()
	{
		if (this instanceof BreakendPair) return (BreakendPair)this;
		StructuralVariant v1 = new StructuralVariant();
		StructuralVariant v2 = new StructuralVariant();
		
		//Positions & Ends
		v1.setChromosome(getChromosome());
		v2.setChromosome(getChromosome());
		v1.setPosition(this.getPosition());
		v2.setPosition(this.getEndPosition());
		
		//Other variant instance properties
		v1.setVariantName(this.getVarID() + "_1");
		v2.setVariantName(this.getVarID() + "_2");
		v1.setRefAllele(this.getRefAllele());
		v2.setRefAllele(this.getRefAllele());
		v1.addAltAllele(this.getAltAllele(0));
		v2.addAltAllele(this.getAltAllele(0));
		v1.setFilterPass(this.passedAllFilters());
		v2.setFilterPass(this.passedAllFilters());
		v1.setQuality(this.getQuality());
		v2.setQuality(this.getQuality());
		String[] filters = this.getFiltersFailed();
		if (filters != null)
		{
			for (String s : filters)
			{
				v1.addFailedFilter(s);
				v2.addFailedFilter(s);
			}
		}
		
		//SVType
		v1.setType(SVType.BND);
		v2.setType(SVType.BND);
		v1.addInfoField(StructuralVariant.INFODEF_INFO_TRUESVTYPE.getKey(), this.getType().toString());
		v2.addInfoField(StructuralVariant.INFODEF_INFO_TRUESVTYPE.getKey(), this.getType().toString());
		
		//SVLen
		v1.setSVLength(this.getSVLength());
		v2.setSVLength(this.getSVLength());
		
		//Imprecise
		v1.setImprecise(this.isImprecise());
		v2.setImprecise(this.isImprecise());
		
		//CIs
		v1.setCIDiff(this.getCIDiff(false, false, false), false, false, false);
		v1.setCIDiff(this.getCIDiff(false, false, true), false, false, true);
		v1.setCIDiff(this.getCIDiff(false, true, false), false, true, false);
		v1.setCIDiff(this.getCIDiff(false, true, true), false, true, true);
		v2.setCIDiff(this.getCIDiff(true, false, false), false, false, false);
		v2.setCIDiff(this.getCIDiff(true, false, true), false, false, true);
		v2.setCIDiff(this.getCIDiff(true, true, false), false, true, false);
		v2.setCIDiff(this.getCIDiff(true, true, true), false, true, true);
		
		//MateID and Event ID
		v1.setEventID(this.getVarID());
		v2.setEventID(this.getVarID());
		v1.addMate(this.getVarID() + "_2");
		v2.addMate(this.getVarID() + "_1");
		
		//Other sv instance properties
		v1.setSecondary(false);
		v2.setSecondary(true);
		v1.evidence_SU = this.evidence_SU;
		v2.evidence_SU = this.evidence_SU;
		v1.evidence_SR = this.evidence_SR;
		v2.evidence_SR = this.evidence_SR;
		v1.evidence_PE = this.evidence_PE;
		v2.evidence_PE = this.evidence_PE;
		v1.evidence_BD = this.evidence_BD;
		v2.evidence_BD = this.evidence_BD;
		v1.prPos = this.prPos;
		v1.highestPRPOS = this.highestPRPOS;
		v2.prPos = this.prEnd;
		v2.highestPRPOS = this.highestPREND;
		//v1.function = this.function;
		//v2.function = this.function;
		
		//Info fields
		Set<String> keys = getAllInfoKeys();
		if (keys != null)
		{
			for(String k : keys)
			{
				v1.addInfoField(k, this.getInfoEntry(k));
				v2.addInfoField(k, this.getInfoEntry(k));
			}	
		}
		
		//Genotypes
		String[] gfields = getOrderedGenotypeFields();
		if (gfields != null)
		{
			for (String g : gfields)
			{
				v1.addGenotypeField(g);
				v2.addGenotypeField(g);
			}	
		}
		Set<String> samps = this.getAllGenotypedSamples();
		if (samps != null)
		{
			for (String s : samps)
			{
				v1.addGenotype(s, this.getSampleGenotype(s));
				v2.addGenotype(s, this.getSampleGenotype(s));
			}	
		}
		
		BreakendPair pair = new BreakendPair(v1, v2);
		
		return pair;
	}
	
	/* --- Other --- */
	
	

}

	