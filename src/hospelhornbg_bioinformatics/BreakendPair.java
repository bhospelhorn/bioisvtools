package hospelhornbg_bioinformatics;

import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import hospelhornbg_bioinformatics.VariantPool.InfoDefinition;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.GeneFunc;

/*
 * UPDATES 
 * 
 * Original Version: 1.0.0 (December 8, 2017)
 * 
 * 1.0.0 -> 1.1.0 | January 12, 2018
 * 	Fixed argument for VCF serializer override
 * 	Added BED and BEDPE serializers
 * 	Fixed construction method to match changes made to Variant and Structural Variant classes
 * 	Overrode LUMPY evidence count getters to return total counts
 * 	Javadoc
 * 
 * 1.1.0 -> 1.1.1 | January 24, 2018
 * 	Added genome browser BED serializer
 * 
 * 1.1.1 -> 1.1.2 | February 7, 2018
 * 	Added ability to directly retrieve internal variants. Use with care.
 * 
 * 1.1.2 -> 1.2.0 | February 20, 2018
 * 	GenomeBuild/ Contig compatibility update
 * 
 * 1.2.0 -> 1.2.1 | March 1, 2018
 * 	Get length of chrom from contig object, not from UCSC static map.
 * 
 * 1.2.1 -> 1.2.2 | March 26, 2018
 * 	Added method to get genotype from individual BNDs. Overrode Variant getSampleGenotype to return var1 genotype.
 *
 * 1.2.2 -> 1.2.3 | May 2, 2018
 * 	Updates super gene func with higher priority gene func from BND pair
 * 
 * 1.2.3 -> 1.2.4 | July 20, 2018
 * 	Added inRegion method
 * 
 * /

/*
 * Possible future improvements:
 * 
 */

/**
 * Structural Variant subtype for representing a pair of breakends that aren't
 * contiguous, such as the breakends in a translocation.
 * Usually, this class is used to hold two BND type SVs.
 * @author Blythe Hospelhorn
 * @version 1.2.4
 * @since July 20, 2018
 *
 */
public class BreakendPair extends StructuralVariant{
	
	public static final int BND_UNKNOWN = 0;
	public static final int BND_INTRACHROMOSOMAL = 1;
	public static final int BND_INTERCHROMOSOMAL = 2;
	public static final int BND_CHROMOCONTIG = 3;
	public static final int BND_DOUBLECONTIG = 4;
	
	public static final int ORIENTATION_UNKNOWN = 0;
	public static final int ORIENTATION_RIGHT = 1;
	public static final int ORIENTATION_LEFT = 2;
	public static final int ORIENTATION_RIGHTREV = 3;
	public static final int ORIENTATION_LEFTREV = 4;

	private StructuralVariant var1;
	private StructuralVariant var2;
	
	private int bndType;
	
	/**
	 * Construct a breakend pair type structural variant from two structural variants.
	 * @param v1 One structural variant in the pair
	 * @param v2 The other structural variant in the pair
	 */
	public BreakendPair(StructuralVariant v1, StructuralVariant v2)
	{
		if (v1 == null || v2 == null) throw new IllegalArgumentException();
		int vcomp = v1.compareTo(v2);
		if (vcomp  > 0)
		{
			var1 = v2;
			var2 = v1;
		}
		else
		{
			var1 = v1;
			var2 = v2;
		}
		constructorCore();
	}
	
	/**
	 * Construct a breakend pair type structural variant from a known structural
	 * variant, and find its mate in a provided collection of structural variants.
	 * @param v Initial structural variant to build breakend pair from.
	 * @param vPool Pool of structural variants to search for mate.
	 * @throws IllegalArgumentException If one or both arguments is null or if a mate
	 * could not be found.
	 */
	public BreakendPair(StructuralVariant v, Collection<StructuralVariant> vPool)
	{
		//Search for a variant with the matching event.
		if (v == null) throw new IllegalArgumentException();
		StructuralVariant mate = null;
		for (StructuralVariant var : vPool)
		{
			if(v.getEventID().equals(var.getEventID()))
			{
				mate = var;
				break;
			}
		}
		if (mate == null) throw new IllegalArgumentException();
		int vcomp = v.compareTo(mate);
		if (vcomp  > 0)
		{
			var1 = mate;
			var2 = v;
		}
		else
		{
			var1 = v;
			var2 = mate;
		}
		//System.out.println("BreakendPair.constructor || " + v.toString() + " paired with " + mate.toString());
		//vPool.remove(mate);
		constructorCore();
	}
	
	private void constructorCore()
	{
		//String doubleChrom = var1.getChromosome() + ":" + var2.getChromosome();
		CompoundChrom doubleChrom = new CompoundChrom(2);
		doubleChrom.addChrom(var1.getChromosome());
		doubleChrom.addChrom(var2.getChromosome());
		super.setChromosome(doubleChrom);
		super.setPosition(var1.getPosition());
		super.setEndPosition(var2.getPosition());
		super.addAltAllele(generateAltString());
		super.setType(SVType.BND);
		super.setVariantName("BNDEvent_" + var1.getVarID() + ":" + var2.getVarID());
		this.bndType = BNDtype();
		super.setSVLength(0);
		//super.setStrandString(var1.getStrandString() + "||" + var2.getStrandString());
		super.setImprecise(var1.isImprecise() || var2.isImprecise());
		super.setCIDiff(var1.getCIDiff(false, false, false), false, false, false);
		super.setCIDiff(var1.getCIDiff(false, false, true), false, false, true);
		super.setCIDiff(var1.getCIDiff(false, true, false), false, true, false);
		super.setCIDiff(var1.getCIDiff(false, true, true), false, true, true);
		super.setCIDiff(var2.getCIDiff(false, false, false), true, false, false);
		super.setCIDiff(var2.getCIDiff(false, false, true), true, false, true);
		super.setCIDiff(var2.getCIDiff(false, true, false), true, true, false);
		super.setCIDiff(var2.getCIDiff(false, true, true), true, true, true);
		
		//super.setMateID("");
		super.setEventID(var1.getEventID());
		super.setSecondary(false);
		//For evidence, set for whichever end has the least total
		//StructuralVariant lessE;
		//if (var1.getTotalEvidenceCount() > var2.getTotalEvidenceCount()) lessE = var2;
		//else lessE = var1;
		
		super.setProbabilityCurve(var1.getProbabilityCurve(false), false);
		super.setProbabilityCurve(var2.getProbabilityCurve(false), false);
		double q = var1.getQuality();
		if (var2.getQuality() < q) q = var2.getQuality();
		super.setQuality(q);
		
		//Position effect!!
		GeneFunc eff = var1.getGeneFunction();
		if (eff != null && var2.getGeneFunction() != null)
		{
			if (var2.getGeneFunction().getPriority() < eff.getPriority()) eff = var2.getGeneFunction();	
		}
		super.setGeneFunction(eff);
	}
	
	private String generateAltString()
	{
		int o1 = getAltOrientation(var1);
		int o2 = getAltOrientation(var2);
		String c1 = var1.getChromosome().getUDPName();
		String c2 = var2.getChromosome().getUDPName();
		switch (o1)
		{
		case ORIENTATION_RIGHT:
			if (o2 == ORIENTATION_LEFT) return c1 + " ][ " + c2;
		case ORIENTATION_LEFT:
			if (o2 == ORIENTATION_RIGHT) return c2 + " ][ " + c1;
		case ORIENTATION_RIGHTREV:
			if (o2 == ORIENTATION_RIGHTREV) return c1 + " ] " + c2 + " ]";
		case ORIENTATION_LEFTREV:
			if (o2 == ORIENTATION_LEFTREV) return "[ " + c1 + " [ " + c2;
		}
		return c1 + " ? " + c2;
	}
	
	/**
	 * Get the orientation of a breakend by examining the alt allele
	 * string in a variant. This will only return a meaningful value
	 * if the alt allele string is formatted according to the VCF 4.2 documentation for BND SV's.
	 * @param v Variant to examine breakpoint orientation of.
	 * @return Breakpoint orientation (see class constants) detected.
	 */
	public static int getAltOrientation(Variant v)
	{
		String alt = v.getAltAllele(0);
		int iR = alt.indexOf(']');
		int iL = alt.indexOf('[');
		if (iR == 0 && iL < 0) return ORIENTATION_LEFT;
		if (iL == 0 && iR < 0) return ORIENTATION_LEFTREV;
		if (iR > 0 && iL < 0) return ORIENTATION_RIGHTREV;
		if (iL > 0 && iR < 0) return ORIENTATION_RIGHT;
		return ORIENTATION_UNKNOWN;
	}
	
	/**
	 * Determine the breakend orientations of two variants, and if
	 * they are valid, determine whether the two are compatible and could
	 * have been generated from the same event. This is an element of the potential
	 * mate check for BNDs.
	 * @param v1 One of the two variants to examine the orientation of.
	 * @param v2 One of the two variants to examine the orientation of.
	 * @return True - If the BND orientation of both variants could be detected
	 * and if they are compatible for forming a breakend pair.
	 * <br>False - Otherwise.
	 */
	public static boolean orientationMatches(Variant v1, Variant v2)
	{
		int o1 = getAltOrientation(v1);
		int o2 = getAltOrientation(v2);
		switch (o1)
		{
		case ORIENTATION_RIGHT:
			if (o2 == ORIENTATION_LEFT) return true;
			return false;
		case ORIENTATION_LEFT:
			if (o2 == ORIENTATION_RIGHT) return true;
			return false;
		case ORIENTATION_RIGHTREV:
			if (o2 == ORIENTATION_RIGHTREV) return true;
			return false;
		case ORIENTATION_LEFTREV:
			if (o2 == ORIENTATION_LEFTREV) return true;
			return false;
		}
		return false;
	}
	
	private int BNDtype()
	{
		int c1t = var1.getChromosome().getType();
		int c2t = var2.getChromosome().getType();
		if (c1t >= Contig.SORTCLASS_AUTOSOME && c1t <= Contig.SORTCLASS_UNKNOWN)
		{
			if (c2t >= Contig.SORTCLASS_AUTOSOME && c2t <= Contig.SORTCLASS_UNKNOWN)
			{
				int chrComp = var1.getChromosome().compareTo(var2.getChromosome());
				if (chrComp == 0) return BND_INTRACHROMOSOMAL;
				else return BND_INTERCHROMOSOMAL;
			}
			return BND_CHROMOCONTIG;
		}
		return BND_DOUBLECONTIG;
	}
	
	/**
	 * Get the sub-type of BND SV for this variant.
	 * @return Int representation of the BND type (see class constants).
	 */
	public int getBNDtype()
	{
		return this.bndType;
	}
	
	public boolean isChromosome()
	{
		return var1.isChromosome() || var2.isChromosome();
	}
	
	public boolean isContig()
	{
		return var1.isContig() || var2.isContig();
	}
	
	public int getTotalEvidenceCount()
	{
		return var1.getTotalEvidenceCount() + var2.getTotalEvidenceCount();
	}

	public int getTotalEvidenceCount(int sampleIndex)
	{
		return var1.getTotalEvidenceCount(sampleIndex) + var2.getTotalEvidenceCount(sampleIndex);
	}
	
	public int getPEEvidenceCount()
	{
		return var1.getPEEvidenceCount() + var2.getPEEvidenceCount();
	}
	
	public int getPEEvidenceCount(int sampleIndex)
	{
		return var1.getPEEvidenceCount(sampleIndex) + var2.getPEEvidenceCount(sampleIndex);
	}
	
	public int getSREvidenceCount()
	{
		return var1.getSREvidenceCount() + var2.getSREvidenceCount();
	}
	
	public int getSREvidenceCount(int sampleIndex)
	{
		return var1.getSREvidenceCount(sampleIndex) + var2.getSREvidenceCount(sampleIndex);
	}
	
	public int getBDEvidenceCount()
	{
		return var1.getBDEvidenceCount() + var2.getBDEvidenceCount();
	}
	
	public int getBDEvidenceCount(int sampleIndex)
	{
		return var1.getBDEvidenceCount(sampleIndex) + var2.getBDEvidenceCount(sampleIndex);
	}
	
	public String toVCFLine(List<String> orderedSamples, List<InfoDefinition> orderedInfoFields)
	{
		String s = var1.toVCFLine(orderedSamples, orderedInfoFields) + "\n" + var2.toVCFLine(orderedSamples, orderedInfoFields);
		return s;
	}
	
	public String toBEDLine()
	{
		//It will have no idea that this is a BND pair. Be careful with BEDs!
		return var1.toBEDLine() + "\n" + var2.toBEDLine();
	}
	
	public String toBEDPELine()
	{
		//chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2	(custom)
		String s = "";
		//chrom1
		s += var1.getChromosome() + "\t";
		//start1
		s += var1.getCIPosition(false, false, false) + "\t";
		//end1
		s += var1.getCIPosition(false, false, true) + "\t";
		//chrom2
		s += var2.getChromosome() + "\t";
		//start2
		s += var2.getCIPosition(true, false, false) + "\t";
		//end2
		s += var2.getCIPosition(true, false, true) + "\t";
		//name
		s += this.getVarID();
		return s;
	}
	
	private String getBEDStrandFromOrientation(boolean var)
	{
		StructuralVariant v = var1;
		if (var) v = var2;
		int o = getAltOrientation(v);
		switch(o)
		{
		case ORIENTATION_RIGHT: return "+";
		case ORIENTATION_LEFT: return "-";
		case ORIENTATION_RIGHTREV: return "+";
		case ORIENTATION_LEFTREV: return "-";
		default: return ".";
		}
	}
	
	public String toViewerBEDLine(String sample)
	{
		String s = "";
		//var1
		Genotype g = var1.getSampleGenotype(sample);
		if (g == null) return null;
		double altpercent = g.getPercentAlt();
		if (altpercent <= 0.0) return null;
		
		//s += "chr" + var1.getChromosome() + "\t"; //chrom
		String cname = var1.getChromosome().getUCSCName();
		if (cname == null) return null;
		s += cname + "\t"; //chrom
		int st = var1.getCIPosition(false, true, false) - 1;
		int ed = var1.getCIPosition(false, true, true);
		if (st < 0) st = 0;
		int clen = (int)var1.getChromosome().getLength();
		if (ed > clen) ed = clen;
		s += st + "\t"; //start
		s += ed + "\t"; //end
		s += var1.getVarID() + "\t"; //name
		
		int score = (int)Math.round(altpercent * 10.0);
		s += score + "\t"; //score
		
		s += getBEDStrandFromOrientation(false) + "\t"; //strand
		s += (var1.getPosition() - 1) + "\t"; //thickstart
		s += var1.getPosition() + "\t"; //thickend
		s += UCSCGVBED.COLOR_BND + "\t"; //color
		s += "1\t"; //blockcount
		s += (ed - st) + "\t"; //blocksize
		s += "0"; //blockstart
		s += "\n";
		
		
		//var2
		g = var2.getSampleGenotype(sample);
		if (g == null) return null;
		altpercent = g.getPercentAlt();
		if (altpercent <= 0.0) return null;
		
		//s += "chr" + var2.getChromosome() + "\t"; //chrom
		cname = var2.getChromosome().getUCSCName();
		if (cname == null) return null;
		s += cname + "\t"; //chrom
		st = var2.getCIPosition(false, true, false) - 1;
		ed = var2.getCIPosition(false, true, true);
		if (st < 0) st = 0;
		clen = (int)var2.getChromosome().getLength();
		if (ed > clen) ed = clen;
		s += st + "\t"; //start
		s += ed + "\t"; //end
		s += var2.getVarID() + "\t"; //name
		
		score = (int)Math.round(altpercent * 10.0);
		s += score + "\t"; //score
		
		s += getBEDStrandFromOrientation(true) + "\t"; //strand
		s += (var2.getPosition() - 1) + "\t"; //thickstart
		s += var2.getPosition() + "\t"; //thickend
		s += UCSCGVBED.COLOR_BND + "\t"; //color
		s += "1\t"; //blockcount
		s += (ed - st) + "\t"; //blocksize
		s += "0"; //blockstart
		
		return s;
	}
	
	public boolean isOnChromosome(Contig chrom)
	{
		return var1.isOnChromosome(chrom) | var2.isOnChromosome(chrom);
	}
	
	public Collection<Contig> getAllChromosomes()
	{
		Set<Contig> cset = new HashSet<Contig>();
		cset.addAll(var1.getAllChromosomes());
		cset.addAll(var2.getAllChromosomes());
		return cset;
	}
	
	public StructuralVariant getBNDVariant(boolean end)
	{
		if (end) return var2;
		else return var1;
	}

	public Genotype getGenotype1(String sample)
	{
		return var1.getSampleGenotype(sample);
	}
	
	public Genotype getGenotype2(String sample)
	{
		return var2.getSampleGenotype(sample);
	}
	
	public Genotype getSampleGenotype(String sampleName)
	{
		return var1.getSampleGenotype(sampleName);
	}
	
	public String getSampleGenotypeString(String sampleName)
	{
		return var1.getSampleGenotypeString(sampleName);
	}
	
	public boolean inRegion(Contig c, int start, int end, boolean anyend)
	{
		boolean c1 = c.equals(var1.getChromosome());
		boolean c2 = c.equals(var2.getChromosome());
		boolean s = c1 && this.getPosition() >= start && this.getPosition() <= end;
		boolean e = c2 && this.getEndPosition() >= start && this.getEndPosition() <= end;
		if (anyend) return s || e;
		else return s && e;
	}
	
}
