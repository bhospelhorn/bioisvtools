package hospelhornbg_bioinformatics;

import hospelhornbg_bioinformatics.VariantPool.InfoDefinition;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.GenomeBuild;

/*
 * UPDATE NOTES
 * 
 * Initial version date: March 1, 2018 (1.0.0)
 * 
 * 1.0.0 -> 1.0.1 | July 20, 2018
 * 	Added inRegion method
 * 
 * 1.0.1 -> 1.0.2 | August 17, 2018
 * 	Added getEndChromosome method
 * 
 * 1.0.2 -> 1.1.0 | January 22, 2019
 * 	Added orientation fields, getters, setters
 * 
 */



/**
 * Specialized StructuralVariant subclass that represents a "translocation" type variant.
 * This is not the same as a BreakendPair.
 * This SV class is represented in VCF form using the "TRA" type.
 * @author Blythe Hospelhorn
 * @version 1.1.0
 * @since January 22, 2019
 *
 */
public class Translocation extends StructuralVariant{
	
	/* --- Constants --- */
	
	/**
	 * Alt allele type representing a translocation in a single variant record.
	 */
	public static final String INFODEF_ALT_TRA = "TRA";
	
	/**
	 * ID = CHR2
	 * <br>Number = 1
	 * <br>Type = String
	 * <br>Description = "Chromosome for END coordinate in case of a translocation"
	 */
	public static final InfoDefinition INFODEF_INFO_CHR2 = new InfoDefinition("CHR2", VariantPool.INFODEF_STRING, "Chromosome for END coordinate in case of a translocation", 1);
	
	/* --- Instance Variables --- */
	
	private Contig chrom2;
	private int orientation1;
	private int orientation2;
	
	/* --- Construction --- */
	
	public Translocation()
	{
		super();
		chrom2 = null;
		orientation1 = BreakendPair.ORIENTATION_UNKNOWN;
		orientation2 = BreakendPair.ORIENTATION_UNKNOWN;
	}
	
	protected Translocation(StructuralVariant sv, GenomeBuild genome)
	{
		super(sv);
		String c2 = super.getSingleInfoEntry(INFODEF_INFO_CHR2.getKey());
		if (c2 == null || c2.isEmpty()) chrom2 = super.getChromosome();
		else
		{
			chrom2 = genome.getContig(c2);
			if (chrom2 == null) chrom2 = super.getChromosome();
		}
		orientation1 = BreakendPair.ORIENTATION_UNKNOWN;
		orientation2 = BreakendPair.ORIENTATION_UNKNOWN;
		
	}
	
	/* --- Getters --- */
	
	public Contig getChromosome()
	{
		CompoundChrom cc = new CompoundChrom(2);
		cc.addChrom(getChromosome1());
		cc.addChrom(chrom2);
		return cc;
	}
	
	public Contig getChromosome1()
	{
		return super.getChromosome();
	}
	
	public Contig getChromosome2()
	{
		return chrom2;
	}
	
	public int getOrientation1()
	{
		return this.orientation1;
	}
	
	public int getOrientation2()
	{
		return this.orientation2;
	}
	
	/* --- Setters --- */
	
	public void setChromosome2(Contig c)
	{
		chrom2 = c;
	}
		
	public void setOrientation1(int o)
	{
		this.orientation1 = o;
	}
	
	public void setOrientation2(int o)
	{
		this.orientation2 = o;
	}
	
	/* --- Serialization --- */
	
	private String getStrand1()
	{
		String strandinfo = super.getSingleStringInfoEntry(StructuralVariant.INFODEF_INFO_STRANDS.getKey());
		if (strandinfo == null || strandinfo.isEmpty()) return ".";
		char c1 = strandinfo.charAt(0);
		if (c1 == '+') return "+";
		else if (c1 == '-') return "-";
		return ".";
	}
	
	private String getStrand2()
	{
		String strandinfo = super.getSingleStringInfoEntry(StructuralVariant.INFODEF_INFO_STRANDS.getKey());
		if (strandinfo == null || strandinfo.isEmpty()) return ".";
		char c2 = strandinfo.charAt(1);
		if (c2 == '+') return "+";
		else if (c2 == '-') return "-";
		return ".";
	}

	public String toBEDLine()
	{
		//chrom	start	end	name	score	strand	(other things)
		String line1 = getChromosome1() + "\t" + this.getCIPosition(false, false, false) + "\t" + this.getCIPosition(false, false, true) + "\t" + this.getVarID() + "_1";
		String line2 = getChromosome2() + "\t" + this.getCIPosition(true, false, false) + "\t" + this.getCIPosition(true, false, true) + "\t" + this.getVarID() + "_2";
		return line1 + "\n" + line2;
	}
	
	public String toBEDPELine()
	{
		//chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2	(custom)
		String s = "";
		//chrom1
		s += getChromosome1() + "\t";
		//start1
		s += this.getCIPosition(false, false, false) + "\t";
		//end1
		s += this.getCIPosition(false, false, true) + "\t";
		//chrom2
		s += getChromosome2() + "\t";
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
		//chrom chrstart chrend name score strand thickstart thickend rgb blockcount blocksizes blockstarts
		
		String s = "";
		//var1
		Genotype g = super.getSampleGenotype(sample);
		if (g == null) return null;
		double altpercent = g.getPercentAlt();
		if (altpercent <= 0.0) return null;
		
		//s += "chr" + var1.getChromosome() + "\t"; //chrom
		String cname = getChromosome1().getUCSCName();
		if (cname == null) return null;
		s += cname + "\t"; //chrom
		int st = getCIPosition(false, true, false) - 1;
		int ed = getCIPosition(false, true, true);
		if (st < 0) st = 0;
		int clen = (int)getChromosome1().getLength();
		if (ed > clen) ed = clen;
		s += st + "\t"; //start
		s += ed + "\t"; //end
		s += getVarID() + "_POS\t"; //name
		
		int score = (int)Math.round(altpercent * 10.0);
		s += score + "\t"; //score
		
		s += getStrand1() + "\t"; //strand
		s += (getPosition() - 1) + "\t"; //thickstart
		s += getPosition() + "\t"; //thickend
		s += UCSCGVBED.COLOR_TRA + "\t"; //color
		s += "1\t"; //blockcount
		s += (ed - st) + "\t"; //blocksize
		s += "0"; //blockstart
		s += "\n";
		
		
		//var2
		
		//s += "chr" + var2.getChromosome() + "\t"; //chrom
		cname = getChromosome2().getUCSCName();
		if (cname == null) return null;
		s += cname + "\t"; //chrom
		st = getCIPosition(true, true, false) - 1;
		ed = getCIPosition(true, true, true);
		if (st < 0) st = 0;
		clen = (int)getChromosome2().getLength();
		if (ed > clen) ed = clen;
		s += st + "\t"; //start
		s += ed + "\t"; //end
		s += getVarID() + "_END\t"; //name
		
		score = (int)Math.round(altpercent * 10.0);
		s += score + "\t"; //score
		
		s += getStrand2() + "\t"; //strand
		s += (getEndPosition() - 1) + "\t"; //thickstart
		s += getEndPosition() + "\t"; //thickend
		s += UCSCGVBED.COLOR_TRA + "\t"; //color
		s += "1\t"; //blockcount
		s += (ed - st) + "\t"; //blocksize
		s += "0"; //blockstart
		
		return s;
	}
	
	public boolean inRegion(Contig c, int start, int end, boolean anyend)
	{
		boolean c1 = c.equals(this.getChromosome1());
		boolean c2 = c.equals(this.getChromosome2());
		boolean s = c1 && this.getPosition() >= start && this.getPosition() <= end;
		boolean e = c2 && this.getEndPosition() >= start && this.getEndPosition() <= end;
		if (anyend) return s || e;
		else return s && e;
	}
	
	public Contig getEndChromosome()
	{
		return chrom2;
	}
	
}
