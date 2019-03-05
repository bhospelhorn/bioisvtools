package hospelhornbg_bioinformatics;

import java.util.Map;

import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.GenomeBuild;

public class LiteSV implements Comparable<LiteSV>{
	
	private Contig iStartChrom;
	private int iStart;
	private Interval cipos;
	
	private Contig iEndChrom;
	private int iEnd;
	private Interval ciend;
	
	private SVType eType;
	private String sName;
	
	private String bndPartner;
	
	public LiteSV()
	{
		iStartChrom = null;
		iStart = -1;
		cipos = null;
		iEndChrom = null;
		iEnd = -1;
		ciend = null;
		eType = null;
		
	}
	
	public static LiteSV readFromVCFLine(String vcfLine, GenomeBuild gb)
	{
		if (vcfLine == null || vcfLine.isEmpty()) return null;
		if (vcfLine.startsWith("#")) return null;
		
		//Chrom	pos	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	(Genotypes)
		String[] fields = vcfLine.split("\t");
		if (fields.length < 8) return null;
		
		LiteSV sv = new LiteSV();
		sv.iStartChrom = gb.getContig(fields[0]);
		
		try
		{
			sv.iStart = Integer.parseInt(fields[1]);
		}
		catch(NumberFormatException e) {return null;}
		
		sv.sName = fields[2];
		
		//Everything else is obtained from the INFO field...
		sv.iEndChrom = sv.iStartChrom; //By default
		Map<String, String[]> infoMap = VCF.mapINFOValues(fields[7]);
		
		//CHR2
		String[] val = infoMap.get("CHR2");
		if (val != null && val.length >= 1)
		{
			sv.iEndChrom = gb.getContig(val[0]);
		}
		
		//END
		val = infoMap.get("END");
		if (val != null && val.length >= 1)
		{
			try {sv.iEnd = Integer.parseInt(val[0]);}
			catch(NumberFormatException e) {e.printStackTrace();}
		}
		
		//SVTYPE
		val = infoMap.get("SVTYPE");
		if (val != null && val.length >= 1)
		{
			sv.eType = SVType.getType(val[0]);
		}
		
		//CIPOS
		val = infoMap.get("CIPOS");
		if (val != null && val.length >= 2)
		{
			try 
			{
				int v1 = Integer.parseInt(val[0]);
				int v2 = Integer.parseInt(val[1]);
				sv.cipos = new Interval(v1, v2);
			}
			catch(NumberFormatException e) {e.printStackTrace();}
		}
		if (sv.cipos == null) sv.cipos = new Interval(0,0);
		
		//CIEND
		val = infoMap.get("CIEND");
		if (val != null && val.length >= 2)
		{
			try 
			{
				int v1 = Integer.parseInt(val[0]);
				int v2 = Integer.parseInt(val[1]);
				sv.ciend = new Interval(v1, v2);
			}
			catch(NumberFormatException e) {e.printStackTrace();}
		}
		if (sv.ciend == null) sv.ciend = new Interval(0,0);
		
		//MATEID
		val = infoMap.get("MATEID");
		if (val != null && val.length >= 1)
		{
			sv.bndPartner = val[0];
		}
		
		return sv;
	}
	
	public void setAsEnd(LiteSV partner)
	{
		iEndChrom = partner.iStartChrom;
		iEnd = partner.iStart;
		ciend = partner.cipos;
	}

	private boolean rangeOverlap(Interval i1, Interval i2)
	{
		//TODO: Write
		if (i1.getEnd() < i2.getStart()) return false;
		if (i2.getEnd() < i1.getStart()) return false;
		
		
		
		return true;
	}
	
	public boolean svIsEquivalent(StructuralVariant sv, int bpLeeway, boolean enforceTypeMatch)
	{
		//First, do type match if requested
		if(enforceTypeMatch)
		{
			if(eType != sv.getType()) return false;
		}
		
		//Check POS
		int v1 = iStart + cipos.getStart() - bpLeeway;
		int v2 = iStart + cipos.getEnd() + bpLeeway;
		Interval myPos = new Interval(v1,v2);
		Interval svPos = sv.getCIPOS90();
		if(!rangeOverlap(myPos, svPos)) return false;
		
		//Check END
		v1 = iEnd + ciend.getStart() - bpLeeway;
		v2 = iEnd + ciend.getEnd() + bpLeeway;
		Interval myEnd = new Interval(v1,v2);
		Interval svEnd = sv.getCIEND90();
		if(!rangeOverlap(myEnd, svEnd)) return false;
		
		return true;
	}
	
	public Contig getChrom1()
	{
		return iStartChrom;
	}
	
	public int getPosition()
	{
		return iStart;
	}
	
	public Interval getCIPOS()
	{
		return cipos;
	}
	
	public Contig getChrom2()
	{
		return iEndChrom;
	}
	
	public int getEndPosition()
	{
		return iEnd;
	}
	
	public Interval getCIEND()
	{
		return ciend;
	}
	
	public SVType getSVType()
	{
		return eType;
	}
	
	public String getVariantID()
	{
		return sName;
	}
	
	public String getMateID()
	{
		return this.bndPartner;
	}

	public boolean equals(Object o)
	{
		if (o == null) return false;
		if (o == this) return true;
		if (!(o instanceof LiteSV)) return false;
		
		return this.compareTo((LiteSV)o) == 0;
	}
	
	public int hashCode()
	{
		return sName.hashCode();
	}
	
	@Override
	public int compareTo(LiteSV other) 
	{
		if (other == null) return 1;
		if (this == other) return 0;
		
		//Compare POS
		Contig tc = this.getChrom1();
		Contig oc = other.getChrom1();
		if (tc != null)
		{
			if (!tc.equals(oc)) return tc.compareTo(oc);	
		}
		else if (tc == null && oc != null) return -1;
		
		int tpos = this.getPosition();
		int opos = other.getPosition();
		
		if (tpos != opos) return tpos - opos;
		
		//Compare END
		tc = this.getChrom2();
		oc = other.getChrom2();
		if (tc != null)
		{
			if (!tc.equals(oc)) return tc.compareTo(oc);	
		}
		else if (tc == null && oc != null) return -1;
		
		tpos = this.getEndPosition();
		opos = other.getEndPosition();
		
		if (tpos != opos) return tpos - opos;
		
		//Compare CIs
		int ciposComp = cipos.compareTo(other.cipos);
		if (ciposComp != 0) return ciposComp;
		
		int ciendComp = ciend.compareTo(other.ciend);
		if (ciendComp != 0) return ciposComp;
		
		//Compare Type
		return this.eType.compareTo(other.eType);
	}
	
	
	
}
