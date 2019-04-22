package hospelhornbg_svdb;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import hospelhornbg_bioinformatics.Interval;
import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.Gene;
import hospelhornbg_genomeBuild.GeneFunc;
import hospelhornbg_genomeBuild.GeneSet;
import hospelhornbg_genomeBuild.GenomeBuild;
import hospelhornbg_segregation.Population;
import waffleoRai_Utils.BinFieldSize;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.SerializedString;

public class DBVariant implements Comparable<DBVariant>{
	
	public static final String ID_INFO_KEY = "SVDBid";
	
	private String sID;
	private long lID;
	private SVType eType;
	private Contig oChrom;
	
	private Interval iStart;
	private Interval iEnd;
	
	private int iCohortTotalCount;
	private int iCohortHomCount;
	private double fCohortPopFreq;
	
	private Map<Population, Integer> mPopTotalCounts;
	private Map<Population, Integer> mPopHomCounts;
	private Map<Population, Double> mPopFreqs;
	
	private List<Gene> lGenes;
	private GeneFunc ePosEff;
	private String sValidationNotes;
	//private String omimNotes;
	
	private Contig oChrom2; //TRA only
	private String sAlt; //INS and INS:ME only
	
	public static class ParsedVariant
	{
		private long size;
		private DBVariant variant;
		
		private ParsedVariant(DBVariant v, long sz)
		{
			size = sz;
			variant = v;
		}
		
		public long getSize()
		{
			return size;
		}
		
		public DBVariant getVariant()
		{
			return variant;
		}
	}
	
	private DBVariant()
	{
		mPopTotalCounts = new TreeMap<Population, Integer>();
		mPopHomCounts = new TreeMap<Population, Integer>();
		mPopFreqs = new TreeMap<Population, Double>();
	}
	
	public static DBVariant getFromDBRecord(String record, GenomeBuild gb, GeneSet gs)
	{
		if (record == null || record.isEmpty()) return null;
		if (gb == null) return null;
		
		DBVariant var = new DBVariant();
		String[] fields = record.split("\t");
		
		//ID	ID(HexInt)	Type	Chrom	Start	End	UDPTotal	UDPHom	UDPPopFreq (Repeat for each pop...)
		//GeneList	PosEff	OMIM?	ValidationNotes	[Type specific]
		
		//Type Specific Fields...
		//TRA
		//	Chrom 2
		//INS & INS:ME
		//	InsertionSeq
		
		try
		{
			var.sID = fields[0];
			var.lID = Long.parseUnsignedLong(fields[1], 16);
			var.eType = SVType.getType(fields[2]);
			var.oChrom = gb.getContig(fields[3]);
			
			String[] st = fields[4].split("-");
			int s1 = Integer.parseInt(st[0]);
			int s2 = Integer.parseInt(st[1]);
			var.iStart = new Interval(s1, s2);
			
			String[] ed = fields[5].split("-");
			int e1 = Integer.parseInt(ed[0]);
			int e2 = Integer.parseInt(ed[1]);
			var.iEnd = new Interval(e1, e2);
			
			var.iCohortTotalCount = Integer.parseInt(fields[6]);
			var.iCohortHomCount = Integer.parseInt(fields[7]);
			var.fCohortPopFreq = Double.parseDouble(fields[8]);
			
			Population[] allpop = Population.values();
			int j = 9;
			for(Population p : allpop)
			{
				int i1 = Integer.parseInt(fields[j]);
				int i2 = Integer.parseInt(fields[j+1]);
				double f1 = Double.parseDouble(fields[j+2]);
				
				var.mPopTotalCounts.put(p, i1);
				var.mPopHomCounts.put(p, i2);
				var.mPopFreqs.put(p, f1);
				
				j+=3;
			}
			
			if (gs == null) j++;
			else
			{
				String[] glist = fields[j].split(";");
				j++;
				var.lGenes = new LinkedList<Gene>();
				for (String gstr : glist)
				{
					List<Gene> g = gs.getGeneByName(gstr);
					var.lGenes.addAll(g);
				}
			}
			
			var.ePosEff = GeneFunc.getFunction(fields[j]);
			j++;
			
			//Skip OMIM
			j++;
			
			if (j < fields.length) var.sValidationNotes = fields[j];
			
			//Read type specific...
			if(var.eType == SVType.TRA)
			{
				j++;
				if(fields.length >= j+1)
				{
					String cstr = fields[j];
					var.oChrom2 = gb.getContig(cstr);
				}
			}
			else if (var.eType == SVType.INS || var.eType == SVType.INSME)
			{
				j++;
				if(fields.length >= j+1)
				{
					var.sAlt = fields[j];
				}
			}
			
		}
		catch(ArrayIndexOutOfBoundsException e)
		{
			System.err.println("DBVariant.getFromDBRecord || Error: Insufficient fields found in variant record! Returning variant...");
			return var;
		}
		catch(NumberFormatException e1)
		{
			System.err.println("DBVariant.getFromDBRecord || Error: An invalid value was found in a number field! Skipping variant...");
			return null;
		}
		
		return var;
	}
	
	public static DBVariant getFromVariant(StructuralVariant sv, String name)
	{
		//Does NOT do anything with population fields!!
		//Does NOT do anything with gene information!
		if (sv == null) return null;
		DBVariant var = new DBVariant();
		
		var.sID = name;
		var.regenerateIntegerID();
		var.eType = sv.getType();
		var.oChrom = sv.getChromosome();
		
		if(sv.isImprecise())
		{
			//Grab the CIs
			var.iStart = sv.getCIPOS90();
			var.iEnd = sv.getCIEND90();
		}
		else
		{
			var.iStart = new Interval(sv.getPosition(), sv.getPosition());
			var.iEnd = new Interval (sv.getEndPosition(), sv.getEndPosition());
		}
		
		if(var.eType == SVType.TRA)
		{
			var.oChrom2 = sv.getEndChromosome();
		}
		else if (var.eType == SVType.INS || var.eType == SVType.INSME)
		{
			var.sAlt = sv.getAltAllele(0);
		}
		
		return var;
	}

	public static ParsedVariant getFromVDBRecord(FileBuffer record, GenomeBuild gb, GeneSet gs, long stoff)
	{
		//VDB Format
		// Magic "VDB_" [4]
		// Version [4]
		//	Record [Variable]
		//		Variant UID [8]
		//		Contig1 [4]
		//		Start1 [4]
		//		Start2 [4]
		//		End1 [4]
		//		End2 [4]
		//		SVType [1]
		//		PosEff [1]
		//		Var Name [2x2 VLS]
		//		Total Count [4]
		//		Total Hom Count [4]
		//		(Population Counts) [4+4]
		//		# Genes[4]
		//		Gene List
		//			Gene Hash [4]...
		//		Flags[2]
		//			0 - Has OMIM String
		//			1 - Has Validation Notes
		//		OMIM String [VLS 2x2, if present]
		//		Validation Notes [VLS 2x2, if present]
		//		Contig2 [4] (TRA Only)
		//		INS Seq [4x4 VLS] (INS Only)
		
		DBVariant var = new DBVariant();
		long cpos = stoff;
		
		var.lID = record.longFromFile(cpos); cpos += 8;
		int cid1 = record.intFromFile(cpos); cpos += 4;
		var.oChrom = gb.getContigByUID(cid1);
		int st1 = record.intFromFile(cpos); cpos += 4;
		int st2 = record.intFromFile(cpos); cpos += 4;
		var.iStart = new Interval(st1, st2);
		int ed1 = record.intFromFile(cpos); cpos += 4;
		int ed2 = record.intFromFile(cpos); cpos += 4;
		var.iEnd = new Interval(ed1, ed2);
		int sve = Byte.toUnsignedInt(record.getByte(cpos)); cpos++;
		var.eType = SVType.getTypeByID(sve);
		int effe = Byte.toUnsignedInt(record.getByte(cpos)); cpos++;
		var.ePosEff = GeneFunc.getByValue(effe);
		SerializedString ss = record.readVariableLengthString(cpos, BinFieldSize.WORD, 2);
		cpos += ss.getSizeOnDisk();
		var.sID = ss.getString();
		
		//Population counts
		var.iCohortTotalCount = record.intFromFile(cpos); cpos += 4;
		var.iCohortHomCount = record.intFromFile(cpos); cpos += 4;
		Population[] pops = Population.values();
		for(Population p : pops)
		{
			var.mPopTotalCounts.put(p, record.intFromFile(cpos)); cpos += 4;
			var.mPopHomCounts.put(p, record.intFromFile(cpos)); cpos += 4;
		}
		
		//Gene List
		int ngenes = record.intFromFile(cpos); cpos += 4;
		var.lGenes = new ArrayList<Gene>(ngenes + 1);
		if (ngenes > 0)
		{
			for(int i = 0; i < ngenes; i++)
			{
				int gid = record.intFromFile(cpos); cpos += 4;
				Gene g = gs.getGeneByTranscriptHashUID(gid);
				if(g != null) var.lGenes.add(g);
			}
		}
		
		//Optional Strings
		int flags = Short.toUnsignedInt(record.shortFromFile(cpos)); cpos += 2;
		//boolean hasomim = (flags & 0x1) != 0;
		boolean hasvalcom = (flags & 0x2) != 0;
		if(hasvalcom)
		{
			ss = record.readVariableLengthString(cpos, BinFieldSize.WORD, 2);
			cpos += ss.getSizeOnDisk();
			var.sValidationNotes = ss.getString();
		}
		
		//Optional fields
		if(var.eType == SVType.TRA)
		{
			//Chrom 2
			int cid2 = record.intFromFile(cpos); cpos += 4;
			var.oChrom2 = gb.getContigByUID(cid2);
		}
		else if (var.eType == SVType.INS || var.eType == SVType.INSME)
		{
			ss = record.readVariableLengthString(cpos, BinFieldSize.DWORD, 4);
			cpos += ss.getSizeOnDisk();
			var.sAlt = ss.getString();
		}
		
		long sz = cpos - stoff;
		
		return new ParsedVariant(var, sz);
	}
	
	public void countIndividual(boolean hom, Collection<Population> popGroups, int total, Map<Population, Integer> groupTotals)
	{
		iCohortTotalCount++;
		if(hom)iCohortHomCount++;
		fCohortPopFreq = (double)iCohortTotalCount/(double)total;
		
		if(popGroups == null) return;
		for(Population p : popGroups)
		{
			Integer i = mPopTotalCounts.get(p);
			if (i == null) mPopTotalCounts.put(p, 1);
			else mPopTotalCounts.put(p, i+1);
			
			if (hom)
			{
				i = mPopHomCounts.get(p);
				if (i == null) mPopHomCounts.put(p, 1);
				else mPopHomCounts.put(p, i+1);
			}
			
			double pfreq = 0.0;
			int tot = groupTotals.get(p);
			pfreq = (double)mPopTotalCounts.get(p)/(double)tot;
			mPopFreqs.put(p, pfreq);
		}
	}
	
	public void removeIndividual(boolean hom, Collection<Population> popGroups, int total, Map<Population, Integer> groupTotals)
	{
		iCohortTotalCount--;
		if(hom)iCohortHomCount--;
		fCohortPopFreq = (double)iCohortTotalCount/(double)total;
		
		if(popGroups == null) return;
		for(Population p : popGroups)
		{
			Integer i = mPopTotalCounts.get(p);
			if (i == null) mPopTotalCounts.put(p, 0);
			else mPopTotalCounts.put(p, i-1);
			
			if (hom)
			{
				i = mPopHomCounts.get(p);
				if (i == null) mPopHomCounts.put(p, 0);
				else mPopHomCounts.put(p, i-1);
			}
			
			double pfreq = 0.0;
			int tot = groupTotals.get(p);
			pfreq = (double)mPopTotalCounts.get(p)/(double)tot;
			mPopFreqs.put(p, pfreq);
		}
	}
	
	public void noteGenes(GeneSet gs)
	{
		this.ePosEff = GeneFunc.INTERGENIC;
		List<Gene> glist = gs.getGenesInRegion(oChrom, iStart.getStart(), iEnd.getEnd());
		lGenes= glist;
		if (glist == null) return;
		
		//Get poseff
		for (Gene g : glist)
		{
			GeneFunc gf = g.getRelativeRegionLocationEffect(iStart.getStart(), iEnd.getEnd());
			if (gf.getPriority() < ePosEff.getPriority()) ePosEff = gf;
			if (ePosEff == GeneFunc.EXONIC) return;
		}
	}
	
	public void setValidationNotes(String notes)
	{
		this.sValidationNotes = notes;
	}
	
	public String toDBRecord()
	{
		StringBuilder sb = new StringBuilder(2048);
		sb.append(sID + "\t");
		sb.append(Long.toHexString(lID) + "\t");
		sb.append(eType.getString() + "\t");
		sb.append(oChrom.getUDPName() + "\t");
		sb.append(iStart.getStart() + "-" + iStart.getEnd() + "\t");
		sb.append(iEnd.getStart() + "-" + iEnd.getEnd() + "\t");
		sb.append(iCohortTotalCount + "\t");
		sb.append(iCohortHomCount + "\t");
		sb.append(fCohortPopFreq + "\t");
		
		for(Population p : Population.values())
		{
			Integer i = mPopTotalCounts.get(p);
			int n = 0;
			if (i != null) n = i;
			sb.append(n + "\t");
			
			i = mPopHomCounts.get(p);
			n = 0;
			if (i != null) n = i;
			sb.append(n + "\t");
			
			Double d = mPopFreqs.get(p);
			double f = 0.0;
			if (d != null) f = d;
			sb.append(f + "\t");
		}
		
		if (lGenes != null)
		{
			boolean first = true;
			Set<String> genenames = new TreeSet<String>();
			for(Gene g : lGenes)
			{
				genenames.add(g.getName());
			}
			for(String s : genenames)
			{
				if(!first) sb.append(";");
				sb.append(s);
				first = false;
			}
		}
		sb.append("\t");
		
		if(ePosEff != null) sb.append(ePosEff.toString());
		sb.append("\t");
		
		sb.append("\t");
		
		if(this.sValidationNotes != null) sb.append(sValidationNotes);
		
		if(eType == SVType.TRA)
		{
			sb.append("\t");
			if(oChrom2 == null) sb.append(oChrom.getUDPName());
			else sb.append(oChrom2.getUDPName());
		}
		else if (eType == SVType.INS || eType == SVType.INSME)
		{
			sb.append("\t");
			sb.append(sAlt);
		}
		
		return sb.toString();
	}

	public FileBuffer toVDBRecord()
	{
		/*VDB Format
		 Magic "VDB_" [4]
		 Version [4]
			Record [Variable]
				Variant UID [8]
				Contig1 [4]
				Start1 [4]
				Start2 [4]
				End1 [4]
				End2 [4]
				SVType [1]
				PosEff [1]
				Var Name [2x2 VLS]
				Total Count [4]
				Total Hom Count [4]
				(Population Counts) [4+4]
				# Genes[4]
				Gene List
					Gene Hash [4]...
				Flags[2]
					0 - Has OMIM String
					1 - Has Validation Notes
				OMIM String [VLS 2x2, if present]
				Validation Notes [VLS 2x2, if present]
				Contig2 [4] (TRA Only)
				INS Seq [4x4 VLS] (INS Only)*/
		
		int minsz = 30;
		minsz += sID.length() + 4;
		minsz += 8 + (8*Population.values().length);
		minsz += 4;
		if(lGenes != null) minsz += 4 * lGenes.size();
		minsz += 2;
		if(sValidationNotes != null) minsz += sValidationNotes.length() + 4;
		minsz += 4;
		if (this.sAlt != null) minsz += sAlt.length() + 8;
		
		FileBuffer record = new FileBuffer(minsz, true);
		record.addToFile(lID);
		if(oChrom != null) record.addToFile(oChrom.getUDPName().hashCode());
		else record.addToFile(-1);
		record.addToFile(iStart.getStart());
		record.addToFile(iStart.getEnd());
		record.addToFile(iEnd.getStart());
		record.addToFile(iEnd.getEnd());
		if(eType != null) record.addToFile((byte)eType.getID());
		else record.addToFile((byte)0xFF);
		if(ePosEff != null) record.addToFile((byte)ePosEff.getPriority());
		else record.addToFile((byte)0xFF);
		record.addVariableLengthString(sID, BinFieldSize.WORD, 2);
		
		//Population Counts
		record.addToFile(iCohortTotalCount);
		record.addToFile(iCohortHomCount);
		Population[] pops = Population.values();
		for(Population p : pops)
		{
			int i = this.mPopTotalCounts.get(p);
			int j = this.mPopHomCounts.get(p);
			record.addToFile(i);
			record.addToFile(j);
		}
		
		//Gene List
		if (lGenes != null)
		{
			record.addToFile(lGenes.size());
			for(Gene g : lGenes) record.addToFile(g.getID().hashCode());
		}
		else record.addToFile(0);
		
		//Optional Fields
		int flags = 0;
		if (sValidationNotes != null && !sValidationNotes.isEmpty()) flags |= 0x2;
		record.addToFile((short)flags);
		if ((flags & 0x2) != 0) record.addVariableLengthString(sValidationNotes, BinFieldSize.WORD, 2);
		
		//Type specific fields
		if(eType == SVType.TRA)
		{
			//Chrom 2
			if(oChrom2 != null) record.addToFile(oChrom2.getUDPName().hashCode());
			else
			{
				if(oChrom != null) record.addToFile(oChrom.getUDPName().hashCode());
				else record.addToFile(-1);
			}
		}
		else if (eType == SVType.INS || eType == SVType.INSME)
		{
			//Insert string
			if(sAlt == null) sAlt = "N";
			record.addVariableLengthString(sAlt, BinFieldSize.DWORD, 4);
		}
		
		return record;
	}
	
	public boolean equals(Object o)
	{
		if (o == null) return false;
		if (o == this) return true;
		if (!(o instanceof DBVariant)) return false;
		
		return this.compareTo((DBVariant)o) == 0;
	}
	
	public int hashCode()
	{
		return sID.hashCode();
	}
	
	@Override
	public int compareTo(DBVariant other) 
	{
		if (other == null) return 1;
		if (this == other) return 0;
		
		//Compare POS
		Contig tc = this.oChrom;
		Contig oc = other.oChrom;
		if (tc != null)
		{
			if (!tc.equals(oc)) return tc.compareTo(oc);	
		}
		else if (tc == null && oc != null) return -1;
		
		int tpos = this.iStart.getStart();
		int opos = other.iStart.getStart();
		
		if (tpos != opos) return tpos - opos;
		
		tpos = this.iStart.getEnd();
		opos = other.iStart.getEnd();
		
		if (tpos != opos) return tpos - opos;
		
		//Compare END
		//tc = this.getChrom2();
		//oc = other.getChrom2();
		//if (tc != null)
		//{
		//	if (!tc.equals(oc)) return tc.compareTo(oc);	
		//}
		//else if (tc == null && oc != null) return -1;
		
		tpos = this.iEnd.getStart();
		opos = other.iEnd.getStart();
		
		if (tpos != opos) return tpos - opos;
		
		tpos = this.iEnd.getEnd();
		opos = other.iEnd.getEnd();
		
		if (tpos != opos) return tpos - opos;
		
		//Compare Type
		return this.eType.compareTo(other.eType);
	}
	
	private boolean rangeOverlap(Interval i1, Interval i2)
	{
		if (i1.getEnd() < i2.getStart()) return false;
		if (i2.getEnd() < i1.getStart()) return false;

		return true;
	}
	
	public boolean svIsEquivalent(StructuralVariant sv, double percLeeway)
	{
		//First, do type match
		if(eType != sv.getType()) return false;
		
		//Check chrom 1
		Contig c1 = sv.getChromosome();
		if (oChrom == null && c1 != null) return false;
		if (oChrom != null && c1 == null) return false;
		if (oChrom != null && c1 != null)
		{
			if (!oChrom.equals(c1)) return false;
		}
		
		int bpLeeway = 0;
		if(eType == SVType.INS || eType == SVType.INSME)
		{
			//Start, end and alt allele must be IDENTICAL!
			//Leeway is 0
			//Check alt allele
			String alt = sv.getAltAllele(0);
			if(sAlt == null && alt != null) return false;
			if(sAlt != null && alt == null) return false;
			if (sAlt != null && alt != null)
			{
				if(!sAlt.equals(alt)) return false;
			}
		}
		else if (eType == SVType.TRA || eType == SVType.BND)
		{
			//Also need to check chrom2
			Contig c2 = sv.getEndChromosome();
			if (oChrom2 == null && c2 != null) return false;
			if (oChrom2 != null && c2 == null) return false;
			if (oChrom2 != null && c2 != null)
			{
				if (!oChrom2.equals(c2)) return false;
			}
			//bpLeeway is constant
			bpLeeway = (int)Math.round(percLeeway * 1000.0);
		}
		else
		{
			int sz = iEnd.getEnd() - iStart.getStart();
			bpLeeway = (int)Math.round(percLeeway * (double)sz);
		}
		
		//Check POS
		int v1 = iStart.getStart() - bpLeeway;
		int v2 = iStart.getEnd() + bpLeeway;
		Interval myPos = new Interval(v1,v2);
		Interval svPos = sv.getCIPOS90();
		if(!rangeOverlap(myPos, svPos)) return false;
		
		//Check END
		v1 = iEnd.getStart() - bpLeeway;
		v2 = iEnd.getEnd() + bpLeeway;
		Interval myEnd = new Interval(v1,v2);
		Interval svEnd = sv.getCIEND90();
		if(!rangeOverlap(myEnd, svEnd)) return false;
		
		return true;
	}
	
	public long getLongID()
	{
		return this.lID;
	}
	
	public int regenerateIntegerID()
	{
		Random r = new Random();
		int iID = 0;
		if(sID == null) iID = r.nextInt();
		else iID = sID.hashCode() ^ r.nextInt();
		return iID;
	}
	
	public void setTotalCount(int total)
	{
		this.iCohortTotalCount = total;
	}
	
	public void setTotalCount(int total, Population p)
	{
		this.mPopTotalCounts.put(p, total);
	}
	
	public void incrementTotalCount()
	{
		iCohortTotalCount++;
	}
	
	public void incrementTotalCount(Population p)
	{
		Integer i = this.mPopTotalCounts.get(p);
		if (i == null) mPopTotalCounts.put(p, 0);
		else this.mPopTotalCounts.put(p, i+1);
	}
	
	public void setHomozygoteCount(int homCount)
	{
		this.iCohortHomCount = homCount;
	}
	
	public void setHomozygoteCount(int homCount, Population p)
	{
		this.mPopHomCounts.put(p, homCount);
	}
	
	public void incrementHomozygoteCount()
	{
		iCohortHomCount++;
	}
	
	public void incrementHomozygoteCount(Population p)
	{
		Integer i = this.mPopHomCounts.get(p);
		if (i == null) mPopHomCounts.put(p, 0);
		else this.mPopHomCounts.put(p, i+1);
	}
	
	public void adjustPopulationFrequency(int totalIndividuals)
	{
		this.fCohortPopFreq = (double)iCohortTotalCount/(double)totalIndividuals;
	}
	
	public void adjustPopulationFrequency(int popIndividuals, Population p)
	{
		Integer nowcount = this.mPopTotalCounts.get(p);
		if (nowcount == null) nowcount = 0;
		double freq = (double)nowcount/(double)popIndividuals;
		this.mPopFreqs.put(p, freq);
	}

	public StructuralVariant toStructuralVariant()
	{
		StructuralVariant sv = new StructuralVariant();
		sv.setChromosome(oChrom);
		int st = iStart.getCenter();
		int ed = iEnd.getCenter();
		sv.setPosition(st);
		sv.setEndPosition(ed);
		sv.setVariantName(sID);
		
		sv.setRefAllele("N");
		sv.addAltAllele("<" + eType.name() + ">");
		sv.setType(eType);
		sv.setSVLength(ed-st);
		
		sv.setCIPOS(iStart);
		sv.setCIEND(iEnd);
		
		sv.addInfoField(ID_INFO_KEY, lID);
		
		return sv;
	}
	
	public Contig getChrom()
	{
		return oChrom;
	}
	
	public boolean inRange(int stPos, int edPos)
	{
		if (stPos >= iEnd.getEnd()) return false;
		if (edPos <= iStart.getStart()) return false;
		return true;
	}
	
	public Interval getStartPosition()
	{
		return iStart;
	}
	
	public Interval getEndPosition()
	{
		return iEnd;
	}
	
	public Contig getEndChrom()
	{
		return oChrom2;
	}
	
	public String getName()
	{
		return sID;
	}
	
	public SVType getType()
	{
		return eType;
	}
	
	public String getAltAlleleString()
	{
		return sAlt;
	}
	
	public int getIndividualCount()
	{
		return this.iCohortTotalCount;
	}
	
	public int getIndividualCount(Population p)
	{
		return this.mPopTotalCounts.get(p);
	}
	
	public int getHomozygoteCount()
	{
		return this.iCohortHomCount;
	}
	
	public int getHomozygoteCount(Population p)
	{
		return this.mPopHomCounts.get(p);
	}
	
	public double getCohortFreq()
	{
		return this.fCohortPopFreq;
	}
	
	public double getCohortFreq(Population p)
	{
		return this.mPopFreqs.get(p);
	}
	
	protected void setLongUID(long uid)
	{
		this.lID = uid;
	}

	public void decrementTotalCount()
	{
		if(iCohortTotalCount > 0) iCohortTotalCount--;
	}
	
	public void decrementTotalCount(Population p)
	{
		Integer i = this.mPopTotalCounts.get(p);
		if (i == null) return;
		this.mPopTotalCounts.put(p, i++);
	}
	
	public void decrementHomozygoteCount()
	{
		if(iCohortHomCount > 0) iCohortHomCount--;
	}
	
	public void decrementHomozygoteCount(Population p)
	{
		Integer i = mPopHomCounts.get(p);
		if (i == null) return;
		mPopHomCounts.put(p, i++);
	}
	
}
