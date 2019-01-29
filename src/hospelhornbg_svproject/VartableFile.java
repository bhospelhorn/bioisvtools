package hospelhornbg_svproject;

import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import hospelhornbg_bioinformatics.BreakendPair;
import hospelhornbg_bioinformatics.Genotype;
import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.Translocation;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_bioinformatics.VariantPool.InfoDefinition;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.Gene;
import hospelhornbg_genomeBuild.GenomeBuild;
import hospelhornbg_segregation.Candidate;
import hospelhornbg_segregation.CandidateFlags;
import hospelhornbg_segregation.Family;
import waffleoRai_Utils.FileBuffer;

public class VartableFile {
	
	public static final String MAGIC = "VARTABLE";
	public static final String MAGIC_HUFF = "VARTHUFF";
	public static final int CURRENT_VERSION = 1;
	
	/** ID = CANDFLAGS
	* <br>Number = Variable
	* <br>Type = String
	* <br>Description = "Flags for transcript specific candidates in form TranscriptID0:FlagBitField0,TranscriptID1:FlagBitField1 (etc)"
	*/
	public static final InfoDefinition INFODEF_INFO_CANDIDATEFLAGS = new InfoDefinition("CANDFLAGS", VariantPool.INFODEF_STRING, "Flags for transcript specific candidates in form TranscriptID0:FlagBitField0,TranscriptID1:FlagBitField1 (etc)", VariantPool.INFODEF_NARGS_VARIABLE);
	
	/** ID = VARIDINT
	* <br>Number = 1
	* <br>Type = Integer
	* <br>Description = "Integer representation of varID"
	*/
	public static final InfoDefinition INFODEF_INFO_VARIDINT = new InfoDefinition("VARIDINT", VariantPool.INFODEF_INT, "Integer representation of varID", 1);
	
	/** ID = MATEIDINT
	* <br>Number = 1
	* <br>Type = Integer
	* <br>Description = "Integer representation of mateID"
	*/
	public static final InfoDefinition INFODEF_INFO_MATEIDINT = new InfoDefinition("MATEIDINT", VariantPool.INFODEF_INT, "Integer representation of mateID", 1);
	
	
	public static void writeVartable(GenomeBuild gb, Collection<Candidate> candidates, Family family, String outpath, boolean huff)
	{
		
	}
	
	public static List<Candidate> readVartable(GenomeBuild gb, Family family, String inpath) throws IOException
	{
		//Header structure
		
		//	Magic [8]
		//	Version [4]
		//	# Indivs [4]
		//	GenomeBuild UID [4]
		
		//	Indiv UID Table (4 x #Indivs)
		
		//Assumption is made that the file is too damn big. Try not to hold all in mem.
		
		//Create a temp file for huff decompression...
		String tpath = FileBuffer.generateTemporaryPath("vartable_reader_temp");
		
		return null;
	}
	
	public static String decompressVarTable(String inpath)
	{
		//Return path of temp file
		return null;
	}
	
	public static SVType getSVType(int i)
	{
		switch(i)
		{
		case 0: return SVType.BED_REGION;
		case 1: return SVType.BND;
		case 2: return SVType.CNV;
		case 3: return SVType.DEL;
		case 4: return SVType.DELME;
		case 5: return SVType.DUP;
		case 6: return SVType.INS;
		case 7: return SVType.INSME;
		case 8: return SVType.INV;
		case 9: return SVType.OTHER;
		case 10: return SVType.TANDEM;
		case 11: return SVType.TRA;
		}
		return null;
	}
	
	public static int getSVTypeInt(SVType svt)
	{	
		switch(svt)
		{
		case BED_REGION:
			return 0;
		case BND:
			return 1;
		case CNV:
			return 2;
		case DEL:
			return 3;
		case DELME:
			return 4;
		case DUP:
			return 5;
		case INS:
			return 6;
		case INSME:
			return 7;
		case INV:
			return 8;
		case OTHER:
			return 9;
		case TANDEM:
			return 10;
		case TRA:
			return 11;
		}
		return -1;
	}
	
	public static StructuralVariant parseVariant(FileBuffer info, long stpos, GenomeBuild gb, List<String> genoSamples, List<String> supportStrings)
	{
		// Big-Endian
		
		// VarID [4]
		// Reserved [4] (Was gonna use for MateID, but BNDs should already be paired!)
		// SVType [2]
		// Info Flags [2]
		// Quality [4]
		
		// Start contig [4]
		// Start pos [4]
		// End contig [4]
		// End pos [4]
		
		// CIs [32]
			//CI90 [16]
			//CI95 [16]
			//(posst posed endst ended)
		
		// Supp Vec [2]
		// BND/TRA Orientation [1] (If applicable)
		// Reserved Flags [1]
		// VarName [44]
		
		//-- 112
		
		// #Candidate Records [4]
		//	Candidate Info (x #Candidates) - It'll repartner etc.
		//		Gene UID [4]
		//		AnnoFlags [4]
		
		//This may need to be updated if more info is needed!
		// Genotypes [4 x #FamMembers]
		//		GenoFlags [2]
		//		Allele 1 [1]
		//		Allele 2 [1]
		
		//Total Size: 116 + (8 x #candidates) + (4 x genotypes)
		
		info.setEndian(true);
		
		long cpos = stpos;
		
		//Skip IDs...
		//int varid = info.intFromFile(cpos); cpos += 4;
		//int mateid = info.intFromFile(cpos); cpos += 4;
		cpos += 8;
		
		//Variant basics
		SVType type = getSVType(Short.toUnsignedInt(info.shortFromFile(cpos))); cpos += 2;
		int rawflags = Short.toUnsignedInt(info.shortFromFile(cpos)); cpos += 2;
		int varqual = info.intFromFile(cpos); cpos += 4;
		
		int stcontigid = info.intFromFile(cpos); cpos += 4;
		int varpos = info.intFromFile(cpos); cpos += 4;
		int edcontigid = info.intFromFile(cpos); cpos += 4;
		int varend = info.intFromFile(cpos); cpos += 4;
		
		//Get contigs
		Contig c1 = gb.getContigByUID(stcontigid);
		Contig c2 = gb.getContigByUID(edcontigid);
		
		//CIs
		int ps90 = info.intFromFile(cpos); cpos += 4;
		int pe90 = info.intFromFile(cpos); cpos += 4;
		int es90 = info.intFromFile(cpos); cpos += 4;
		int ee90 = info.intFromFile(cpos); cpos += 4;
		
		int ps95 = info.intFromFile(cpos); cpos += 4;
		int pe95 = info.intFromFile(cpos); cpos += 4;
		int es95 = info.intFromFile(cpos); cpos += 4;
		int ee95 = info.intFromFile(cpos); cpos += 4;
		
		//Support vector
		short sVec = info.shortFromFile(cpos); cpos += 2;
		//TRA/BND orientation
		byte trao = info.getByte(cpos); cpos+=2; //Skip next byte
		//Variant name
		String vname = info.getASCII_string(cpos, 44); cpos += 44;
		
		//Instantiate SV based upon the type...
		StructuralVariant sv = null;
		BreakendPair bp = null; //Clumsy, but I'm too tired to do something more clever
		if (type == SVType.BND) {
			bp = new BreakendPair(c1, varpos, c2, varend);
			sv = bp.getBNDVariant(false);
		}
		else if (type == SVType.TRA) {
			Translocation tra = new Translocation();
			tra.setChromosome(c1);
			tra.setPosition(varpos);
			tra.setChromosome2(c2);
			tra.setEndPosition(varend);
			sv = tra;
		}
		else {
			sv = new StructuralVariant();
			sv.setChromosome(c1);
			sv.setPosition(varpos);
			sv.setEndPosition(varend);
		}
		
		//Candidate Flags
		int ccount = info.intFromFile(cpos); cpos += 4;
		Map<Integer, Integer> cflagmap = new HashMap<Integer,Integer>();
		for (int c = 0; c < ccount; c++)
		{
			int cid = info.intFromFile(cpos); cpos += 4;
			int flags = info.intFromFile(cpos); cpos += 4;
			cflagmap.put(cid, flags);
		}
		
		//Genotypes
		//Load directly into sv
		for(String s : genoSamples)
		{
			int gflags = Short.toUnsignedInt(info.shortFromFile(cpos)); cpos += 2;
			int allele1 = Byte.toUnsignedInt(info.getByte(cpos)); cpos++;
			int allele2 = Byte.toUnsignedInt(info.getByte(cpos)); cpos++;
			//Geno flags:
			//	[0] Is phased?
			
			boolean phased = (gflags & 0x01) != 0;
			Genotype g = new Genotype();
			g.setPhased(phased);
			int[] alleles = {allele1, allele2};
			g.setAlleles(alleles);
			
			sv.addGenotype(s, g);
		}
		
		//Load the remaining data into sv
		//Raw flags:
		//	[0] Is imprecise
		if ((rawflags & 0x1) != 0) sv.setImprecise(true);
		sv.setCIDiff(ps90, false, false, false);
		sv.setCIDiff(pe90, false, false, true);
		sv.setCIDiff(es90, true, false, false);
		sv.setCIDiff(ee90, true, false, true);
		sv.setCIDiff(ps95, false, true, false);
		sv.setCIDiff(pe95, false, true, true);
		sv.setCIDiff(es95, true, true, false);
		sv.setCIDiff(ee95, true, true, true);
		sv.setQuality((double)varqual);
		
		//BND specific stuff
		if (type == SVType.BND)
		{
			StructuralVariant var1 = bp.getBNDVariant(false);
			StructuralVariant var2 = bp.getBNDVariant(true);
			sv = bp;
			var1.setVariantName(vname + "_1");
			var2.setVariantName(vname + "_2");
		}
		
		//Load var name
		sv.setVariantName(vname);
		
		//Load support vector
		int suppvec = Short.toUnsignedInt(sVec);
		int mask = 0x1;
		for (String s : supportStrings)
		{
			if ((suppvec & mask) != 0) sv.addSupportMark(s);
			mask = mask << 1;
		}
		
		//Load raw candidate flags
		Set<Integer> keyset = cflagmap.keySet();
		for (Integer k : keyset)
		{
			int rawcflags = cflagmap.get(k);
			sv.addCandidateFlagsNotation(k, CandidateFlags.readVarTableField(rawcflags));
		}
		
		
		//If bp or tra, load orientation
		if (type == SVType.BND || type == SVType.TRA)
		{
			int trai = Byte.toUnsignedInt(trao);
			int o1 = (trai >>> 4) & 0xF;
			int o2 = trai & 0xF;
			
			if (type == SVType.BND)
			{
				bp.setOrientation1(o1);
				bp.setOrientation2(o2);
			}
			else if (type == SVType.TRA)
			{
				Translocation tra = (Translocation)sv;
				tra.setOrientation1(o1);
				tra.setOrientation2(o2);
			}	
		}
		
		sv.setRefAllele("N");
		sv.addAltAllele("<" + sv.getType().name() + ">");
		

		return sv;
	}
	
	public static FileBuffer serializeVariant(Collection<Candidate> varCandidates, List<String> genoSamples, List<String> supportStrings)
	{
		//Extracts the variant from the first Candidate, then tosses any other Candidates that don't have the variant!
		if (varCandidates == null || varCandidates.isEmpty()) return null;
		
		int cnum = varCandidates.size();
		int ngeno = genoSamples.size();
		
		int maxrsz = 116 + (cnum * 8) + (ngeno * 4);
		FileBuffer svar = new FileBuffer(maxrsz, true);
		List<int[]> candinfo = new LinkedList<int[]>();
		
		Variant v = null;
		for (Candidate c : varCandidates)
		{
			if (v == null) v = c.getVariant();
			if (v != c.getVariant()) continue; //Variant must be SAME REFERENCE!
			int tid = -1; //If intergenic
			Gene g = c.getGene();
			if (g != null) tid = g.getID().hashCode();
			CandidateFlags cf = c.getFlags();
			int cfi = 0;
			if (cf != null) cfi = cf.getVarTableField();
			int[] cfield = {tid, cfi};
			candinfo.add(cfield);
		}
		if (v == null) return null;
		
		//**Serialize variant info
		
		svar.addToFile(v.getVarID().hashCode()); //VarID
		svar.addToFile(0); //MateID
		if (v instanceof StructuralVariant)
		{
			StructuralVariant sv = (StructuralVariant)v;
			svar.addToFile((short)getSVTypeInt(sv.getType())); //SVType
			int flags = 0;
			if (sv.isImprecise()) flags |= 0x1;
			svar.addToFile((short)flags); //Info flags
		}
		else
		{
			svar.addToFile((short)-1); //SVType
			svar.addToFile((short)0); //Info flags
		}
		svar.addToFile((int)Math.round(v.getQuality())); //Quality
		Contig c = v.getChromosome();
		int pos = v.getPosition();
		svar.addToFile(c.getUDPName().hashCode()); //Start Contig
		svar.addToFile(pos); //Start pos
		if (v instanceof StructuralVariant)
		{
			StructuralVariant sv = (StructuralVariant)v;
			c = sv.getEndChromosome();
			pos = sv.getEndPosition();
			svar.addToFile(c.getUDPName().hashCode()); //End Contig
			svar.addToFile(pos); //End pos
			svar.addToFile(sv.getCIDiff(false, false, false)); //CIPos90 Start
			svar.addToFile(sv.getCIDiff(false, false, true)); //CIPos90 End
			svar.addToFile(sv.getCIDiff(true, false, false)); //CIEnd90 Start
			svar.addToFile(sv.getCIDiff(true, false, true)); //CIEnd90 End
			svar.addToFile(sv.getCIDiff(false, true, false)); //CIPos95 Start
			svar.addToFile(sv.getCIDiff(false, true, true)); //CIPos95 End
			svar.addToFile(sv.getCIDiff(true, true, false)); //CIEnd95 Start
			svar.addToFile(sv.getCIDiff(true, true, true)); //CIEnd95 End
		}
		else
		{
			svar.addToFile(c.getUDPName().hashCode()); //End Contig
			svar.addToFile(pos); //End pos
			svar.addToFile(0); //CIPos90 Start
			svar.addToFile(0); //CIPos90 End
			svar.addToFile(0); //CIEnd90 Start
			svar.addToFile(0); //CIEnd90 End
			svar.addToFile(0); //CIPos95 Start
			svar.addToFile(0); //CIPos95 End
			svar.addToFile(0); //CIEnd95 Start
			svar.addToFile(0); //CIEnd95 End
		}
		
		//Support Vector
		int suppvec = 0;
		int mask = 0x1;
		for (String suppstr : supportStrings)
		{
			boolean ev = v.supportMarked(suppstr);
			if (ev) suppvec |= mask;
			mask = mask << 1;
		}
		svar.addToFile((short)suppvec);
		
		if (v instanceof BreakendPair)
		{
			BreakendPair bp = (BreakendPair)v;
			int o = bp.getOrientation1();
			o = o << 4;
			o |= bp.getOrientation2();
			svar.addToFile((byte)o); //Translocation orientation
		}
		else if (v instanceof Translocation)
		{
			Translocation tra = (Translocation)v;
			int o = tra.getOrientation1();
			o = o << 4;
			o |= tra.getOrientation2();
			svar.addToFile((byte)o); //Translocation orientation
		}
		else
		{
			svar.addToFile((byte)0);
		}
		
		svar.addToFile((byte)0); //Reserved
		
		//Variant name
		String vname = v.getVarID();
		if (vname.length() > 44) vname = vname.substring(0, 44);
		int vlen = vname.length();
		svar.printASCIIToFile(vname);
		if (vlen < 44)
		{
			for (int i = vlen; i < 44; i++) svar.addToFile((byte)0);
		}
		
		//**Serialize candidate info
		int cn = candinfo.size();
		svar.addToFile(cn); //# Candidates
		for (int[] cinfo : candinfo)
		{
			if (cinfo.length != 2) continue;
			svar.addToFile(cinfo[0]); //Transcript ID
			svar.addToFile(cinfo[1]); //Candidate flags
		}
		
		//**Serialize genotype info
		for (String sname : genoSamples)
		{
			Genotype geno = v.getSampleGenotype(sname);
			if (geno != null)
			{
				int[] alleles = geno.getAlleles();
				int flags = 0;
				if (geno.isPhased()) flags |= 0x1;
				int a1 = 0;
				int a2 = 0;
				if (alleles.length >= 1) a1 = alleles[0];
				if (alleles.length >= 2) a2 = alleles[1];
				svar.addToFile((short)flags);
				svar.addToFile((byte)a1);
				svar.addToFile((byte)a2);
			}
			else
			{
				svar.addToFile((short)0);
				svar.addToFile((byte)0);
				svar.addToFile((byte)0);
			}
		}
		
		return svar;
	}

	public static void indexVarTable(String path)
	{
		//Two indices - one by transcript, one by chrom
		
		//Transcript Index .gidx
		//	Transcript UID [4]
		//	File offset [4]
		//	#Variants [4]
		//	
		
		//Chrom Index .cidx
		//	Contig UID [4]
		//	File offset [4]
		//	#Variants [4]
		//	
		
		
	}
	
	
}
