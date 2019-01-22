package hospelhornbg_svproject;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import hospelhornbg_bioinformatics.BreakendPair;
import hospelhornbg_bioinformatics.Genotype;
import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.Translocation;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_bioinformatics.VariantPool.InfoDefinition;
import hospelhornbg_genomeBuild.Contig;
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
	
	public static List<Candidate> readVartable(GenomeBuild gb, Family family, String inpath)
	{
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
		// VarName [46]
		
		// #Candidate Records [4]
		//	Candidate Info (x #Candidates) - It'll repartner etc.
		//		Gene UID [4]
		//		AnnoFlags [4]
		
		//This may need to be updated if more info is needed!
		// Genotypes [4 x #FamMembers]
		//		GenoFlags [2]
		//		Allele 1 [1]
		//		Allele 2 [1]
		
		
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
		String vname = info.getASCII_string(cpos, 46); cpos += 46;
		
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
		

		return sv;
	}
	
	public static FileBuffer serializeVariant(Collection<Candidate> varCandidates, List<String> genoSamples, List<String> supportStrings)
	{
		//Extracts the variant from the first Candidate, then tosses any other Candidates that don't have the variant!
		
		
		return null;
	}

}
