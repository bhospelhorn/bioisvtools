package hospelhornbg_svtools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import hospelhornbg_bioinformatics.Genotype;
import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.UCSCGVBED;
import hospelhornbg_bioinformatics.VCF;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_bioinformatics.VariantPool.InfoDefinition;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.GenomeBuild;

public class PennToVCF {

	public static final String OP_SAMPLENAME = "-n"; 
	public static final String OP_MINSNP = "-i"; 
	public static final String OP_TABLE = "-t"; 
	public static final String OP_OUTPATH = "-o"; 
	public static final String OP_BEDTRACK = "-b"; 
	
	//Read a set of a minsnp file (SNP chip SOP/ PennCNV output) and Genome Studio SNP table
	//Convert minsnp calls to VCF format!
	
	public static class SNP implements Comparable<SNP>
	{
		private String name;
		private Contig chrom;
		private int position;
		
		public SNP(String n, Contig c, int p)
		{
			name = n;
			chrom = c;
			position = p;
		}

		public String getName()
		{
			return name;
		}
		
		public Contig getChrom()
		{
			return chrom;
		}
		
		public int getPos()
		{
			return position;
		}
		
		public int hashCode()
		{
			return name.hashCode();
		}
		
		public boolean equals(Object o)
		{
			if (o == null) return false;
			if (this == o) return true;
			if (!(o instanceof SNP)) return false;	
			SNP s = (SNP)o;
			if (!this.chrom.equals(s.chrom)) return false;
			if (this.position != s.position) return false;
			
			return this.name.equals(s.name);
		}
		
		public int compareTo(SNP o) 
		{
			if (o == null) return 1;
			//Compare chrom
			int comp = this.chrom.compareTo(o.chrom);
			if (comp != 0) return comp;
			
			//Compare pos
			comp = this.position - o.position;
			if (comp != 0) return comp;
			
			return this.name.compareTo(o.name);
		}
	
		public boolean matchesName(String n)
		{
			return name.equals(n);
		}
		
	}
	
	public static class SNPCNV
	{
		private int state;
		private int CN;
		
		private String startSNP;
		private String endSNP;
		
		public SNPCNV(String minsnpline)
		{
			if (minsnpline == null || minsnpline.isEmpty())
			{
				state = -1;
				CN = -1;
				startSNP = null;
				endSNP = null;
				throw new IllegalArgumentException();
			}
			//Replace spaces with tabs...
			minsnpline = minsnpline.replaceAll(" +", "\t");
			
			String[] fields = minsnpline.split("\t");
			if (fields.length != 7)
			{
				System.err.println("PennToVCF.SNPCNV.<init> || Illegal line found: " + minsnpline); 
				System.err.println("\tFields from split: " + fields.length);
				for (String s : fields) System.err.println("\t" + s);
				throw new IllegalArgumentException();
			}
			String[] f4 = fields[3].split(",");
			state = Integer.parseInt(f4[0].replaceAll("state", ""));
			int eq = f4[1].indexOf('=');
			CN = Integer.parseInt(f4[1].substring(eq + 1));
			
			eq = fields[5].indexOf('=');
			if (eq < 0) throw new IllegalArgumentException();
			startSNP = fields[5].substring(eq + 1);
			
			eq = fields[6].indexOf('=');
			if (eq < 0) throw new IllegalArgumentException();
			endSNP = fields[6].substring(eq + 1);
		}
		
		public int getState()
		{
			return state;
		}
		
		public int getCopynumber()
		{
			return CN;
		}
		
		public String getStartSNPName()
		{
			return startSNP;
		}
		
		public String getEndSNPName()
		{
			return endSNP;
		}
		
	}
	
	public static void printUsage()
	{
		System.out.println("--------------------------------------------------------------------------------");
		System.out.println("BioisvTools || PennToVCF");
		System.out.println();
		System.out.println("Purpose: For UDP specifically - converts a PennCNV minsnp output to SV standard VCF");
		System.out.println("by using a Genome Studio generated SNP table for reference.");
		System.out.println();
		System.out.println("Input Formats:");
		System.out.println("\tInput callset must be in [minsnp#] format");
		System.out.println("\tInput table is a [txt] file usually. See -t flag for columns.");
		System.out.println();
		System.out.println("Output Formats:");
		System.out.println("\tPlain text Variant Call Format [vcf]");
		System.out.println();
		System.out.println("Flags:");
		System.out.println("\t-n\tSTRING\t[Required]\t\tSample Name");
		System.out.println("\t-i\tFILE\t[Required]\t\tInput callset in the form of a minsnp file.");
		System.out.println("\t-t\tFILE\t[Required]\t\tGS SNP table. Tab delimited, first three columns must be: SNPID, chrom, pos");
		System.out.println("\t-o\tFILE\t[Required]\t\tOutput file path");
		System.out.println("\t-b\tFILE\t[Optional]\t\tOutput path for optional BED track for UCSC genome browser.");
		System.out.println();
		System.out.println("Sample Usage:");
		System.out.println("java -jar bioisvtools.jar pennconvert -g GRCh37 -n NA12878 -i NA12878.minsnp5 -t mychip_cohort.txt -o NA12878_PennCNV.vcf");
		System.out.println("java -jar bioisvtools.jar pennconvert -g hg38 -n NA12878 -i NA12878.minsnp10 -t mychip_cohort.txt -o NA12878_PennCNV_minsnp10.vcf -b NA12878_PennCNV_track.bed");
		System.out.println();
		System.out.println("--------------------------------------------------------------------------------");
	}
	
	public static int getTableIndex(List<SNP> SNPtable, String SNPname)
	{
		int tlen = SNPtable.size();
		for (int i = 0; i < tlen; i++)
		{
			if (SNPtable.get(i).matchesName(SNPname)) return i;
		}
		return -1;
	}
	
	public static VariantPool spawnPool(GenomeBuild g, String samplename)
	{
		//Populates it with header information
		VariantPool pool = new VariantPool(1);
		pool.addSample(samplename);
		pool.setGenomeBuild(g);
		
		StructuralVariant.addStandardDefs(pool, false);
		
		InfoDefinition def = StructuralVariant.INFODEF_INFO_CIPOS95;
		pool.addInfoFieldDefinition(def.getKey(), def);
		pool.addInfoKeyToActiveList(def.getKey());
		def = StructuralVariant.INFODEF_INFO_CIEND95;
		pool.addInfoFieldDefinition(def.getKey(), def);
		pool.addInfoKeyToActiveList(def.getKey());
		
		
		def = Genotype.INFODEF_CN;
		pool.addFormatFieldDefinition(def.getKey(), def);
		pool.addFormatKeyToActiveList(def.getKey());
		
		return pool;
	}
	
	public static VariantPool generatePool(List<SNP> SNPtable, Set<SNPCNV> myVars, GenomeBuild g, String samplename)
	{
		VariantPool pool = spawnPool(g, samplename);
		int tSize = SNPtable.size();
		
		for (SNPCNV v : myVars)
		{
			//Find start SNP (and flanking)
			String sn = v.getStartSNPName();
			int si = getTableIndex(SNPtable, sn);
			if (si < 0) continue;
			SNP sP = SNPtable.get(si);
			SNP sB = null;
			SNP sA = null;
			if (si - 1 >= 0) sB = SNPtable.get(si - 1);
			if (si + 1 < tSize) sA = SNPtable.get(si + 1);
				//Make sure they are on the same contig
			Contig cs = sP.getChrom();
			if (sB != null && !sB.getChrom().equals(cs)) sB = null;
			if (sA != null && !sA.getChrom().equals(cs)) sA = null;
			
			int stPos = sP.getPos();
			int stCILo = stPos;
			int stCIHi = stPos;
			
			if (sB != null) stCILo = sB.getPos();
			if (sA != null) stCIHi = sA.getPos();
			
			//Find end SNP (and flanking)
			String en = v.getEndSNPName();
			int ei = getTableIndex(SNPtable, en);
			if (ei < 0) continue;
			SNP eP = SNPtable.get(ei);
			SNP eB = null;
			SNP eA = null;
			if (ei - 1 >= 0) eB = SNPtable.get(ei - 1);
			if (ei + 1 < tSize) eA = SNPtable.get(ei + 1);
				//Make sure they are on the same contig
			Contig ce = sP.getChrom();
			if (eB != null && !eB.getChrom().equals(ce)) eB = null;
			if (eA != null && !eA.getChrom().equals(ce)) eA = null;
			
			int edPos = eP.getPos();
			int edCILo = edPos;
			int edCIHi = edPos;
			
			if (eB != null) edCILo = eB.getPos();
			if (eA != null) edCIHi = eA.getPos();
			
			//Generate genotype
			Genotype gt = new Genotype();
			int cn = v.getCopynumber();
			gt.setCopyNumber(cn);
			if (cn == 2) {
				int[] arr = {0,0};
				gt.setAlleles(arr);
			}
			else if (cn != 2 && (cn % 2 == 0))
			{
				int[] arr = {1,1};
				gt.setAlleles(arr);
			}
			else
			{
				int[] arr = {0,1};
				gt.setAlleles(arr);
			}
			
			//Set up Variant and add information
				//Remember, if the ends are on different chroms, make it a BND pair...
			if (cs.equals(ce))
			{
				StructuralVariant sv = new StructuralVariant();
				
				sv.setChromosome(cs);
				sv.setPosition(stPos);
				sv.setEndPosition(edPos);
				if (v.getCopynumber() > 2) sv.setType(SVType.DUP);
				else if (v.getCopynumber() < 2) sv.setType(SVType.DEL);
				else sv.setType(SVType.CNV);
				sv.setRefAllele("N");
				sv.addAltAllele("<" + sv.getType().toString() + ">");
				if (sv.getType() == SVType.DEL) sv.setSVLength(stPos - edPos);
				else sv.setSVLength(edPos - stPos);
				sv.setCIDiff(stCILo - stPos, false, false, false);
				sv.setCIDiff(stCIHi - stPos, false, false, true);
				sv.setCIDiff(edCILo - edPos, true, false, false);
				sv.setCIDiff(edCIHi - edPos, true, false, true);
				sv.setCIDiff(stCILo - stPos, false, true, false);
				sv.setCIDiff(stCIHi - stPos, false, true, true);
				sv.setCIDiff(edCILo - edPos, true, true, false);
				sv.setCIDiff(edCIHi - edPos, true, true, true);
				sv.addGenotype(samplename, gt);
				
				//Add to pool
				pool.addVariant(sv);
			}
			else
			{
				//No support at the moment
				continue;
			}
		
		}
		
		pool.sortVariants();
		return pool;
	}
	
	public static List<SNP> readTable(String tablePath, GenomeBuild g) throws IOException
	{
		if (tablePath == null || tablePath.isEmpty()) throw new IllegalArgumentException();
		
		FileReader fr = new FileReader(tablePath);
		BufferedReader br = new BufferedReader(fr);
		
		//Array Size
		final int arraySize = 1000000;
		List<SNP> table = new ArrayList<SNP>(arraySize);
		
		//Eat the first line
		br.readLine();
		
		String line = null;
		while ((line = br.readLine()) != null)
		{
			String[] fields = line.split("\t");
			if (fields.length < 3){
				br.close();
				fr.close();
				throw new IllegalArgumentException();
			}
			
			String name = fields[0];
			String chromStr = fields[1];
			String posStr = fields[2];
			
			Contig c = g.getContig(chromStr);
			if (c == null) System.err.println("PennToVCF.pennToVCF || WARNING: Contig not found for chrom \"" + chromStr + "\". SNP will be tossed.");
			int pos = Integer.parseInt(posStr);
			
			if (c != null)
			{
				SNP s = new SNP(name, c, pos);
				table.add(s);	
			}
		}
		
		br.close();
		fr.close();
		return table;
	}
	
	public static Set<SNPCNV> readMinsnp(String minsnpPath) throws IOException
	{
		if (minsnpPath == null || minsnpPath.isEmpty()) throw new IllegalArgumentException();
		
		FileReader fr = new FileReader(minsnpPath);
		BufferedReader br = new BufferedReader(fr);
		
		Set<SNPCNV> myset = new HashSet<SNPCNV>();
		String line = null;
		while ((line = br.readLine()) != null)
		{
			if (line.isEmpty()) continue;
			SNPCNV cnv = new SNPCNV(line);
			myset.add(cnv);
		}
		
		br.close();
		fr.close();
		
		return myset;
	}
	
	public static void pennToVCF(String[] args, GenomeBuild g, boolean verbose)
	{
		String minsnpPath = null;
		String tablePath = null;
		String outPath = null;
		String sampleName = null;
		String bedpath = null;
		
		for (int i = 0; i < args.length; i++)
		{
			String s = args[i];
			if (s.equals(OP_SAMPLENAME))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: -n flag MUST be followed by sample name!");
					printUsage();
					System.exit(1);
				}
				sampleName = args[i+1];
			}
			if (s.equals(OP_MINSNP))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: -i flag MUST be followed by path to minsnp file!");
					printUsage();
					System.exit(1);
				}
				minsnpPath = args[i+1];
			}
			if (s.equals(OP_TABLE))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: -t flag MUST be followed by path to GS table file!");
					printUsage();
					System.exit(1);
				}
				tablePath = args[i+1];
			}
			if (s.equals(OP_OUTPATH))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: -o flag MUST be followed by output file path!");
					printUsage();
					System.exit(1);
				}
				outPath = args[i+1];
			}
			if (s.equals(OP_BEDTRACK))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_BEDTRACK + " flag MUST be followed by output path for BED track file!");
					printUsage();
					System.exit(1);
				}
				bedpath = args[i+1];
			}
		}
		
		if (sampleName == null || sampleName.isEmpty())
		{
			System.err.println("ERROR: Sample name is required!");
			printUsage();
			System.exit(1);
		}
		if (minsnpPath == null || minsnpPath.isEmpty())
		{
			System.err.println("ERROR: minsnp file path is required!");
			printUsage();
			System.exit(1);
		}
		if (tablePath == null || tablePath.isEmpty())
		{
			System.err.println("ERROR: GS table file path is required!");
			printUsage();
			System.exit(1);
		}
		if (outPath == null || outPath.isEmpty())
		{
			System.err.println("ERROR: Output file path is required!");
			printUsage();
			System.exit(1);
		}
		
		// Now we can assume all arguments are good.
		
		//Read table
		List<SNP> SNPtable = null;
		if (verbose) System.err.println("PennToVCF.pennToVCF || Loading SNP table...");
		try 
		{
			SNPtable = readTable(tablePath, g);
		} 
		catch (IOException e) {
			System.err.println("ERROR: Table " + tablePath + " is not a valid GS table file!");
			e.printStackTrace();
			System.exit(1);
		}
		if (verbose) System.err.println("PennToVCF.pennToVCF || SNP Table loaded.");
		if (verbose) System.err.println("PennToVCF.pennToVCF || Table size: " + SNPtable.size());
		if (verbose) System.err.println("PennToVCF.pennToVCF || Sorting SNP table...");
		Collections.sort(SNPtable); //Order is super important
		if (verbose) System.err.println("PennToVCF.pennToVCF || SNP Table sorted.");
		
		//Read callset
		Set<SNPCNV> callset = null;
		if (verbose) System.err.println("PennToVCF.pennToVCF || Loading callset...");
		try 
		{
			callset = readMinsnp(minsnpPath);
		} 
		catch (IOException e) {
			System.err.println("ERROR: " + tablePath + " is not a valid minsnp file!");
			e.printStackTrace();
			System.exit(1);
		}
		if (verbose) System.err.println("PennToVCF.pennToVCF || Callset loaded.");
		if (verbose) System.err.println("PennToVCF.pennToVCF || Size: " + callset.size());
		
		
		VariantPool pool = generatePool(SNPtable, callset, g, sampleName);
		try 
		{
			VCF.writeVCF(pool, "BioisvTools", outPath);
		} 
		catch (IOException e) {
			System.err.println("ERROR: " + outPath + " could not be written!");
			e.printStackTrace();
			System.exit(1);
		}
		
		try
		{
			if (bedpath != null && !bedpath.isEmpty())
			{
				int sep = minsnpPath.lastIndexOf(File.separator);
				String minsnpName = minsnpPath;
				if (sep >= 0) minsnpName = minsnpPath.substring(sep + 1);
				sep = tablePath.lastIndexOf(File.separator);
				String tablename = tablePath;
				if (sep >= 0) tablename = tablePath.substring(sep + 1);
				String trackname = "PennCNV " + sampleName;
				String trackDesc = "PennCNV converted callset. minsnp=" + minsnpName + " SNPtable=" + tablename + " genome=" + g.getBuildName();
				UCSCGVBED track = new UCSCGVBED(trackname, pool.getVariants());
				track.setDescription(trackDesc);
				
				track.write(bedpath, sampleName);
			}	
		}
		catch (IOException e)
		{
			System.err.println("ERROR: " + bedpath + " could not be written!");
			e.printStackTrace();
			System.exit(1);
		}
		
		
	}
	
}
