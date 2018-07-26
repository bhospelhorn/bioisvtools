package hospelhornbg_svtools;

import java.io.BufferedWriter;
//import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.VCF;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_genomeBuild.GeneFunc;
import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class VarCounter {
	
//Also make it so that it splits everything into different files?
//Don't just output counts to stdout, also to table files
	
	//Also split all by...
		//Non-BND
		//Confirmed (if possible)
	
	//By type	
	//By size
	//By location effect
	//Imprecise?
	//Chrom or contig?
	
	/* --- Constants --- */
	
	public static final String OP_VCFIN = "-i"; 
	public static final String OP_OUTPUT = "-o";
	public static final String OP_SPLIT = "-s";
	public static final String OP_HARDSPLIT = "-S";
	//public static final String OP_THREADS = "-t"; 
	
	public static final int SZCLASS_0_MIN = 0;
	public static final int SZCLASS_1_MIN = 50;
	public static final int SZCLASS_2_MIN = 500;
	public static final int SZCLASS_3_MIN = 1000;
	public static final int SZCLASS_4_MIN = 5000;
	public static final int SZCLASS_5_MIN = 10000;
	public static final int SZCLASS_6_MIN = 50000;
	public static final int SZCLASS_7_MIN = 100000;
	public static final int SZCLASS_8_MIN = 1000000;
	public static final int SZCLASS_9_MIN = 5000000;
	
	public static final int SZCLASS_COUNT = 10;
	public static final int[] SZMINS = {SZCLASS_0_MIN, SZCLASS_1_MIN, SZCLASS_2_MIN,
										SZCLASS_3_MIN, SZCLASS_4_MIN, SZCLASS_5_MIN,
										SZCLASS_6_MIN, SZCLASS_7_MIN, SZCLASS_8_MIN,
										SZCLASS_9_MIN};
	
	/* --- Instance Variables --- */
	
	private Map<SVType, TypePool> typeMap;
	private VariantPool skeletonReference;
	private int nonSV_count;
	private int SV_count;
	
	/* --- Inner Structures --- */
	
	private class TypePool
	{
		private Map<GeneFunc, EffectPool> epools;
		private EffectPool unkPool;
		
		public TypePool()
		{
			epools = new HashMap<GeneFunc, EffectPool>();
		}
		
		public void addVariant(StructuralVariant sv)
		{
			if (sv == null) return;
			GeneFunc eff = sv.getGeneFunction();
			if (eff == null)
			{
				if (unkPool == null) unkPool = new EffectPool();
				unkPool.addVariant(sv);
			}
			else
			{
				EffectPool e = epools.get(eff);
				if (e == null)
				{
					e = new EffectPool();
					epools.put(eff, e);
				}
				e.addVariant(sv);
			}
		}
		
		public int getTotal()
		{
			int tot = 0;
			if (unkPool != null) tot += unkPool.getTotal();
			Set<GeneFunc> allkeys = epools.keySet();
			for (GeneFunc eff : allkeys)
			{
				tot += epools.get(eff).getTotal();
			}
			return tot;
		}
		
		public int getTotal(GeneFunc eff)
		{
			if (eff == null)
			{
				if (unkPool != null) return unkPool.getTotal();
				return 0;
			}
			EffectPool ep = epools.get(eff);
			if (ep == null) return 0;
			return ep.getTotal();
		}
		
		public int getTotal(GeneFunc eff, int sizeClass)
		{
			if (eff == null)
			{
				if (unkPool != null) return unkPool.getTotal(sizeClass);
				return 0;
			}
			EffectPool ep = epools.get(eff);
			if (ep == null) return 0;
			return ep.getTotal(sizeClass);
		}
		
		public int getTotal(int sizeClass)
		{
			int tot = 0;
			if (unkPool != null) tot += unkPool.getTotal(sizeClass);
			Set<GeneFunc> allkeys = epools.keySet();
			for (GeneFunc eff : allkeys)
			{
				tot += epools.get(eff).getTotal(sizeClass);
			}
			return tot;
		}
		
		public VariantPool makePool()
		{
			VariantPool pool = skeletonReference.makeSkeletonCopy();
			Set<GeneFunc> keys = epools.keySet();
			for (GeneFunc eff : keys)
			{
				EffectPool ep = epools.get(eff);
				if (ep == null) continue;
				pool.addVariants(ep.getVariants());
			}
			if (unkPool != null) pool.addVariants(unkPool.getVariants());
			
			return pool;
		}
		
		public VariantPool makePool(GeneFunc eff)
		{
			if (eff == null)
			{
				if (unkPool == null) return null;
				else return unkPool.makePool();
			}
			else
			{
				EffectPool ep = epools.get(eff);
				if (ep == null) return null;
				return ep.makePool();
			}
		}
		
	}
	
	private class EffectPool
	{
		private int[] sizeTally;
		private List<StructuralVariant> variants;
		private boolean plzupdate;
	
		public EffectPool()
		{
			sizeTally = new int[SZCLASS_COUNT];
			variants = new LinkedList<StructuralVariant>();
			plzupdate = true;
		}
		
		public void addVariant(StructuralVariant sv)
		{
			variants.add(sv);
			plzupdate = true;
		}
		
		public void clearSizeTally()
		{
			for (int i = 0; i < SZCLASS_COUNT; i++) sizeTally[i] = 0;
		}
		
		public void updateSizeTally()
		{
			clearSizeTally();
			for (StructuralVariant sv : variants)
			{
				int sz = sv.getAbsoluteSVLength();
				for (int i = 0; i < SZCLASS_COUNT; i++)
				{
					if (sz >= SZMINS[i])
					{
						if ((i + 1) < SZCLASS_COUNT)
						{
							if (sz < SZMINS[i+1])
							{
								sizeTally[i]++;
								break;
							}
						}
						else{
							sizeTally[i]++;
							break;
						}
					}
				}
			}
			plzupdate = false;
		}
		
		public int getTotal()
		{
			return variants.size();
		}
		
		public int getTotal(int sizeclass)
		{
			if (sizeclass < 0 || sizeclass >= SZCLASS_COUNT) return 0;
			if (plzupdate) updateSizeTally();
			return sizeTally[sizeclass];
		}
		
		public List<Variant> getVariants()
		{
			List<Variant> vlist = new ArrayList<Variant>(variants.size());
			vlist.addAll(variants);
			return vlist;
		}
		
		public VariantPool makePool()
		{
			VariantPool pool = skeletonReference.makeSkeletonCopy();
			pool.addVariants(getVariants());
			return pool;
		}
		
	}
	
	/* --- Construction --- */
	
	public VarCounter(VariantPool pool)
	{
		typeMap = new HashMap<SVType, TypePool>();
		nonSV_count = 0;
		SV_count = 0;
		if (pool == null) return;
		skeletonReference = pool.makeSkeletonCopy();
		dumpVariantPool(pool);
	}
	
	protected void dumpVariantPool(VariantPool pool)
	{
		List<Variant> vlist = pool.getVariants();
		for (Variant v : vlist)
		{
			if (v instanceof StructuralVariant)
			{
				StructuralVariant sv = (StructuralVariant)v;
				SVType t = sv.getType();
				if (t == null)
				{
					nonSV_count++;
				}
				else
				{
					SV_count++;
					TypePool tp = typeMap.get(t);
					if (tp == null)
					{
						tp = new TypePool();
						typeMap.put(t, tp);
					}
					tp.addVariant(sv);
				}
			}
			else
			{
				nonSV_count++;
			}
		}
	}
	
	/* --- Totals --- */
	
	public Map<GeneFunc, Integer> tallyByEffect()
	{
		Collection<GeneFunc> allfunc = GeneFunc.getAll();
		Map<GeneFunc, Integer> tallyMap = new HashMap<GeneFunc, Integer>();
		Set<SVType> typeset = typeMap.keySet();
		for(GeneFunc eff : allfunc) tallyMap.put(eff, 0);
		for(SVType t : typeset)
		{
			TypePool tp = typeMap.get(t);
			if (tp != null)
			{
				for(GeneFunc eff : allfunc)
				{
					int count = tp.getTotal(eff);
					tallyMap.put(eff, tallyMap.get(eff) + count);
				}
			}
		}
		return tallyMap;
	}
	
	public Map<Integer, Integer> tallyBySize()
	{
		//Collection<GeneFunc> allfunc = GeneFunc.getAll();
		Map<Integer, Integer> tallyMap = new HashMap<Integer, Integer>();
		Set<SVType> typeset = typeMap.keySet();
		for(int i = 0; i < SZCLASS_COUNT; i++) tallyMap.put(i, 0);
		for(SVType t : typeset)
		{
			TypePool tp = typeMap.get(t);
			if (tp != null)
			{
				for(int i = 0; i < SZCLASS_COUNT; i++)
				{
					int count = tp.getTotal(i);
					tallyMap.put(i, tallyMap.get(i) + count);
				}
			}
		}
		return tallyMap;
	}
	
	/* --- Output --- */
	
	public String sizeClassToString(int sizeclass)
	{
		if (sizeclass < 0) return null;
		if (sizeclass >= SZCLASS_COUNT) return null;
		String s = "";
		s += SZMINS[sizeclass] + "bp";
		if (sizeclass < SZCLASS_COUNT - 1)
		{
			s += " - " + (SZMINS[sizeclass + 1] - 1) + "bp";
		}
		else s += " +";
		
		return s;
	}

	public void printReport()
	{
		//Print some stats to stdout
		System.out.println("Total Variants: " + (nonSV_count + SV_count));
		System.out.println("Structural Variants: " + SV_count);
		System.out.println("Non SVs: " + nonSV_count);
		System.out.println();
		
		System.out.println("--- Type Breakdown ---");
		List<SVType> allTypes = new ArrayList<SVType>(typeMap.size());
		allTypes.addAll(typeMap.keySet());
		Collections.sort(allTypes);
		for (SVType t : allTypes)
		{
			TypePool tp = typeMap.get(t);
			int tot = 0;
			if (tp != null) tot = tp.getTotal();
			System.out.println("\t" + t.getString() + "\t" + tot);
		}
		System.out.println();
		
		System.out.println("--- Effect Breakdown ---");
		Map<GeneFunc, Integer> efftally = tallyByEffect();
		Collection<GeneFunc> allfunc = GeneFunc.getAll();
		for (GeneFunc eff : allfunc)
		{
			System.out.println("\t" + eff.toString() + "\t" + efftally.get(eff));
		}
		System.out.println();
		
		System.out.println("--- Size Breakdown ---");
		Map<Integer, Integer> sztally = tallyBySize();
		for (int i = 0; i < SZCLASS_COUNT; i++)
		{
			System.out.println("\t" + sizeClassToString(i) + "\t" + sztally.get(i));
		}
		System.out.println();
		
	}
	
	public void generateTables(String prefix)
	{
		//A tsv file for each type.
		//Columns: Effect
		//Rows: Size class
		
		List<SVType> allTypes = new ArrayList<SVType>(typeMap.size());
		allTypes.addAll(typeMap.keySet());
		Collection<GeneFunc> allfunc = GeneFunc.getAll();
		
		for (SVType t : allTypes)
		{
			String outfile = prefix + "_" + t.getString() + ".tsv";
			TypePool tp = typeMap.get(t);
			if (tp == null) continue;
			try 
			{
				FileWriter fw = new FileWriter(outfile);
				BufferedWriter bw = new BufferedWriter(fw);
				
				bw.write("SizeClass\t");
				for (GeneFunc eff : allfunc)
				{
					bw.write(eff.toString() + "\t");
				}
				bw.write("UNKNOWN");
				bw.newLine();
				
				for (int i = 0; i < SZCLASS_COUNT; i++)
				{
					bw.write(sizeClassToString(i) + "\t");
					for (GeneFunc eff : allfunc)
					{
						bw.write(tp.getTotal(eff, i) + "\t");
					}
					bw.write(Integer.toString(tp.getTotal(null, i)));
					if (i < SZCLASS_COUNT - 1) bw.newLine();
				}
				
				bw.close();
				fw.close();
			} 
			catch (IOException e) 
			{
				System.err.println("ERROR: Table " + outfile + " could not be written!");
				e.printStackTrace();
			}
			
		}
	}
	
	public void printVCFs(String prefix, boolean splitByEffect)
	{
		List<SVType> allTypes = new ArrayList<SVType>(typeMap.size());
		allTypes.addAll(typeMap.keySet());
		Collection<GeneFunc> allfunc = GeneFunc.getAll();
		
		for (SVType t : allTypes)
		{
			TypePool tp = typeMap.get(t);
			if (tp == null) continue;
			if (splitByEffect)
			{
				for (GeneFunc eff : allfunc)
				{
					String outfile = prefix + "_" + t.getString() + "_" + eff.toString() + ".vcf";
					VariantPool pool = tp.makePool(eff);
					if (pool == null) continue;
					pool.sortVariants();
					try 
					{
						VCF.writeVCF(pool, "bioisvtools", outfile);
					} 
					catch (IOException e) 
					{
						System.err.println("ERROR: Set " + outfile + " could not be written!");
						e.printStackTrace();
					}
				}
				String outfile = prefix + "_" + t.getString() + "_UNK.vcf";
				VariantPool pool = tp.makePool(null);
				if (pool == null) continue;
				pool.sortVariants();
				try 
				{
					VCF.writeVCF(pool, "bioisvtools", outfile);
				} 
				catch (IOException e) 
				{
					System.err.println("ERROR: Set " + outfile + " could not be written!");
					e.printStackTrace();
				}
			}
			else
			{
				String outfile = prefix + "_" + t.getString() + ".vcf";
				VariantPool pool = tp.makePool();
				if (pool == null) continue;
				pool.sortVariants();
				try 
				{
					VCF.writeVCF(pool, "bioisvtools", outfile);
				} 
				catch (IOException e) 
				{
					System.err.println("ERROR: Set " + outfile + " could not be written!");
					e.printStackTrace();
				}
			}
		}
		
	}
	
	/* --- bioisvtools Front End --- */
	
	public static void printUsage()
	{
		System.out.println("--------------------------------------------------------------------------------");
		System.out.println("BioisvTools || svqstats");
		System.out.println();
		System.out.println("Purpose: Obtaining some simple tallies of structural variant characteristics");
		System.out.println("broken down by type, location effect (if annotated by svanno), and size range.");
		System.out.println("This tool can also split the input vcf into smaller vcfs by sv type and location effect.");
		System.out.println();
		System.out.println("Input Formats:");
		System.out.println("\tInput callset must be in [vcf] format");
		System.out.println();
		System.out.println("Output Formats:");
		System.out.println("\tTab delimited table [tsv]");
		System.out.println("\tVariant call format [vcf]");
		System.out.println();
		System.out.println("Flags:");
		System.out.println("\t-i\tFILE\t[Required]\t\tInput vcf path.");
		System.out.println("\t-o\tPATH\t[Optional]\t\tOutput tsv tables file prefix.");
		System.out.println("\t-s\tPATH\t[Optional]\t\tOutput split vcf file prefix.");
		System.out.println("\t-S\t    \t[Optional]\t\tHard split - if splitting vcf files, split not only by type, but by effect");
		//System.out.println("\t-t\tINT \t[Optional]\t\tNumber of threads. Defaults to 1.");
		System.out.println();
		System.out.println("Sample Usage:");
		System.out.println("java -jar bioisvtools.jar svqstats -g GRCh37 -v -i NA12878svcalls.vcf");
		System.out.println("java -jar bioisvtools.jar svqstats -g GRCh38 -i NA12878svcalls.vcf -o my/dir/svstattables/NA12878svcalls");
		System.out.println("java -jar bioisvtools.jar svqstats -g hg19 -i NA12878svcalls.vcf -o my/dir/svstattables_hg19/NA12878svcalls -s my/dir/NA12878svcalls_split");
		System.out.println("java -jar bioisvtools.jar svqstats -g hg38 -i NA12878svcalls.vcf -s my/dir/NA12878svcalls_split -S");
		System.out.println();
		System.out.println("--------------------------------------------------------------------------------");
	}
	
	public static void svQuickStats(String[] args, GenomeBuild g)
	{
		String inFile = null;
		String tblDir = null;
		String vcfDir = null;
		//int threads = 1;
		boolean hardsplit = false;
		
		for (int i = 0; i < args.length; i++)
		{
			String s = args[i];
			if (s.equals(OP_VCFIN))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_VCFIN + " flag MUST be followed by input VCF path!");
					printUsage();
					System.exit(1);
				}
				inFile = args[i+1];
			}
			else if (s.equals(OP_OUTPUT))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_OUTPUT + " flag MUST be followed by table output directory!");
					printUsage();
					System.exit(1);
				}
				tblDir = args[i+1];
			}
			else if (s.equals(OP_SPLIT))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_SPLIT + " flag MUST be followed by vcf output directory!");
					printUsage();
					System.exit(1);
				}
				vcfDir = args[i+1];
			}
			else if (s.equals(OP_HARDSPLIT))
			{
				hardsplit = true;
			}
		}
		
		if (inFile == null || inFile.isEmpty())
		{
			System.err.println("ERROR: Input path is required!");
			printUsage();
			System.exit(1);
		}
		
		/*if(tblDir != null)
		{
			int lastSlash = tblDir.lastIndexOf(File.separator);
			if (lastSlash == (tblDir.length() - 1)) tblDir = tblDir.substring(0, (tblDir.length() - 1));
		}
		
		if(vcfDir != null)
		{
			int lastSlash = vcfDir.lastIndexOf(File.separator);
			if (lastSlash == (vcfDir.length() - 1)) vcfDir = vcfDir.substring(0, (vcfDir.length() - 1));
		}*/
		
		//Open VCF
		VariantPool myVar = null;
		try 
		{
			myVar = VCF.readVCF(inFile, g, true);
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: " + inFile + " could not be opened!");
			e.printStackTrace();
			printUsage();
			System.exit(1);
		} 
		catch (UnsupportedFileTypeException e) 
		{
			System.err.println("ERROR: " + inFile + " could not be read!");
			e.printStackTrace();
			printUsage();
			System.exit(1);
		}
		
		//Count
		VarCounter counter = new VarCounter(myVar);
		
		//Print
		counter.printReport();
		
		//Print tables, if requested.
		if(tblDir != null && !tblDir.isEmpty()) counter.generateTables(tblDir);
		
		//Print split VCFs, if requested.
		if(vcfDir != null && !vcfDir.isEmpty()) counter.printVCFs(vcfDir, hardsplit);
		
	}
	
}
