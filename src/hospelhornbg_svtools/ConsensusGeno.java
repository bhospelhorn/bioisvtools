package hospelhornbg_svtools;

import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import hospelhornbg_bioinformatics.Genotype;
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.VCF;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_bioinformatics.VariantPool.InfoDefinition;
import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class ConsensusGeno {
	
	public static final String OP_INPUT = "-i";
	public static final String OP_OUTPUT = "-o";
	public static final String OP_SAMPLE = "-s";
	//public static final String OP_BEDPATH = "-b";
	public static final String OP_NEW = "-n";
	public static final String OP_SUPPKEY = "-k";
	
	//NOTE: DR - Number of supporting reads (REF,ALT).
	
	public static final String SUPPDES = "Number of sample VCFs containing evidence for SURVIVOR call.";
	public static final String SUPPDES_VEC = "Vector of sample VCFs containing evidence for SURVIVOR call.";
	
	public static void printUsage()
	{
		System.out.println("--------------------------------------------------------------------------------");
		System.out.println("BioisvTools || SURVIVOR Re-Genotyper");
		System.out.println();
		System.out.println("Purpose: For \"cleaning\" SURVIVOR output VCFs so that genotype calls are");
		System.out.println("chosen by consensus and placed into a single column for the specified sample.");
		System.out.println();
		System.out.println("Input Formats:");
		System.out.println("\tInput callset must be in [vcf] format");
		System.out.println();
		System.out.println("Output Formats:");
		System.out.println("\tPlain text Variant Call Format [vcf]");
		System.out.println();
		System.out.println("Flags:");
		System.out.println("\t-i\tFILE\t[Required]\t\tInput file path");
		System.out.println("\t-o\tFILE\t[Required]\t\tOutput file path");
		System.out.println("\t-s\tFILE\t[Required]\t\tSample name");
		System.out.println("\t-k\tFILE\t[Optional]\t\tINFO key to move SURVIVOR SUPP/SUPP_VEC fields to. These will simply be deleted if this option isn't given.");
		System.out.println("\t-n\t\t[Required]\t\tRefresh the file (\"new file\") ");
		System.out.println();
		System.out.println("Sample Usage:");
		System.out.println("java -jar bioisvtools.jar survivorgeno -g GRCh37 -i mergedset.vcf -o mergedset_consensus.vcf -s NA12878 -n -k MERGESUPP");
		System.out.println();
		System.out.println("--------------------------------------------------------------------------------");
	}
	
	private static class Genocall
	{
		private String gt;
		private int votes;
		
		private int leastR;
		private int mostA;
		
		public Genocall(String genotype)
		{
			gt = genotype;
			votes = 0;
			leastR = Integer.MAX_VALUE;
			mostA = 0;
		}
		
		public void addVote(String DR)
		{
			try
			{
				String fields[] = DR.split(",");
				int R = Integer.parseInt(fields[0]);
				int A = Integer.parseInt(fields[1]);
				int diff = A-R;
				int bestdiff = mostA - leastR;
				if (diff > bestdiff)
				{
					leastR = R;
					mostA = A;
				}
			}
			catch (ArrayIndexOutOfBoundsException | NumberFormatException e)
			{
				//Eat
				System.err.println("Genocall.addVote | Parsing error: Input = " + DR);
			}
			votes++;
		}
		
		public int countVotes()
		{
			return votes;
		}
		
		public String getGenotypeString()
		{
			return gt;
		}
		
		public int getARDiff()
		{
			return mostA-leastR;
		}
		
		public int getARTotal()
		{
			return mostA+leastR;
		}
		
		public String getDRField()
		{
			return leastR + "," + mostA;
		}
		
		
	}
	
	public static void condenseGenotype(VariantPool pool, String samplename)
	{
		System.err.println("survivorgeno condenseGenotype || Function called...");
		System.err.println("\tsamplename = " + samplename);
		List<String> sIndexes = pool.getAllSamples();
		Map<String, Genocall> votemap = new HashMap<String, Genocall>();
		List<Variant> varList = pool.getVariants();
		
		//int ecount = 0;
		//int counter = 0;
		//boolean debuginfo = false;
		for (Variant v : varList)
		{
			//debuginfo = (counter % 500 == 0);
			//if (debuginfo) System.err.println("survivorgeno condenseGenotype || DEBUG Annotation info for variant " + (counter+1));
			String suppvec = v.getSingleInfoEntry(StructuralVariant.INFODEF_INFO_SUPP_VEC.getKey());
			int unkgeno = 0;
			//if (debuginfo) System.err.println("\tsuppvec = " + suppvec);
			for (int i = 0; i < sIndexes.size(); i++)
			{
				//if (debuginfo) System.err.println("\tChecking support for sample " + i + "(" + sIndexes.get(i) + ")");
				if (suppvec.length() <= i) continue; 
				//if (debuginfo) System.err.println("\tSupport bool: " + suppvec.charAt(i));
				if (suppvec.charAt(i) == '1')
				{
					//Get the gt string
					//if (debuginfo) System.err.println("\tSupport found from this sample.");
					Genotype g = v.getSampleGenotype(sIndexes.get(i));
					if (g.getPercentAlt() < 0.01) {
						if (g.isGenotypeUnknown()) unkgeno++;
						continue; //homref
					}
					//if (debuginfo) System.err.println("\tSample genotype call is not homref. (Percent Alt = " + String.format("%.3f", g.getPercentAlt()));
					String gt = g.getField(Genotype.INFODEF_GT.getKey());
					//if (debuginfo) System.err.println("\tGT = " + gt);
					String dr = g.getField("DR");
					//if (debuginfo) System.err.println("\tDR = " + dr);
					Genocall gc = votemap.get(gt);
					if (gc == null) gc = new Genocall(gt);
					gc.addVote(dr);
					//if (debuginfo) System.err.println("\tVote added.");
					votemap.put(gt, gc);
					//if (debuginfo) System.err.println("\tVote added to map.");
				}
				else
				{
					//This sample did not support the call. Replace genotype with 0/0
					//Genotype string will be ./. - We want to replace this with 0/0
					Genotype g = v.getSampleGenotype(sIndexes.get(i));
					final int[] refref = {0,0};
					g.setAlleles(refref);
					//System.err.println("DEBUG survivorgeno.condenseGenotype || Genotype set to: " + g.getField("GT"));
				}
			}
			//Tally votes
			Set<String> allcalls = votemap.keySet();
			Genocall best = null;
			if (allcalls.isEmpty() && unkgeno > 0)
			{
				best = new Genocall("./.");
				best.votes = 0;
				best.leastR = 0;
				best.mostA = 0;
			}
			for (String k : allcalls)
			{
				//if (debuginfo) System.err.println("\tVote Key: " + k);
				Genocall c = votemap.get(k);
				if (best == null) {
					//if (debuginfo) System.err.println("\tFirst genotype... Setting as best.");
					best = c;
				}
				else if (c.countVotes() > best.countVotes()) {
					//if (debuginfo) System.err.println("\tCurrent has more votes. Setting as best...");
					//if (debuginfo) System.err.println("\t\tOld Best: " + best.countVotes());
					//if (debuginfo) System.err.println("\t\tNew Best: " + c.countVotes());
					best = c;
				}
				else if (c.countVotes() == best.countVotes())
				{
					//if (debuginfo) System.err.println("\tSame number of votes. Comparing supporting read counts...");
					int diff1 = best.getARDiff();
					int diff2 = c.getARDiff();
					//if (debuginfo) System.err.println("\t\tBest Diff: " + diff1);
					//if (debuginfo) System.err.println("\t\tThis Diff: " + diff2);
					if (diff1 == diff2)
					{
						//if (debuginfo) System.err.println("\tSame supporting read split. Comparing totals...");
						int tot1 = best.getARTotal();
						int tot2 = c.getARTotal();
						//if (debuginfo) System.err.println("\t\tBest Total: " + tot1);
						//if (debuginfo) System.err.println("\t\tThis Total: " + tot2);
						if (tot2 > tot1) best = c;
					}
					else if (diff2 > diff1) best = c;
				}
				
			}
			//Clear map
			votemap.clear();
			//if (debuginfo) System.err.println("\tVote map cleared!");
			//Save genotype
			if (best != null)
			{
				//if (debuginfo) System.err.println("\tBest call pointer is non-null.");
				Genotype ng = new Genotype();
				ng.setAlleles(best.getGenotypeString());
				//if (debuginfo) System.err.println("\tGenotype set: " + ng.getField("GT"));
				ng.setField("DR", best.getDRField());
				//if (debuginfo) System.err.println("\tDR set: " + ng.getField("DR"));
				v.addGenotype(samplename, ng);
				//if (debuginfo) System.err.println("\tGenotype added? Genotype retrieval is null: " + (v.getSampleGenotype(samplename) == null));
			}
			//counter++;
		}
		
		pool.addSample(samplename);
		
	}
	
	public static void cleanPool(VariantPool pool, boolean clearOldSamples, String samplename, String suppkey)
	{
		//Remove all FORMAT fields that aren't used in consensus call
		//System.err.println("survivorgeno cleanPool || Function called...");
		//System.err.println("\tclearOldSamples = " + clearOldSamples);
		//System.err.println("\tsamplename = " + samplename);
		//System.err.println("\tsuppkey = " + suppkey);
		pool.clearActiveFormatKeys();
		//System.err.println("survivorgeno cleanPool || Format keys cleared.");
		pool.addFormatKeyToActiveList(Genotype.INFODEF_GT.getKey());
		//System.err.println("survivorgeno cleanPool || Added GT key back...");
		pool.addFormatField("DR", new InfoDefinition("DR", VariantPool.INFODEF_INT, "# supporting reference,variant reads in that order", 2));
		//pool.addFormatKeyToActiveList("DR");
		//System.err.println("survivorgeno cleanPool || Added DR key back...");
		//List<String> gkeys = pool.getOrderedGenotypeKeys();
		//System.err.println("survivorgeno cleanPool || DEBUG: Checking geno FORMAT keys...");
		//for (String k : gkeys) System.err.println("\t" + k);
		
		//Remove all INFO fields that could conflict with later SURVIVOR runs
		pool.removeInfoKeyFromActiveList(StructuralVariant.INFODEF_INFO_SUPP.getKey());
		pool.removeInfoKeyFromActiveList(StructuralVariant.INFODEF_INFO_SUPP_VEC.getKey());
		pool.removeInfoKeyFromActiveList(StructuralVariant.INFODEF_INFO_AVGLEN.getKey());
		pool.removeInfoKeyFromActiveList(StructuralVariant.INFODEF_INFO_SVMETHOD.getKey());
		//System.err.println("survivorgeno cleanPool || Certain SURVIVOR INFO keys removed...");
		
		InfoDefinition newsupp = null;
		InfoDefinition newsuppvec = null;
		if (suppkey != null && !suppkey.isEmpty())
		{
			newsupp = new InfoDefinition(suppkey, VariantPool.INFODEF_INT, SUPPDES, 1);
			newsuppvec = new InfoDefinition(suppkey + "_VEC", VariantPool.INFODEF_STRING, SUPPDES_VEC, 1);
			
			pool.addInfoField(newsupp.getKey(), newsupp);
			pool.addInfoField(newsuppvec.getKey(), newsuppvec);	
		}
		//System.err.println("survivorgeno cleanPool || New SUPP vector infodef added...");
		
		//Remove all variants with no consensus genotype
		//Transfer support information to a different INFO key (if specified)
		System.err.println("survivorgeno cleanPool || Getting variants from pool. Variant Count: " + pool.countVariants());
		List<Variant> varList = pool.getVariants();
		pool.clearVariants();
		System.err.println("survivorgeno cleanPool || Variant list retrieved. Pool Variant Count: " + pool.countVariants() + " | Local Variant Count: " + varList.size());
		List<Variant> passed = new LinkedList<Variant>();
		int rej_nosamplename = 0;
		for (Variant v : varList)
		{
			if (v.getSampleGenotype(samplename) == null) {
				rej_nosamplename++;
				continue;
			}
			if (suppkey != null && !suppkey.isEmpty())
			{
				String supp = v.getSingleInfoEntry(StructuralVariant.INFODEF_INFO_SUPP.getKey());
				String vec = v.getSingleInfoEntry(StructuralVariant.INFODEF_INFO_SUPP_VEC.getKey());
				if (supp != null) v.addInfoField(newsupp.getKey(), supp);
				if (vec != null) v.addInfoField(newsuppvec.getKey(), vec);	
			}
			passed.add(v);
		}
		System.err.println("survivorgeno cleanPool || Variants Scanned... Passed Variant Count: " + passed.size());
		System.err.println("survivorgeno cleanPool || Variants Rejected (No sample record): " + rej_nosamplename);
		pool.addVariants(passed);
		
		System.err.println("survivorgeno cleanPool || Variants added back to pool: " + pool.countVariants());
		//Clear out original sample genotypes if requested
		if (clearOldSamples)
		{
			pool.clearSamples();
			pool.addSample(samplename);
		}
		
		pool.sortVariants();
		//System.err.println("survivorgeno cleanPool || Function returning...");
		
	}
	
	public static void scgeno(String[] args, GenomeBuild g)
	{
		String inPath = null;
		String outPath = null;
		boolean makenew = false;
		String samplename = null;
		String suppkey = null;
		
		for (int i = 0; i < args.length; i++)
		{
			String s = args[i];
			if (s.equals(OP_INPUT))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_INPUT + " flag MUST be followed by input VCF path!");
					printUsage();
					System.exit(1);
				}
				inPath = args[i+1];
			}
			if (s.equals(OP_OUTPUT))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_OUTPUT + " flag MUST be followed by output VCF path!");
					printUsage();
					System.exit(1);
				}
				outPath = args[i+1];
			}
			if (s.equals(OP_SAMPLE))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_SAMPLE + " flag MUST be followed by sample name!");
					printUsage();
					System.exit(1);
				}
				samplename = args[i+1];
			}
			if (s.equals(OP_SUPPKEY))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_SUPPKEY + " flag MUST be followed by INFO key stem to transfer SURVIVOR information to!");
					printUsage();
					System.exit(1);
				}
				suppkey = args[i+1];
			}
			if (s.equals(OP_NEW))
			{
				makenew = true;
			}
			
		}
		
		//Check required arguments
		boolean pass = true;
		if (inPath == null || inPath.isEmpty())
		{
			System.err.println("ERROR: Input path is required!");
			pass = false;
		}
		if (outPath == null || outPath.isEmpty())
		{
			System.err.println("ERROR: Output path is required!");
			pass = false;
		}
		if (samplename == null || samplename.isEmpty())
		{
			System.err.println("ERROR: Sample name is required!");
			pass = false;
		}
		
		if (!pass)
		{
			printUsage();
			System.exit(1);
		}
		if (g == null)
		{
			System.err.println("ERROR: Reference genome build is required!");
			printUsage();
			System.exit(1);
		}
		System.err.println("survivorgeno || Genome Selected: " + g.getBuildName());
		
		System.err.println("survivorgeno || Argument Summary: ");
		System.err.println("\tInput: " + inPath);
		System.err.println("\tOutput: " + outPath);
		System.err.println("\tSample Name: " + samplename);
		System.err.println("\tSuppkey: " + suppkey);
		System.err.println("\tMakenew: " + makenew);
		
		//Import VCF
		
		VariantPool vars = null;
		try 
		{
			VCF reader = new VCF(inPath, true, g);
			vars = reader.getVariants();
		} 
		catch (UnsupportedFileTypeException e) 
		{
			System.err.println("ERROR: Input VCF " + inPath + " could not be read!");
			e.printStackTrace();
			System.exit(1);
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: Input VCF " + inPath + " could not be opened!");
			e.printStackTrace();
			System.exit(1);
		}
		
		System.err.println("survivorgeno || VCF Read! ");
		
		//Create sample index array (for quicker checking against SUPP_VEC)
		//Run genotype evaluator
		condenseGenotype(vars, samplename);
		System.err.println("survivorgeno || Genotype condensed! Variant count: " + vars.getVariants().size());
		cleanPool(vars, makenew, samplename, suppkey);
		System.err.println("survivorgeno || Pool Cleaned! Variant count: " + vars.getVariants().size());
		
		//Handle write - if n flag is used, then write only consensus to new file (taking out certain SURVIVOR fields)
		//If not n flag, simply add the new sample column to the end of the new file.
		try 
		{
			VCF.writeVCF(vars, "BioisvTools", outPath);
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: Output " + outPath + " could not be written!");
			e.printStackTrace();
			System.exit(1);
		}
		System.err.println("survivorgeno || VCF Written!");
		
		
	}

}
