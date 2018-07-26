package hospelhornbg_svtools;

import java.io.IOException;
import java.util.Collection;
import java.util.List;

import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.VCF;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_bioinformatics.VariantPool.InfoDefinition;
import hospelhornbg_genomeBuild.GeneSet;
import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class SVAnno {
	
	public static final String OP_VCFIN = "-i"; 
	public static final String OP_VCFOUT = "-o"; 
	public static final String OP_THREADS = "-t"; 
	
	public static class RunnerThread extends Thread
	{
		private List<Variant> queue;
		private int counter;
		private GeneSet genes;
		
		public RunnerThread(List<Variant> vlist, GeneSet gs, int number)
		{
			queue = vlist;
			counter = 0;
			genes = gs;
			this.setName("SVANNO_workerThread_" + number);
			this.setDaemon(true);
		}
		
		private synchronized void incrementCounter()
		{
			counter++;
		}
		
		public void run()
		{
			for (Variant v : queue)
			{
				if (v instanceof StructuralVariant)
				{
					StructuralVariant sv = (StructuralVariant)v;
					genes.annotateStructuralVariant(sv, false);
					//incrementCounter();
				}
				else
				{
					System.err.println("Variant " + v.getVarID() + " was not recognized as a structural variant. Skipping...");
				}
				incrementCounter();
			}
		}
		
		public synchronized int countProcessed()
		{
			return counter;
		}
		
	}
	
	public static void printUsage()
	{
		System.out.println("--------------------------------------------------------------------------------");
		System.out.println("BioisvTools || svanno");
		System.out.println();
		System.out.println("Purpose: For annotating structural variants with refGene data...");
		System.out.println("Note which (if any) genes are affected by each variant, and what region of the genes are affected.");
		System.out.println("Flanking gene information included for intergenic variants.");
		System.out.println();
		System.out.println("Input Formats:");
		System.out.println("\tInput callset must be in [vcf] format");
		System.out.println();
		System.out.println("Output Formats:");
		System.out.println("\tPlain text Variant Call Format [vcf]");
		System.out.println();
		System.out.println("Flags:");
		System.out.println("\t-i\tFILE\t[Required]\t\tInput vcf path.");
		System.out.println("\t-o\tFILE\t[Required]\t\tOutput vcf path.");
		System.out.println("\t-t\tINT\t[Optional]\t\tNumber of threads. Defaults to 1.");
		System.out.println();
		System.out.println("Sample Usage:");
		System.out.println("java -jar bioisvtools.jar svanno -g GRCh37 -v -i NA12878_svset.vcf -o NA12878_svset_refGene.vcf");
		System.out.println("java -jar bioisvtools.jar svanno -g hg38 -i NA12878_svset.vcf -o NA12878_svset_refGene.vcf -t 24");
		System.out.println();
		System.out.println("--------------------------------------------------------------------------------");
	}
	
	public static boolean checkForAlive(Thread[] threads)
	{
		if (threads == null) return false;
		if (threads.length < 1) return false;
		for (int i = 0; i < threads.length; i++)
		{
			if (threads[i] != null)
			{
				if (threads[i].isAlive()) return true;
			}
		}
		return false;
	}
	
	public static void annotatePool(VariantPool pool, GeneSet genes, boolean verbose)
	{
		//Annotate variants
		List<Variant> vlist = pool.getVariants();
		int counter = 0;
		for (Variant v : vlist)
		{
			if (v instanceof StructuralVariant)
			{
				StructuralVariant sv = (StructuralVariant)v;
				genes.annotateStructuralVariant(sv, false);
				counter++;
				if (verbose && counter % 1000 == 0) System.err.println(counter + " variants processed!");
			}
			else
			{
				System.err.println("Variant " + v.getVarID() + " was not recognized as a structural variant. Skipping...");
			}
		}
		
		//Add infodefs
		Collection<InfoDefinition> icoll = GeneSet.getInfoDefinitions();
		for (InfoDefinition def : icoll) pool.addInfoField(def.getKey(), def);
		
	}
	
	public static void annotatePool(VariantPool pool, GeneSet genes, int threads, boolean verbose)
	{
		//Get variants
		List<Variant> vlist = pool.getVariants();
		int vCount = vlist.size();
		if (verbose) System.err.println(vCount + " variants retrieved! Now distributing to worker threads...");
		
		//Distribute variants
		RunnerThread[] tarr = new RunnerThread[threads - 1];
		int nvar = (vCount / tarr.length) + 1; //Ceiling
		int ed = 0;
		for (int i = 0; i < tarr.length; i++)
		{
			int st = ed;
			ed = (nvar * (i + 1));
			if (ed > vCount) ed = vCount;
			List<Variant> sublist = vlist.subList(st, ed);
			tarr[i] = new RunnerThread(sublist, genes, i);
			tarr[i].start();
		}
		
		//Main thread monitors progress!
		while(checkForAlive(tarr))
		{
			try 
			{
				Thread.sleep(100);
			} 
			catch (InterruptedException e) 
			{
				System.err.println("Main thread sleep interrupted unexpectedly! Going back to sleep...");
				e.printStackTrace();
			}
			if(verbose)
			{
				int tot = 0;
				for (int i = 0; i < tarr.length; i++) tot += tarr[i].countProcessed();
				double percent = ((double)tot / (double)vCount) * 100.0;
				System.err.println(tot + " of " + vCount + " variants processed (" + String.format("%.2f", percent) + "%)...");
			}
		}
		if (verbose) System.err.println("All structural variants have been processed!");
		
		//Add infodefs
		Collection<InfoDefinition> icoll = GeneSet.getInfoDefinitions();
		for (InfoDefinition def : icoll) pool.addInfoField(def.getKey(), def);
		if (verbose) System.err.println("InfoDefinitions added to VCF...");
		
		//Variants are all references, so everything in pool should be up to date!! ...?
		
	}
	
	public static void svanno(String[] args, GenomeBuild genome, GeneSet mygenes, boolean verbose)
	{
		String inPath = null;
		String outPath = null;
		int threads = 1;
		
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
				inPath = args[i+1];
			}
			else if (s.equals(OP_VCFOUT))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_VCFOUT + " flag MUST be followed by output VCF path!");
					printUsage();
					System.exit(1);
				}
				outPath = args[i+1];
			}
			else if (s.equals(OP_THREADS))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_THREADS + " flag MUST be followed by an integer!");
					printUsage();
					System.exit(1);
				}
				String t = args[i+1];
				try{
					threads = Integer.parseInt(t);
				}
				catch (NumberFormatException e)
				{
					System.err.println("ERROR: " + OP_THREADS + " flag MUST be followed by an integer!");
					printUsage();
					System.exit(1);
				}
			}
		}
		
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
		
		if (!pass)
		{
			printUsage();
			System.exit(1);
		}
		
		if (threads < 1) threads = 1;
		
		VariantPool pool = null;
		try 
		{
			pool = VCF.readVCF(inPath, genome, true);
		} 
		catch (IOException e) 
		{
			System.err.println("IO ERROR: Input VCF " + inPath + " could not be opened!");
			e.printStackTrace();
			System.exit(1);
		} 
		catch (UnsupportedFileTypeException e) 
		{
			System.err.println("PARSE ERROR: Input VCF " + inPath + " could not be read!");
			e.printStackTrace();
			System.exit(1);
		}
		
		if (pool == null)
		{
			System.err.println("UNKNOWN ERROR: Input VCF " + inPath + " was not read!");
			System.exit(1);
		}
		
		if (threads == 1) annotatePool(pool, mygenes, verbose);
		else annotatePool(pool, mygenes, threads, verbose);
		
		try 
		{
			VCF.writeVCF(pool, "bioisvtools", outPath);
		} 
		catch (IOException e) 
		{
			System.err.println("IO ERROR: Output VCF " + outPath + " could not be written!");
			e.printStackTrace();
		}
		
	}

}
