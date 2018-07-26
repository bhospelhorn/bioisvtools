package hospelhornbg_svtools;

import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.VCF;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class SizeFilter {
	
	public static final String OP_VCFIN = "-i"; 
	public static final String OP_VCFOUT = "-o"; 
	public static final String OP_MIN = "-m"; 
	public static final String OP_MAX = "-M"; 
	
	public static void printUsage()
	{
		System.out.println("--------------------------------------------------------------------------------");
		System.out.println("BioisvTools || sizefilter");
		System.out.println();
		System.out.println("Purpose: To filter out variants larger or smaller than the given size threshold(s).");
		System.out.println();
		System.out.println("Input Formats:");
		System.out.println("\tInput callset must be in [vcf] format");
		System.out.println();
		System.out.println("Output Formats:");
		System.out.println("\tVariant call format [vcf]");
		System.out.println();
		System.out.println("Flags:");
		System.out.println("\t-i\tFILE\t[Required]\t\tInput vcf path.");
		System.out.println("\t-o\tFILE\t[Required]\t\tOutput vcf path.");
		System.out.println("\t-m\tINT\t[Optional]\t\tMinimum size of variants to pass.");
		System.out.println("\t-M\tINT\t[Optional]\t\tMaximum size of variants to pass.");
		System.out.println("Note: You can omit both -m and -M, but this tool will just output the unfiltered set.");
		System.out.println();
		System.out.println("Sample Usage:");
		System.out.println("java -jar bioisvtools.jar sizefilter -g hg19 -i NA12878_svset.vcf -o NA12878_svset_large.vcf -m 1000000");
		System.out.println("java -jar bioisvtools.jar sizefilter -g grch37 -i NA12878_svset.vcf -o NA12878_svset_small.vcf -M 1000");
		System.out.println("java -jar bioisvtools.jar sizefilter -g hg38 -i NA12878_svset.vcf -o NA12878_svset_medium.vcf -m 1000000 -M 1000");
		System.out.println();
		System.out.println("--------------------------------------------------------------------------------");
	}
	
	public static void filterAllAbove(VariantPool pool, int maxSize)
	{
		List<Variant> vlist = pool.getVariants();
		pool.clearVariants();
		
		Set<Variant> keepset = new HashSet<Variant>();
		for (Variant v : vlist)
		{
			if (v instanceof StructuralVariant)
			{
				StructuralVariant sv = (StructuralVariant)v;
				if (Math.abs(sv.getSVLength()) <= maxSize) keepset.add(v);
			}
			else
			{
				if (v.getLargestAbsoluteLength() <= maxSize) keepset.add(v);
			}
		}
		
		pool.addVariants(keepset);
		pool.sortVariants();
		
	}
	
	public static void filterAllBelow(VariantPool pool, int minSize)
	{
		List<Variant> vlist = pool.getVariants();
		pool.clearVariants();
		
		Set<Variant> keepset = new HashSet<Variant>();
		for (Variant v : vlist)
		{
			if (v instanceof StructuralVariant)
			{
				StructuralVariant sv = (StructuralVariant)v;
				if (Math.abs(sv.getSVLength()) >= minSize) keepset.add(v);
			}
			else
			{
				if (v.getSmallestAbsoluteLength() >= minSize) keepset.add(v);
			}
		}
		
		pool.addVariants(keepset);
		pool.sortVariants();

	}

	public static void filterBySize(String[] args, GenomeBuild gb)
	{
		String inFile = null;
		String outFile = null;
		int minsize = -1;
		int maxsize = -1;
		
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
			else if (s.equals(OP_VCFOUT))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_VCFOUT + " flag MUST be followed by output VCF path!");
					printUsage();
					System.exit(1);
				}
				outFile = args[i+1];
			}
			else if (s.equals(OP_MIN))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_MIN + " flag MUST be followed by a valid positive integer!");
					printUsage();
					System.exit(1);
				}
				String minStr = args[i+1];
				try
				{
					minsize = Integer.parseInt(minStr);
				}
				catch(NumberFormatException e)
				{
					System.err.println("ERROR: " + OP_MIN + " flag MUST be followed by a valid positive integer!");
					printUsage();
					System.exit(1);
				}
				if (minsize < 0)
				{
					System.err.println("ERROR: " + OP_MIN + " flag MUST be followed by a valid positive integer!");
					printUsage();
					System.exit(1);
				}
			}
			else if (s.equals(OP_MAX))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_MAX + " flag MUST be followed by a valid positive integer!");
					printUsage();
					System.exit(1);
				}
				String maxStr = args[i+1];
				try
				{
					maxsize = Integer.parseInt(maxStr);
				}
				catch(NumberFormatException e)
				{
					System.err.println("ERROR: " + OP_MAX + " flag MUST be followed by a valid positive integer!");
					printUsage();
					System.exit(1);
				}
				if (maxsize < 0)
				{
					System.err.println("ERROR: " + OP_MAX + " flag MUST be followed by a valid positive integer!");
					printUsage();
					System.exit(1);
				}
			}
		}
		
		boolean pass = true;
		if (inFile == null || inFile.isEmpty())
		{
			System.err.println("ERROR: Input path is required!");
			pass = false;
		}
		if (outFile == null || outFile.isEmpty())
		{
			System.err.println("ERROR: Output path is required!");
			pass = false;
		}
		if (!pass)
		{
			printUsage();
			System.exit(1);
		}
		
		VariantPool mypool = null;
		try 
		{
			mypool = VCF.readVCF(inFile, gb, true);
		} 
		catch (IOException e) {
			System.err.println("IO ERROR: Input file " + inFile + " could not be opened!");
			e.printStackTrace();
			System.exit(1);
		} 
		catch (UnsupportedFileTypeException e) {
			System.err.println("PARSE ERROR: Input file " + inFile + " could not be read!");
			e.printStackTrace();
			System.exit(1);
		}
		
		if(minsize >= 0) filterAllBelow(mypool, minsize);
		if(maxsize >= 0) filterAllAbove(mypool, maxsize);
		
		try 
		{
			VCF.writeVCF(mypool, "bioisvtools", outFile);
		} 
		catch (IOException e) {
			System.err.println("IO ERROR: Output file " + outFile + " could not be written!");
			e.printStackTrace();
		}
		
	}
	
}
