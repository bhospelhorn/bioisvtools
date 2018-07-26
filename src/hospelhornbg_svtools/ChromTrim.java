package hospelhornbg_svtools;

import java.io.IOException;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;

import hospelhornbg_bioinformatics.VCF;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

//For taking out any SV calls associated with any non-autosomes (or sex chromosomes that shouldn't be present)

public class ChromTrim {
	
	public static final String OP_INPUT = "-i"; 
	public static final String OP_OUTPUT = "-o"; 
	public static final String OP_ALLFEMALE = "-f"; 
	
	public static void printUsage()
	{
		System.out.println("--------------------------------------------------------------------------------");
		System.out.println("BioisvTools || TrimChrom");
		System.out.println();
		System.out.println("Purpose: To remove any variants supposedly associated with any unmapped contigs, mitochrondrial");
		System.out.println("contigs, or optionally, the Y chromosome.");
		System.out.println();
		System.out.println("Input Formats:");
		System.out.println("\tInput callset must be in [vcf] format.");
		System.out.println();
		System.out.println("Output Formats:");
		System.out.println("\tPlain text Variant Call Format [vcf]");
		System.out.println();
		System.out.println("Flags:");
		System.out.println("\t-i\tFILE\t[Required]\t\tInput callset in the form of a vcf file.");
		System.out.println("\t-o\tFILE\t[Required]\t\tOutput file path");
		System.out.println("\t-f\tFLAG\t[Optional]\t\tAll samples in vcf are standard mammalian female (ie. remove all Y chromosome variants)");
		System.out.println();
		System.out.println("Sample Usage:");
		System.out.println("java -jar bioisvtools.jar trimchrom -g hg38 -i NA24385_svset.vcf -o NA24385_svset_noGLM.vcf");
		System.out.println("java -jar bioisvtools.jar trimchrom -g hg19 -i NA12878_svset.vcf -o NA12878_svset_noGLMY.vcf -f");
		System.out.println("java -jar bioisvtools.jar trimchrom -g grch37 -i NA12878_trio_svset.vcf -o NA12878_trio_svset_noGLM.vcf");
		System.out.println();
		System.out.println("--------------------------------------------------------------------------------");
	}
	
	public static void removeChromosomes_mammal(VariantPool pool, boolean isFemale, boolean verbose)
	{
		List<Variant> vcoll = pool.getVariants();
		
		List<Variant> passed = new LinkedList<Variant>();
		
		int toss_contig = 0;
		int toss_mito = 0;
		int toss_y = 0;
		int toss_unk = 0;
		
		for (Variant v : vcoll)
		{
			boolean toss = false;
			Collection<Contig> chroms = v.getAllChromosomes();
			for (Contig c : chroms)
			{
				if (c.getType() == Contig.SORTCLASS_CONTIG)
				{
					toss = true;
					toss_contig++;
					break;
				}
				else if (c.getType() == Contig.SORTCLASS_MITO)
				{
					toss = true;
					toss_mito++;
					break;
				}
				else if (isFemale && c.getUDPName().equals("Y"))
				{
					toss = true;
					toss_y++;
					break;
				}
				else if (c.getType() == Contig.SORTCLASS_UNKNOWN)
				{
					toss = true;
					toss_unk++;
					break;
				}
			}
			if(!toss) passed.add(v);
		}
		if (verbose)
		{
			System.out.println("ChromTrim Summary -----------");
			System.out.println("Initial Variants: " + vcoll.size());
			System.out.println("Variants Rejected for Unmapped Contig Association: " + toss_contig);
			System.out.println("Variants Rejected for Mitochondrial Association: " + toss_mito);
			System.out.println("Variants Rejected for Female Y Association: " + toss_y);
			System.out.println("Variants Rejected for Unknown Contig Association: " + toss_unk);
			System.out.println("Final Variants: " + passed.size());
		}
		
		pool.clearVariants();
		pool.addVariants(passed);
		
	}
	
	public static void runTool(String[] args, GenomeBuild g, boolean verbose)
	{
		String inPath = null;
		String outPath = null;
		boolean isfem = false;
		
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
			if (s.equals(OP_ALLFEMALE))
			{
				isfem = true;
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
		
		VariantPool pool = null;
		
		try 
		{
			pool = VCF.readVCF(inPath, g, false);
		} 
		catch (IOException e) 
		{
			System.err.println("ChromTrim.runTool || ERROR: Input VCF could not be opened!");
			e.printStackTrace();
			System.exit(1);
		} 
		catch (UnsupportedFileTypeException e) 
		{
			System.err.println("ChromTrim.runTool || ERROR: Input VCF could not be read!");
			e.printStackTrace();
			System.exit(1);
		}
		
		removeChromosomes_mammal(pool, isfem, verbose);
		
		try 
		{
			VCF.writeVCF(pool, "bioisvtools", outPath);
		} 
		catch (IOException e) {
			System.err.println("ChromTrim.runTool || ERROR: Output VCF could not be written!");
			e.printStackTrace();
		}
		
	}

}
