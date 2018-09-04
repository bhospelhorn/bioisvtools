package hospelhornbg_svtools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import hospelhornbg_bioinformatics.BreakendPair;
import hospelhornbg_bioinformatics.Genotype;
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.Translocation;
import hospelhornbg_bioinformatics.VCF;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_genomeBuild.Gene;
import hospelhornbg_genomeBuild.GeneFunc;
import hospelhornbg_genomeBuild.GeneSet;
import hospelhornbg_genomeBuild.GenomeBuild;
import hospelhornbg_segregation.Candidate;
import hospelhornbg_segregation.Individual;
import hospelhornbg_segregation.Inheritance;
import hospelhornbg_segregation.Inheritor;
import hospelhornbg_segregation.Pedigree;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class SVFam {
	
	public static final double DEFO_MINMAP = 0.5;
	public static final int DEFO_MAXSIZE = 50000;

	public static final String OP_VCFIN = "-i"; 
	public static final String OP_OUTDIR = "-o"; 
	
	public static final String OP_PEDIN = "-p"; 
	
	public static final String OP_MAXSIZE = "-s"; 
	public static final String OP_MINMAP = "-m"; 
	
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
	
	public static void writeFullTable(List<Candidate> clist, Pedigree fam, String outpath) throws IOException
	{
		//int ncand = clist.size();
		String header = "VARID\tPOS\tTYPE\tSIZE\tPOSEFF\tGENE\tTID\tALLELE";
		List<Individual> famlist = fam.getAllMembers();
		List<Individual> affected = fam.getAllAffected();
		for (Individual aff : affected){
			header += "\t" + "SEG_" + aff.getName();
			header += "\t" + "SEGPARTNERS_" + aff.getName();
		}
		for (Individual i : famlist) header += "\t" + i.getName();
		
		FileWriter fw = new FileWriter(outpath);
		BufferedWriter bw = new BufferedWriter(fw);
		
		bw.write("##" + header + "\n");
		
		for (Candidate c : clist)
		{
			Variant v = c.getVariant();
			String rec = "";
			
			//VARID
			rec += v.getVarID() + "\t";
			
			//CHROM/POS
			//TYPE
			//SIZE
			//This field depends upon what v is...
			if (v instanceof StructuralVariant)
			{
				StructuralVariant sv = (StructuralVariant) v;
				if (sv instanceof BreakendPair)
				{
					rec += sv.getChromosome().getUDPName() + ":";
					rec += sv.getPosition() + "|";
					rec += sv.getEndChromosome().getUDPName() + ":";
					rec += sv.getEndPosition() + "\t";
				}
				else if (sv instanceof Translocation)
				{
					rec += ((Translocation) sv).getChromosome1().getUDPName() + ":";
					rec += sv.getPosition() + "|";
					rec += sv.getEndChromosome().getUDPName() + ":";
					rec += sv.getEndPosition() + "\t";
				}
				else
				{
					rec += sv.getChromosome().getUDPName() + ":";
					rec += sv.getPosition() + "-" + sv.getEndPosition() + "\t";
				}
				rec += sv.getType().getString() + "\t";
				rec += sv.getAbsoluteSVLength() + "bp\t";
			}
			else
			{
				rec += v.getChromosome().getUDPName() + ":" + v.getPosition() + "\t";
				rec += "[None]\t";
				rec += "1bp\t";
			}
			
			//POSEFF
			GeneFunc eff = c.getPositionEffect();
			if (eff == null) rec += "[UNKNOWN]\t";
			else rec += eff.toString() + "\t";
			
			//GENE, TID
			Gene g = c.getGene();
			if (g != null)
			{
				rec += g.getName() + "\t";
				rec += g.getID() + "\t";
			}
			else
			{
				rec += "[N/A]\t";
				rec += "[N/A]\t";	
			}
			
			//ALLELE
			rec += c.getAllele() + "\t";
			
			//SEG & SEGPARTNERS
			for (Individual aff : affected)
			{
				Inheritance ip = c.getInheritancePattern(aff);
				if (ip != null) rec += ip.toString() + "\t";
				else rec += "[N/A]\t";
				
				List<Candidate> plist = c.getAllPartners(aff);
				if (plist != null && !plist.isEmpty())
				{
					boolean first = true;
					for (Candidate p : plist)
					{
						if(!first) rec += ";";
						rec += p.getVariant().getVarID() + "," + c.getAllele();
						first = false;
					}	
				}
				else rec += "[None]";
				rec += "\t";
			}
			
			//GENOTYPES
			boolean first = true;
			for (Individual i : famlist)
			{
				Genotype gt = v.getSampleGenotype(i.getName());
				if (!first) rec += "\t";
				if (gt == null)
				{
					rec += "null";
				}
				else
				{
					rec += gt.getField(Genotype.INFODEF_GT.getKey());
				}
				first = false;
			}
			
			//Print
			bw.write(rec + "\n");
		}
		
		bw.close();
		fw.close();
	}
	
	public static void RunSVFamily(String[] args, GenomeBuild g)
	{
		String inPath = null;
		String outPath = null;
		String pedPath = null;
		int maxsize = DEFO_MAXSIZE;
		double minmap = DEFO_MINMAP;
		
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
			if (s.equals(OP_OUTDIR))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_OUTDIR + " flag MUST be followed by output directory path!");
					printUsage();
					System.exit(1);
				}
				outPath = args[i+1];
			}
			if (s.equals(OP_PEDIN))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_PEDIN + " flag MUST be followed by input PED path!");
					printUsage();
					System.exit(1);
				}
				pedPath = args[i+1];
			}
			if (s.equals(OP_MAXSIZE))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_MAXSIZE + " flag MUST be followed by integer!");
					printUsage();
					System.exit(1);
				}
				String rarg1 = args[i+1];
				try
				{
					maxsize = Integer.parseInt(rarg1);
				}
				catch (NumberFormatException e)
				{
					System.err.println("ERROR: " + OP_MAXSIZE + " flag MUST be followed by integer!");
					printUsage();
					System.exit(1);
				}
			}
			if (s.equals(OP_MINMAP))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_MINMAP + " flag MUST be followed by float!");
					printUsage();
					System.exit(1);
				}
				String rarg2 = args[i+1];
				try
				{
					minmap = Double.parseDouble(rarg2);
				}
				catch (NumberFormatException e)
				{
					System.err.println("ERROR: " + OP_MINMAP + " flag MUST be followed by float!");
					printUsage();
					System.exit(1);
				}
			}
		}
		
		
		//Check args
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
		if (pedPath == null || pedPath.isEmpty())
		{
			System.err.println("ERROR: Pedigree file path is required!");
			pass = false;
		}
		if (maxsize < 0)
		{
			System.err.println("ERROR: Max variant size is invalid!");
			pass = false;
		}
		if (minmap < 0.00)
		{
			System.err.println("ERROR: Min mapability value is invalid!");
			pass = false;
		}
		
		if (!pass)
		{
			printUsage();
			System.exit(1);
		}
		
		//Read files
		
		Pedigree fam = null;
		try 
		{
			fam = new Pedigree(pedPath);
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: PED file \"" + pedPath + "\" could not be opened!");
			e.printStackTrace();
			System.exit(1);
		} 
		catch (UnsupportedFileTypeException e) 
		{
			
			System.err.println("ERROR: PED file \"" + pedPath + "\" could not be read!");
			e.printStackTrace();
			System.exit(1);
		}
		
		VariantPool vcf = null;
		try 
		{
			vcf = VCF.readVCF(inPath, g, true);
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: VCF file \"" + inPath + "\" could not be opened!");
			e.printStackTrace();
			System.exit(1);
		} 
		catch (UnsupportedFileTypeException e) 
		{
			System.err.println("ERROR: VCF file \"" + inPath + "\" could not be read!");
			e.printStackTrace();
			System.exit(1);
		}
		
		//Load gene set
		GeneSet gs = GeneSet.loadRefGene(g);
		if (g == null)
		{
			System.err.println("ERROR: Gene set could not be opened!");
			System.exit(1);
		}
		System.err.println("Gene info loaded!");
		
		//Filter size - toss anything above max size
		List<Variant> vlist = vcf.getVariants();
		List<Variant> smallv = new LinkedList<Variant>();
		for (Variant v : vlist)
		{
			if (v.getSmallestAbsoluteLength() <= maxsize)
			{
				smallv.add(v);
			}
		}
		vcf.clearVariants();
		vcf.addVariants(smallv);
		
		//Filter mapability - not functional at the moment
		
		//Get candidates
		List<Candidate> clist = Inheritor.getCandidates(vcf, fam, gs);
		
		//Sort and write
		//Print full table
		String fulltablepath = outPath + File.separator + "fulltable.tsv";
		try 
		{
			writeFullTable(clist, fam, fulltablepath);
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: Table could not be written to " + fulltablepath);
			e.printStackTrace();
		}
		//Filter out ref alleles
		//Go by type
		
		
	}
	
}
