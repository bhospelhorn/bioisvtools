package hospelhornbg_svtools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import hospelhornbg_bioinformatics.BreakendPair;
import hospelhornbg_bioinformatics.Genotype;
import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.Translocation;
import hospelhornbg_bioinformatics.UCSCGVBED;
import hospelhornbg_bioinformatics.VCF;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_bioinformatics.VariantPool.InfoDefinition;
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
		System.out.println("BioisvTools || svreport");
		System.out.println();
		System.out.println("Purpose: For analyzing a family SV callset VCF, annotating it, and reprinting information into");
		System.out.println("new, separated files so that the data can be more easily analyzed in pieces.");
		System.out.println();
		System.out.println("Input Formats:");
		System.out.println("\tInput callset must be in [vcf] format");
		System.out.println("\tInput family pedigree must be in [ped] format");
		System.out.println();
		System.out.println("Output Formats:");
		System.out.println("\tVariant Call Format [vcf]");
		System.out.println("\tBrowser Extensible Data [bed]");
		System.out.println("\tTab-Separated Value table [tsv]");
		System.out.println("\tPlain text [txt]");
		System.out.println();
		System.out.println("Flags:");
		System.out.println("\t-i\tFILE\t[Required]\t\tInput vcf path.");
		System.out.println("\t-o\tFILE\t[Required]\t\tOutput directory path.");
		System.out.println("\t-p\tFILE\t[Required]\t\tInput PED (pedigree) file path.");
		System.out.println("\t-s\tINT\t[Optional]\t\tMaximum size of CNVs kept in analysis. [Default: 50000bp]");
		//System.out.println("\t-m\tFLOAT\t[Optional]\t\tMinimum mappability score of regions where SV endpoints fall for SV to be kept. [Default: 0.5]");
		System.out.println();
		System.out.println("Sample Usage:");
		System.out.println("java -jar bioisvtools.jar svreport -g GRCh37 -v -i NA12878_svset_fullfam.vcf -o ./NA12878/SVanalysis -p NA12878.ped");
		System.out.println("java -jar bioisvtools.jar svreport -g hg38 -i NA12878_svset_fullfam.vcf -o ./NA12878/SVanalysis -p NA12878.ped -s 1000000");
		System.out.println();
		System.out.println("--------------------------------------------------------------------------------");
	}
	
	public static void writeTable(List<Candidate> clist, Pedigree fam, String outpath) throws IOException
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
	
	public static void writeGeneList(Collection<Gene> genes, String outpath)
	{
		List<Gene> sorted = new LinkedList<Gene>();
		sorted.addAll(genes);
		Collections.sort(sorted);
		
		FileWriter fw = null;
		BufferedWriter bw = null;
		
		try
		{
			fw = new FileWriter(outpath);
			bw = new BufferedWriter(fw);
			
			for (Gene g : sorted)
			{
				bw.write(g.getName() + "\t" + g.getID() + "\n");
			}
			
			bw.close();
			fw.close();	
		}
		catch (IOException e)
		{
			System.err.println("ERROR: Gene list could not be written to " + outpath);
			e.printStackTrace();
		}
	}
	
	public static void writeVCF(List<Candidate> clist, Pedigree fam, String outpath, GenomeBuild gb)
	{
		List<Individual> famlist = fam.getAllMembers();
		VariantPool pool = new VariantPool(famlist.size());
		pool.setGenomeBuild(gb);
		//Add keys
		StructuralVariant.addStandardDefs(pool, true);
		Collection<InfoDefinition> geneinfos = GeneSet.getInfoDefinitions();
		for (InfoDefinition infodef : geneinfos) pool.addInfoField(infodef.getKey(), infodef);
		//Add samples
		for (Individual i : famlist) pool.addSample(i.getName());
		//Put together variant list
		Set<Variant> vset = new HashSet<Variant>();
		for (Candidate c : clist)
		{
			vset.add(c.getVariant());
			vset.addAll(c.getAllPartnerVariants());
		}
		pool.addVariants(vset);
		pool.sortVariants();
		try 
		{
			VCF.writeVCF(pool, "bioisvtools", outpath);
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: VCF could not be written to " + outpath);
			e.printStackTrace();
		}
	}
	
	public static void writeUCSCBED(List<Candidate> clist, String outpath, String trackname, String trackdesc, String pbname)
	{
		Set<Variant> vset = new HashSet<Variant>();
		for (Candidate c : clist)
		{
			vset.add(c.getVariant());
			vset.addAll(c.getAllPartnerVariants());
		}
		List<Variant> vlist = new LinkedList<Variant>();
		vlist.addAll(vset);
		UCSCGVBED bed = new UCSCGVBED(trackname, vlist);
		bed.setDescription(trackdesc);
		try 
		{
			bed.write(outpath, pbname);
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: BED could not be written to " + outpath);
			e.printStackTrace();
		}
		
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
			if (v instanceof StructuralVariant)
			{
				StructuralVariant sv = (StructuralVariant)v;
				if (sv.getType() == SVType.INV) continue;
				if (sv.getType() == SVType.TRA) continue;
				if (sv.getType() == SVType.BND) continue;
				if (sv.getType() == SVType.INS) continue;
				if (sv.getType() == SVType.INSME) continue;
			}
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
			writeTable(clist, fam, fulltablepath);
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: Table could not be written to " + fulltablepath);
			e.printStackTrace();
		}
		//Filter out ref alleles (They will only appear as halfhet partners)
		List<Candidate> altlist = new LinkedList<Candidate>();
		for (Candidate c : clist)
		{
			if (c.getAllele() != 0) altlist.add(c);
		}
		
		//Write by type
		SVType[] types = SVType.values();
		for (SVType t : types)
		{
			String tstr = t.getString();
			
			//Make directory
			String tdir = outPath + File.separator + tstr;
			try 
			{
				Files.createDirectories(Paths.get(tdir));
			}
			catch (IOException e) 
			{
				System.err.println("ERROR: Could not create directory " + tdir);
				e.printStackTrace();
				continue;
			}
			
			//Get subset
			List<Candidate> tlist = new LinkedList<Candidate>();
			for (Candidate c : altlist)
			{
				Variant v = c.getVariant();
				if (v instanceof StructuralVariant)
				{
					StructuralVariant sv = (StructuralVariant)v;
					if (sv.getType() == t) tlist.add(c);
				}
			}
			
			//Sort by effect
			GeneFunc[] effects = GeneFunc.values();
			for (GeneFunc e : effects)
			{
				String estr = effects.toString();
				List<Candidate> elist = new LinkedList<Candidate>();
				for (Candidate c : tlist)
				{
					if (c.getPositionEffect() == e) elist.add(c);
				}
				
				String pathbase = tdir + File.separator + tstr + "_" + estr;
				
				//Separate out DNS and UPHH candidates
				List<Candidate> finalc = new LinkedList<Candidate>();
				List<Candidate> uphh = new LinkedList<Candidate>();
				List<Candidate> dns = new LinkedList<Candidate>();
				for (Candidate c : elist)
				{
					if (!c.segregatesWithAny())
					{
						dns.add(c);
						continue;
					}
					if (c.isUnpairedHalfHet())
					{
						uphh.add(c);
						continue;
					}
					finalc.add(c);
				}
				
				//Write tsv table
				String tsvpath_std = pathbase + "_table.tsv";
				String tsvpath_uphh = pathbase + "_table_uphh.tsv";
				String tsvpath_dns = pathbase + "_table_dns.tsv";
				
				try 
				{
					writeTable(finalc, fam, tsvpath_std);
				} 
				catch (IOException ex) 
				{
					System.err.println("ERROR: Table could not be written to " + tsvpath_std);
					ex.printStackTrace();
				}
				try 
				{
					writeTable(uphh, fam, fulltablepath);
				} 
				catch (IOException ex) 
				{
					System.err.println("ERROR: Table could not be written to " + tsvpath_uphh);
					ex.printStackTrace();
				}
				try 
				{
					writeTable(dns, fam, tsvpath_dns);
				} 
				catch (IOException ex) 
				{
					System.err.println("ERROR: Table could not be written to " + tsvpath_dns);
					ex.printStackTrace();
				}
				
				
				//Write gene lists
				Set<Gene> genes_std = new HashSet<Gene>();
				Set<Gene> genes_uphh = new HashSet<Gene>();
				Set<Gene> genes_dns = new HashSet<Gene>();
				
				String glpath_std = pathbase + "_genes.txt";
				String glpath_uphh = pathbase + "_genes_uphh.txt";
				String glpath_dns = pathbase + "_genes_dns.txt";
				
				for (Candidate c : finalc) genes_std.add(c.getGene());
				for (Candidate c : uphh) genes_uphh.add(c.getGene());
				for (Candidate c : dns) genes_dns.add(c.getGene());
				
				writeGeneList(genes_std, glpath_std);
				writeGeneList(genes_uphh, glpath_uphh);
				writeGeneList(genes_dns, glpath_dns);
				
				//Write vcf
				String vcfpath_std = pathbase + ".vcf";
				String vcfpath_uphh = pathbase + "_uphh.vcf";
				String vcfpath_dns = pathbase + "_dns.vcf";
				
				writeVCF(finalc, fam, vcfpath_std, g);
				writeVCF(uphh, fam, vcfpath_uphh, g);
				writeVCF(dns, fam, vcfpath_dns, g);
				
				//Write BED track
				String bedsuf_std = ".bed";
				String bedsuf_uphh = "_uphh.bed";
				String bedsuf_dns = "_dns.bed";
				
				List<Individual> aff = fam.getAllAffected();
				for (Individual i : aff)
				{
					String tname_base = i.getName() + " " + tstr + " " + estr;
					String tdesc_base = "Structural Variant Calls for " + i.getName() + " (" + estr + " " + tstr + ") MaxSize = " + maxsize + "bp";
					writeUCSCBED(finalc, pathbase + "_" + i.getName() + "_" + bedsuf_std, tname_base, tdesc_base, i.getName());
					writeUCSCBED(uphh, pathbase + "_" + i.getName() + "_" + bedsuf_uphh, tname_base + " [UPHH]", tdesc_base + " [Unpaired HalfHets]", i.getName());
					writeUCSCBED(dns, pathbase + "_" + i.getName() + "_" + bedsuf_dns, tname_base + " [DNS]", tdesc_base + " [Called Non-Segregating]", i.getName());
				}
				
			}
			
			
		}
		
		
	}
	
}
