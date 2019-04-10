package hospelhornbg_svtools;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import hospelhornbg_bioinformatics.LiteSV;
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.VCF;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_bioinformatics.VariantPool.InfoDefinition;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class LargeMerger {
	
	public static final String OP_INPUT_MAIN = "-V"; 
	public static final String OP_INPUT_SNP_VCF = "-i";
	public static final String OP_RESCUE = "-r"; //Comma delimited
	public static final String OP_OUTPUT = "-o"; 
	
	/**
	 * ID = SNPARR_VISABLE
	 * <br>Number = 0
	 * <br>Type = Flag
	 * <br>Description = "Structural variant was visible on the SNP array"
	 */
	public static final InfoDefinition INFODEF_INFO_SNPARRVIS = new InfoDefinition("SNPARR_VISABLE", VariantPool.INFODEF_FLAG, "Structural variant was visible on the SNP array", 0);
	
	public static void printUsage()
	{
		System.out.println("--------------------------------------------------------------------------------------------------------------");
		System.out.println("BioisvTools || lrgmrg");
		System.out.println();
		System.out.println("Purpose: Merges a set of highly imprecise variants (such as those called from a SNP array) into a");
		System.out.println("pre-existing callset.");
		System.out.println();
		System.out.println("Input Formats:");
		System.out.println("\tVCF - Variant call format (all inputs)");
		System.out.println();
		System.out.println("Flags:");
		System.out.println("\t-V\tFILE\t[Required]\t\tPath to VCF that serves as merge target.");
		System.out.println("\t-i\tFILE\t[Required]\t\tPath to input VCF with imprecise variants to merge into target VCF.");
		System.out.println("\t-o\tFILE\t[Required]\t\tDesired output VCF path.");
		System.out.println("\t-r\tFILES\t[Optional]\t\tComma-delimited list of VCF paths to use as rescue sets.");
		System.out.println();
		System.out.println("Sample Usage:");
		System.out.println("java -jar bioisvtools.jar lrgmrg -g grch37 -V indiv_ngscalls.vcf -i indiv_snparr.vcf -o indiv_merged.vcf");
		System.out.println("java -jar bioisvtools.jar lrgmrg -g grch37 -V indiv_ngscalls.vcf -i indiv_snparr.vcf -o indiv_merged.vcf -r dellyset.vcf,mantaset.vcf,cnvnaterset.vcf");
		System.out.println();
		System.out.println("--------------------------------------------------------------------------------------------------------------");
	}
	
	public static List<LiteSV> readLargeVarVCF(String vcfpath, GenomeBuild gb) throws IOException
	{
		if(gb == null) return null;
		List<LiteSV> svlist = new LinkedList<LiteSV>();
		
		BufferedReader br = new BufferedReader(new FileReader(vcfpath));
		String line = null;
		while((line = br.readLine()) != null)
		{
			if (line.isEmpty()) continue;
			if (line.startsWith("#")) continue;
			LiteSV sv = LiteSV.readFromVCFLine(line, gb);
			if(sv != null) svlist.add(sv);
		}
		br.close();
		
		
		return svlist;
	}
	
	public static List<LiteSV> mergeLargeVariants(VariantPool variants, List<LiteSV> largeVars)
	{
		variants.addInfoField(INFODEF_INFO_SNPARRVIS.getKey(), INFODEF_INFO_SNPARRVIS);
		List<LiteSV> missed = new LinkedList<LiteSV>();
		List<Variant> vlist = variants.getVariants();
		
		//Map by chrom
		Map<Contig, List<StructuralVariant>> map = new HashMap<Contig, List<StructuralVariant>>();
		for(Variant v : vlist)
		{
			if (v instanceof StructuralVariant)
			{
				StructuralVariant sv = (StructuralVariant)v;
				Contig c = sv.getChromosome();
				if (c == null) continue;
				List<StructuralVariant> svlist = map.get(c);
				if (svlist == null)
				{
					svlist = new LinkedList<StructuralVariant>();
					map.put(c, svlist);
				}
				svlist.add(sv);
			}
		}
		
		for(LiteSV lsv : largeVars)
		{
			//Look for a match
			Contig c1 = lsv.getChrom1();
			if (c1 == null)
			{
				missed.add(lsv);
				continue;
			}
			List<StructuralVariant> svlist = map.get(c1);
			if(svlist == null)
			{
				System.err.println("LargeMerge.mergeLargeVariants || No variants mapped to contig " + c1.getUDPName() + " found in source set...");
				missed.add(lsv);
				continue;
			}
			boolean matched = false;
			for(StructuralVariant sv: svlist)
			{
				//Look for a match
				matched = lsv.svIsEquivalent(sv, 0, true);
				if(matched)
				{
					sv.addInfoFlag(INFODEF_INFO_SNPARRVIS.getKey());
					System.err.println("LargeMerge.mergeLargeVariants || -DEBUG- Variants matched: " + lsv.getVariantID() + " to " + sv.getVarID());
					break;
				}
			}
			if(!matched) missed.add(lsv);
		}
		
		return missed; //Returns the list of unmatched large variants. Can be compared against other VCFs
	}

	public static List<LiteSV> rescueFromVCF(VariantPool mainPool, List<LiteSV> rejects, Collection<Variant> rescueVars)
	{
		List<LiteSV> missed = new LinkedList<LiteSV>();
		
		Map<Contig, List<StructuralVariant>> map = new HashMap<Contig, List<StructuralVariant>>();
		for(Variant v : rescueVars)
		{
			if (v instanceof StructuralVariant)
			{
				StructuralVariant sv = (StructuralVariant)v;
				Contig c = sv.getChromosome();
				if (c == null) continue;
				List<StructuralVariant> svlist = map.get(c);
				if (svlist == null)
				{
					svlist = new LinkedList<StructuralVariant>();
					map.put(c, svlist);
				}
				svlist.add(sv);
			}
		}
		
		for(LiteSV lsv : rejects)
		{
			//Look for a match
			Contig c1 = lsv.getChrom1();
			if (c1 == null)
			{
				missed.add(lsv);
				continue;
			}
			List<StructuralVariant> svlist = map.get(c1);
			if(svlist == null)
			{
				System.err.println("LargeMerge.rescueFromVCF || No variants mapped to contig " + c1.getUDPName() + " found in rescue set...");
				missed.add(lsv);
				continue;
			}
			boolean matched = false;
			for(StructuralVariant sv: svlist)
			{
				//Look for a match
				matched = lsv.svIsEquivalent(sv, 0, true);
				if(matched)
				{
					sv.addInfoFlag(INFODEF_INFO_SNPARRVIS.getKey());
					mainPool.addVariant(sv);
					System.err.println("LargeMerge.rescueFromVCF || -DEBUG- Variants matched: " + lsv.getVariantID() + " to " + sv.getVarID());
					break;
				}
			}
			if(!matched) missed.add(lsv);
		}
		
		return missed; //Again, returns rejects
	}
	
	public static void runMerger(String[] args, GenomeBuild gb)
	{
		String vcfPath = null;
		String snpvcfPath = null;
		String outPath = null;
		String rescueList = null;
		
		//Get args
		for (int i = 0; i < args.length; i++)
		{
			String s = args[i];
			if (s.equals(OP_INPUT_MAIN))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_INPUT_MAIN + " flag MUST be followed by input VCF path!");
					printUsage();
					System.exit(1);
				}
				vcfPath = args[i+1];
			}
			else if (s.equals(OP_INPUT_SNP_VCF))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_INPUT_SNP_VCF + " flag MUST be followed by input VCF path!");
					printUsage();
					System.exit(1);
				}
				snpvcfPath = args[i+1];
			}
			else if (s.equals(OP_RESCUE))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_RESCUE + " flag MUST be followed by comma delimited list of VCF files!");
					printUsage();
					System.exit(1);
				}
				rescueList = args[i+1];
			}
			else if (s.equals(OP_OUTPUT))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_OUTPUT + " flag MUST be followed by output VCF path!");
					printUsage();
					System.exit(1);
				}
				outPath = args[i+1];
			}
		}
		
		//Check args
		if(vcfPath == null || vcfPath.isEmpty())
		{
			System.err.println("ERROR! Input variant set is required!");
			printUsage();
			System.exit(1);
		}
		if(snpvcfPath == null || snpvcfPath.isEmpty())
		{
			System.err.println("ERROR! Input large variant set is required!");
			printUsage();
			System.exit(1);
		}
		if(outPath == null || outPath.isEmpty())
		{
			System.err.println("ERROR! Output path is required!");
			printUsage();
			System.exit(1);
		}
		
		//Parse rescue set
		List<String> rescuePaths = new LinkedList<String>();
		if(rescueList != null && !rescueList.isEmpty())
		{
			String[] list = rescueList.split(",");
			for(String s : list) rescuePaths.add(s);
		}
		
		//Open main files
		VariantPool main = null;
		List<LiteSV> largeList = null;
		
		try {main = VCF.readVCF(vcfPath, gb, true);} 
		catch (IOException e) 
		{
			System.err.println("ERROR! Input vcf \"" + vcfPath + "\" could not be opened!");
			e.printStackTrace();
			System.exit(1);
		} 
		catch (UnsupportedFileTypeException e) 
		{
			System.err.println("ERROR! Input vcf \"" + vcfPath + "\" could not be read!");
			e.printStackTrace();
			System.exit(1);
		}
		
		try 
		{
			largeList = readLargeVarVCF(snpvcfPath, gb);
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR! Input vcf \"" + snpvcfPath + "\" could not be opened!");
			e.printStackTrace();
			System.exit(1);
		}
		System.err.println("LargeMerge.runMerger || -DEBUG- Input Variants (Large): " + largeList.size());
		
		//Run initial merge
		List<LiteSV> rejects = mergeLargeVariants(main, largeList);
		System.err.println("LargeMerge.runMerger || -DEBUG- Unmatched variants (First pass): " + rejects.size());
		
		//Run rescues
		for(String p : rescuePaths)
		{
			if (rejects == null || rejects.isEmpty()) break;
			VariantPool rpool = null;
			try {rpool = VCF.readVCF(p, gb, true);}
			catch (IOException e) 
			{
				System.err.println("ERROR! Rescue vcf \"" + p + "\" could not be opened!");
				e.printStackTrace();
				continue;
			} 
			catch (UnsupportedFileTypeException e) 
			{
				System.err.println("ERROR! Rescue vcf \"" + p + "\" could not be read!");
				e.printStackTrace();
				continue;
			}
			
			Collection<Variant> rescues = rpool.getVariants();
			rejects = rescueFromVCF(main, rejects, rescues);
			System.err.println("LargeMerge.runMerger || -DEBUG- Unmatched variants (Rescue): " + rejects.size());
		}
		
		//If any rejects remain, print to stdout
		if(rejects != null && !rejects.isEmpty())
		{
			System.out.println("No matches were found for the following variants:");
			for(LiteSV sv: rejects)
			{
				System.out.println("\t" + sv.getVariantID());
			}
			System.out.println("These variants will not appear in the output.");
		}
		
		//Write output VCF
		main.sortVariants();
		try 
		{
			VCF.writeVCF(main, "bioisvtools", outPath);
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR! Output vcf \"" + outPath + "\" could not be written!");
			e.printStackTrace();
			System.exit(1);
		}
		
	}

}
