package hospelhornbg_svtools;

import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import hospelhornbg_bioinformatics.BreakendPair;
import hospelhornbg_bioinformatics.Genotype;
import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.VCF;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class LumpyFilter {
	
	public static final String OP_INPUT = "-i"; 
	public static final String OP_OUTPUT = "-o"; 
	public static final String OP_XMODE = "-x"; 
	public static final String OP_SOFT = "-s"; //Instead of throwing out variants that don't pass, just mark in FILTER column. 
	public static final String OP_POSTGENO_THRESHOLD = "-t"; //Defo [0.01] (default mode), TBD for x-mode 
	
	public static final double DEFO_THRESHOLD_STANDARD = 0.01;
	public static final double DEFO_THRESHOLD_XMODE = 0.01; //TODO: Update when ready
	
	public static final String FILTER_S_KEY = "svtQual";
	public static final String FILTER_Z_KEY = "size";
	public static final String FILTER_G_KEY = "homref";
	//public static final String FILTER_S_DEF = "Called with a quality";
	public static final String FILTER_X_KEY = "inEv";
	//public static final String FILTER_X_DEF = "";
	
	public static final int MIN_SIZE = 50;
	public static final int MAX_SIZE = 5000000;
	
	//Filter a LUMPY output VCF, either directly out of LUMPY (mode x),
	// or after svtyper genotyping.
	
	//The post-geno filter very simply removes any variants with a quality of 0
	//	or that are homref.
	//The x-mode filter actually attempts to score based on evidence and prob
	//	curves and mark anything that doesn't pass.
	
	public static void printUsage()
	{
		System.out.println("--------------------------------------------------------------------------------");
		System.out.println("BioisvTools || LumpyFilter");
		System.out.println();
		System.out.println("Purpose: To filter LUMPY sv caller and svtyper genotyper outputs.");
		System.out.println();
		System.out.println("Input Formats:");
		System.out.println("\tInput callset must be in [vcf] format with standard SV annotations.");
		System.out.println("\t\tIf using x-mode, LUMPY annotations are also required.");
		System.out.println();
		System.out.println("Output Formats:");
		System.out.println("\tPlain text Variant Call Format [vcf]");
		System.out.println();
		System.out.println("Flags:");
		System.out.println("\t-i\tFILE\t[Required]\t\tInput callset in the form of a vcf file.");
		System.out.println("\t-o\tFILE\t[Required]\t\tOutput file path");
		System.out.println("\t-x\t\t[Optional]\t\tRun in x-mode. X-mode evaluates LUMPY evidence instead of variant quality score.");
		System.out.println("\t-s\t\t[Optional]\t\tSoft mode. Instead of throwing away failed variants, mark as failed in the FILTER field.");
		System.out.println("\t-t\tFLOAT\t[Optional]\t\tThreshold value for filtering. Default [0.01] for standard mode, [0.01] for x-mode");
		System.out.println();
		System.out.println("NOTE ~~~");
		System.out.println("Keep in mind that the default filter operates on svtyper output - it fails any variants");
		System.out.println("given a quality score (variant score from QUAL VCF field) below the threshold.");
		System.out.println("If run on a set with no QUAL scores, it will output a set identical to the original!");
		System.out.println();
		System.out.println("Sample Usage:");
		System.out.println("java -jar bioisvtools.jar lumpyfilter -g GRCh37 -i NA12878_geno.vcf -o NA12878_genofiltered.vcf");
		System.out.println("java -jar bioisvtools.jar lumpyfilter -g GRCh37 -v -i NA12878_geno.vcf -o NA12878_genofiltered.vcf -s -t 0.1");
		System.out.println("java -jar bioisvtools.jar lumpyfilter -g GRCh37 -i NA12878_lumpy.vcf -o NA12878_lumpyfiltered.vcf -x");
		System.out.println();
		System.out.println("--------------------------------------------------------------------------------");
	}
	
	public static VariantPool pregeno_filter(VariantPool fullset, double threshold, boolean soft)
	{
		//TODO: Write this when filter is designed.
		return fullset;
	}
	
	public static VariantPool postgeno_filter(VariantPool fullset, double threshold, boolean soft, String pathstem)
	{
		Set<Variant> passedVar = new HashSet<Variant>();
		Set<Variant> failedVar = new HashSet<Variant>();
		Set<Variant> smallVar = new HashSet<Variant>();
		final String FILTER_DEF = "Genotyped by svtyper with a variant quality below " + threshold;
		final String SIZEFILTER_DEF = "Called size falls outside of reliably detectable range (" + MIN_SIZE + " bp - " + MAX_SIZE + " bp)";
		final String GENOFILTER_DEF = "Variant was genotyped as homozygous reference in sample of interest";
		
		Collection<Variant> allVar = fullset.getVariants();
		String sample = fullset.getAllSamples().get(0);
		if (sample == null) return null;
		
		int c_qual_bnd = 0;
		int c_qual = 0;
		int c_sz_s = 0;
		int c_sz_l = 0;
		int c_homref_bnd = 0;
		int c_homref = 0;
		
		int c = 0;
		int max = 0;
		for (Variant v : allVar)
		{
			//c++;
			//Check quality
			boolean failed = false;
			double qual = v.getQuality();
			if (v instanceof BreakendPair)
			{
				if (qual < threshold && qual != -1.0) {
					v.addFailedFilter(FILTER_S_KEY);
					failedVar.add(v);
					failed = true;
					c_qual_bnd++;
				}
				BreakendPair bp = (BreakendPair)v;
				Genotype geno1 = bp.getGenotype1(sample);
				Genotype geno2 = bp.getGenotype2(sample);
				boolean failg = false;
				if (geno1 != null)
				{
					if(geno1.getPercentAlt() < 0.01)
					{
						v.addFailedFilter(FILTER_G_KEY);
						failedVar.add(v);
						failed = true;
						failg = true;
						c_homref_bnd++;
					}
				}
				else c++;
				if (geno2 != null)
				{
					if(!failg && geno2.getPercentAlt() < 0.01)
					{
						v.addFailedFilter(FILTER_G_KEY);
						failedVar.add(v);
						failed = true;
						c_homref_bnd++;
					}
				}
				else c++;
			}
			else
			{
				if (qual < threshold && qual != -1.0) {
					v.addFailedFilter(FILTER_S_KEY);
					failedVar.add(v);
					failed = true;
					c_qual++;
				}
				if(v.getSmallestAbsoluteLength() < MIN_SIZE)
				{
					v.addFailedFilter(FILTER_Z_KEY);
					failedVar.add(v);
					failed = true;
					smallVar.add(v);
					c_sz_s++;
				}
				int len = v.getLargestAbsoluteLength();
				if (len > max) max = len;
				if(v.getLargestAbsoluteLength() > MAX_SIZE)
				{ 
					if (v instanceof StructuralVariant)
					{
						StructuralVariant sv = (StructuralVariant) v;
						if (sv.getType() != SVType.INV)
						{
							v.addFailedFilter(FILTER_Z_KEY);
							failedVar.add(v);
							failed = true;
							c_sz_l++;
						}
					}
					else
					{
						v.addFailedFilter(FILTER_Z_KEY);
						failedVar.add(v);
						failed = true;
						c_sz_l++;
					}
				}
				Genotype geno = v.getSampleGenotype(sample);
				if (geno != null)
				{
					if(geno.getPercentAlt() < 0.01)
					{
						v.addFailedFilter(FILTER_G_KEY);
						failedVar.add(v);
						failed = true;
						c_homref++;
					}
				}
				else
				{
					//System.err.println("LumpyFilter.postgeno_filter || WARNING: Genotype for sample " + sample + " was not found for variant " + v.getVarID());
					c++;
					//System.exit(1);
					//v.addFailedFilter(FILTER_G_KEY);
					//failedVar.add(v);
					//failed = true;
				}
			}
			
			if (!failed) passedVar.add(v);
		}
		//Now we check 
		if (c > 0) System.err.println("LumpyFilter.postgeno_filter || WARNING: Genotype not found for " + c + " of " + allVar.size() + " variants!");
		int tCount = allVar.size();
		int pCount = passedVar.size();
		int fCount = failedVar.size();
		System.err.println("LumpyFilter.postgeno_filter || Input variants: " + tCount);
		System.err.println("LumpyFilter.postgeno_filter || Passed variants: " + pCount);
		System.err.println("LumpyFilter.postgeno_filter || Failed variants: " + fCount);
		System.err.println("LumpyFilter.postgeno_filter || \t\tQual below " + threshold + " : " + c_qual);
		System.err.println("LumpyFilter.postgeno_filter || \t\tQual below " + threshold + " (BND): " + c_qual_bnd);
		System.err.println("LumpyFilter.postgeno_filter || \t\tCalled Homozygous Ref: " + c_homref);
		System.err.println("LumpyFilter.postgeno_filter || \t\tCalled Homozygous Ref (BND): " + c_homref_bnd);
		System.err.println("LumpyFilter.postgeno_filter || \t\tBelow minimum size (50 bp): " + c_sz_s);
		System.err.println("LumpyFilter.postgeno_filter || \t\tExceeding maximum size (5 Mbp): " + c_sz_l);
		System.err.println("LumpyFilter.postgeno_filter || Largest variant found: " + max + " bp");
		
		if (soft)
		{
			fullset.addFilter(FILTER_S_KEY, FILTER_DEF);
			fullset.addFilter(FILTER_Z_KEY, SIZEFILTER_DEF);
			fullset.addFilter(FILTER_G_KEY, GENOFILTER_DEF);
			for (Variant v : passedVar) v.setFilterPass(true);
			for (Variant v : failedVar)
			{
				v.setFilterPass(false);
				//v.addFailedFilter(FILTER_S_KEY);
			}
			/*Redundant, as variants are referenced
			fullset.clearVariants();
			fullset.addVariants(passedVar);
			fullset.addVariants(failedVar);
			 * */
			fullset.sortVariants();
		}
		else
		{
			
			fullset.clearVariants();
			for (Variant v : passedVar) v.setFilterPass(true);
			fullset.addVariants(passedVar);
			fullset.sortVariants();
			
			VariantPool smallset = fullset.makeSkeletonCopy();
			smallset.addFilter(FILTER_S_KEY, FILTER_DEF);
			smallset.addFilter(FILTER_Z_KEY, SIZEFILTER_DEF);
			smallset.addFilter(FILTER_G_KEY, GENOFILTER_DEF);
			for (Variant v : smallVar) v.setFilterPass(false);
			smallset.addVariants(smallVar);
			String outPath = pathstem + "_small.vcf";
			try 
			{
				VCF.writeVCF(smallset, "BioisvTools", outPath);
			} 
			catch (IOException e) 
			{
				System.err.println("ERROR: File \"" + outPath + "\" could not be written!");
				e.printStackTrace();
			}
		}
	
		
		return fullset;
	}
	
	public static void lumpyFilter(String[] args, GenomeBuild g)
	{
		String inPath = null;
		String outPath = null;
		boolean xMode = false; //Optional
		boolean softMode = false; //Optional
		double threshold = DEFO_THRESHOLD_STANDARD; //Optional
		boolean tSpecified = false;
		
		for (int i = 0; i < args.length; i++)
		{
			String s = args[i];
			if (s.equals(OP_INPUT))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: -i flag MUST be followed by input VCF path!");
					printUsage();
					System.exit(1);
				}
				inPath = args[i+1];
			}
			if (s.equals(OP_OUTPUT))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: -o flag MUST be followed by output VCF path!");
					printUsage();
					System.exit(1);
				}
				outPath = args[i+1];
			}
			if (s.equals(OP_POSTGENO_THRESHOLD))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: -t flag MUST be followed by valid double float value!");
					printUsage();
					System.exit(1);
				}
				try
				{
					threshold = Double.parseDouble(args[i+1]);
				}
				catch (NumberFormatException e)
				{
					System.err.println("ERROR: Threshold (-t) argument must be a valid floating point value!");
					printUsage();
					System.exit(1);
				}
				tSpecified = true;
			}
			if (s.equals(OP_XMODE))
			{
				xMode = true;
			}
			if (s.equals(OP_SOFT))
			{
				softMode = true;
			}
		}
		//End for
		
		//Check for required arguments
		if (inPath == null || inPath.isEmpty())
		{
			System.err.println("ERROR: Input path is required!");
			printUsage();
			System.exit(1);
		}
		if (outPath == null || outPath.isEmpty())
		{
			System.err.println("ERROR: Output path is required!");
			printUsage();
			System.exit(1);
		}
		
		//Adjust threshold
		if (!tSpecified && xMode)
		{
			//Change threshold to defo xmode value
			threshold = DEFO_THRESHOLD_XMODE;
		}
		
		//Attempt to read input VCF
		VariantPool input = null;
		try 
		{
			VCF vcfreader = new VCF(inPath, true, g);
			input = vcfreader.getVariants();
		} 
		catch (UnsupportedFileTypeException e) 
		{
			System.err.println("ERROR: File \"" + inPath + "\" does not appear to be a valid VCF file!");
			e.printStackTrace();
			System.exit(1);
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: File \"" + inPath + "\" could not be read!");
			e.printStackTrace();
			System.exit(1);
		}
		
		if (input == null)
		{
			System.err.println("ERROR: There was a problem parsing the input variants!");
			System.exit(1);
		}
		input.setGenomeBuild(g);
		
		//Filter the set
		VariantPool output = null;
		if (xMode)
		{
			System.err.println("x mode is not currently available!");
			System.exit(1);
			output = pregeno_filter(input, threshold, softMode);
		}
		else
		{
			int vcfi = outPath.lastIndexOf(".vcf");
			String pathstem = outPath;
			if (vcfi >= 0) pathstem = outPath.substring(0, vcfi);
			output = postgeno_filter(input, threshold, softMode, pathstem);
		}
		
		//Write the set.
		if (output == null)
		{
			System.err.println("ERROR: There was a problem filtering the set!");
			System.exit(1);
		}
		try 
		{
			VCF.writeVCF(output, "BioisvTools", outPath);
			System.err.println("Filtered VCF was successfully written to " + outPath);
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: File \"" + outPath + "\" could not be written!");
			e.printStackTrace();
			System.exit(1);
		}
		
	}//End main
	
	
}//End class
