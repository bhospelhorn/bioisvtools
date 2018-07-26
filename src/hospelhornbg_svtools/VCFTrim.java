package hospelhornbg_svtools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import hospelhornbg_bioinformatics.Genotype;
import hospelhornbg_bioinformatics.Variant;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

//Not multithreaded because that makes it easier to keep in order

public class VCFTrim {
	
	public static final String OP_VCFIN = "-i"; 
	public static final String OP_VCFOUT = "-o"; 
	public static final String OP_SAMPLES = "-s"; 
	
	public static void printUsage()
	{
		System.out.println("--------------------------------------------------------------------------------");
		System.out.println("BioisvTools || vcftrimsamp");
		System.out.println();
		System.out.println("Purpose: To trim samples out of a multi-sample VCF.");
		System.out.println("Removes both unwanted genotype columns and variants for which all desired");
		System.out.println("samples are homozygous reference.");
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
		System.out.println("\t-s\tSTRING\t[Required]\t\tComma delimited list of samples to keep.");
		System.out.println();
		System.out.println("Sample Usage:");
		System.out.println("java -jar bioisvtools.jar vcftrimsamp -i mycohort.vcf -o myfam.vcf -s PB123,UM456,UF789,AS098,US765");
		System.out.println("java -jar bioisvtools.jar vcftrimsamp -v -i mycohort.vcf -o PB123.vcf -s PB123");
		System.out.println();
		System.out.println("--------------------------------------------------------------------------------");
	}
	
	public static boolean keepVar(Variant v, Collection<String> validSamples)
	{
		for (String s: validSamples)
		{
			Genotype g = v.getSampleGenotype(s);
			if (g != null)
			{
				int z = g.getZygosity();
				switch(z)
				{
				case Genotype.ZYGOSITY_HETEROAA: return true;
				case Genotype.ZYGOSITY_HETEROAA_CNV: return true;
				case Genotype.ZYGOSITY_HETERORA: return true;
				case Genotype.ZYGOSITY_HETERORA_CNV: return true;
				case Genotype.ZYGOSITY_HETERORAA_CNV: return true;
				case Genotype.ZYGOSITY_HOMOALT: return true;
				case Genotype.ZYGOSITY_HOMOALT_CNV: return true;
				default: continue;
				}
			}
		}
		return false;
	}

	public static String trimLine(String varLine, List<String> allSamples, Collection<String> validSamples)
	{
		if (varLine == null) return null;
		if (varLine.isEmpty()) return null;
		
		String[] fields = varLine.split("\t");
		String outline = "";
		for (int i = 0; i < 9; i++)
		{
			if (i < fields.length)
			{
				if (i < fields.length - 1 && i < 8) outline += fields[i] + "\t";
				else outline += fields[i];
			}
		}
		
		if (fields.length > 9)
		{
			for (int i = 9; i < fields.length; i++)
			{
				String sample = allSamples.get(i-9);
				if (validSamples.contains(sample)) {
					outline += "\t" + fields[i];
				}
			}
		}
		else return varLine;
		
		return outline;
	}
	
	public static void trimFile(String inPath, String outPath, List<String> validSamples) throws IOException, UnsupportedFileTypeException
	{
		FileReader fr = new FileReader(inPath);
		BufferedReader br = new BufferedReader(fr);
		
		FileWriter fw = new FileWriter(outPath);
		BufferedWriter bw = new BufferedWriter(fw);
		
		String line = null;
		List<String> samplelist = new ArrayList<String>();
		while ((line = br.readLine()) != null)
		{
			if (line.isEmpty()) continue;
			if (line.charAt(0) == '#')
			{
				if (line.startsWith("#CHROM"))
				{
					String[] fields = line.split("\t");
					String outline = "";
					for (int i = 0; i < 9; i++)
					{
						if (i < fields.length)
						{
							if (i < fields.length - 1 && i < 8) outline += fields[i] + "\t";
							else outline += fields[i];
						}
					}
					if (fields.length > 9)
					{
						for (int i = 9; i < fields.length; i++)
						{
							samplelist.add(fields[i]);
							if (validSamples.contains(fields[i])) {
								outline += "\t" + fields[i];
							}
						}
					}
					bw.write(outline);
					//System.out.println(outline);
				}
				else bw.write(line + "\n");
			}
			else
			{
				try
				{
					Variant v = new Variant(line, samplelist, null, true);
					if (keepVar(v, validSamples))
					{
						bw.write("\n" + trimLine(line, samplelist, validSamples));
					}
				}
				catch (UnsupportedFileTypeException e)
				{
					System.err.println("PARSING ERROR! Variant Line: ");
					System.err.println(line);
					e.printStackTrace();
					fr.close();
					br.close();
					//fw.close();
					bw.close();
					throw e;
				}
				

			}
		}
		
		fr.close();
		br.close();
		
		//fw.close();
		bw.close();
		
	}
	
	public static void trimVCF(String[] args)
	{
		String inFile = null;
		String outFile = null;
		List<String> slist = new ArrayList<String>();
		
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
			else if (s.equals(OP_SAMPLES))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_SAMPLES + " flag MUST be followed by comma delimited list of samples!");
					printUsage();
					System.exit(1);
				}
				String samps = args[i+1];
				if (samps == null || samps.isEmpty())
				{
					System.err.println("ERROR: " + OP_SAMPLES + " flag MUST be followed by comma delimited list of samples!");
					printUsage();
					System.exit(1);
				}
				String[] sarr = samps.split(",");
				for (String x : sarr) slist.add(x);
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
		if (slist.isEmpty())
		{
			System.err.println("ERROR: Sample list is required!");
			pass = false;
		}
		if (!pass)
		{
			printUsage();
			System.exit(1);
		}
		
		System.err.println("Samples Requested: ");
		for (String s : slist)
		{
			System.err.println("\t" + s);	
		}
		
		try 
		{
			trimFile(inFile, outFile, slist);
		} 
		catch (IOException e) 
		{
			System.err.println("IO ERROR: One or more files could not be read/written!");
			e.printStackTrace();
		} 
		catch (UnsupportedFileTypeException e)
		{
			System.err.println("Parsing ERROR: Variant record could not be parsed!");
			e.printStackTrace();
		}
		
	}
	
}
