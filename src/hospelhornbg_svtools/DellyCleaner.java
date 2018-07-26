package hospelhornbg_svtools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class DellyCleaner {
	
	public static final String OP_INPUT = "-i"; 
	public static final String OP_OUTPUT = "-o"; 
	
	public static void printUsage()
	{
		System.out.println("--------------------------------------------------------------------------------");
		System.out.println("BioisvTools || CleanDelly");
		System.out.println();
		System.out.println("Purpose: To filter DELLY sv caller output variants - removes variants marked as");
		System.out.println("low quality and variants that are over 5000000bp in size");
		System.out.println();
		System.out.println("Input Formats:");
		System.out.println("\tInput callset must be in [vcf] format with standard SV annotations.");
		System.out.println();
		System.out.println("Output Formats:");
		System.out.println("\tPlain text Variant Call Format [vcf]");
		System.out.println();
		System.out.println("Flags:");
		System.out.println("\t-i\tFILE\t[Required]\t\tInput callset in the form of a vcf file.");
		System.out.println("\t-o\tFILE\t[Required]\t\tOutput file path");
		System.out.println();
		System.out.println("Sample Usage:");
		System.out.println("java -jar bioisvtools.jar cleandelly -i NA12878_delly.vcf -o NA12878_dellyfiltered.vcf");
		System.out.println("java -jar bioisvtools.jar cleandelly -v -i NA12878_delly.vcf -o NA12878_dellyfiltered.vcf");
		System.out.println();
		System.out.println("--------------------------------------------------------------------------------");
	}
	
	public static void cleanDelly(String inpath, String outpath, boolean tossLQ, boolean verbose) throws IOException, UnsupportedFileTypeException
	{
		FileReader fr = new FileReader(inpath);
		BufferedReader br = new BufferedReader(fr);
		
		FileWriter fw = new FileWriter(outpath);
		BufferedWriter bw = new BufferedWriter(fw);
		
		String line = null;
		boolean first = true;
		int tosscount = 0;
		int vcount = 0;
		int passcount = 0;
		while ((line = br.readLine()) != null)
		{
			//If header line, just copy
			if (line.charAt(0) == '#') {
				if (!first) bw.write("\n");
				first = false;
				bw.write(line);
				continue;
			}
			vcount++;
			//Split
			String[] fields = line.split("\t");
			//Check filter
			if (fields.length < 8) continue;
			if (tossLQ && ((fields[6]==null) || !(fields[6].equals("PASS")))) {
				if(verbose) System.err.println("DellyCleaner.cleanDelly || Variant tossed: " + fields[2] + " | Reason: Low quality");
				tosscount++;
				//TODO: DELLY appears to write BNDs like translocations instead of as pairs... Update SV parser?
				continue;
			}
			//Check type. If it's BND or TRA, then write and continue (don't check size)
			if(fields[2].startsWith("BND"))
			{
				if (!first) bw.write("\n");
				first = false;
				bw.write(line);
				passcount++;
				continue;
			}
			//Check SV size - if it's too big, then toss
			final int MAXSIZE = 5000000;
			final int MINSIZE = 300;
			//Get position
			int pos = -1;
			try {pos = Integer.parseInt(fields[1]);}
			catch (NumberFormatException e) {br.close(); bw.close(); fw.close(); fr.close(); throw new FileBuffer.UnsupportedFileTypeException();}
			//Get end
			int end = -1;
			String info = fields[7];
			int ECHAR = info.indexOf("END=");
			if (ECHAR < 0)
			{
				bw.close();
				fw.close();
				br.close();
				fr.close();
				throw new FileBuffer.UnsupportedFileTypeException();
			}
			int nextsemi = info.indexOf(';', ECHAR);
			try {end = Integer.parseInt(info.substring(ECHAR + 4, nextsemi));}
			catch (NumberFormatException e) {br.close(); bw.close(); fw.close(); fr.close(); throw new FileBuffer.UnsupportedFileTypeException();}
			int svlen = end - pos;
			if (Math.abs(svlen) > MAXSIZE) {
				if(verbose) System.err.println("DellyCleaner.cleanDelly || Variant tossed: " + fields[2] + " | Reason: Too large (Size = " + svlen + " bp)");
				tosscount++;
				continue;
			}
			if (Math.abs(svlen) < MINSIZE) {
				if(verbose) System.err.println("DellyCleaner.cleanDelly || Variant tossed: " + fields[2] + " | Reason: Too small (Size = " + svlen + " bp)");
				tosscount++;
				continue;
			}
			if (!first) bw.write("\n");
			first = false;
			bw.write(line);
			passcount++;
		}
		
		
		bw.close();
		fw.close();
		
		br.close();
		fr.close();
		
		if(verbose)
		{
			System.out.println("DellyCleaner.cleanDelly || Total Variants Processed: " + vcount);
			System.out.println("DellyCleaner.cleanDelly || Variants Passed: " + passcount);
			System.out.println("DellyCleaner.cleanDelly || Variants Rejected: " + tosscount);
		}
	}
	
	public static void runTool(String[] args, boolean verbose)
	{
		String inPath = null;
		String outPath = null;
		
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
		
		try 
		{
			cleanDelly(inPath, outPath, true, verbose);
		} 
		catch (IOException e) 
		{
			System.err.println("DellyCleaner.runTool || ERROR: Input or output file path was not valid!");
			e.printStackTrace();
			System.exit(1);
		} 
		catch (UnsupportedFileTypeException e) 
		{
			System.err.println("DellyCleaner.runTool || ERROR: Input file could not be processed!");
			e.printStackTrace();
			System.exit(1);
		}
		
	}

}
