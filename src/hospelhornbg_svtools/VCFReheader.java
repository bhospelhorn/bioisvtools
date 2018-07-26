package hospelhornbg_svtools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * For changing the names of samples in a VCF header.
 * @author Blythe Hospelhorn
 * @version 1.0.0
 * @since June 21, 2018
 *
 */
public class VCFReheader {
	
	public static final String OP_INPUT = "-i";
	public static final String OP_OUTPUT = "-o";
	public static final String OP_TABLE = "-t";

	public static void printUsage()
	{
		System.out.println("--------------------------------------------------------------------------------");
		System.out.println("BioisvTools || VCF Sample Name Swapper");
		System.out.println();
		System.out.println("Purpose: For substituting the sample names in the header column of a VCF file.");
		System.out.println();
		System.out.println("Input Formats:");
		System.out.println("\tInput callset must be in [vcf] format");
		System.out.println("\tName replacement table must be tab-delimited plain text [tsv]");
		System.out.println("\t\tEach record should be formatted like so: [OLDNAME]\\t[NEWNAME]");
		System.out.println();
		System.out.println("Output Formats:");
		System.out.println("\tPlain text Variant Call Format [vcf]");
		System.out.println();
		System.out.println("Flags:");
		System.out.println("\t-i\tFILE\t[Required]\t\tInput file path");
		System.out.println("\t-o\tFILE\t[Required]\t\tOutput file path");
		System.out.println("\t-t\tFILE\t[Required]\t\tReplacement table name path");
		System.out.println();
		System.out.println("Sample Usage:");
		System.out.println("java -jar bioisvtools.jar vcfsns -i mycalls.vcf -o mycalls_clean.vcf -t nametable.tsv");
		System.out.println();
		System.out.println("--------------------------------------------------------------------------------");
	}
	
	public static Map<String, String> readNameTable(String tablePath) throws IOException
	{
		Map<String, String> map = new HashMap<String, String>();
		FileReader fr = new FileReader(tablePath);
		BufferedReader br = new BufferedReader(fr);
		
		String line = null;
		while ((line = br.readLine()) != null)
		{
			//Skip any lines that are less than two tab separated fields
			//If there is bad data, then too bad.
			String[] fields = line.split("\t");
			if (fields.length < 2) continue;
			map.put(fields[0], fields[1]);
		}
		
		br.close();
		fr.close();
		
		return map;
	}
	
	public static void changeSampleNames(String myVCF, String nametable, String outVCF)
	{
		Map<String, String> namemap = null;
		try
		{
			namemap = readNameTable(nametable);
		}
		catch (IOException e)
		{
			System.err.println("VCFReheader.changeSampleNames || ERROR: Sample table could not be read!");
			e.printStackTrace();
			System.exit(1);
		}
		
		if (namemap == null)
		{
			System.err.println("VCFReheader.changeSampleNames || ERROR: Sample table could not be read!");
			System.exit(1);
		}
		
		try
		{
			FileReader fr = new FileReader(myVCF);
			BufferedReader br = new BufferedReader(fr);
			
			FileWriter fw = new FileWriter(outVCF);
			BufferedWriter bw = new BufferedWriter(fw);
			
			boolean nl = false;
			String line = null;
			while ((line = br.readLine()) != null)
			{
				if (nl) bw.write("\n");
				nl = true;
				
				if (line.startsWith("#CHROM"))
				{
					String[] fields = line.split("\t");
					if (fields.length > 9)
					{
						for (int i = 9; i < fields.length; i++)
						{
							String sname = namemap.get(fields[i]);
							if (sname != null && !sname.isEmpty()) fields[i] = sname;
						}
					}
					//Rewrite line and output again
					String newline = "";
					for (int i = 0; i < fields.length; i++)
					{
						if (i > 0) newline += "\t";
						newline += fields[i];
					}
					bw.write(newline);
				}
				else bw.write(line);
			}
			
			bw.close();
			fw.close();
			br.close();
			fr.close();
		}
		catch (IOException e)
		{
			System.err.println("VCFReheader.changeSampleNames || ERROR: Reading or writing VCF");
			e.printStackTrace();
			System.exit(1);
		}
		
	}
	
	public static void reheader(String[] args)
	{
		String inPath = null;
		String outPath = null;
		String tblPath = null;
		
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
					System.err.println("ERROR: " + OP_OUTPUT + " flag MUST be followed by input VCF path!");
					printUsage();
					System.exit(1);
				}
				outPath = args[i+1];
			}
			if (s.equals(OP_TABLE))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_TABLE + " flag MUST be followed by input VCF path!");
					printUsage();
					System.exit(1);
				}
				tblPath = args[i+1];
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
		if (tblPath == null || tblPath.isEmpty())
		{
			System.err.println("ERROR: Name table (tab delimited) is required!");
			pass = false;
		}
		
		if (!pass)
		{
			printUsage();
			System.exit(1);
		}
		
		changeSampleNames(inPath, tblPath, outPath);
		
	}
	
}
