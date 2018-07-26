package hospelhornbg_svtools;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;

import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class TrackMaker {
	
	public static final String OP_INPUT = "-i";
	public static final String OP_OUTPUT = "-o";
	public static final String OP_TRACKNAME = "-n";
	public static final String OP_DESCRIPTION = "-d";
	public static final String OP_TEMPLATE = "-t";
	public static final String OP_SAMPLE = "-s";

	public static final String TEMPLATE_DEFAULT = "default";
	public static final String TEMPLATE_STANDARDSV = "svdefo";
	
	public static Map<String, TrackTemplate> templateMap;
	
	public static void printUsage()
	{
		System.out.println("--------------------------------------------------------------------------------");
		System.out.println("BioisvTools || TrackMaker");
		System.out.println();
		System.out.println("Purpose: To generate a UCSC genome browser ready BED file with the variants from");
		System.out.println("a provided VCF file represented as blocks on a track.");
		System.out.println();
		System.out.println("Input Formats:");
		System.out.println("\tInput callset must be in [vcf] format");
		System.out.println("\tTemplate name must be know. See list below.");
		System.out.println();
		System.out.println("Output Formats:");
		System.out.println("\tBED12 format [bed]");
		System.out.println();
		System.out.println("Flags:");
		System.out.println("\t-i\tFILE\t[Required]\t\tInput VCF");
		System.out.println("\t-o\tFILE\t[Required]\t\tOutput file path.");
		System.out.println("\t\t\t\t\tIf sample is not specified, the path will be truncated to its last separator character (/ on UNIX, \\ on Windows)");
		System.out.println("\t\t\t\t\tto determine target directory, and tracks for all samples in the VCF will be saved to files named sample_trackname.bed ");
		System.out.println("\t-n\tSTRING\t[Required]\t\tName of BED track to appear in genome browser.");
		System.out.println("\t-d\tSTRING\t[Required]\t\tDescription of BED track to appear in genome browser.");
		System.out.println("\t-t\tSTRING\t[Required]\t\tName of template to use when generating track.");
		System.out.println("\t-s\tSTRING\t[Optional]\t\tName of single sample to generate track for.");
		System.out.println();
		System.out.println("Valid Templates:");
		System.out.println("\t" + TEMPLATE_DEFAULT);
		System.out.println("\t" + TEMPLATE_STANDARDSV);
		System.out.println();
		System.out.println("Sample Usage:");
		System.out.println("java -jar bioisvtools.jar trackmaker -g GRCh37 -i NA12878_svset.vcf -o /somedirs/NA12878/ucscbed/NA12878sv.bed + -n \"NA12878 SV\"");
		System.out.println("-d \"A structural variant callset for NA12878\" -t svdefo -s NA12878");
		System.out.println("java -jar bioisvtools.jar trackmaker -g hg38 -i somesvs.vcf -o /somedirs/Feb2018Cohort/ucscbed/ + -n \"Some SVs\"");
		System.out.println("-d \"A structural variant callset.\" -t svdefo");
		System.out.println();
		System.out.println("--------------------------------------------------------------------------------");
	}
	
	public static void populateTemplateMap(GenomeBuild g)
	{
		templateMap = new HashMap<String, TrackTemplate>();
		templateMap.put(TEMPLATE_DEFAULT, new STDefoTemplate(g));
		templateMap.put(TEMPLATE_STANDARDSV, new SVDefoTemplate(g));
	}
	
	public static void trackMaker(String[] args, GenomeBuild g)
	{
		String inPath = null;
		String outPath = null;
		String tName = null;
		String tDesc = null;
		String template = null;
		String samplename = null;
		
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
					System.err.println("ERROR: -o flag MUST be followed by output directory path!");
					printUsage();
					System.exit(1);
				}
				outPath = args[i+1];
			}
			if (s.equals(OP_TRACKNAME))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: -n flag MUST be followed by desired track name!");
					printUsage();
					System.exit(1);
				}
				tName = args[i+1];
			}
			if (s.equals(OP_DESCRIPTION))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: -d flag MUST be followed by desired track description!");
					printUsage();
					System.exit(1);
				}
				tDesc = args[i+1];
			}
			if (s.equals(OP_TEMPLATE))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: -t flag MUST be followed by valid template name!");
					printUsage();
					System.exit(1);
				}
				template = args[i+1];
			}
			if (s.equals(OP_SAMPLE))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_SAMPLE + " flag MUST be followed by sample name!");
					printUsage();
					System.exit(1);
				}
				samplename = args[i+1];
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
		if (tName == null || tName.isEmpty())
		{
			System.err.println("ERROR: Track name is required!");
			pass = false;
		}
		if (tDesc == null || tDesc.isEmpty())
		{
			System.err.println("ERROR: Track description is required!");
			pass = false;
		}
		if (template == null || template.isEmpty())
		{
			System.err.println("ERROR: Template name is required!");
			pass = false;
		}
		
		if (!pass)
		{
			printUsage();
			System.exit(1);
		}
		
		//Map templates to template names, and check to see if template requested exists.
		populateTemplateMap(g);
		TrackTemplate tt = templateMap.get(template);
		if (tt == null)
		{
			System.err.println("ERROR: Template \"" + template + "\" is not a valid template!");
			printUsage();
			System.exit(1);
		}
		
		//Check for the sampleName. If not there, need to generate an outdir
		if (samplename == null || samplename.isEmpty())
		{
			//Check for outDir if that will be needed...
			String outDir = null;
			int li = outPath.length() - 1;
			if (outPath.charAt(li) == File.separatorChar) outDir = outPath.substring(0,li);
			if (!FileBuffer.directoryExists(outDir))
			{
				System.err.println("Output directory " + outDir + " does not exist. Creating...");
				try
				{
					Files.createDirectories(Paths.get(outDir));
				} 
				catch (IOException e) {
					System.err.println("ERROR: Output directory " + outDir + " could not be created!");
					e.printStackTrace();
					System.exit(1);
				}
			}
			//Now actually do the track generation.
			try 
			{
				tt.setTrackName(tName);
				tt.setTrackDescription(tDesc);
				tt.generateGroup(inPath, outDir);
			} 
			catch (UnsupportedFileTypeException e) 
			{
				System.err.println("ERROR: Input file could not be parsed as a VCF!");
				e.printStackTrace();
				System.exit(1);
			} 
			catch (IOException e) 
			{
				System.err.println("ERROR: I/O error reading VCF or writing track!");
				e.printStackTrace();
				System.exit(1);
			}
		}
		else
		{
			try 
			{
				tt.setTrackName(tName);
				tt.setTrackDescription(tDesc);
				tt.generateSingle(inPath, outPath, samplename);
			} 
			catch (UnsupportedFileTypeException e) 
			{
				System.err.println("ERROR: Input file could not be parsed as a VCF!");
				e.printStackTrace();
				System.exit(1);
			} 
			catch (IOException e) 
			{
				System.err.println("ERROR: I/O error reading VCF or writing track!");
				e.printStackTrace();
				System.exit(1);
			}
		}
	
	
	}
	
}
