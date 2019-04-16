package hospelhornbg_svtools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import hospelhornbg_genomeBuild.GeneSet;
import hospelhornbg_genomeBuild.GenomeBuild;
import hospelhornbg_svdb.ConsoleFrontEnd;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class ConsoleMain {
	
	public static final String TOOL_INSTALLME = "install";
	public static final String TOOL_CHECKINSTALL = "checkinstall";
	public static final String TOOL_TRACKMAKER = "trackmaker";
	public static final String TOOL_PENNTOVCF = "pennconvert";
	public static final String TOOL_FILTERLUMPY = "lumpyfilter";
	public static final String TOOL_GBMANAGER = "gbmanager";
	public static final String TOOL_MERGEGENO = "survivorgeno";
	public static final String TOOL_STDCHROM = "stdchrom";
	public static final String TOOL_SVANNO = "svanno";
	public static final String TOOL_SVQSTATS = "svqstats";
	public static final String TOOL_SZTALLY = "varsztally";
	public static final String TOOL_VCFTRIM = "vcftrimsamp";
	public static final String TOOL_SZFILTER = "sizefilter";
	public static final String TOOL_SPLITMERGE = "smsv";
	public static final String TOOL_VCFSAMPLENAMESWAP = "vcfsns";
	public static final String TOOL_DELLYCLEANER = "cleandelly";
	public static final String TOOL_TRIMCHROM = "trimchrom";
	public static final String TOOL_SVANALYZE = "svreport";
	public static final String TOOL_SVGENEHITTALLY = "svght";
	public static final String TOOL_SCANSAM = "scansam";
	public static final String TOOL_FIXSAM = "fixsam";
	public static final String TOOL_PEDTOFAMI = "ped2fami";
	public static final String TOOL_VIEWFAMI = "viewfami";
	public static final String TOOL_SVDB = "svdb";
	public static final String TOOL_LRGMRG = "lrgmrg";

	public static final String OP_GENOMEBUILD = "-g";
	public static final String OP_VERBOSE = "-v";
	
	public static final int OS_WINDOWS = 0;
	public static final int OS_LINUX = 1;
	public static final int OS_MAC = 2;
	
	public static final String DIR_CONFIG = "bioisvtools";
	public static final String FILE_GBDAT = "gb.dat";
	public static final String FILE_INIT = "init.cfg";

	public static void printUsage()
	{
		System.out.println("--------------------------------------------------------------------------------");
		System.out.println(" == BioisvTools ==");
		System.out.println();
		System.out.println("Purpose: A collection of various tools originally written for the NHGRI UDP SV detection pipeline.");
		System.out.println();
		System.out.println("Tools:");
		System.out.println("\tinstall\t\tGenerate an init file in the user's home directory that points to a directory containing genome build metadata.");
		System.out.println("\tcheckinstall\tCheck whether an init file has been setup for this user. Prints NOINSTALL if not. Prints reference paths init file is present.");
		System.out.println("\tgbmanager\tManage the genome build metadata table for user.");
		System.out.println("\ttrackmaker\tGenerate a BED file formatted for the UCSC genome browser from a VCF file.");
		System.out.println("\tpennconvert\tConvert PennCNV SNP array CNV call information to VCF.");
		System.out.println("\tlumpyfilter\tChoose one of two filters for quick filtering of lumpy/svtyper callsets.");
		System.out.println("\tsurvivorgeno\tAdd a consensus genotype column from a SURVIVOR output VCF.");
		System.out.println("\tstdchrom\tStandardize the name of the chromosomes in a SAM file.");
		System.out.println("\tsvanno\tAnnotate structural variants with refGene. (Note: This tool does not work with hg18 by default)");
		System.out.println("\tvcfsns\tChange sample names in a VCF file (VCF Sample Name Swapper)");
		System.out.println("\tsvreport\tPrint files separated by SV type and position effect containing candidate information from a family merged callset.");
		System.out.println("\tsvght\tTally number of times gene hits are found in different families from svreport output.");
		System.out.println("\tscansam\tScan a SAM file to look for syntax errors.");
		System.out.println("\tfixsam\tScan a SAM file to look for contig name errors, and output a fixed version.");
		System.out.println();
		System.out.println("Flags:");
		System.out.println("\t-g\tSTRING\t[Usually Required]\t\tName (case insensitive) of genome build to use with input.");
		System.out.println("\t-v\t\t[Optional]\t\t\tVerbose");
		System.out.println();
		System.out.println("Sample Usage:");
		System.out.println("java -jar bioisvtools.jar install [...]");
		System.out.println("java -jar bioisvtools.jar checkinstall [...]");
		System.out.println("java -jar bioisvtools.jar gbmanager -v [...]");
		System.out.println("java -jar bioisvtools.jar trackmaker -g GRCh37 -v [...]");
		System.out.println("java -jar bioisvtools.jar pennconvert -g hg18 [...]");
		System.out.println("java -jar bioisvtools.jar lumpyfilter -g hg38 -v [...]");
		System.out.println("java -jar bioisvtools.jar survivorgeno -g hg19 -v [...]");
		System.out.println("java -jar bioisvtools.jar stdchrom -g GRCh38 -v [...]");
		System.out.println("java -jar bioisvtools.jar svanno -g grch37 -v [...]");
		System.out.println("java -jar bioisvtools.jar svreport -g hg19 -v [...]");
		System.out.println("java -jar bioisvtools.jar svght -v [...]");
		System.out.println("java -jar bioisvtools.jar samscan (-g grch37) -v [...]");
		System.out.println("java -jar bioisvtools.jar fixsam -g grch37 -v [...]");
		System.out.println();
		System.out.println("--------------------------------------------------------------------------------");
	}
	
	public static void cleanDeadRefPaths(String homedir)
	{
		String initPath = homedir + File.separator + DIR_CONFIG + File.separator + FILE_INIT;
		List<String> paths = getRefPaths(homedir);
		
		try 
		{
			Files.deleteIfExists(Paths.get(initPath));
			
			FileWriter fw = new FileWriter(initPath);
			BufferedWriter bw = new BufferedWriter(fw);
			
			boolean first = true;
			for (String p : paths)
			{
				if (FileBuffer.fileExists(p))
				{
					if (!first) bw.write("\n");
					bw.write(p);
					first = false;
				}
			}
			
			bw.close();
			fw.close();
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: There was an error updating the reference paths config file!");
			e.printStackTrace();
			System.exit(1);
		}
		
	}
	
	public static List<String> getRefPaths(String homedir)
	{
		String initPath = homedir + File.separator + DIR_CONFIG + File.separator + FILE_INIT;
		//System.err.println("DEBUG || Checking for init file at: " + initPath);
		if (!FileBuffer.fileExists(initPath)) return null;
		
		List<String> paths = new LinkedList<String>();
		try
		{
			FileReader fr = new FileReader(initPath);
			BufferedReader br = new BufferedReader(fr);
			
			String line = null;
			while ((line = br.readLine()) != null)
			{
				paths.add(line);
			}
			
			br.close();
			fr.close();
			
		}
		catch (IOException e)
		{
			System.err.println("ERROR: There was an error reading the init file!");
			e.printStackTrace();
			return null;
		}		
		
		return paths;
	}
	
	public static Map<String, String> readTable(String refdir)
	{
		String buildTablePath = refdir + File.separator + ConsoleMain.FILE_GBDAT;
		
		if (!FileBuffer.fileExists(buildTablePath)){
			System.err.println("WARNING: Genome build metadata file does not exist in " + refdir + " !");
			return new HashMap<String,String>();
		}
		
		Map<String, String> buildmap = new HashMap<String, String>();
		
		try
		{
			FileReader fr = new FileReader(buildTablePath);
			BufferedReader br = new BufferedReader(fr);
			
			String line = null;
			while ((line = br.readLine()) != null)
			{
				String[] fields = line.split("\t");
				if (fields.length < 2) continue;
				buildmap.put(fields[0], fields[1]);
			}
			
			br.close();
			fr.close();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		
		return buildmap;
	}
	
	public static Map<String, String> readRefSeqTable(String refdir)
	{
		String buildTablePath = refdir + File.separator + ConsoleMain.FILE_GBDAT;
		
		if (!FileBuffer.fileExists(buildTablePath)){
			System.err.println("WARNING: Genome build metadata file does not exist in " + refdir + " !");
			return new HashMap<String,String>();
		}
		
		Map<String, String> buildmap = new HashMap<String, String>();
		
		try
		{
			FileReader fr = new FileReader(buildTablePath);
			BufferedReader br = new BufferedReader(fr);
			
			String line = null;
			while ((line = br.readLine()) != null)
			{
				String[] fields = line.split("\t");
				if (fields.length != 3) continue;
				buildmap.put(fields[0], fields[2]);
			}
			
			br.close();
			fr.close();
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		
		return buildmap;
	}
	
	public static int detectOS(boolean verbose)
	{
		String osname = System.getProperty("os.name");
		if (verbose) System.err.println("Operating System Detected: " + osname);
		osname = osname.toLowerCase();
		if (osname.indexOf("win") >= 0) return OS_WINDOWS;
		if (osname.indexOf("mac") >= 0) return OS_MAC;
		if (osname.indexOf("nix") >= 0 || osname.indexOf("nux") >= 0) return OS_LINUX;
		return -1;
	}
	
	public static String getUserHome(boolean verbose)
	{
		int os = detectOS(verbose);
		String username = System.getProperty("user.name");
		if (username == null || username.isEmpty()) return null;
		switch(os)
		{
		case OS_WINDOWS:
			if (verbose) System.err.println("Operating System Class: Windows");
			//I really hate putting stuff in My Documents. That's for -documents-, yo.
			//And I don't trust Windows to not make Documents the home directory.
			String dir = "C:\\Users\\" + username;
			if (verbose) System.err.println("Trying directory " + dir + " ...");
			if (FileBuffer.directoryExists(dir)) return dir;
			if (verbose) System.err.println("Default Windows user directory not found.");
			dir = System.getProperty("user.home");
			if (verbose) System.err.println("Trying directory " + dir + " ...");
			if (FileBuffer.directoryExists(dir)) return dir;
			else
			{
				if (verbose) System.err.println("ERROR: Home directory could not be found!");
				return null;
			}
		case OS_MAC:
			if (verbose) System.err.println("Operating System Class: Macintosh");
			String mdir = System.getProperty("user.home");
			if (verbose) System.err.println("Trying directory " + mdir + " ...");
			if (FileBuffer.directoryExists(mdir)) return mdir;
			else
			{
				if (verbose) System.err.println("ERROR: Home directory could not be found!");
				return null;
			}
		case OS_LINUX:
			if (verbose) System.err.println("Operating System Class: Unix");
			String udir = System.getProperty("user.home");
			if (verbose) System.err.println("Trying directory " + udir + " ...");
			if (FileBuffer.directoryExists(udir)) return udir;
			else
			{
				if (verbose) System.err.println("ERROR: Home directory could not be found!");
				return null;
			}
		}
		if (verbose) System.err.println("Operating system not recognized! Will attempt to retrieve home directory path anyway...");
		String dir = System.getProperty("user.home");
		if (verbose) System.err.println("Trying directory " + dir + " ...");
		if (FileBuffer.directoryExists(dir)) return dir;
		else
		{
			if (verbose) System.err.println("ERROR: Home directory could not be found!");
			return null;
		}
	}
	
	public static void checkInstall(String homedir)
	{
		List<String> rpaths = getRefPaths(homedir);
		if (rpaths == null)
		{
			System.out.println("NOINSTALL");
			System.exit(0);
		}
		if (rpaths.isEmpty())
		{
			System.out.println("NOINSTALL");
			System.exit(0);
		}
		//System.out.println("YES");
		for (String p : rpaths)
		{
			System.out.println(p);
		}
		System.exit(0);
				
	}
	
	public static GenomeBuild loadBuild(String homedir, String name, boolean verbose)
	{
		if (name == null || name.isEmpty())
		{
			System.err.println("ERROR: Name of genome build is required for this tool!");
			printUsage();
			System.exit(1);
		}
		
		if (verbose) System.err.println("Attempting to load information for genome: \"" + name + "\" ...");
		List<String> refpaths = getRefPaths(homedir);
		GenomeBuild gb = null;
		String path = null;
		
		for (String p : refpaths)
		{
			Map<String, String> gbtable = readTable(p);
			Set<String> keys = gbtable.keySet();
			String file = null;
			for (String k : keys)
			{
				if (name.equalsIgnoreCase(k)) file = gbtable.get(k);
			}
			if (file != null){
				path = p + File.separator + file;
				break;
			}
		}
		
		if (path == null)
		{
			System.err.println("ERROR: Genome build \"" + name + "\" could not be found!");
			printUsage();
			System.exit(1);
		}
		
		try 
		{
			gb = new GenomeBuild(path);
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: Genome build metadata file \"" + path + "\" could not be opened!");
			e.printStackTrace();
			System.exit(1);
		} 
		catch (UnsupportedFileTypeException e) 
		{
			System.err.println("ERROR: Genome build metadata file \"" + path + "\" could not be parsed!");
			e.printStackTrace();
			System.exit(1);
		}
		
		return gb;
	}
	
	public static GeneSet loadRefSeq(String homedir, GenomeBuild build, boolean verbose)
	{
		if (build == null)
		{
			System.err.println("ERROR: Genome build is required to load refGene!");
			printUsage();
			System.exit(1);
		}
		
		String buildname = build.getBuildName();
		
		if (verbose) System.err.println("Attempting to load refGene database for genome: \"" + buildname + "\" ...");
		List<String> refpaths = getRefPaths(homedir);
		GeneSet gs = null;
		String path = null;
		
		for (String p : refpaths)
		{
			Map<String, String> gstable = readRefSeqTable(p);
			Set<String> keys = gstable.keySet();
			String file = null;
			for (String k : keys)
			{
				if (buildname.equalsIgnoreCase(k)) file = gstable.get(k);
			}
			if (file != null){
				path = p + File.separator + file;
				break;
			}
		}
		
		if (path == null)
		{
			System.err.println("ERROR: RefSeq data for \"" + buildname + "\" could not be found!");
			printUsage();
			System.exit(1);
		}
		
		try 
		{
			gs = new GeneSet(path, build, true);
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: RefSeq dataset file \"" + path + "\" could not be opened!");
			e.printStackTrace();
			System.exit(1);
		} 
		catch (UnsupportedFileTypeException e) 
		{
			System.err.println("ERROR: RefSeq dataset file \"" + path + "\" could not be parsed!");
			e.printStackTrace();
			System.exit(1);
		}
		
		return gs;
	}
	
	public static void main(String[] args) 
	{
		String program = null;
		String genome = null;
		boolean verbose = false;
		
		//Get requested program name
		if (args.length < 1)
		{
			System.err.println("ERROR: No program name was provided!");
			printUsage();
			System.exit(1);
		}
		program = args[0];
		
		//Get other arguments: genome build and verbose flag
		for (int i = 0; i < args.length; i++)
		{
			String s = args[i];
			if (s.equals(OP_GENOMEBUILD))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_GENOMEBUILD + " flag MUST be followed by name of reference genome build!");
					printUsage();
					System.exit(1);
				}
				genome = args[i+1];
			}
			if (s.equals(OP_VERBOSE)) verbose = true;
		}
		
		//Get homedir
		String homedir = getUserHome(verbose);
		
		//Try to load genome build if needed
		if (program.equals(TOOL_CHECKINSTALL))
		{
			if (verbose){
				System.err.println("Tool Selected: Check Install...");
				System.err.println();
			}
			checkInstall(homedir);
		}
		else if (program.equals(TOOL_FILTERLUMPY))
		{
			if (verbose){
				System.err.println("Tool Selected: Filter Lumpy...");
				System.err.println();
			}
			GenomeBuild gb = loadBuild(homedir, genome, verbose);
			if(gb == null)
			{
				System.err.println("Genome build for \"" + genome + "\" could not be loaded!");
				printUsage();
				System.exit(1);
			}
			LumpyFilter.lumpyFilter(args, gb);
		}
		else if (program.equals(TOOL_GBMANAGER))
		{
			if (verbose){
				System.err.println("Tool Selected: Genome Build Manager");
				System.err.println();
			}
			GenomeManager.genomeManager(args, homedir);
		}
		else if (program.equals(TOOL_INSTALLME))
		{
			if (verbose){
				System.err.println("Tool Selected: Install...");
				System.err.println();
			}
			Install.installBIOISVTOOLS(args, homedir);
		}
		else if (program.equals(TOOL_PENNTOVCF))
		{
			if (verbose){
				System.err.println("Tool Selected: PennToVCF...");
				System.err.println();
			}
			GenomeBuild gb = loadBuild(homedir, genome, verbose);
			if(gb == null)
			{
				System.err.println("Genome build for \"" + genome + "\" could not be loaded!");
				printUsage();
				System.exit(1);
			}
			PennToVCF.pennToVCF(args, gb, verbose);
		}
		else if (program.equals(TOOL_TRACKMAKER))
		{
			if (verbose){
				System.err.println("Tool Selected: TrackMaker...");
				System.err.println();
			}
			GenomeBuild gb = loadBuild(homedir, genome, verbose);
			if(gb == null)
			{
				System.err.println("Genome build for \"" + genome + "\" could not be loaded!");
				printUsage();
				System.exit(1);
			}
			TrackMaker.trackMaker(args, gb);
		}
		else if (program.equals(TOOL_MERGEGENO))
		{
			if (verbose){
				System.err.println("Tool Selected: SURVIVOR genotype merger...");
				System.err.println();
			}
			GenomeBuild gb = loadBuild(homedir, genome, verbose);
			if(gb == null)
			{
				System.err.println("Genome build for \"" + genome + "\" could not be loaded!");
				printUsage();
				System.exit(1);
			}
			ConsensusGeno.scgeno(args, gb);
		}
		else if (program.equals(TOOL_STDCHROM))
		{
			if (verbose){
				System.err.println("Tool Selected: stdchrom");
				System.err.println();
			}
			GenomeBuild gb = loadBuild(homedir, genome, verbose);
			if(gb == null)
			{
				System.err.println("Genome build for \"" + genome + "\" could not be loaded!");
				printUsage();
				System.exit(1);
			}
			Stdchrom.standardizeChrom(args, gb);
		}
		else if (program.equals(TOOL_SVANNO))
		{
			if (verbose){
				System.err.println("Tool Selected: svanno");
				System.err.println();
			}
			GenomeBuild gb = loadBuild(homedir, genome, verbose);
			if(gb == null)
			{
				System.err.println("Genome build for \"" + genome + "\" could not be loaded!");
				printUsage();
				System.exit(1);
			}
			GeneSet gs = loadRefSeq(homedir, gb, verbose);
			if (gs == null)
			{
				System.err.println("RefGene database could not be loaded for " + genome + "!");
				System.exit(1);
			}
			SVAnno.svanno(args, gb, gs, verbose);
		}
		else if (program.equals(TOOL_SVQSTATS))
		{
			if (verbose){
				System.err.println("Tool Selected: " + TOOL_SVQSTATS);
				System.err.println();
			}
			GenomeBuild gb = loadBuild(homedir, genome, verbose);
			if(gb == null)
			{
				System.err.println("Genome build for \"" + genome + "\" could not be loaded!");
				printUsage();
				System.exit(1);
			}
			VarCounter.svQuickStats(args, gb);
		}
		else if (program.equals(TOOL_SZTALLY))
		{
			if (verbose){
				System.err.println("Tool Selected: " + TOOL_SZTALLY);
				System.err.println();
			}
			VarSizes.varsizes(args);
		}
		else if (program.equals(TOOL_VCFTRIM))
		{
			if (verbose){
				System.err.println("Tool Selected: " + TOOL_VCFTRIM);
				System.err.println();
			}
			VCFTrim.trimVCF(args);
		}
		else if (program.equals(TOOL_VCFSAMPLENAMESWAP))
		{
			if (verbose){
				System.err.println("Tool Selected: " + TOOL_VCFSAMPLENAMESWAP);
				System.err.println();
			}
			VCFReheader.reheader(args);
		}
		else if (program.equals(TOOL_SZFILTER))
		{
			if (verbose){
				System.err.println("Tool Selected: " + TOOL_SZFILTER);
				System.err.println();
			}
			GenomeBuild gb = loadBuild(homedir, genome, verbose);
			if(gb == null)
			{
				System.err.println("Genome build for \"" + genome + "\" could not be loaded!");
				printUsage();
				System.exit(1);
			}
			SizeFilter.filterBySize(args, gb);
		}
		else if (program.equals(TOOL_SPLITMERGE))
		{
			if (verbose){
				System.err.println("Tool Selected: " + TOOL_SPLITMERGE);
				System.err.println();
			}
			GenomeBuild gb = loadBuild(homedir, genome, verbose);
			if(gb == null)
			{
				System.err.println("Genome build for \"" + genome + "\" could not be loaded!");
				printUsage();
				System.exit(1);
			}
			SplitMerge.splitmerge(args, gb);
		}
		else if (program.equals(TOOL_DELLYCLEANER))
		{
			if (verbose){
				System.err.println("Tool Selected: " + TOOL_DELLYCLEANER);
				System.err.println();
			}
			DellyCleaner.runTool(args, verbose);
		}
		else if (program.equals(TOOL_TRIMCHROM))
		{
			if (verbose){
				System.err.println("Tool Selected: " + TOOL_TRIMCHROM);
				System.err.println();
			}
			GenomeBuild gb = loadBuild(homedir, genome, verbose);
			if(gb == null)
			{
				System.err.println("Genome build for \"" + genome + "\" could not be loaded!");
				printUsage();
				System.exit(1);
			}
			ChromTrim.runTool(args, gb, verbose);
		}
		else if (program.equals(TOOL_SVANALYZE))
		{
			if (verbose){
				System.err.println("Tool Selected: " + TOOL_SVANALYZE);
				System.err.println();
			}
			GenomeBuild gb = loadBuild(homedir, genome, verbose);
			if(gb == null)
			{
				System.err.println("Genome build for \"" + genome + "\" could not be loaded!");
				printUsage();
				System.exit(1);
			}
			SVFam.RunSVFamily(args, gb);
		}
		else if (program.equals(TOOL_SVGENEHITTALLY))
		{
			if (verbose){
				System.err.println("Tool Selected: " + TOOL_SVGENEHITTALLY);
				System.err.println();
			}
			GeneTally.runGeneTally(args);
		}
		else if (program.equals(TOOL_SCANSAM))
		{
			if (verbose){
				System.err.println("Tool Selected: " + TOOL_SCANSAM);
				System.err.println();
			}
			GenomeBuild gb = loadBuild(homedir, genome, verbose);
			SamScanner.runSamScanner(args, gb, verbose);
		}
		else if (program.equals(TOOL_FIXSAM))
		{
			if (verbose){
				System.err.println("Tool Selected: " + TOOL_FIXSAM);
				System.err.println();
			}
			GenomeBuild gb = loadBuild(homedir, genome, verbose);
			if(gb == null)
			{
				System.err.println("Genome build for \"" + genome + "\" could not be loaded!");
				printUsage();
				System.exit(1);
			}
			SamFixer.runSamFixer(args, gb, verbose);
		}
		else if (program.equals(TOOL_PEDTOFAMI))
		{
			if (verbose){
				System.err.println("Tool Selected: " + TOOL_PEDTOFAMI);
				System.err.println();
			}
			Ped2Fami.runPed2Fami(args);
		}
		else if (program.equals(TOOL_VIEWFAMI))
		{
			if (verbose){
				System.err.println("Tool Selected: " + TOOL_VIEWFAMI);
				System.err.println();
			}
			ViewFami.runViewFami(args);
		}
		else if (program.equals(TOOL_SVDB))
		{
			if (verbose){
				System.err.println("Tool Selected: " + TOOL_SVDB);
				System.err.println();
			}
			GenomeBuild gb = loadBuild(homedir, genome, verbose);
			if(gb == null)
			{
				System.err.println("Genome build for \"" + genome + "\" could not be loaded!");
				printUsage();
				System.exit(1);
			}
			ConsoleFrontEnd.runSVDB(args, gb, verbose);
		}
		else if (program.equals(TOOL_LRGMRG))
		{
			if (verbose){
				System.err.println("Tool Selected: " + TOOL_LRGMRG);
				System.err.println();
			}
			GenomeBuild gb = loadBuild(homedir, genome, verbose);
			if(gb == null)
			{
				System.err.println("Genome build for \"" + genome + "\" could not be loaded!");
				printUsage();
				System.exit(1);
			}
			LargeMerger.runMerger(args, gb);
		}
		else
		{
			System.err.println("ERROR: \"" + program + "\" is not a valid tool!");
			printUsage();
			System.exit(1);
		}
		
	}

}
