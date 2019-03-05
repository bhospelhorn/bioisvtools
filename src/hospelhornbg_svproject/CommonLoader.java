package hospelhornbg_svproject;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import hospelhornbg_genomeBuild.Gene;
import hospelhornbg_genomeBuild.GeneSet;
import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class CommonLoader {
	
	/* --- Constants --- */
	
	public static final int OS_WINDOWS = 0;
	public static final int OS_LINUX = 1;
	public static final int OS_MAC = 2;
	
	public static final String HOME_INSTALL_DIR = "bioisvtools";
	public static final String HOME_INSTALL_FILE = "svproject.cfg";
	
	public static final String OMIM_FILE_NAME = "omim_whitelist";
	public static final String CALLER_TABLE_NAME = "suppcallers.tsv";
	public static final String GENOMES_DIR = "genomes";
	public static final String GB_PATH_TABLE = "genomes.tsv";
	public static final String GS_PATH_TABLE = "refseqs.tsv";
	
	public static final String PROJECT_DIR = "projects";
	public static final String COHORT_TABLE = "cohortNames.tsv";
	public static final String COHORT_PROJECT_LIST_STEM = "cohort_";
	
	public static final String BLACKLIST_FILE_LC = "BL_LowComplexity";
	public static final String BLACKLIST_FILE_VT = "BL_HighlyVariable";
	public static final String BLACKLIST_FILE_SG = "BL_LargeFamily";
	public static final String BLACKLIST_FILE_PG = "BL_Pseudogenes";
	public static final String WHITELIST_FILE = "WL_UserWhitelist";
	public static final String GREYLIST_FILE_EXT = ".txt";
	
	public static final String COMMON_DIR_FIELD_KEY = "COMMONDIR";
	public static final String USER_DIR_FIELD_KEY = "USERDIR";
	
	public static final String OMIM_URL = "https://www.omim.org/static/omim/data/mim2gene.txt";
	
	/* --- Static Variables --- */
	
	private static String common_dir;
	private static String user_dir;
	
	private static Map<Integer, String> genomebuild_paths;
	private static Map<Integer, String> geneset_paths;
	
	private static Map<Integer, String> cohort_UID_map;
	
	
	/* --- Home Paths --- */
	
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
			if (verbose) {
				System.err.println("Operating System Class: Macintosh");
				System.err.println("You should be ashamed of yourself.");
			}
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

	/* --- Read Program Files --- */
	
	public static boolean loadProgramData(boolean verbose) throws IOException
	{
		if(!loadMainDirNames(verbose)) return false;
		loadGenomePaths(verbose);
		loadCohortIDs(verbose);
		return true;
	}
	
	private static boolean loadMainDirNames(boolean verbose) throws IOException
	{
		String home = getUserHome(verbose);
		String hdir = home + File.separator + HOME_INSTALL_DIR;
		if (!FileBuffer.directoryExists(hdir))
		{
			if(verbose) System.err.println("CommonLoader.loadMainDirNames || ERROR: Installation path record does not exist!");
			return false;
		}
		
		String cfgPath = hdir + File.separator + HOME_INSTALL_FILE;
		if (!FileBuffer.fileExists(cfgPath))
		{
			if(verbose) System.err.println("CommonLoader.loadMainDirNames || ERROR: Installation path record does not exist!");
			return false;
		}
		
		//Read
		Map<String, String> cfgmap = new HashMap<String, String>();
		
		BufferedReader br = new BufferedReader(new FileReader(cfgPath));
		String line = null;
		while((line = br.readLine()) != null)
		{
			String[] split = line.split("=");
			if (split.length != 2) continue;
			cfgmap.put(split[0], split[1]);
		}
		br.close();
		
		//Get the fields of interest
		common_dir = cfgmap.get(COMMON_DIR_FIELD_KEY);
		user_dir = cfgmap.get(USER_DIR_FIELD_KEY);
		return true;
	}

	private static boolean loadGenomePaths(boolean verbose) throws IOException
	{
		if (common_dir == null || common_dir.isEmpty())
		{
			if(verbose) System.err.println("CommonLoader.loadGenomePaths || ERROR: No common data directory found!");
			return false;
		}
		String genomedir = common_dir + File.separator + GENOMES_DIR;
		if (!FileBuffer.directoryExists(genomedir))
		{
			if(verbose) System.err.println("CommonLoader.loadGenomePaths || ERROR: Genome data directory not found!");
			return false;
		}
		
		//Read tables...
		String tblpath = genomedir + File.separator + GB_PATH_TABLE;
		if (!FileBuffer.fileExists(tblpath))
		{
			if(verbose) System.err.println("CommonLoader.loadGenomePaths || ERROR: Genome Build path table not found!");
			return false;
		}
		genomebuild_paths = new HashMap<Integer, String>();
		
		BufferedReader br = new BufferedReader(new FileReader(tblpath));
		String line = null;
		while((line = br.readLine()) != null)
		{
			String[] split = line.split("\t");
			if (split.length != 2) continue;
			//Parse split[0]
			try
			{
				int uid = Integer.parseInt(split[0]);
				genomebuild_paths.put(uid, split[1]);
			}
			catch (NumberFormatException e)
			{
				if(verbose) System.err.println("CommonLoader.loadGenomePaths || ERROR: Genome Build path table could not be read!");
				br.close();
				return false;
			}
		}
		br.close();
		
		//Now the other table...
		tblpath = genomedir + File.separator + GS_PATH_TABLE;
		if (!FileBuffer.fileExists(tblpath))
		{
			if(verbose) System.err.println("CommonLoader.loadGenomePaths || ERROR: Gene Set path table not found!");
			return false;
		}
		geneset_paths = new HashMap<Integer, String>();
		
		br = new BufferedReader(new FileReader(tblpath));
		line = null;
		while((line = br.readLine()) != null)
		{
			String[] split = line.split("\t");
			if (split.length != 2) continue;
			try
			{
				int uid = Integer.parseInt(split[0]);
				geneset_paths.put(uid, split[1]);
			}
			catch (NumberFormatException e)
			{
				if(verbose) System.err.println("CommonLoader.loadGenomePaths || ERROR: Genome Set path table could not be read!");
				br.close();
				return false;
			}
		}
		br.close();
		
		return true;
	}
	
	private static void loadCohortIDs(boolean verbose) throws IOException
	{
		cohort_UID_map = new HashMap<Integer, String>();
		if (user_dir == null || user_dir.isEmpty())
		{
			if(verbose) System.err.println("CommonLoader.loadCohortIDs || ERROR: No user data directory found!");
			return;
		}
		
		String ctable = user_dir + File.separator + COHORT_TABLE;
		if (!FileBuffer.fileExists(ctable))
		{
			if(verbose) System.err.println("CommonLoader.loadCohortIDs || ERROR: Cohort ID table not found!");
			return;
		}
		
		BufferedReader br = new BufferedReader(new FileReader(ctable));
		String line = null;
		while((line = br.readLine()) != null)
		{
			String[] split = line.split("\t");
			if (split.length != 2) continue;
			//Parse split[0]
			try
			{
				int uid = Integer.parseInt(split[0]);
				cohort_UID_map.put(uid, split[1]);
			}
			catch (NumberFormatException e)
			{
				if(verbose) System.err.println("CommonLoader.loadCohortIDs || ERROR: Cohort ID path table could not be read!");
				br.close();
				return;
			}
		}
		br.close();
		
	}
	
	/* --- Static Getters --- */
	
	public static String getCommonDirPath()
	{
		return common_dir;
	}
	
	public static String getUserDirPath()
	{
		return user_dir;
	}
	
	public static String generateGreylistFilePath(String fncore, int gbUID)
	{
		return user_dir + File.separator + fncore + "_" + Integer.toHexString(gbUID) + GREYLIST_FILE_EXT;
	}
	
	public static String generateOMIMWhitelistFilePath(int gbUID)
	{
		return common_dir + File.separator + OMIM_FILE_NAME + "_" + Integer.toHexString(gbUID) + GREYLIST_FILE_EXT;
	}
	
	public static String getUserProjectDirPath()
	{
		return user_dir + File.separator + PROJECT_DIR;
	}
	
	/* --- GenomeBuild & GeneSet Loaders --- */
	
	public static GenomeBuild loadGenomeBuild(int gbUID) throws IOException, UnsupportedFileTypeException
	{
		String gbPath = genomebuild_paths.get(gbUID);
		if (gbPath == null) return null;
		
		GenomeBuild gb = new GenomeBuild(gbPath);
		
		return gb;
	}
	
	public static GeneSet loadGeneSet(GenomeBuild build, boolean verbose) throws UnsupportedFileTypeException, IOException
	{
		if (build == null) return null;
		int uid = -1;
		if(build.getUIDEnum() == null) uid = build.getBuildName().hashCode();
		else uid = build.getUIDEnum().getUID();
		String gsPath = geneset_paths.get(uid);
		if (gsPath == null) return null;
		
		GeneSet gs = new GeneSet(gsPath, build, true);
		
		//Load lists...
		String listpath = generateGreylistFilePath(BLACKLIST_FILE_LC, uid);
		if(!FileBuffer.fileExists(listpath))
		{
			if(verbose) System.err.println("CommonLoader.loadGeneSet || WARNING: Low Complexity blacklist not found!");
		}
		else gs.importLowComplexityGreylist(listpath, verbose);
		
		listpath = generateGreylistFilePath(BLACKLIST_FILE_VT, uid);
		if(!FileBuffer.fileExists(listpath))
		{
			if(verbose) System.err.println("CommonLoader.loadGeneSet || WARNING: Variability Tolerant blacklist not found!");
		}
		else gs.importHighlyVariableGreylist(listpath, verbose);
		
		listpath = generateGreylistFilePath(BLACKLIST_FILE_SG, uid);
		if(!FileBuffer.fileExists(listpath))
		{
			if(verbose) System.err.println("CommonLoader.loadGeneSet || WARNING: Large Family blacklist not found!");
		}
		else gs.importManySimilarGreylist(listpath, verbose);
		
		listpath = generateGreylistFilePath(BLACKLIST_FILE_PG, uid);
		if(!FileBuffer.fileExists(listpath))
		{
			if(verbose) System.err.println("CommonLoader.loadGeneSet || WARNING: Pseudogene blacklist not found!");
		}
		else gs.importPseudogeneGreylist(listpath, verbose);
		
		listpath = generateGreylistFilePath(WHITELIST_FILE, uid);
		if(!FileBuffer.fileExists(listpath))
		{
			if(verbose) System.err.println("CommonLoader.loadGeneSet || WARNING: User whitelist not found!");
		}
		else gs.importWhitelist(listpath, verbose);
		
		listpath = generateOMIMWhitelistFilePath(uid);
		if(!FileBuffer.fileExists(listpath))
		{
			if(verbose) System.err.println("CommonLoader.loadGeneSet || WARNING: OMIM whitelist not found!");
		}
		else gs.importOMIMWhitelist(listpath, verbose);
		
		return gs;
	}
	
	/* --- Other Loaders/Getters --- */
	
	public static List<String> getCohortProjectPaths(int cohort_uid, boolean verbose) throws IOException
	{
		List<String> list = new LinkedList<String>();
		//int cuid = cohort_UID_map.get(cohort_name);
		
		String cohortListFile = user_dir + File.separator +  COHORT_PROJECT_LIST_STEM + Integer.toHexString(cohort_uid) + ".txt";
		if (!FileBuffer.fileExists(cohortListFile))
		{
			System.err.println("CommonLoader.getCohortProjectPaths || No previous record for cohort 0x" + Integer.toHexString(cohort_uid) + " (" + cohort_UID_map.get(cohort_uid) + ") found!");
			return list;
		}
		
		BufferedReader br = new BufferedReader(new FileReader(cohortListFile));
		String line = null;
		while((line = br.readLine()) != null)
		{
			list.add(line);
		}
		br.close();
		
		return list;
	}
	
	public static List<String> getSupportStrings_ordered()
	{
		//TODO: Write
		return null;
	}
	
	public static String getProjectDirectory(String projectName)
	{
		if (user_dir == null) return null;
		return user_dir + File.separator + PROJECT_DIR + File.separator + projectName;
	}
	
	/* --- Install --- */
	
	public static void installSVProject(String commonDir, String userDir, boolean verbose) throws IOException
	{
		String homedir = getUserHome(verbose);
		common_dir = commonDir;
		user_dir = userDir;
		
		//Make cfg file and directory...
		String cfgPath = homedir + File.separator + HOME_INSTALL_DIR;
		if (!FileBuffer.directoryExists(cfgPath))
		{
			Files.createDirectories(Paths.get(cfgPath));
		}
		
		//Overwrite any other cfg that happens to be there...
		cfgPath += File.separator + HOME_INSTALL_FILE;
		BufferedWriter bw = new BufferedWriter(new FileWriter(cfgPath));
		bw.write(COMMON_DIR_FIELD_KEY + "=" + common_dir + "\n");
		bw.write(USER_DIR_FIELD_KEY + "=" + user_dir + "\n");
		bw.close();
		
		//Check if the common and user dirs exist. If not, make.
		String subdir = common_dir + File.separator + GENOMES_DIR;
		if (!FileBuffer.directoryExists(subdir))
		{
			Files.createDirectories(Paths.get(subdir));
		}
		
		subdir = user_dir + File.separator + PROJECT_DIR;
		if (!FileBuffer.directoryExists(subdir))
		{
			Files.createDirectories(Paths.get(subdir));
		}
		
	}
	
	public static void installGenome(GenomeBuild gb, GeneSet gs, boolean verbose)
	{
		//TODO: Write
		//This will overwrite existing files in the common dir. Use with caution.
	}
	
	private static void downloadOMIMList(String targetPath, OMIMUpdateListener l) throws IOException
	{
		URL txturl = new URL(OMIM_URL);
		System.setProperty("http.agent", "Chrome");
		InputStream is = txturl.openStream();
		BufferedWriter bw = new BufferedWriter(new FileWriter(targetPath));
		
		int b = -1;
		int counter = 0;
		int sz = is.available();
		while ((b = is.read()) >= 0)
		{
			bw.write(b);
			counter++;
			if (counter % 1000 == 0 && l != null)
			{
				l.onDownloadProgressUpdate(counter, sz);
			}
		}
		
		is.close();
		bw.close();
	}
	
	public static void updateOMIMList(int gbUID, OMIMUpdateListener l)
	{
		//TODO: Write
	}
	
	/* --- Write --- */
	
	public static void newCohort(List<String> projectpaths)
	{
		//TODO: Write
	}
	
	public static void writeLowComplexity_Blacklist(GeneSet genes, int gbUID)
	{
		//TODO: Write
	}
	
	public static void writeHighlyVariable_Blacklist(List<Gene> flaggedGenes, int gbUID)
	{
		//TODO: Write
	}
	
	public static void writeManySimilar_Blacklist(List<Gene> flaggedGenes, int gbUID)
	{
		//TODO: Write
	}
	
	public static void writePseudogene_Blacklist(List<Gene> flaggedGenes, int gbUID)
	{
		//TODO: Write
	}
	
	public static void writeUserWhitelist(List<Gene> flaggedGenes, int gbUID)
	{
		//TODO: Write
	}
	
}
