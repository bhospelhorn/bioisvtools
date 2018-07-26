package hospelhornbg_svtools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;


public class GenomeManager {
	
	public static final String MODE_ADD = "add";
	public static final String MODE_ALIAS = "alias";
	public static final String MODE_REMOVE = "remove";
	public static final String MODE_LIST = "list";
	
	public static final String OP_MODE = "-m";
	public static final String OP_KEY = "-n";
	public static final String OP_REF = "-r";
	public static final String OP_PATH = "-p";
	
	//Automatically cleans init file of directories that don't work
	
	public static void printUsage()
	{
		System.out.println("--------------------------------------------------------------------------------");
		System.out.println("BioisvTools || Genome Build Metadata Manager");
		System.out.println();
		System.out.println("Purpose: To manage paths and keys to genome build metadata files without directly editing the table.");
		System.out.println();
		System.out.println("Modes:");
		System.out.println("\tadd\tAdd a new genome build information file to available references.");
		System.out.println("\tremove\tRemove a genome build key/path pair from list of available references.");
		System.out.println("\talias\tAdd a new alias key for an existing genome reference.");
		System.out.println("\tlist\tList available genome build name keys.");
		System.out.println();
		System.out.println("Flags:");
		System.out.println("\t-m\tSTRING\t[Required]\t\tMode to call.");
		System.out.println("\t-n\tSTRING\t[Required]\t\tName associated with path link to edit, or new alias.");
		System.out.println("\t-r\tSTRING\t[Conditional]\t\tFor \"alias\" mode: known name associated with path to add alias to.");
		System.out.println("\t-p\tFILE\t[Conditional]\t\tPath of genome build metadata file to copy to references and link with key.");
		System.out.println();
		System.out.println("Sample Usage:");
		System.out.println("java -jar bioisvtools.jar gbmanager -m add -n GRCh37 -p /somedirs/refs/GRCh37.gbdh");
		System.out.println("java -jar bioisvtools.jar gbmanager -m alias -n hg19 -r GRCh37");
		System.out.println("java -jar bioisvtools.jar gbmanager -m remove -n hg19");
		System.out.println("java -jar bioisvtools.jar gbmanager -m list");
		System.out.println();
		System.out.println("--------------------------------------------------------------------------------");
	}
	
	public static void add(String refdir, String key, String path)
	{
		Map <String, String> nowtable = ConsoleMain.readTable(refdir);
		
		//Determine file's name
		int sep = path.lastIndexOf(File.separator);
		String fname = "";
		if (sep >= 0) fname = path.substring(sep + 1);
		else fname = path;
		
		//Copy the file at path to the ref dir
		String newpath = refdir + File.separator + fname;
		try 
		{
			Files.copy(Paths.get(path), Paths.get(newpath));
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: Could not copy " + path + " to reference directory!");
			e.printStackTrace();
			System.exit(1);
		}
		
		//Place in map
		nowtable.put(key, fname);
		
		//Delete current gb dat file, if it exists.
		//Write out everything in the map to the gb dat file.
		String gbdat = refdir + File.separator + ConsoleMain.FILE_GBDAT;
		try 
		{
			Files.deleteIfExists(Paths.get(gbdat));
			
			FileWriter fw = new FileWriter(gbdat);
			BufferedWriter bw = new BufferedWriter(fw);
			
			Set<String> keyset = nowtable.keySet();
			int c = 1;
			for (String s : keyset)
			{
				bw.write(s + "\t" + nowtable.get(s));
				if (c < keyset.size()) bw.write("\n");
				c++;
			}
				
			bw.close();
			fw.close();
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: There was an error updating the metadata table!");
			e.printStackTrace();
			System.exit(1);
		}
	
		
	}
	
	public static void alias(String refdir, String alias, String refkey)
	{
		Map <String, String> nowtable = ConsoleMain.readTable(refdir);
		
		//Get original file name
		String fname = nowtable.get(refkey);
		
		if (fname == null || fname.isEmpty())
		{
			System.err.println("ERROR: Provided reference key does not refer to any existing resource!");
			System.exit(1);
		}
		
		//Place in map with new alias
		nowtable.put(alias, fname);
		
		//Delete current gb dat file, if it exists.
		//Write out everything in the map to the gb dat file.
		String gbdat = refdir + File.separator + ConsoleMain.FILE_GBDAT;
		try 
		{
			Files.deleteIfExists(Paths.get(gbdat));
			
			FileWriter fw = new FileWriter(gbdat);
			BufferedWriter bw = new BufferedWriter(fw);
			
			Set<String> keyset = nowtable.keySet();
			int c = 1;
			for (String s : keyset)
			{
				bw.write(s + "\t" + nowtable.get(s));
				if (c < keyset.size()) bw.write("\n");
				c++;
			}
				
			bw.close();
			fw.close();
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: There was an error updating the metadata table!");
			e.printStackTrace();
			System.exit(1);
		}
	
	}
	
	public static void remove(String refdir, String key)
	{
		Map <String, String> nowtable = ConsoleMain.readTable(refdir);
		
		//Remove targeted entry from map
		nowtable.remove(key);
		
		//Delete current gb dat file, if it exists.
		//Write out everything in the map to the gb dat file.
		String gbdat = refdir + File.separator + ConsoleMain.FILE_GBDAT;
		try 
		{
			Files.deleteIfExists(Paths.get(gbdat));
			
			FileWriter fw = new FileWriter(gbdat);
			BufferedWriter bw = new BufferedWriter(fw);
			
			Set<String> keyset = nowtable.keySet();
			int c = 1;
			for (String s : keyset)
			{
				bw.write(s + "\t" + nowtable.get(s));
				if (c < keyset.size()) bw.write("\n");
				c++;
			}
				
			bw.close();
			fw.close();
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: There was an error updating the metadata table!");
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	public static void list(String homedir)
	{
		ConsoleMain.cleanDeadRefPaths(homedir);
		List<String> reflist = ConsoleMain.getRefPaths(homedir);
		
		Set<String> keys = new HashSet<String>();
		for (String refpath : reflist)
		{
			Map<String, String> reftbl = ConsoleMain.readTable(refpath);
			keys.addAll(reftbl.keySet());
		}
		
		List<String> sortedkeys = new ArrayList<String>(keys.size());
		sortedkeys.addAll(keys);
		Collections.sort(sortedkeys);
		
		for (String k : sortedkeys) System.out.println(k);
	}
	
	public static void genomeManager(String[] args, String homedir)
	{
		String mode = null;
		String key = null;
		String refname = null;
		String path = null;
		
		for (int i = 0; i < args.length; i++)
		{
			String s = args[i];
			if (s.equals(OP_MODE))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_MODE + " flag MUST be followed by valid mode name!");
					printUsage();
					System.exit(1);
				}
				mode = args[i+1];
			}
			if (s.equals(OP_KEY))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_KEY + " flag MUST be followed by name to set for genome build!");
					printUsage();
					System.exit(1);
				}
				key = args[i+1];
			}
			if (s.equals(OP_REF))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_REF + " flag MUST be followed by name of genome build to make alias for!");
					printUsage();
					System.exit(1);
				}
				refname = args[i+1];
			}
			if (s.equals(OP_PATH))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_PATH + " flag MUST be followed by path of build table to add to references!");
					printUsage();
					System.exit(1);
				}
				path = args[i+1];
			}
		}
		
		//Check mode
		//Check args required for that mode
		
		if (mode == null || mode.isEmpty())
		{
			System.err.println("ERROR: Valid mode name is required!");
			printUsage();
			System.exit(1);
		}
		
		List<String> refpaths = ConsoleMain.getRefPaths(homedir);
		if (refpaths == null)
		{
			System.err.println("ERROR: INIT file could not be found! Please check to ensure that a reference path is installed for this user!");
			System.exit(1);
		}
		if (refpaths.isEmpty())
		{
			System.err.println("ERROR: INIT file contains no valid reference paths! Please reinstall for this user!");
			System.exit(1);
		}
		String mainpath = refpaths.get(0);
		
		if (mode.equals(MODE_ADD))
		{
			// Need key and path
			if (key == null || key.isEmpty())
			{
				System.err.println("ERROR: Desired key is required for add mode.");
				printUsage();
				System.exit(1);
			}
			if (path == null || path.isEmpty())
			{
				System.err.println("ERROR: Path to reference file is required for add mode.");
				printUsage();
				System.exit(1);
			}
			add(mainpath, key, path);
		}
		else if (mode.equals(MODE_ALIAS))
		{
			// Need key and refname
			if (key == null || key.isEmpty())
			{
				System.err.println("ERROR: Desired new key is required for alias mode.");
				printUsage();
				System.exit(1);
			}
			if (refname == null || refname.isEmpty())
			{
				System.err.println("ERROR: Exisiting reference key is required for alias mode.");
				printUsage();
				System.exit(1);
			}
			alias(mainpath, key, refname);
		}
		else if (mode.equals(MODE_REMOVE))
		{
			// Need key
			if (key == null || key.isEmpty())
			{
				System.err.println("ERROR: Existing key is required for remove mode.");
				printUsage();
				System.exit(1);
			}
			remove(mainpath, key);
		}
		else if (mode.equals(MODE_LIST))
		{
			list(homedir);
		}
		else
		{
			System.err.println("\"" + mode + "\" is not a recognized mode.");
			printUsage();
			System.exit(1);
		}
				
	}

}
