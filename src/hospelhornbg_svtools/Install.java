package hospelhornbg_svtools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import hospelhornbg_genomeBuild.GeneSet;
import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.FileBuffer;

public class Install {
	
	public static final String OP_REFPATH = "-d";
	
	//Look for init file in homedir/bioisvtools
	//		(or C:\Users\(You)\bioisvtools on Windows)
	//If there, tack this new path onto the end (look through list until hit valid directory
	//	or find desired reference)
	//If not, make new init file and write this path to it.
	
	//Check the directory specified by d for ref files contained in the jar.
	//Any ref files not there should be copied from the jar to this directory.
	//If there isn't permission to do so, program will exit with error.
	
	private static class BuildRecord
	{
		public String gbPath;
		public String refSeqPath;
	}
	
	public static void printUsage()
	{
		System.out.println("--------------------------------------------------------------------------------");
		System.out.println("BioisvTools || Install");
		System.out.println();
		System.out.println("Purpose: Locally \"install\" a reference file directory path for usage by BioisvTools.");
		System.out.println();
		System.out.println("Flags:");
		System.out.println("\t-d\tPATH\t[Required]\t\tPath to directory where program should look for reference materials.");
		System.out.println();
		System.out.println("Sample Usage:");
		System.out.println("java -jar bioisvtools.jar install -d /somedir/bioisvtools/reference");
		System.out.println();
		System.out.println("--------------------------------------------------------------------------------");
	}
	
	public static boolean copyReferences(String targetDir)
	{
		String buildTablePath = targetDir + File.separator + ConsoleMain.FILE_GBDAT;
		Map<String, BuildRecord> existing = new HashMap<String, BuildRecord>();
		
		//Check table and copy existing records
		if (FileBuffer.fileExists(buildTablePath))
		{
			try
			{
				FileReader fr = new FileReader(buildTablePath);	
				BufferedReader br = new BufferedReader(fr);
				
				String line = null;
				while ((line = br.readLine()) != null)
				{
					if (line.isEmpty()) continue;
					String[] fields = line.split("\t");
					if (fields.length != 3){
						//br.close();
						//fr.close();
						System.err.println("Build table path not formatted correctly! It's likely that you are using an old version... Now reinstalling...");
						//return false;
						break;
					}
					else
					{
						BuildRecord rec = new BuildRecord();
						rec.gbPath = fields[1];
						rec.refSeqPath = fields[2];
						existing.put(fields[0], rec);	
					}
				}
				
				br.close();
				fr.close();
				
				Files.deleteIfExists(Paths.get(buildTablePath));
			}
			catch (IOException e)
			{
				System.err.println("ERROR: Genome build table registered as existing, but could not be read or modified!");
				e.printStackTrace();
				return false;
			}
		}
		
		//Compile set of JAR resources
		Set<String> jarkeys = GenomeBuild.getAllStandardBuildKeys();
		
		//Remove records that
			// a. Refer to a resource in the JAR
			// b. Are missing a field
			// c. Refer to a file that doesn't exist
		Set<String> existingKeys = existing.keySet();
		Set<String> keepKeys = new HashSet<String>();
		for (String k : existingKeys)
		{
			if (jarkeys.contains(k)) continue;
			BuildRecord rec = existing.get(k);
			String bpath = rec.gbPath;
			if (bpath == null || bpath.isEmpty()) continue;
			String fullpath = targetDir + File.separator + bpath;
			if (!FileBuffer.fileExists(fullpath)) continue;
			String rpath = rec.refSeqPath;
			if (rpath == null || rpath.isEmpty()) continue;
			fullpath = targetDir + File.separator + rpath;
			if (!FileBuffer.fileExists(fullpath)) continue;
			keepKeys.add(k);
		}
		
		//Replace all standard references in the dir with the JAR current
		for (String k : jarkeys)
		{
			try
			{
				String packPath = GenomeBuild.getStandardBuildJARPath(k);
				InputStream is = GenomeBuild.class.getResourceAsStream(packPath);
				//InputStream is = GenomeBuild.class.getClassLoader().getResourceAsStream("/hospelhornbg_genomeBuild." + packPath);
				
				if (is == null){
					System.err.println("ERROR: Resource \"" + packPath + "\" could not be opened...");
					return false;
				}
				System.err.println("Reading resource " + packPath);
				FileBuffer buffer = new FileBuffer(512 * 1024); //512 KB buffer
				int i = -1;
				while ((i = is.read()) != -1)
				{
					buffer.addToFile((byte)i);
				}
				
				String packName = packPath.substring(packPath.lastIndexOf('/') + 1);
				String tPath = targetDir + File.separator + packName;
				System.err.println("Copying " + packPath + " to " + tPath);
				buffer.writeFile(tPath);
				
				is.close();	
				
				packPath = GeneSet.getStandardDB_packagePath(k);
				if (packPath == null || packPath.isEmpty()) continue;
				is = GeneSet.class.getResourceAsStream(packPath);
				if (is == null){
					System.err.println("ERROR: Resource \"" + packPath + "\" could not be opened...");
					return false;
				}
				System.err.println("Reading resource " + packPath);
				buffer = new FileBuffer(512 * 1024); //512 KB buffer
				i = -1;
				while ((i = is.read()) != -1)
				{
					buffer.addToFile((byte)i);
				}
				
				packName = packPath.substring(packPath.lastIndexOf('/') + 1);
				tPath = targetDir + File.separator + packName;
				System.err.println("Copying " + packPath + " to " + tPath);
				buffer.writeFile(tPath);
				
				is.close();	
				
			}
			catch (IOException e)
			{
				System.err.println("ERROR: One or more resources could not be copied from JAR");
				e.printStackTrace();
				return false;
			}

		}
		
		//Rewrite the table
		try
		{
			FileWriter fw = new FileWriter(buildTablePath);
			BufferedWriter bw = new BufferedWriter(fw);	
			
			//Standard
			for (String k : jarkeys){
				String packPath = GenomeBuild.getStandardBuildJARPath(k);
				String packName = packPath.substring(packPath.lastIndexOf('/') + 1);
				String rgPackPath = GeneSet.getStandardDB_packagePath(k);
				if (rgPackPath != null)
				{
					String rgPackName = rgPackPath.substring(rgPackPath.lastIndexOf('/') + 1);
					bw.write(k + "\t" + packName + "\t" + rgPackName + "\n");
				}
				else bw.write(k + "\t" + packName + "\tNONE\n");
			}

			//Existing
			for (String k : keepKeys)
			{
				BuildRecord rec = existing.get(k);
				bw.write(k + "\t" + rec.gbPath + "\t" + rec.refSeqPath + "\n");
			}
			
			bw.close();
			fw.close();
		}
		catch (IOException e)
		{
			System.err.println("ERROR: Genome build lookup table could not be written.");
			e.printStackTrace();
			return false;
		}
	
		return true;
	}
	
	public static void installBIOISVTOOLS(String[] args, String homeDir)
	{
		String refDir = null;
		for (int i = 0; i < args.length; i++)
		{
			String s = args[i];
			if (s.equals(OP_REFPATH))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: -d flag MUST be followed by directory path!");
					printUsage();
					System.exit(1);
				}
				refDir = args[i+1];
			}
		}
		
		if (refDir == null || refDir.isEmpty())
		{
			System.err.println("ERROR: Target reference directory path is required!");
			printUsage();
			System.exit(1);
		}
		
		//Remove separator char on path, if present.
		if (refDir.charAt(refDir.length() - 1) == File.separatorChar) {
			refDir = refDir.substring(0, refDir.length() - 1);
		}
		
		//Set up init file
		String initDir = homeDir + File.separator + "bioisvtools";
		String initFile = initDir + File.separator + "init.cfg";
		
		if (!FileBuffer.directoryExists(initDir))
		{
			try 
			{
				Files.createDirectory(Paths.get(initDir));
			} 
			catch (IOException e) 
			{
				System.err.println("ERROR: Directory " + initDir + " could not be created! (I/O error)");
				e.printStackTrace();
				System.exit(1);
			}
		}
		
		if (FileBuffer.fileExists(initFile))
		{
			System.err.println("Init file " + initFile + " already exists!");
			System.err.println("Adding new installation path to existing paths...");
			List<String> paths = new LinkedList<String>();
			
			try
			{
				FileReader fr = new FileReader(initFile);
				BufferedReader br = new BufferedReader(fr);
				
				String line = null;
				while ((line = br.readLine()) != null) paths.add(line);
				
				br.close();
				fr.close();
				
				Files.deleteIfExists(Paths.get(initFile));
				
				FileWriter fw = new FileWriter(initFile);
				BufferedWriter bw = new BufferedWriter(fw);
				bw.write(refDir + "\n");
				for (String s : paths) {
					//bw.write(s + "\n");
					if (!s.equals(refDir)) bw.write(s + "\n");
				}
				
				bw.close();
				fw.close();
			}
			catch (IOException e) {
				System.err.println("ERROR: Init file could not be read or written!");
				e.printStackTrace();
				System.exit(1);
			}

		}
		else
		{
			try
			{
				FileWriter fw = new FileWriter(initFile);
				BufferedWriter bw = new BufferedWriter(fw);
				bw.write(refDir);
				bw.close();
				fw.close();
			}
			catch (IOException e) {
				System.err.println("ERROR: Init file could not be written!");
				e.printStackTrace();
				System.exit(1);
			}
		}
		
		//Now, check the refdir for references.
		//See if the references contained in this jar are there.
		//If not, try to copy them into the refdir folder.
		String buildTablePath = refDir + File.separator + ConsoleMain.FILE_GBDAT;
		
		System.err.println("Checking for existence of requested reference directory...");
		if (!FileBuffer.directoryExists(refDir))
		{
			System.err.println(refDir + " was not found. Attempting to create directory...");
			try 
			{
				Files.createDirectories(Paths.get(refDir));
			} 
			catch (IOException e) {
				System.err.println(refDir + " could not be created! Exiting...");
				e.printStackTrace();
				System.exit(1);
			}
			if (!copyReferences(refDir))
			{
				System.err.println("There was an error updating reference materials!");
				System.exit(1);
			}
		}
		else
		{
			System.err.println(refDir + " was already found! Checking existing references...");
			System.err.println("Checking for genome build reference table...");
			if (FileBuffer.fileExists(buildTablePath))
			{
				System.err.println("Pre-existing Genome Build reference table found.");
				System.err.println("Updating...");
				if (!copyReferences(refDir))
				{
					System.err.println("There was an error updating reference materials!");
					System.exit(1);
				}
			}
			else
			{
				System.err.println("Genome Build reference table was not found!");
				if (!copyReferences(refDir))
				{
					System.err.println("There was an error updating reference materials!");
					System.exit(1);
				}
			}
		}
		
		
	}

}
