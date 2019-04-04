package hospelhornbg_svdb;

import java.io.IOException;

import hospelhornbg_genomeBuild.GenomeBuild;
import hospelhornbg_segregation.Family;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class ConsoleFrontEnd {
	
	/* ----- Constants ----- */
	
	public static final String PROG_NEWDB = "newdb";
	public static final String PROG_ADDFAM = "addfam";
	public static final String PROG_UPDATEFAM = "updatefam";
	public static final String PROG_REMOVEFAM = "removefam";
	public static final String PROG_DUMPFAM = "famvardump";
	
	public static final String PROG_VARINFO = "varinfo"; //Dumps all info on a single variant, including all genotypes!
	
	public static final String OP_GENOME = "-g";
	public static final String OP_FAMILY = "-f";
	
	public static final String OP_DBNAME = "-n";
	public static final String OP_DBPATH = "-d";
	
	public static final String OP_INPUTFAMI = "-F";
	public static final String OP_INPUTVCF = "-V";
	
	public static final String OP_OUTPUTPATH = "-o";
	
	public static final String OP_LEEWAY = "-l";
	
	public static final String OP_UID = "-I"; //Can be used for sample or variant IDs

	/* ----- Usage Message ----- */
	
	public void printUsage()
	{
		//TODO
	}
	
	/* ----- Helper Methods ----- */
	
	/* ----- Program Primary Methods ----- */
	
	public static void newDB(String dbName, String dbDir, GenomeBuild gb, int leeway)
	{
		if (gb == null)
		{
			System.err.println("ConsoleFrontEnd.newDB || ERROR: Genome build is invalid! Database could not be created!");
			return;
		}
		try 
		{
			System.err.println("Creating new SV database...");
			System.err.println("DB Directory: " + dbDir);
			System.err.println("DB Name: " + dbName);
			System.err.println("Genome Build: " + gb.toString());
			double ld = (double)leeway/ 1000.0;
			System.err.println("Leeway: " + leeway + " (" + ld + ")");
			SVDatabase db = SVDatabase.newDatabase(dbDir, dbName, gb, leeway);
			db.saveDatabase();
			System.err.println("Empty database has been created!");
		} 
		catch (IOException e) 
		{
			System.err.println("ConsoleFrontEnd.newDB || ERROR: Database could not be created!");
			e.printStackTrace();
		}
	}
	
	public static void addFam(String dbName, String dbDir, String famiPath, String vcfPath)
	{
		//Open FAMI file
		Family fam = null;
		System.err.println("Reading family file: " + famiPath);
		try 
		{
			fam = Family.readFromFAMI(famiPath);
		} 
		catch (IOException e1) 
		{
			System.err.println("ConsoleFrontEnd.addFam || ERROR: Family FAMI file could not be opened!");
			System.err.println("Family File Path: " + famiPath);
			e1.printStackTrace();
			return;
		} 
		catch (UnsupportedFileTypeException e1) 
		{
			System.err.println("ConsoleFrontEnd.addFam || ERROR: Family FAMI file could not be read!");
			System.err.println("Family File Path: " + famiPath);
			e1.printStackTrace();
			return;
		}
		
		//Open DB
		System.err.println("Opening database...");
		SVDatabase db = null;
		try {db = SVDatabase.readDatabase(dbDir, dbName);} 
		catch (IOException e) 
		{
			System.err.println("ConsoleFrontEnd.addFam || ERROR: Database could not be opened!");
			System.err.println("DB Name: " + dbName);
			System.err.println("DB Location: " + dbDir);
			e.printStackTrace();
			return;
		}
		System.err.println("Database opened! Now reading family information/ variants...");
		
		//Read into db
		try 
		{
			db.addFamily(vcfPath, fam);
		} 
		catch (IOException e) 
		{
			System.err.println("ConsoleFrontEnd.addFam || ERROR: An I/O error occurred when reading VCF or writing to database!");
			e.printStackTrace();
			return;
		} 
		catch (UnsupportedFileTypeException e) 
		{
			System.err.println("ConsoleFrontEnd.addFam || ERROR: One or more files could not be read!");
			e.printStackTrace();
			return;
		}
	}

	public static void updateFam(String dbName, String dbDir, String famName, String newVCF)
	{
		//Open DB
		System.err.println("Opening database...");
		SVDatabase db = null;
		try {db = SVDatabase.readDatabase(dbDir, dbName);} 
		catch (IOException e) 
		{
			System.err.println("ConsoleFrontEnd.updateFam || ERROR: Database could not be opened!");
			System.err.println("DB Name: " + dbName);
			System.err.println("DB Location: " + dbDir);
			e.printStackTrace();
			System.exit(1);
		}
		
		//Get family
		Family f = db.getFamily(famName);
		if (f == null)
		{
			System.err.println("ConsoleFrontEnd.updateFam || ERROR: Family not found in database!");
			System.err.println("Family Name: " + famName);
			System.exit(1);
		}
		
		try 
		{
			boolean b = db.updateFamily(f, newVCF);
			if(b) System.err.println("ConsoleFrontEnd.updateFam || Update Succeeded!");
			else System.err.println("ConsoleFrontEnd.updateFam || Update Failed!");
		} 
		catch (IOException e) 
		{
			System.err.println("ConsoleFrontEnd.updateFam || ERROR: One or more db files could not be read/written to!");
			System.err.println("ConsoleFrontEnd.updateFam || Update Failed!");
			e.printStackTrace();
			System.exit(1);
		} 
		catch (UnsupportedFileTypeException e) 
		{
			System.err.println("ConsoleFrontEnd.updateFam || ERROR: Input VCF could not be read!");
			System.err.println("ConsoleFrontEnd.updateFam || Update Failed!");
			e.printStackTrace();
			System.exit(1);
		}
		
	}
	
	public static void removeFam(String dbName, String dbDir, String famName)
	{
		//Open DB
		System.err.println("Opening database...");
		SVDatabase db = null;
		try {db = SVDatabase.readDatabase(dbDir, dbName);} 
		catch (IOException e) 
		{
			System.err.println("ConsoleFrontEnd.removeFam || ERROR: Database could not be opened!");
			System.err.println("DB Name: " + dbName);
			System.err.println("DB Location: " + dbDir);
			e.printStackTrace();
			System.exit(1);
		}
		
		//Get family
		Family f = db.getFamily(famName);
		if (f == null)
		{
			System.err.println("ConsoleFrontEnd.removeFam || ERROR: Family not found in database!");
			System.err.println("Family Name: " + famName);
			System.exit(1);
		}
		
		try 
		{
			boolean b = db.removeFamily(f);
			if(b) System.err.println("ConsoleFrontEnd.removeFam || Update Succeeded!");
			else System.err.println("ConsoleFrontEnd.removeFam || Update Failed!");
		} 
		catch (IOException e) 
		{
			System.err.println("ConsoleFrontEnd.removeFam || ERROR: One or more db files could not be read/written to!");
			System.err.println("ConsoleFrontEnd.removeFam || Update Failed!");
			e.printStackTrace();
			System.exit(1);
		} 

	}
	
	public static void dumpFam(String dbName, String dbDir, String famName, String outpath)
	{
		//TODO
	}
	
	/* ----- Main Method ----- */
	
	public static void runSVDB(String[] args, GenomeBuild gb, boolean verbose)
	{
		//TODO
		
		//Expects svdb mode name to be after program name
		// ie. should be args[1]
		
		
		
	}
	
}
