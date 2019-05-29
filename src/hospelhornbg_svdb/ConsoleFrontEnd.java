package hospelhornbg_svdb;

import java.io.IOException;

import hospelhornbg_genomeBuild.GenomeBuild;
import hospelhornbg_segregation.Family;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class ConsoleFrontEnd {
	
	//TODO:
	//	1. When making new db, allow for setting omim path 
	//	2. Note which samples have had vars imported.
	//		Reason: For telling the difference between 0/0 and ./. genotypes
	
	/* ----- Constants ----- */
	
	public static final String PROG_NEWDB = "newdb";
	public static final String PROG_ADDFAM = "addfam";
	public static final String PROG_UPDATEFAM = "updatefam";
	public static final String PROG_REMOVEFAM = "removefam";
	public static final String PROG_DUMPFAM = "famvardump";
	public static final String PROG_ADDVARS = "addvars";
	
	public static final String PROG_VARINFO = "varinfo"; //Dumps all info on a single variant, including all genotypes!
	
	public static final String OP_GENOME = "-g";
	public static final String OP_FAMILY = "-f";
	
	public static final String OP_DBNAME = "-n";
	public static final String OP_DBPATH = "-d";
	
	public static final String OP_INPUTFAMI = "-F";
	public static final String OP_INPUTVCF = "-V";
	public static final String OP_OMIMPATH = "-O";
	
	public static final String OP_OUTPUTPATH = "-o";
	
	public static final String OP_LEEWAY = "-l";
	
	public static final String OP_UID = "-I"; //Can be used for sample or variant IDs

	public static final int DEFO_MERGE_FACTOR = 50;
	
	/* ----- Usage Message ----- */
	
	public static void printUsage()
	{
		//TODO
	}
	
	/* ----- Helper Methods ----- */
	
	/* ----- Program Primary Methods ----- */
	
	public static void newDB(String dbName, String dbDir, GenomeBuild gb, int leeway, String omimpath) throws IOException
	{
		SVDatabase db = SVDatabase.newDatabase(dbDir, dbName, leeway, gb);
		if(omimpath != null) db.setOMIMTablePath(omimpath);
		db.saveDatabase();
	}
	
	public static void addOrUpdateFam(String dbDir, String famiPath) throws IOException, UnsupportedFileTypeException
	{
		Family fam = Family.readFromFAMI(famiPath);
		SVDatabase db = SVDatabase.loadDatabase(dbDir);
		db.addFamily(fam);
		db.saveDatabase();
	}
	
	public static void removeFam(String dbDir, String famName) throws IOException
	{
		SVDatabase db = SVDatabase.loadDatabase(dbDir);
		db.removeFamily(famName);
	}
	
	public static void dumpFam(String dbDir, String famName, String outpath, boolean verbose) throws IOException
	{
		SVDatabase db = SVDatabase.loadDatabase(dbDir);
		Family fam = db.getFamily(famName);
		if(fam == null) return;
		db.dumpFamily(fam, outpath, verbose);
	}
	
	public static void addVCF(String dbDir, String vcfPath, boolean verbose) throws IOException
	{
		SVDatabase db = SVDatabase.loadDatabase(dbDir);
		db.addVCF(vcfPath, verbose);
		db.saveDatabase();
	}
	
	/* ----- Main Method ----- */
	
	public static void runSVDB(String[] args, GenomeBuild gb, boolean verbose)
	{
		//TODO
		
		//Expects svdb mode name to be after program name
		// ie. should be args[1]
		if(args.length < 2)
		{
			printUsage();
			System.exit(1);
		}
		
		String mode = args[1];
		String famname = null;
		String dbname = null;
		String dbdir = null;
		String famiPath = null;
		String vcfPath = null;
		String outPath = null;
		String omimPath = null;
		int mFactor = DEFO_MERGE_FACTOR;
		
		for (int i = 0; i < args.length; i++)
		{
			String s = args[i];
			if (s.equals(OP_DBPATH))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_DBPATH + " flag MUST be followed by path to a DB directory!");
					printUsage();
					System.exit(1);
				}
				dbdir = args[i+1];
			}
			else if (s.equals(OP_DBNAME))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_DBNAME + " flag MUST be followed by a DB name!");
					printUsage();
					System.exit(1);
				}
				dbname = args[i+1];
			}
			else if (s.equals(OP_FAMILY))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_FAMILY + " flag MUST be followed by the ID of a family!");
					printUsage();
					System.exit(1);
				}
				famname = args[i+1];
			}
			else if (s.equals(OP_INPUTFAMI))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_INPUTFAMI + " flag MUST be followed by a path to a FAMI (.fam) file!");
					printUsage();
					System.exit(1);
				}
				famiPath = args[i+1];
			}
			else if (s.equals(OP_INPUTVCF))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_INPUTVCF + " flag MUST be followed by a path to a VCF file!");
					printUsage();
					System.exit(1);
				}
				vcfPath = args[i+1];
			}
			else if (s.equals(OP_OMIMPATH))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_OMIMPATH + " flag MUST be followed by a path to OMIM genemap2 file!");
					printUsage();
					System.exit(1);
				}
				omimPath = args[i+1];
			}
		}
		
		
		
		if(mode.equals(PROG_NEWDB))
		{
			//Need dbname, dbdir
			//mfactor optional - defaults to 50
			if(dbname == null || dbname.isEmpty())
			{
				System.err.println(PROG_NEWDB + " ERROR | Database name is required!");
				printUsage();
				System.exit(1);
			}
			if(dbdir == null || dbdir.isEmpty())
			{
				System.err.println(PROG_NEWDB + " ERROR | Database directory is required!");
				printUsage();
				System.exit(1);
			}
			
			try 
			{
				newDB(dbname, dbdir, gb, mFactor, omimPath);
			} 
			catch (IOException e) 
			{
				System.err.println(PROG_NEWDB + " ERROR | Database creation failed!");
				e.printStackTrace();
				System.exit(1);
			}
		}
		else if(mode.equals(PROG_ADDFAM))
		{
			//Need dbdir, famipath
			if(dbdir == null || dbdir.isEmpty())
			{
				System.err.println(PROG_ADDFAM + " ERROR | Database directory is required!");
				printUsage();
				System.exit(1);
			}
			if(famiPath == null || famiPath.isEmpty())
			{
				System.err.println(PROG_ADDFAM + " ERROR | FAMI path is required!");
				printUsage();
				System.exit(1);
			}
			
			try 
			{
				addOrUpdateFam(dbdir, famiPath);
			} 
			catch (IOException e) 
			{
				System.err.println(PROG_ADDFAM + " ERROR | Family could not be added (IO Error)!");
				e.printStackTrace();
				System.exit(1);
			} 
			catch (UnsupportedFileTypeException e) 
			{
				System.err.println(PROG_ADDFAM + " ERROR | Family could not be added (File Error)!");
				e.printStackTrace();
				System.exit(1);
			}
			
		}
		else if(mode.equals(PROG_REMOVEFAM))
		{
			//Need dbdir, family name
			if(dbdir == null || dbdir.isEmpty())
			{
				System.err.println(PROG_REMOVEFAM + " ERROR | Database directory is required!");
				printUsage();
				System.exit(1);
			}
			if(famname == null || famname.isEmpty())
			{
				System.err.println(PROG_REMOVEFAM + " ERROR | Family name is required!");
				printUsage();
				System.exit(1);
			}
			
			try 
			{
				removeFam(dbdir, famname);
			} 
			catch (IOException e) 
			{
				System.err.println(PROG_REMOVEFAM + " ERROR | Family could not be removed!");
				printUsage();
				System.exit(1);
			}
		}
		else if(mode.equals(PROG_DUMPFAM))
		{
			//Need dbdir, famname, outdir, verbose
			if(dbdir == null || dbdir.isEmpty())
			{
				System.err.println(PROG_DUMPFAM + " ERROR | Database directory is required!");
				printUsage();
				System.exit(1);
			}
			if(famname == null || famname.isEmpty())
			{
				System.err.println(PROG_DUMPFAM + " ERROR | Family name is required!");
				printUsage();
				System.exit(1);
			}
			if(outPath == null || outPath.isEmpty())
			{
				System.err.println(PROG_DUMPFAM + " ERROR | Output dir is required!");
				printUsage();
				System.exit(1);
			}
			
			try 
			{
				dumpFam(dbdir, famname, outPath, verbose);
			} 
			catch (IOException e) 
			{
				System.err.println(PROG_DUMPFAM + " ERROR | Dump failed!");
				e.printStackTrace();
				System.exit(1);
			}
			
		}
		else if(mode.equals(PROG_ADDVARS))
		{
			//Need dbdir and vcfpath
			if(dbdir == null || dbdir.isEmpty())
			{
				System.err.println(PROG_ADDVARS + " ERROR | Database directory is required!");
				printUsage();
				System.exit(1);
			}
			if(vcfPath == null || vcfPath.isEmpty())
			{
				System.err.println(PROG_ADDVARS + " ERROR | VCF path is required!");
				printUsage();
				System.exit(1);
			}
			
			try 
			{
				addVCF(dbdir, vcfPath, verbose);
			} 
			catch (IOException e) 
			{
				System.err.println(PROG_ADDVARS + " ERROR | VCF could not be added!");
				e.printStackTrace();
				System.exit(1);
			}
			
		}
		else
		{
			System.err.println("ERROR: Mode \"" + mode + "\" not recognized!");
			printUsage();
			System.exit(1);	
		}
		
		
		
	}
	
}
