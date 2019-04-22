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
		//TODO
	}
	
	public static void addFam(String dbName, String dbDir, String famiPath, String vcfPath)
	{
		//TODO
	}

	public static void updateFam(String dbName, String dbDir, String famName, String newVCF)
	{
		//TODO
	}
	
	public static void removeFam(String dbName, String dbDir, String famName)
	{
		//TODO

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
