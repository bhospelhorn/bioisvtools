package hospelhornbg_svdb;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_genomeBuild.GeneSet;
import hospelhornbg_genomeBuild.GenomeBuild;
import hospelhornbg_genomeBuild.GenomeBuildUID;
import hospelhornbg_segregation.Family;
import hospelhornbg_segregation.FamilyMember;
import waffleoRai_Utils.FileBuffer;

public class SVDatabase {
	
	/* ----- Constants ----- */
	
	public static final double POPFREQ_CUTOFF_G2 = 0.05;
	public static final double POPFREQ_CUTOFF_G3 = 0.02;
	public static final double POPFREQ_CUTOFF_G4 = 0.01;
	
	public static final String SETTINGS_FILE = "SVDB_settings.bin";
	
	/* ----- Instance Variables ----- */
	
	private DBSampleTable sampleTable;
	private DBVariantTable variantTable;
	
	private String directory;
	private String dbName;
	
	private int mergeFactor;
	private GenomeBuild genome;
	private GeneSet genes;
	
	private String omim_table_path;
	
	/* ----- Private Constructors ----- */
	
	private SVDatabase(String dir, String name)
	{
		//TODO
		//Load existing
	}
	
	private SVDatabase(String dir, String name, int leewayValue, GenomeBuild gb)
	{
		//TODO
		//New DB
	}
	
	private void loadSettings()
	{
		//TODO
		
	}
	
	public void saveDatabase()
	{
		//TODO
	}
	
	/* ----- Sample Management ----- */
	
	public void addFamily(Family fam)
	{
		sampleTable.addOrReplaceFamily(fam);
	}

	public boolean removeFamily(int famUID)
	{
		//Don't forget to remove from variant table...
		Family f = sampleTable.removeFamily(famUID);
		if(f == null) return false;
		
		return variantTable.removeFamily(f);
	}
	
	public boolean removeFamily(String famName)
	{
		Family f = sampleTable.removeFamily(famName);
		if(f == null) return false;
		
		return variantTable.removeFamily(f);
	}
	
	public boolean removeSample(int sampleUID)
	{
		FamilyMember s = sampleTable.removeSample(sampleUID);
		if(s == null) return false;
		return variantTable.removeSample(s);
	}
	
	public Family getFamily(int famUID)
	{
		return sampleTable.getFamily(famUID);
	}
	
	public Family getFamily(String famName)
	{
		return sampleTable.getFamily(famName);
	}
	
	public FamilyMember getSample(int sampleUID)
	{
		return sampleTable.getSample(sampleUID);
	}
	
	public FamilyMember getSample(String sampleName)
	{
		return sampleTable.getSample(sampleName);
	}
	
	/* ----- Variant Management ----- */
	
	private List<String> getVCFSampleNames(String vcfpath) throws IOException
	{
		List<String> slist = new LinkedList<String>();
		BufferedReader br = new BufferedReader(new FileReader(vcfpath));
		String line = null;
		while((line = br.readLine()) != null)
		{
			if(line.isEmpty()) continue;
			if(line.startsWith("#CHROM"))
			{
				String[] fields = line.split("\t");
				if(fields.length > 9)
				{
					for(int i = 9; i < fields.length; i++)
					{
						slist.add(fields[i]);
					}
				}
				break;
			}
		}
		
		br.close();
		
		return slist;
	}
	
	public boolean addVCF(String vcfpath, boolean verbose) throws IOException
	{
		//See if VCF even exists
		if(!FileBuffer.fileExists(vcfpath)) return false;
		
		//Get sample names and pull sample info
		List<String> snames = getVCFSampleNames(vcfpath);
		Map<String, FamilyMember> smap = new HashMap<String, FamilyMember>();
		for(String s : snames)
		{
			FamilyMember fm = sampleTable.getSample(s);
			if (fm == null) 
			{
				if(verbose) System.err.println("SVDatabase.addVCF || WARNING: Sample \"" + s + "\" was not found in database! Skipping...");
				continue;
			}
			smap.put(s, fm);
		}
		
		//Add to variant table...
		return variantTable.addVCF(vcfpath, smap, mergeFactor);
	}
	
	/* ----- Family Data Dumping ----- */
	
	public void dumpFamily(int famUID, String dumpDir)
	{
		//TODO
		/*
		 * Table Fields ----- (CSV)
		 * 		Var Name
		 * 		Var ID (Hex)
		 * 		Chrom 1
		 * 		Chrom 2 (Empty if not TRA)
		 * 		Pos (Range)
		 * 		End (Range)
		 * 		Insertion Sequence (Empty if not INS/INS:ME)
		 * 		SVType
		 * 		Size
		 * 		PosEff
		 * 		Gene Name
		 * 		Transcript ID
		 * 		OMIM Info (Blank if none found)
		 * 		DECIPHER Score (Blank if none found)
		 * 		UDP Cohort Total Indivs with variant
		 * 		UDP Cohort Hom Count
		 * 		UDP Cohort Pop Freq
		 * 		UDP Cohort Gene Hit Count (Variant)
		 * 		UDP Cohort Gene Hit Count (Indivs)
		 * 		UDP Cohort Exon Hit Count (Variant)
		 * 		UDP Cohort Exon Hit Count (Indivs)
		 * 		[Population] Allele Count
		 * 		[Population] Hom Count
		 * 		[Population] Pop Freq
		 * 			...
		 * 		[Seg Patterns]
		 * 			Per affected...
		 * 		CompHet Partners (String IDs(HexID))
		 * 		[Genotypes] #Occurrences of allele
		 * 		[Genotypes]	Start Position
		 * 		[Genotypes]	End Position
		 * 			...
		 */
	}
	
	public void setOMIMTablePath(String path)
	{
		//TODO
	}
	
	public void dumpToBED(int sampleUID, SVType type, String outpath)
	{
		//TODO
	}
	
	/* ----- Database Loaders ----- */
	
	public static SVDatabase loadDatabase(String dir, String name)
	{
		//TODO
		return null;
	}
	
	public static SVDatabase newDatabase(String dir, String name, int mergeFactor, GenomeBuildUID gbid)
	{
		//TODO
		return null;
	}
	
}
