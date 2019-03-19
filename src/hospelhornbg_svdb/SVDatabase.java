package hospelhornbg_svdb;

import java.io.File;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_genomeBuild.GenomeBuild;
import hospelhornbg_segregation.Population;

public class SVDatabase {
	
	/* --- Constants --- */
	
	public static final String TABLE_EXTENSION = ".tsv";
	public static final String INDEX_EXTENSION = ".idx";
	
	/* --- Instance Variables --- */
	
	private String dirPath;
	private String dbName;
	private GenomeBuild genomeBuild;
	
	private int indivCount;
	private int famCount;
	
	private Map<Population, Integer> popIndivCount;
	
	/* --- Construction --- */
	
	public SVDatabase(String directory, String name, GenomeBuild genome)
	{
		dirPath = directory;
		dbName = name;
		genomeBuild = genome;
	}
	
	/* --- Paths --- */
	
	public String getVariantTablePath(SVType svtype)
	{
		return dirPath + File.separator + "svdb_" + dbName + "_" + svtype.getString() + TABLE_EXTENSION;
	}
	
	public String getGenoTablePath(SVType svtype)
	{
		return dirPath + File.separator + "svdb_" + dbName + "_genoTable_" + svtype.getString() + TABLE_EXTENSION;
	}
	
	/* --- Read --- */
	
	/* --- Write --- */
	
	/* --- Variant Addition --- */
	
	public boolean addFamilyVCF(String vcfpath, Collection<Population> pGroups)
	{
		return false;
	}
	
	/* --- Variant Query --- */
	
	public List<String> queryForRecords(Collection<QueryCondition> conditions)
	{
		return null;
	}

}
