package hospelhornbg_segregation;

import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

public enum Population {
	
	NFE(0, "Non-Finnish European", "NFE"),
	AFR(1, "African/African American", "AFR"),
	AMR(2, "Latino/Native American", "AMR"),
	FIN(3, "Finnish European", "FIN"),
	EAS(4, "East Asian", "EAS"),
	SAS(5, "South Asian", "SAS"),
	ASJ(6, "Ashkenazi Jewish", "ASJ"),
	OTH(7, "Other", "OTH");

	private int eid;
	private String desc;
	private String shrt;
	
	private Population(int i, String str, String shortStr)
	{
		desc = str;
	}
	
	public int getIDNumber()
	{
		return eid;
	}
	
	public String toString()
	{
		return desc;
	}
	
	private static Map<Integer, Population> idmap;
	
	private static void populateMap()
	{
		idmap = new TreeMap<Integer, Population>();
		Population[] all = Population.values();
		for (Population p : all) idmap.put(p.getIDNumber(), p);
	}
	
	public static Population getPopulation(int id)
	{
		if(idmap == null) populateMap();
		return idmap.get(id);
	}
	
	private static Map<String, Population> codemap;
	
	private static void populateCodeMap()
	{
		codemap = new HashMap<String, Population>();
		Population[] all = Population.values();
		for (Population p : all) codemap.put(p.getShortString(), p);
	}
	
	public static Population getPopulation(String shortCode)
	{
		if(codemap == null) populateCodeMap();
		return codemap.get(shortCode);
	}
	
	public String getShortString()
	{
		return shrt;
	}

}
