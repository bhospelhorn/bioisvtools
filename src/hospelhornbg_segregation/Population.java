package hospelhornbg_segregation;

import java.util.Map;
import java.util.TreeMap;

public enum Population {
	
	NFE(0, "Non-Finnish European"),
	AFR(1, "African/African American"),
	AMR(2, "Latino/Native American"),
	FIN(3, "Finnish European"),
	EAS(4, "East Asian"),
	SAS(5, "South Asian"),
	ASJ(6, "Ashkenazi Jewish"),
	OTH(7, "Other");

	private int eid;
	private String desc;
	
	private Population(int i, String str)
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

}
