package hospelhornbg_genomeBuild;

public enum GenomeBuildUID {
	
	GRCh37(2110496020, "GRCh37", "Homo sapiens"),
	GRCh38(-1237832459, "GRCh38", "Homo sapiens"),
	NCBI36(3198822, "hg18", "Homo sapiens");
	
	private int UID;
	private String lookup_name;
	private String species;
	
	private GenomeBuildUID(int intID, String buildname, String speciesName)
	{
		UID = intID;
		lookup_name = buildname;
	}
	
	public int getUID()
	{
		return UID;
	}
	
	public String getName()
	{
		return lookup_name;
	}
	
	public static GenomeBuildUID getByID(int ID)
	{
		GenomeBuildUID[] all = GenomeBuildUID.values();
		for (GenomeBuildUID i : all)
		{
			if (i.UID == ID) return i;
		}
		return null;
	}
	
	public static GenomeBuildUID getByName(String name)
	{
		GenomeBuildUID[] all = GenomeBuildUID.values();
		for (GenomeBuildUID i : all)
		{
			if (i.lookup_name.equalsIgnoreCase(name)) return i;
		}
		return null;
	}

	public String toString()
	{
		return getName() + " (" + species + ")";
	}
	
}
