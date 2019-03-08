package hospelhornbg_segregation;

public enum Population {
	
	NFE("Non-Finnish European"),
	AFR("African/African American"),
	AMR("Latino/Native American"),
	FIN("Finnish European"),
	EAS("East Asian"),
	SAS("South Asian"),
	ASJ("Ashkenazi Jewish"),
	OTH("Other");

	private String desc;
	
	private Population(String str)
	{
		desc = str;
	}
	
	public String toString()
	{
		return desc;
	}

}
