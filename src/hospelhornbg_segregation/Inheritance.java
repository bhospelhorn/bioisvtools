package hospelhornbg_segregation;

public enum Inheritance {
	
	HOM_REC("Homozygous Recessive"),
	DOMINANT("Dominant"),
	HALF_HET("Half Het"),
	HALF_HET_SV("Half Het (CNV)"),
	COMP_HET("Compound Heterozygote"),
	COMP_HET_DEL("Compound Het w/ Deletion"),
	COMP_HET_SV("Compound Het w/ CNV"),
	DENOVO_DOM("DeNovo Dominant"),
	DENOVO_HET("DeNovo Het"),
	DENOVO_HET_SV("DeNovo Het (CNV)"),
	DENOVO_REC("DeNovo Recessive"),
	UNRESOLVED("Inheritance Unresolved"),
	MVIOL("Medelian Violation");
	
	public static final int HALFHET_PHASE_NONE = 0;
	public static final int HALFHET_PHASE_PARENT1 = 1;
	public static final int HALFHET_PHASE_PARENT2 = 2;
	public static final int HALFHET_PHASE_UNDETERMINED = -1;
	
	private String str;
	
	private Inheritance(String name)
	{
		str = name;
	}
	
	public String toString()
	{
		return str;
	}
	
	public static boolean isHalfhet(Inheritance i)
	{
		switch (i)
		{
		case COMP_HET:
			return false;
		case COMP_HET_DEL:
			return false;
		case COMP_HET_SV:
			return false;
		case DENOVO_DOM:
			return false;
		case DENOVO_HET:
			return true;
		case DENOVO_HET_SV:
			return true;
		case DENOVO_REC:
			return false;
		case DOMINANT:
			return false;
		case HALF_HET:
			return true;
		case HALF_HET_SV:
			return true;
		case HOM_REC:
			return false;
		case MVIOL:
			return false;
		case UNRESOLVED:
			return false;
		default:
			return false;
		}
	}

}
