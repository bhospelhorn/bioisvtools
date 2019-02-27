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
	DENOVO_REC("DeNovo Recessive (MV)"),
	UNRESOLVED("Inheritance Unresolved"),
	MVIOL("Mendelian Violation"),
	X_LINKED_REC("X-Linked Recessive"),
	X_LINKED_DOM("X-Linked Dominant"),
	X_LINKED_DN("X-Linked De Novo"),
	X_LINKED_MV("X-Linked Mendelian Violation"),
	X_PATIMPRINT_HH("X-Linked Half-Het (Imprint Rescue)"),
	X_PATIMPRINT_DOM("X-Linked Dominant (Imprint Rescue)"),
	Y_LINKED_DOM("Y-Linked Dominant"),
	Y_LINKED_DN("Y-Linked De Novo"),
	Y_LINKED_SPECIAL("Y-Linked Special Case"),;
	
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
		case X_PATIMPRINT_HH:
			return true;
		default:
			return false;
		}
	}

}
