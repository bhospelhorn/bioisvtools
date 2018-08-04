package hospelhornbg_segregation;

public enum Inheritance {
	
	HOM_REC,
	DOMINANT,
	HALF_HET,
	HALF_HET_SV,
	COMP_HET,
	COMP_HET_DEL,
	COMP_HET_SV,
	DENOVO_DOM,
	DENOVO_HET,
	DENOVO_HET_SV,
	DENOVO_REC,
	UNRESOLVED,
	MVIOL;
	
	public static final int HALFHET_PHASE_NONE = 0;
	public static final int HALFHET_PHASE_PARENT1 = 1;
	public static final int HALFHET_PHASE_PARENT2 = 2;
	public static final int HALFHET_PHASE_UNDETERMINED = -1;

}
