package hospelhornbg_svdb;

public interface QueryCondition {
	
	//This one is only for the variants themselves, not genotype conditions!
	
	public boolean passes(String varRecord);

}
