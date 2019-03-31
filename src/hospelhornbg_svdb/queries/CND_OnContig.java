package hospelhornbg_svdb.queries;

import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_svdb.DBVariant;
import hospelhornbg_svdb.QueryCondition;

public class CND_OnContig implements QueryCondition{
	
	private Contig c;
	
	public CND_OnContig(Contig contig)
	{
		c = contig;
	}

	public boolean passes(String varRecord)
	{
		if (c == null) return false;
		if(varRecord == null || varRecord.isEmpty()) return false;
		String[] fields = varRecord.split("\t");
		if (fields.length < 4) return false;
		String chr = fields[3];
		if(c.hasName(chr)) return true;
		
		return false;
	}
	
	public boolean passes(DBVariant varRecord)
	{
		if (c == null) return false;
		if (varRecord == null) return false;
		Contig vc = varRecord.getChrom();
		if (vc == null) return false;
		return c.equals(vc);
	}
	
}
