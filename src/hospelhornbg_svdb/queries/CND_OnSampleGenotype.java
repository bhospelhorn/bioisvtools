package hospelhornbg_svdb.queries;

import java.util.Collection;
import java.util.List;

import hospelhornbg_svdb.DBVariant;
import hospelhornbg_svdb.QueryCondition;
import hospelhornbg_svdb.SVDBGenotype;
import hospelhornbg_svdb.SVDBGenotype.SVDBAllele;
import hospelhornbg_svdb.SVDatabase;

public class CND_OnSampleGenotype implements QueryCondition{
	
	public static final int COMPARE_EQUALS = 0;
	public static final int COMPARE_ATLEAST = 1;
	public static final int COMPARE_ATMOST = 2;
	public static final int COMPARE_LESSTHAN = 3;
	public static final int COMPARE_MORETHAN = 4;
	
	private SVDatabase db;
	private int sampleUID;
	private int alleleCount;
	private int compareMode;
	
	public CND_OnSampleGenotype(SVDatabase dbLink, int sampleID, int alleleCopies, int compareMode)
	{
		db = dbLink;
		sampleUID = sampleID;
		alleleCount = alleleCopies;
		this.compareMode = compareMode;
	}
	
	public boolean passes(String varRecord)
	{
		if (db == null) return false;
		if (varRecord == null || varRecord.isEmpty()) return false;
		String[] fields = varRecord.split("\t");
		if (fields.length < 2) return false;
		
		//Nab sample ID from the record
		//Field 2
		int sid = 0;
		try {sid = Integer.parseUnsignedInt(fields[1], 16);}
		catch(NumberFormatException e) {return false;}
		
		List<SVDBGenotype> glist = db.getGenotypesForVariant(sid);
		if(glist == null || glist.isEmpty()) return false;
		//See if there is a geno for this sample at all...
		SVDBGenotype mygeno = null;
		for(SVDBGenotype g : glist)
		{
			if (g.getIndividualUID() == sampleUID)
			{
				mygeno = g;
				break;
			}
		}
		if (mygeno == null) return false;
		Collection<SVDBAllele> alist = mygeno.getAlleles();
		for(SVDBAllele a : alist)
		{
			switch(compareMode)
			{
			case COMPARE_EQUALS:
				if(a.getAlleleCount() == alleleCount) return true;
				break;
			case COMPARE_ATLEAST:
				if(a.getAlleleCount() >= alleleCount) return true;
				break;
			case COMPARE_ATMOST:
				if(a.getAlleleCount() <= alleleCount) return true;
				break;
			case COMPARE_LESSTHAN:
				if(a.getAlleleCount() < alleleCount) return true;
				break;
			case COMPARE_MORETHAN:
				if(a.getAlleleCount() > alleleCount) return true;
				break;
			}
		}
		
		return false;
	}
	
	public boolean passes(DBVariant varRecord)
	{
		if (db == null) return false;
		if (varRecord == null) return false;
		
		List<SVDBGenotype> glist = db.getGenotypesForVariant(varRecord.getIntegerID());
		if(glist == null || glist.isEmpty()) return false;
		//See if there is a geno for this sample at all...
		SVDBGenotype mygeno = null;
		for(SVDBGenotype g : glist)
		{
			if (g.getIndividualUID() == sampleUID)
			{
				mygeno = g;
				break;
			}
		}
		if (mygeno == null) return false;
		Collection<SVDBAllele> alist = mygeno.getAlleles();
		for(SVDBAllele a : alist)
		{
			switch(compareMode)
			{
			case COMPARE_EQUALS:
				if(a.getAlleleCount() == alleleCount) return true;
				break;
			case COMPARE_ATLEAST:
				if(a.getAlleleCount() >= alleleCount) return true;
				break;
			case COMPARE_ATMOST:
				if(a.getAlleleCount() <= alleleCount) return true;
				break;
			case COMPARE_LESSTHAN:
				if(a.getAlleleCount() < alleleCount) return true;
				break;
			case COMPARE_MORETHAN:
				if(a.getAlleleCount() > alleleCount) return true;
				break;
			}
		}
		
		return false;
	}
	

}
