package hospelhornbg_segregation;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import hospelhornbg_bioinformatics.AffectedStatus;
import hospelhornbg_bioinformatics.Sex;

public class FamilyMember extends Individual{
	
	private int uid;
	
	private String firstName;
	private String lastName;
	private List<String> middleNames;
	
	private int birthYear;
	private int deathYear;
	
	private Sex phenotypicSex;
	
	private Map<String, AffectedStatus> affected;
	
	public FamilyMember(String samplename)
	{
		super(samplename);
		uid = samplename.hashCode();
		affected = new HashMap<String, AffectedStatus>();
		birthYear = -1;
		deathYear = -1;
	}
	
	public int getUID()
	{
		return uid;
	}
	
	public String getFirstName()
	{
		return firstName;
	}
	
	public String getLastName()
	{
		return lastName;
	}
	
	public List<String> getMiddleNames()
	{
		int ncount = 1;
		if (middleNames != null) ncount += middleNames.size();
		List<String> list = new ArrayList<String>(ncount);
		if (middleNames != null) list.addAll(middleNames);
		return list;
	}
	
	public int getBirthYear()
	{
		return birthYear;
	}
	
	public int getDeathYear()
	{
		return deathYear;
	}
	
	public Sex getPhenotypicSex()
	{
		return phenotypicSex;
	}
	
	public AffectedStatus getAffectedStatus(String pheno)
	{
		return affected.get(pheno);
	}

}
