package hospelhornbg_segregation;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import hospelhornbg_bioinformatics.AffectedStatus;
import hospelhornbg_bioinformatics.Sex;

public class Individual {

	private Individual iParent1;
	private Individual iParent2;
	
	private Map<String, Individual> iChildren;
	
	private String sName;
	private Sex eSex;
	private AffectedStatus eAffected;
	
	public Individual(String name)
	{
		sName = name;
		iParent1 = null;
		iParent2 = null;
		//iChildren = new HashSet<Individual>();
		iChildren = new HashMap<String, Individual>();
		eSex = Sex.UNKNOWN;
		eAffected = AffectedStatus.UNKNOWN;
	}
	
	public String getName()
	{
		return sName;
	}
	
	public Sex getSex()
	{
		return eSex;
	}
	
	public AffectedStatus getAffectedStatus()
	{
		return eAffected;
	}
	
	public boolean isAffected()
	{
		return eAffected == AffectedStatus.AFFECTED || eAffected == AffectedStatus.PARTIALLY_AFFECTED;
	}
	
	public Individual getParent1()
	{
		return iParent1;
	}
	
	public Individual getParent2()
	{
		return iParent2;
	}
	
	public Individual getChild(String childID)
	{
		return iChildren.get(childID);
	}
	
	public Set<Individual> getChildren()
	{
		Set<Individual> set = new HashSet<Individual>();
		set.addAll(iChildren.values());
		return set;
	}
	
	public void setSex(Sex s)
	{
		eSex = s;
	}
	
	public void setAffectedStatus(AffectedStatus a)
	{
		eAffected = a;
	}
	
	public void setParent1(Individual p1)
	{
		Individual oldp1 = iParent1;
		iParent1 = p1;
		if (oldp1 != null) oldp1.removeChild(this);
		if (p1 != null) p1.addChild(this);
	}
	
	public void setParent2(Individual p2)
	{
		Individual oldp2 = iParent2;
		iParent2 = p2;
		if (oldp2 != null) oldp2.removeChild(this);
		if (p2 != null) p2.addChild(this);
	}
	
	protected void addChild(Individual child)
	{
		if (child == null) return;
		iChildren.put(child.getName(), child);
	}
	
	protected void removeChild(Individual child)
	{
		if (child == null) return;
		iChildren.remove(child.getName());
	}
	
	public String getENGString_sex()
	{
		switch(eSex)
		{
		case FEMALE:
			return "F";
		case MALE:
			return "M";
		case OTHER:
			return "O";
		case UNKNOWN:
			return "U";
		default:
			return "U";
		
		}
	}
	
	public String getENGString_affected()
	{
		switch(eAffected)
		{
		case AFFECTED:
			return "AFF";
		case PARTIALLY_AFFECTED:
			return "PAF";
		case POSSIBLY_AFFECTED:
			return "MAF";
		case UNAFFECTED:
			return "UAF";
		case UNKNOWN:
			return "UNK";
		default:
			return "UNK";
		}
	}
	
	public Individual getMother()
	{
		if (iParent1 != null && iParent1.eSex == Sex.FEMALE) return iParent1;
		if (iParent2 != null && iParent2.eSex == Sex.FEMALE) return iParent2;
		return null;
	}
	
	public Individual getFather()
	{
		if (iParent1 != null && iParent1.eSex == Sex.MALE) return iParent1;
		if (iParent2 != null && iParent2.eSex == Sex.MALE) return iParent2;
		return null;
	}
	
}
