package hospelhornbg_segregation;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

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
	private Set<Population> population_tags;
	
	public FamilyMember(String samplename)
	{
		super(samplename);
		Random r = new Random();
		uid = samplename.hashCode() ^ r.nextInt();
		affected = new HashMap<String, AffectedStatus>();
		birthYear = -1;
		deathYear = -1;
		population_tags = new TreeSet<Population>();
	}
	
	public FamilyMember(Individual indiv)
	{
		super(indiv);
		uid = super.getName().hashCode();
		affected = new HashMap<String, AffectedStatus>();
		birthYear = -1;
		deathYear = -1;
		phenotypicSex = super.getSex();
		population_tags = new TreeSet<Population>();
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

	public FamilyMember getParent1AsFamilyMember()
	{
		//See if parent 1 is an instance of FamilyMember.
		//If not, return null
		Individual p1 = super.getParent1();
		if (p1 instanceof FamilyMember) return (FamilyMember)p1;
		return null;
	}
	
	public FamilyMember getParent2AsFamilyMember()
	{
		//See if parent 2 is an instance of FamilyMember.
		//If not, return null
		Individual p2 = super.getParent2();
		if (p2 instanceof FamilyMember) return (FamilyMember)p2;
		return null;
	}
	
	public Collection<FamilyMember> getChildrenAsFamilyMembers()
	{
		//Get any children that are FamilyMembers
		Set<Individual> children = super.getChildren();
		Set<FamilyMember> set = new HashSet<FamilyMember>();
		if (children != null && !children.isEmpty())
		{
			for (Individual c : children)
			{
				if (c instanceof FamilyMember)
				{
					set.add((FamilyMember)c);
				}
			}
		}
		return set;
	}
	
	private void dumpRelativesIntoSet(Set<FamilyMember> set)
	{
		set.add(this);
		FamilyMember p1 = this.getParent1AsFamilyMember();
		if (p1 != null) p1.dumpRelativesIntoSet(set);
		FamilyMember p2 = this.getParent2AsFamilyMember();
		if (p2 != null) p2.dumpRelativesIntoSet(set);
		Collection<FamilyMember> children = this.getChildrenAsFamilyMembers();
		if (!children.isEmpty())
		{
			for (FamilyMember c : children) c.dumpRelativesIntoSet(set);
		}
	}
	
	public Collection<FamilyMember> getAllRelativesAsFamilyMembers()
	{
		Set<FamilyMember> set = new HashSet<FamilyMember>();
		this.dumpRelativesIntoSet(set);
		return set;
	}
	
	public void setAffectedCondition(String pheno)
	{
		super.setAffectedStatus(this.getAffectedStatus(pheno));
	}
	
	public void setAffectedStatus(String pheno, AffectedStatus status)
	{
		this.affected.put(pheno, status);
		
	}
	
	public void removeAffectedStatus(String pheno)
	{
		this.affected.remove(pheno);
	}
	
	public void setPhenotypicSex(Sex s)
	{
		this.phenotypicSex = s;
	}
	
	public void setUID(int UID)
	{
		this.uid = UID;
	}
	
	public void setFirstName(String s)
	{
		this.firstName = s;
	}
	
	public void setLastName(String s)
	{
		this.lastName = s;
	}
	
	public void setBirthYear(int year)
	{
		this.birthYear = year;
	}
	
	public void setDeathYear(int year)
	{
		this.deathYear = year;
	}
	
	public void addMiddleName(String s)
	{
		if (this.middleNames == null)
		{
			this.middleNames = new LinkedList<String>();
		}
		middleNames.add(s);
	}
	
	public void clearMiddleNames()
	{
		this.middleNames = null;
	}
	
	public Collection<Population> getPopulationTags()
	{
		List<Population> list = new ArrayList<Population>(population_tags.size() + 1);
		list.addAll(population_tags);
		return list;
	}
	
	public void addPopulationTag(Population p)
	{
		population_tags.add(p);
	}
	
	public void clearPopulationTags()
	{
		population_tags.clear();
	}
	
	public boolean hasPopulationTag(Population p)
	{
		return population_tags.contains(p);
	}
	
}
