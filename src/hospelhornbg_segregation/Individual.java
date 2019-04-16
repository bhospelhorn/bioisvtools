package hospelhornbg_segregation;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import hospelhornbg_bioinformatics.AffectedStatus;
import hospelhornbg_bioinformatics.Sex;

public class Individual implements Comparable<Individual>{

	private Individual iParent1;
	private Individual iParent2;
	
	private Map<String, Individual> iChildren;
	
	private String sName;
	private Sex eSex;
	private AffectedStatus eAffected;
	
	private int custom_x;
	private int custom_y;
	
	private Map<Individual, Integer> iAncestors;
	
	public Individual(String name)
	{
		sName = name;
		iParent1 = null;
		iParent2 = null;
		//iChildren = new HashSet<Individual>();
		iChildren = new HashMap<String, Individual>();
		eSex = Sex.UNKNOWN;
		eAffected = AffectedStatus.UNKNOWN;
		custom_x = 2;
		custom_y = 1;
	}
	
	protected Individual (Individual source)
	{
		iParent1 = source.iParent1;
		iParent2 = source.iParent2;
		sName = source.sName;
		eSex = source.eSex;
		eAffected = source.eAffected;
		iChildren = new HashMap<String, Individual>();
		Set<String> keyset = source.iChildren.keySet();
		for (String k : keyset) iChildren.put(k, source.iChildren.get(k));
		custom_x = source.custom_x;
		custom_y = source.custom_y;
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
	
	private void addParentsToAncestorTree(Individual i, int level)
	{
		Individual p1 = i.getParent1();
		Individual p2 = i.getParent2();
		if (p1 != null)
		{
			addParentsToAncestorTree(p1, level+1);
			iAncestors.put(p1, level);
		}
		if (p2 != null)
		{
			addParentsToAncestorTree(p2, level+1);
			iAncestors.put(p2, level);
		}
	}
	
	public Map<Individual, Integer> getAncestorTree()
	{
		if (iAncestors != null) return iAncestors;
		iAncestors = new HashMap<Individual, Integer>();
		iAncestors.put(this, 0);
		addParentsToAncestorTree(this, 1); //Recursive
		return iAncestors;
	}
	
	public void disposeOfAncestorTree()
	{
		iAncestors = null;
	}
	
	private void checkAncestor(List<Individual> list, Individual i, Map<Individual, Integer> othertree)
	{
		if (othertree.containsKey(i))
		{
			list.add(i);
			return;
		}
		
		Individual p1 = i.getParent1();
		Individual p2 = i.getParent2();
		
		if (p1 != null) checkAncestor(list, p1, othertree);
		if (p2 != null) checkAncestor(list, p2, othertree);
	}
	
	public List<Individual> getRootCommonAncestors(Individual other)
	{
		List<Individual> alist = new LinkedList<Individual>();
		if (other == null) return alist;
		
		Map<Individual, Integer> othertree = other.getAncestorTree();
		
		checkAncestor(alist, this, othertree);
		
		return alist;
	}

	public int getGenerationDifference(Individual ancestor)
	{
		Map<Individual, Integer> mytree = getAncestorTree();
		return mytree.get(ancestor);
	}
	
	public Relationship getRelationship(Individual other)
	{
		if (other == null) return null;
		if (other == this) return new SelfRelationship(this);
		List<Individual> alist = getRootCommonAncestors(other);
		if (alist == null) return null;
		if (alist.isEmpty()) return null;
		
		//See if this or other are in alist (direct descendant/ancestor)
		if (alist.contains(this))
		{
			//this is an ancestor of other
			int g = other.getGenerationDifference(this) * -1;
			return new DescendantRelationship(this, other, g);
		}
		else if (alist.contains(other))
		{
			//other is an ancestor of this
			int g = this.getGenerationDifference(other);
			return new DescendantRelationship(this, other, g);
		}
		
		//Otherwise, see if share one or two ancestors with the same gen offsets
		//Use ancestor(s) with smallest offsets
		Individual a1 = null;
		Individual a2 = null;
		int gmin = Integer.MAX_VALUE;
		for (Individual a : alist)
		{
			int gdiff = getGenerationDifference(a);
			if (gdiff < gmin)
			{
				a1 = a;
				a2 = null;
				gmin = gdiff;
			}
			else if (gdiff == gmin)
			{
				if (a2 == null) a2 = a;
				else System.err.println("Individual.getRelationship || Warning: More than 2 common ancestors of same generational distance found...");
			}
		}
		
		if (a1 == null) return null;
		if (a2 != null)
		{
			return new PairRelationship(this, other, a1, a2);
		}
		else
		{
			return new SingleRelationship(this, other, a1);
		}
		
	}
	
	public boolean equals(Object o)
	{
		return (this == o);
	}

	@Override
	public int compareTo(Individual o) 
	{
		if (o == null) return 1;
		if (o == this) return 0;
		return this.sName.compareTo(o.sName);
	}
	
	private void dumpRelativesIntoSet(Set<Individual> set)
	{
		//Dump self, parents, and children
		set.add(this);
		if (this.iParent1 != null) {
			iParent1.dumpRelativesIntoSet(set);
		}
		if (this.iParent2 != null) {
			iParent2.dumpRelativesIntoSet(set);
		}
		//Children
		if (iChildren != null && !iChildren.isEmpty())
		{
			Collection<Individual> children = iChildren.values();
			for (Individual c : children)
			{
				c.dumpRelativesIntoSet(set);
			}
		}
		
	}
	
	public Collection<Individual> getAllRelatives()
	{
		//Set includes this individual (recursive)
		Set<Individual> rset = new HashSet<Individual>();
		this.dumpRelativesIntoSet(rset);
		return rset;
	}
	
	public void setSexChromCount(int X, int Y)
	{
		this.custom_x = X;
		this.custom_y = Y;
	}
	
	public int getExpectedXCount()
	{
		if (eSex == Sex.FEMALE) return 2;
		if (eSex == Sex.MALE) return 1;
		return custom_x;
	}
	
	public int getExpectedYCount()
	{
		if (eSex == Sex.FEMALE) return 0;
		if (eSex == Sex.MALE) return 1;
		return custom_y;
	}
	
	protected void setSampleName(String newname)
	{
		this.sName = newname;
	}

}
