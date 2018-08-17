package hospelhornbg_segregation;

import hospelhornbg_bioinformatics.Sex;

public class DescendantRelationship implements Relationship{
	
	private Individual iTarget;
	private Individual iRelative;
	
	private int nGenerations;
	
	public DescendantRelationship(Individual t, Individual r, int generations)
	{
		iTarget = t;
		iRelative = r;
		nGenerations = generations;
	}
	
	public Individual getTarget()
	{
		return iTarget;
	}
	
	public Individual getRelative()
	{
		return iRelative;
	}
	
	public Individual getCommonAncestor()
	{
		if (nGenerations >= 0) return iRelative;
		else return iTarget;
	}
	
	public Individual[] getCommonAncestors()
	{
		Individual[] iarr = new Individual[1];
		iarr[0] = getCommonAncestor();
		return iarr;
	}
	
	public boolean isAncestor()
	{
		//If relative is ancestor of target
		return nGenerations > 0;
	}
	
	public boolean isDescendant()
	{
		return nGenerations < 0;
	}
	
	public boolean isParent()
	{
		return nGenerations == 1;
	}
	
	public boolean isChild()
	{
		return nGenerations == -1;
	}
	
	public boolean isSibling()
	{
		return false;
	}
	
	public boolean isHalfSibling()
	{
		return false;
	}
	
	public boolean isSelf()
	{
		return nGenerations == 0;
	}
	
	public boolean commonAncestorCouple()
	{
		return false;
	}
	
	public int targetGenerationsToCommonAncestor()
	{
		if (nGenerations > 0) return nGenerations;
		else return 0;
	}
	
	public int relativeGenerationsToCommonAncestor()
	{
		if (nGenerations < 0) return nGenerations * -1;
		else return 0;
	}
	
	public String toString_English()
	{
		if (nGenerations == 0) return "Self";
		if (nGenerations > 0)
		{
			//Relative is ancestor of target
			if (nGenerations == 1)
			{
				if (iRelative.getSex() == Sex.FEMALE) return "Mother";
				else if (iRelative.getSex() == Sex.MALE) return "Father";
				else return "Parent";
			}
			else
			{
				String s = "Grand";
				if (iRelative.getSex() == Sex.FEMALE) s += "mother";
				else if (iRelative.getSex() == Sex.MALE) s += "father";
				else s += "parent";
				//2 grandparent, 3 great grandparent, 4 great great grandparent
				if (nGenerations >= 3)
				{
					int greats = nGenerations - 2;
					if (greats < 3)
					{
						for (int i = 0; i < greats; i++) s = "Great " + s;
					}
					else
					{
						s = "Great[" + greats + "] " + s;
					}
				}
				return s;
			}
		}
		else
		{
			//Target is ancestor of relative
			int offset = nGenerations * -1;
			if (offset == 1)
			{
				if (iRelative.getSex() == Sex.FEMALE) return "Daughter";
				else if (iRelative.getSex() == Sex.MALE) return "Son";
				else return "Child";
			}
			else
			{
				String s = "Grand";
				if (iRelative.getSex() == Sex.FEMALE) s += "daughter";
				else if (iRelative.getSex() == Sex.MALE) s += "son";
				else s += "child";
				//2 grandparent, 3 great grandparent, 4 great great grandparent
				if (offset >= 3)
				{
					int greats = offset - 2;
					if (greats < 3)
					{
						for (int i = 0; i < greats; i++) s = "Great " + s;
					}
					else
					{
						s = "Great[" + greats + "] " + s;
					}
				}
				return s;
			}
		}
	}
	

}
