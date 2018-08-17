package hospelhornbg_segregation;

import hospelhornbg_bioinformatics.Sex;

public class SingleRelationship implements Relationship{

	private Individual iTarget;
	private Individual iRelative;
	
	private Individual iCommonAncestor;
	
	private int nGenerations_target;
	private int nGenerations_relative;
	
	public SingleRelationship(Individual t, Individual r, Individual ancestor)
	{
		iTarget = t;
		iRelative = r;
		iCommonAncestor = ancestor;
		nGenerations_target = iTarget.getGenerationDifference(iCommonAncestor);
		nGenerations_relative = iRelative.getGenerationDifference(iCommonAncestor);
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
		return iCommonAncestor;
	}
	
	public Individual[] getCommonAncestors()
	{
		Individual[] iarr = new Individual[1];
		iarr[0] = getCommonAncestor();
		return iarr;
	}
	
	public boolean isAncestor()
	{
		return false;
	}
	
	public boolean isDescendant()
	{
		return false;
	}
	
	public boolean isParent()
	{
		return false;
	}
	
	public boolean isChild()
	{
		return false;
	}
	
	public boolean isSibling()
	{
		return false;
	}
	
	public boolean isHalfSibling()
	{
		return (nGenerations_target == 1) && (nGenerations_relative == 1);
	}
	
	public boolean isSelf()
	{
		return false;
	}
	
	public boolean commonAncestorCouple()
	{
		return false;
	}
	
	public int targetGenerationsToCommonAncestor()
	{
		return this.nGenerations_target;
	}
	
	public int relativeGenerationsToCommonAncestor()
	{
		return this.nGenerations_relative;
	}
	
	public String toString_English()
	{
		if (nGenerations_target == nGenerations_relative)
		{
			int gen = nGenerations_target;
			if (gen == 1)
			{
				//Half-siblings
				if (iRelative.getSex() == Sex.FEMALE) return "Half-Sister";
				else if (iRelative.getSex() == Sex.MALE) return "Half-Brother";
				else return "Half-Sibling";
			}
			else
			{
				//Some kind of cousin
				String s = "Half-Cousin"; //Made up this term, but WTH
				int dist = gen - 2;
				if (dist > 0)
				{
					if (dist < Relationship.DISTANCE_ENGLISH.length)
					{
						s = Relationship.DISTANCE_ENGLISH[dist] + " " + s;
					}
					else s = "Distant " + s;
				}
				return s;
			}
		}
		else
		{
			//See who's closer to common ancestor...
			if (nGenerations_target >= nGenerations_relative)
			{
				//Relative is in upper generation
				//See if target is descended from relative's half-sib
				if (nGenerations_relative == 1)
				{
					String s = "Half-Nt";
					if (iRelative.getSex() == Sex.FEMALE) s = "Half-Aunt";
					else if (iRelative.getSex() == Sex.MALE) s = "Half-Uncle";
					int gdiff = nGenerations_target - nGenerations_relative;
					if (gdiff > 1)
					{
						int greats = gdiff - 1;
						if (greats < 3)
						{
							for (int i = 0; i < greats; i++) s = "Great " + s;
						}
						else s = "Great[" + greats + "] " + s;
					}
					return s;
				}
				else
				{
					//Otherwise target is descended from relative's cousin to some degree
					String s = "Half-Cousin";
					int dist = nGenerations_relative - 2;
					if (dist > 0)
					{
						if (dist < Relationship.DISTANCE_ENGLISH.length)
						{
							s = Relationship.DISTANCE_ENGLISH[dist] + " " + s;
						}
						else s = "Distant " + s;
					}
					//Add removal
					int removal = nGenerations_target - nGenerations_relative;
					String rstr = "Many Times";
					if (removal <= 20)
					{
						rstr = Relationship.REMOVAL_ENGLISH[removal-1];
					}
					s += " " + rstr + " Removed";
					return s;
				}
			}
			else
			{
				//Target is in upper generation
				//See if relative is descended from target's half-sib
				if (nGenerations_relative == 1)
				{
					String s = "Half-Niece/Nephew";
					if (iRelative.getSex() == Sex.FEMALE) s = "Half-Niece";
					else if (iRelative.getSex() == Sex.MALE) s = "Half-Nephew";
					int gdiff = nGenerations_relative - nGenerations_target;
					if (gdiff > 1)
					{
						int greats = gdiff - 1;
						if (greats < 3)
						{
							for (int i = 0; i < greats; i++) s = "Great " + s;
						}
						else s = "Great[" + greats + "] " + s;
					}
					return s;
				}
				else
				{
					//Otherwise relative is descended from target's cousin to some degree
					String s = "Half-Cousin";
					int dist = nGenerations_target - 2;
					if (dist > 0)
					{
						if (dist < Relationship.DISTANCE_ENGLISH.length)
						{
							s = Relationship.DISTANCE_ENGLISH[dist] + " " + s;
						}
						else s = "Distant " + s;
					}
					//Add removal
					int removal = nGenerations_relative - nGenerations_target;
					String rstr = "Many Times";
					if (removal <= 20)
					{
						rstr = Relationship.REMOVAL_ENGLISH[removal-1];
					}
					s += " " + rstr + " Removed";
					return s;
				}
			}
		}
	}
	

}
