package hospelhornbg_segregation;

import hospelhornbg_bioinformatics.Sex;

public class PairRelationship implements Relationship{
	
	private Individual iTarget;
	private Individual iRelative;
	
	private Individual iCommonAncestor_1;
	private Individual iCommonAncestor_2;
	
	private int nGenerations_target;
	private int nGenerations_relative;
	
	public PairRelationship(Individual t, Individual r, Individual ancestor1, Individual ancestor2)
	{
		iTarget = t;
		iRelative = r;
		iCommonAncestor_1 = ancestor1;
		iCommonAncestor_2 = ancestor2;
		nGenerations_target = iTarget.getGenerationDifference(iCommonAncestor_1);
		nGenerations_relative = iRelative.getGenerationDifference(iCommonAncestor_1);
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
		return iCommonAncestor_1;
	}
	
	public Individual[] getCommonAncestors()
	{
		Individual[] iarr = new Individual[2];
		iarr[1] = iCommonAncestor_1;
		iarr[2] = iCommonAncestor_2;
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
		return (nGenerations_target == 1) && (nGenerations_relative == 1);
	}
	
	public boolean isHalfSibling()
	{
		return false;
	}
	
	public boolean isSelf()
	{
		return false;
	}
	
	public boolean commonAncestorCouple()
	{
		return true;
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
				if (iRelative.getSex() == Sex.FEMALE) return "Sister";
				else if (iRelative.getSex() == Sex.MALE) return "Brother";
				else return "Sibling";
			}
			else
			{
				//Some kind of cousin
				String s = "Cousin";
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
					String s = "Nt";
					if (iRelative.getSex() == Sex.FEMALE) s = "Aunt";
					else if (iRelative.getSex() == Sex.MALE) s = "Uncle";
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
					String s = "Cousin";
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
					String s = "Niece/Nephew";
					if (iRelative.getSex() == Sex.FEMALE) s = "Niece";
					else if (iRelative.getSex() == Sex.MALE) s = "Nephew";
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
					String s = "Cousin";
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
