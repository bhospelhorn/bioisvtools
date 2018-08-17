package hospelhornbg_segregation;

public class SelfRelationship implements Relationship{
	
	private Individual iIndiv;
	
	public SelfRelationship(Individual i)
	{
		iIndiv = i;
	}
	
	public Individual getTarget()
	{
		return iIndiv;
	}
	
	public Individual getRelative()
	{
		return iIndiv;
	}
	
	public Individual getCommonAncestor()
	{
		return iIndiv;
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
		return false;
	}
	
	public boolean isSelf()
	{
		return true;
	}
	
	public boolean commonAncestorCouple()
	{
		return true;
	}
	
	public int targetGenerationsToCommonAncestor()
	{
		return 0;
	}
	
	public int relativeGenerationsToCommonAncestor()
	{
		return 0;
	}
	
	public String toString_English()
	{
		return "Self";
	}
	

}
