package hospelhornbg_bioinformatics;

public class Interval implements Comparable<Interval>{
	
	private int iStart;
	private int iEnd;
	
	public Interval(int st, int ed)
	{
		iStart = st;
		iEnd = ed;
	}
	
	public int getStart()
	{
		return iStart;
	}
	
	public int getEnd()
	{
		return iEnd;
	}

	@Override
	public int compareTo(Interval o) 
	{
		if (o == null) return 1;
		if (o == this) return 0;
		
		if(this.iStart != o.iStart) return this.iStart - o.iStart;
		if(this.iEnd != o.iEnd) return this.iEnd - o.iEnd;
		
		return 0;
	}
	
	public boolean equals(Object o)
	{
		if (o == null) return false;
		if (o == this) return true;
		if (!(o instanceof Interval)) return false;
		
		Interval i = (Interval)o;
		if (this.iStart != i.iStart) return false;
		if (this.iEnd != i.iEnd) return false;
		
		return true;
	}
	
	public int hashCode()
	{
		return iStart ^ iEnd;
	}

	public int getCenter()
	{
		return (iEnd + iStart)/2;
	}
	
}
