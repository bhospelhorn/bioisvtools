package hospelhornbg_bioinformatics;

public class FlagInfoFilter implements VariantFilter{

	private String fieldKey;
	
	public FlagInfoFilter(String field)
	{
		if (field == null) throw new IllegalArgumentException();
		if (field.isEmpty()) throw new IllegalArgumentException();
		fieldKey = field;
	}
	
	public boolean passes(Variant v)
	{
		return (v.getInfoFlag(fieldKey));
	}
	
	public int getType()
	{
		return VariantFilter.FILTERTYPE_FLAG;
	}
	
	public String toString()
	{
		String s = "";
		s += fieldKey;
		return s;
	}
	
	public boolean equals(Object o)
	{
		if (o == null) return false;
		if (this == o) return true;
		if (!(o instanceof FlagInfoFilter)) return false;
		FlagInfoFilter f = (FlagInfoFilter)o;
		return this.fieldKey.equals(f.fieldKey);
	}
	
	public int compareTo(VariantFilter o)
	{
		if (o instanceof FlagInfoFilter)
		{
			FlagInfoFilter f = (FlagInfoFilter)o;
			return this.fieldKey.compareTo(f.fieldKey);
		}
		else return this.getType() - o.getType();
	}

	
}
