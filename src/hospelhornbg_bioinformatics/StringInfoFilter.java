package hospelhornbg_bioinformatics;

public class StringInfoFilter implements VariantFilter{

	private int mode;
	private int multi;
	private String match;
	private String fieldKey;
	
	public StringInfoFilter(String field, int condition, String reference)
	{
		if (field == null) throw new IllegalArgumentException();
		if (field.isEmpty()) throw new IllegalArgumentException();
		if (!conditionValid(condition)) throw new IllegalArgumentException();
		fieldKey = field;
		mode = condition;
		this.match = reference;
		multi = VariantFilter.MULTIFIELD_ANY;
	}
	
	public void set_multi_AND()
	{
		multi = VariantFilter.MULTIFIELD_ALL;
	}
	
	public void set_multi_OR()
	{
		multi = VariantFilter.MULTIFIELD_ANY;
	}
	
	public void set_multi_NOT()
	{
		multi = VariantFilter.MULTIFIELD_NONE;
	}
	
	public boolean conditionValid(int condition)
	{
		if (condition < 0 || condition > CONDITION_LESSOREQUAL) return false;
		return true;
	}
	
	private boolean valuePasses(String val)
	{
		switch(mode)
		{
		case CONDITION_EQUAL: return (match.equals(val));
		case CONDITION_NOTEQUAL:  return (!match.equals(val));
		case CONDITION_GREATERTHAN: return (val.compareTo(match) > 0);
		case CONDITION_GREATEROREQUAL: return (val.compareTo(match) >= 0);
		case CONDITION_LESSTHAN: return (val.compareTo(match) < 0);
		case CONDITION_LESSOREQUAL: return (val.compareTo(match) <= 0);
		}
		return false;
	}
	
	public boolean passes(Variant v)
	{
		String[] val = v.getInfoEntry(fieldKey);
		if (val == null) return false;
		if (val.length < 1) return false;
		switch(multi)
		{
		case MULTIFIELD_ALL:
			for (String s : val)
			{
				if (!valuePasses(s)) return false;
			}
			return true;
		case MULTIFIELD_ANY:
			for (String s : val)
			{
				if (valuePasses(s)) return true;
			}
			return false;
		case MULTIFIELD_NONE: 
			for (String s : val)
			{
				if (valuePasses(s)) return false;
			}
			return true;
		}
		return false;
	}
	
	public int getType()
	{
		return VariantFilter.FILTERTYPE_STRING;
	}
	
	public String toString()
	{
		String s = "";
		s += VariantFilter.multifieldStringENG(multi) + " ";
		s += fieldKey;
		s += " " + VariantFilter.conditionStringSymbol(mode) + " ";
		s += match;
		return s;
	}
	
	public boolean equals(Object o)
	{
		if (o == null) return false;
		if (this == o) return true;
		if (!(o instanceof StringInfoFilter)) return false;
		StringInfoFilter f = (StringInfoFilter)o;
		if (this.multi != f.multi) return false;
		if (this.mode != f.mode) return false;
		if (!this.match.equals(f.match)) return false;
		if (!this.fieldKey.equals(f.fieldKey)) return false;
		return true;
	}
	
	public int compareTo(VariantFilter o)
	{
		if (o instanceof StringInfoFilter)
		{
			StringInfoFilter f = (StringInfoFilter)o;
			if (this.fieldKey.equals(f.fieldKey))
			{
				if (this.match.equals(f.match))
				{
					if (this.multi == f.multi)
					{
						if (this.mode == f.mode) return 0;
						else return this.mode - f.mode;
					}
					else return this.multi - f.multi;
				}
				else return this.match.compareTo(f.match);
			}
			else return this.fieldKey.compareTo(f.fieldKey);
		}
		else return this.getType() - o.getType();
	}

	
}
