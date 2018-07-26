package hospelhornbg_bioinformatics;

public class FloatInfoFilter implements VariantFilter{

	private int mode;
	private int multi;
	private double threshold;
	private String fieldKey;
	
	public FloatInfoFilter(String field, int condition, double threshold)
	{
		if (field == null) throw new IllegalArgumentException();
		if (field.isEmpty()) throw new IllegalArgumentException();
		if (!conditionValid(condition)) throw new IllegalArgumentException();
		fieldKey = field;
		mode = condition;
		this.threshold = threshold;
		multi = VariantFilter.MULTIFIELD_ALL;
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
	
	private boolean valuePasses(double val)
	{
		switch(mode)
		{
		case CONDITION_EQUAL: return (val == threshold);
		case CONDITION_NOTEQUAL:  return (val != threshold);
		case CONDITION_GREATERTHAN: return (val > threshold);
		case CONDITION_GREATEROREQUAL: return (val >= threshold);
		case CONDITION_LESSTHAN: return (val < threshold);
		case CONDITION_LESSOREQUAL: return (val <= threshold);
		}
		return false;
	}
	
	public boolean passes(Variant v)
	{
		double[] val = v.getFloatInfoEntry(fieldKey);
		if (val == null) return false;
		if (val.length < 1) return false;
		switch(multi)
		{
		case MULTIFIELD_ALL:
			for (double i : val)
			{
				if (!valuePasses(i)) return false;
			}
			return true;
		case MULTIFIELD_ANY:
			for (double i : val)
			{
				if (valuePasses(i)) return true;
			}
			return false;
		case MULTIFIELD_NONE: 
			for (double i : val)
			{
				if (valuePasses(i)) return false;
			}
			return true;
		}
		return false;
	}
	
	public int getType()
	{
		return VariantFilter.FILTERTYPE_FLOAT;
	}
	
	public String toString()
	{
		String s = "";
		s += VariantFilter.multifieldStringENG(multi) + " ";
		s += fieldKey;
		s += " " + VariantFilter.conditionStringSymbol(mode) + " ";
		s += Double.toString(threshold);
		return s;
	}

	public boolean equals(Object o)
	{
		if (o == null) return false;
		if (this == o) return true;
		if (!(o instanceof FloatInfoFilter)) return false;
		FloatInfoFilter f = (FloatInfoFilter)o;
		if (this.multi != f.multi) return false;
		if (this.mode != f.mode) return false;
		if (this.threshold != f.threshold) return false;
		if (!this.fieldKey.equals(f.fieldKey)) return false;
		return true;
	}
	
	public int compareTo(VariantFilter o)
	{
		if (o instanceof FloatInfoFilter)
		{
			FloatInfoFilter f = (FloatInfoFilter)o;
			if (this.fieldKey.equals(f.fieldKey))
			{
				if (this.threshold == f.threshold)
				{
					if (this.multi == f.multi)
					{
						if (this.mode == f.mode) return 0;
						else return this.mode - f.mode;
					}
					else return this.multi - f.multi;
				}
				else
				{
					if (this.threshold > f.threshold) return 1;
					else return -1;
				}
			}
			else return this.fieldKey.compareTo(f.fieldKey);
		}
		else return this.getType() - o.getType();
	}
	
}
