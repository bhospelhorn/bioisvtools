package hospelhornbg_bioinformatics;

public interface VariantFilter extends Comparable<VariantFilter>{
	
	public static final int CONDITION_EQUAL = 0;
	public static final int CONDITION_NOTEQUAL = 1;
	public static final int CONDITION_GREATERTHAN = 2;
	public static final int CONDITION_GREATEROREQUAL = 3;
	public static final int CONDITION_LESSTHAN = 4;
	public static final int CONDITION_LESSOREQUAL = 5;
	
	public static final int MULTIFIELD_ALL = 0;
	public static final int MULTIFIELD_ANY = 1;
	public static final int MULTIFIELD_NONE = 2;
	
	public static final int FILTERTYPE_FLAG = 0;
	public static final int FILTERTYPE_INT = 1;
	public static final int FILTERTYPE_FLOAT = 2;
	public static final int FILTERTYPE_STRING = 3;

	public boolean passes(Variant v);
	public int getType();

 	public static String conditionStringSymbol(int condition)
	{
		switch(condition)
		{
		case CONDITION_EQUAL: return "=";
		case CONDITION_NOTEQUAL: return "!=";
		case CONDITION_GREATERTHAN: return ">";
		case CONDITION_GREATEROREQUAL: return ">=";
		case CONDITION_LESSTHAN: return "<";
		case CONDITION_LESSOREQUAL: return "<=";
		default: return "";
		}
	}
	
	public static String conditionStringENG(int condition)
	{
		switch(condition)
		{
		case CONDITION_EQUAL: return "is equal to";
		case CONDITION_NOTEQUAL: return "is not equal to";
		case CONDITION_GREATERTHAN: return "is greater than";
		case CONDITION_GREATEROREQUAL: return "is greater than or equal to";
		case CONDITION_LESSTHAN: return "is less than";
		case CONDITION_LESSOREQUAL: return "is less than or equal to";
		default: return "";
		}
	}
	
	public static String multifieldStringENG(int condition)
	{
		switch(condition)
		{
		case MULTIFIELD_ALL: return "all";
		case MULTIFIELD_ANY: return "any";
		case MULTIFIELD_NONE: return "no";
		default: return "";
		}
	}
	
}
