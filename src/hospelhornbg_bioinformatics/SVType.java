package hospelhornbg_bioinformatics;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

public enum SVType{
	
	BND("BND"),
	DEL("DEL"),
	DUP("DUP"),
	INS("INS"),
	INV("INV"),
	CNV("CNV"),
	TANDEM("DUP:TANDEM"),
	OTHER("OTHER"),
	BED_REGION("REG"),
	DELME("DEL:ME"),
	INSME("INS:ME"),
	TRA("TRA");
	
	private String abrv;
	
	private SVType (String tWritten)
	{
		abrv = tWritten;
	}
	
	public String getString()
	{
		return abrv;
	}
	
	public String toString()
	{
		return abrv;
	}
	
	public static SVType getType(String typeString)
	{
		if (typeString == null) return null;
		if (typeString.isEmpty()) return null;
		if (typeString.equals("BND")) return BND;
		else if (typeString.equals("DEL")) return DEL;
		else if (typeString.equals("DUP")) return DUP;
		else if (typeString.equals("INS")) return INS;
		else if (typeString.equals("INV")) return INV;
		else if (typeString.equals("CNV")) return CNV;
		else if (typeString.equals("DUP:TANDEM")) return TANDEM;
		else if (typeString.equals("DEL:ME")) return DELME;
		else if (typeString.equals("INS:ME")) return INSME;
		else if (typeString.equals("TRA")) return TRA;
		else if (typeString.equals("REG")) return BED_REGION;
		return OTHER;
	}
	
	public static Collection<SVType> allTypes()
	{
		Set<SVType> allTypes = new HashSet<SVType>();
		allTypes.add(DEL);
		allTypes.add(DELME);
		allTypes.add(BND);
		allTypes.add(DUP);
		allTypes.add(INS);
		allTypes.add(INSME);
		allTypes.add(INV);
		allTypes.add(CNV);
		allTypes.add(TANDEM);
		allTypes.add(OTHER);
		allTypes.add(BED_REGION);
		allTypes.add(TRA);
		return allTypes;
	}
	
}
