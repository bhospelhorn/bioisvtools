package hospelhornbg_bioinformatics;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public enum SVType{
	
	BND("BND", 0),
	DEL("DEL", 1),
	DUP("DUP", 2),
	INS("INS", 3),
	INV("INV", 4),
	CNV("CNV", 5),
	TANDEM("DUP:TANDEM", 6),
	OTHER("OTHER", 7),
	BED_REGION("REG", 8),
	DELME("DEL:ME", 9),
	INSME("INS:ME", 10),
	TRA("TRA", 11);
	
	private String abrv;
	private int eID;
	
	private SVType (String tWritten, int id)
	{
		abrv = tWritten;
		eID = id;
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
	
	public int getID()
	{
		return eID;
	}
	
	private static Map<Integer, SVType> idmap;
	
	private static void populateIDMap()
	{
		idmap = new HashMap<Integer, SVType>();
		for(SVType t : SVType.values())
		{
			idmap.put(t.getID(), t);
		}
	}
	
	public static SVType getTypeByID(int id)
	{
		if (idmap == null) populateIDMap();
		return idmap.get(id);
	}
	
}
