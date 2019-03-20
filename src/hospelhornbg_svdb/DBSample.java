package hospelhornbg_svdb;

import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

import hospelhornbg_segregation.Population;
import waffleoRai_Utils.BinFieldSize;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.SerializedString;

public class DBSample {
	
	private String name;
	private int uid;
	private Set<Population> popFlags;
	private String famName;
	
	private DBSample()
	{
		popFlags = new TreeSet<Population>();
	}
	
	public static DBSample readFromSampleTable(FileBuffer file, long offset)
	{
		//UID[4]
		//Name [2x2 VLS]
		//Pop Flags [8]
		//FamName [2x2 VLS]
		DBSample s = new DBSample();
		
		s.uid = file.intFromFile(offset); offset += 4;
		SerializedString ss = file.readVariableLengthString(offset, BinFieldSize.WORD, 2);
		s.name = ss.getString();
		offset += ss.getSizeOnDisk();
		
		//Pop flags
		int flags = 0;
		for(int i = 0; i < 8; i++)
		{
			byte b = file.getByte(offset); offset++;
			if (b == 0xFF) break;
			int p = Byte.toUnsignedInt(b);
			Population pop = Population.getPopulation(p);
			if (pop != null) s.popFlags.add(pop);
			flags++;
		}
		while(flags < 8) {offset++; flags++;}
		
		ss = file.readVariableLengthString(offset, BinFieldSize.WORD, 2);
		s.famName = ss.getString();
		//offset += ss.getSizeOnDisk();
		
		return s;
	}
	
	public static DBSample createNewSample(String name)
	{
		DBSample s = new DBSample();
		s.name = name;
		s.regenerateUID();
		return s;
	}
	
	public int getUID()
	{
		return uid;
	}
	
	public void regenerateUID()
	{
		Random r = new Random();
		if(name == null) uid = r.nextInt();
		uid = name.hashCode() ^ r.nextInt();
	}

	public String getName()
	{
		return name;
	}
	
	public String toString()
	{
		return name;
	}
	
	public void addPopulationFlag(Population p)
	{
		popFlags.add(p);
	}
	
	public void clearPopulationFlags()
	{
		popFlags.clear();
	}
	
	public List<Population> getPopulationFlags()
	{
		List<Population> plist = new LinkedList<Population>();
		plist.addAll(popFlags);
		return plist;
	}

	public int calculateSerializedSize()
	{
		int sz = 4+8;
		sz += name.length() + 2;
		if (sz % 2 != 0) sz++;
		sz += famName.length() + 2;
		if (sz % 2 != 0) sz++;
		return sz;
	}
	
	public FileBuffer serializeSample()
	{
		FileBuffer sdata = new FileBuffer(calculateSerializedSize(), true);
		
		sdata.addToFile(uid);
		sdata.addVariableLengthString(name, BinFieldSize.WORD, 2);
		int fadded = 0;
		for (Population p : popFlags)
		{
			if (fadded >= 8) break;
			sdata.addToFile((byte)p.getIDNumber());
			fadded++;
		}
		while(fadded < 8)
		{
			sdata.addToFile((byte)0xFF);
			fadded++;
		}
		sdata.addVariableLengthString(famName, BinFieldSize.WORD, 2);
		
		return sdata;
	}
	
	public void setFamilyName(String name)
	{
		famName = name;
	}
	
	public void setSampleName(String newname)
	{
		name = newname;
	}
	
	public String getFamilyName()
	{
		return famName;
	}
	
	public boolean isInPopulation(Population p)
	{
		if (popFlags == null) return false;
		return popFlags.contains(p);
	}
	
}
