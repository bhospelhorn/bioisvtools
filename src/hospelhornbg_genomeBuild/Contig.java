package hospelhornbg_genomeBuild;

import java.awt.Point;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import waffleoRai_Utils.FileBuffer;

public class Contig implements Comparable<Contig>{

	public static final int SORTCLASS_AUTOSOME = 0;
	public static final int SORTCLASS_SEXCHROM = 1;
	public static final int SORTCLASS_MITO = 2;
	public static final int SORTCLASS_CONTIG = 3;
	public static final int SORTCLASS_COMPOUND = 4;
	public static final int SORTCLASS_UNKNOWN = 5;
	
	public static final int SORTMODE_NORMAL = 0;
	public static final int SORTMODE_NORMAL_CTGLEN = 1;
	private static int SORTMODE = SORTMODE_NORMAL;
	
	private Set<String> names;
	private String UCSC_name;
	private String UDP_name;
	
	private int sortClass;
	
	private int guid;
	private long length;
	
	public Contig()
	{
		names = new HashSet<String>();
		length = -1;
		UCSC_name = null;
		UDP_name = null;
		sortClass = SORTCLASS_UNKNOWN;
	}
	
	public Collection<String> getAllNames()
	{
		Collection<String> copy = new HashSet<String>();
		copy.addAll(names);
		copy.add(UCSC_name);
		copy.add(UDP_name);
		return copy;
	}
	
	public String getUDPName()
	{
		return UDP_name;
	}
	
	public String getUCSCName()
	{
		return UCSC_name;
	}
	
	public int getUID()
	{
		return guid;
	}
	
	public long getLength()
	{
		return length;
	}
	
	public int getType()
	{
		return sortClass;
	}
	
	public void addName(String name)
	{
		if (name == null) return;
		names.add(name);
	}
	
	public void removeName(String name)
	{
		names.remove(name);
	}
	
	public void setUID(int uid)
	{
		guid = uid;
	}
	
	public void setUDPName(String name)
	{
		UDP_name = name;
	}
	
	public void setUCSCName(String name)
	{
		UCSC_name = name;
	}
	
	public void setLength(long len)
	{
		length = len;
	}
	
	public void setType(int type)
	{
		if (type > 3) return;
		sortClass = type;
	}
	
	public FileBuffer serialize()
	{
		return serializeWithPARs(null);
	}
	
	public FileBuffer serializeWithPARs(Point[] pars)
	{
		FileBuffer spars = null;
		int psz = 0;
		if (pars != null)
		{
			psz = 4 + (8 * pars.length);
			spars = new FileBuffer(psz, true);
			spars.addToFile(pars.length);
			for (int i = 0; i < pars.length; i++)
			{
				spars.addToFile(pars[i].x);
				spars.addToFile(pars[i].y);
			}
		}
		
		FileBuffer altNames = new FileBuffer(names.size() * 32, true);
		
		for (String n : names)
		{
			short len = (short)n.length();
			altNames.addToFile(len);
			altNames.printASCIIToFile(n);
			if (len % 2 != 0) altNames.addToFile((byte)0x00);
		}
		
		int namesSize = (int)altNames.getFileSize();
		FileBuffer myContig = new FileBuffer(100 + namesSize, true);
		
		int bsz = 100 + namesSize + psz;
		myContig.addToFile(bsz); //Contig block size excluding actual size record
		myContig.addToFile(guid);
		myContig.addToFile(length); //Contig length
		myContig.addToFile(sortClass); //Type
		//myContig.addToFile(localUID); //Contig ID
		if (spars != null) myContig.addToFile(spars);
		
		String name = UCSC_name;
		if (name != null)
		{
			if (name.length() > 64) name = name.substring(0, 63);
			myContig.printASCIIToFile(name);
			if (name.length() < 64)
			{
				for (int i = name.length(); i < 64; i++) myContig.addToFile((byte)0x00);
			}
		}
		else
		{
			for (int i = 0; i < 64; i++) myContig.addToFile((byte)0x00);
		}
		
		name = UDP_name;
		if (name != null)
		{
			if (name.length() > 16) name = name.substring(0, 15);
			myContig.printASCIIToFile(name);
			if (name.length() < 16)
			{
				for (int i = name.length(); i < 16; i++) myContig.addToFile((byte)0x00);
			}	
		}
		else
		{
			for (int i = 0; i < 16; i++) myContig.addToFile((byte)0x00);
		}
		
		myContig.addToFile(names.size());
		
		for (int i = 0; i < namesSize; i++){
			myContig.addToFile(altNames.getByte(i));
		}
		
		return myContig;
	}
	
	public boolean equals(Object o)
	{
		if (o == null) return false;
		if (this == o) return true;
		if (!(o instanceof Contig)) return false;
		
		Contig c = (Contig)o;
		
		if (this.length != c.length) return false;
		if (this.sortClass != c.sortClass) return false;
		if (this.UCSC_name == null)
		{
			System.err.println("Contig.equals || UCSC name field is null!!");
			System.err.println("Contig.equals || Contig info:");
			printInfo();
			throw new NullPointerException();
		}
		if (!(this.UCSC_name.equals(c.UCSC_name))) return false;
		if (!(this.UDP_name.equals(c.UDP_name))) return false;
		if (this.names.size() != c.names.size()) return false;
		if (!this.names.containsAll(c.names)) return false;
		if (!c.names.containsAll(this.names)) return false;
		
		return true;
	}
	
	public int hashCode()
	{
		return UCSC_name.hashCode();
	}

	public int compareTo(Contig o) 
	{
		if (o == null) return 1;
		//First, compare sort type
		if (this.sortClass != o.sortClass) return this.sortClass - o.sortClass;
		if (this.sortClass == SORTCLASS_AUTOSOME)
		{
				//Try to parse it as a number first
			int num1 = -1;
			try {num1 = Integer.parseInt(this.UDP_name);}
			catch (NumberFormatException e){num1 = -1;}
			
			int num2 = -1;
			try {num2 = Integer.parseInt(o.UDP_name);}
			catch (NumberFormatException e){num2 = -1;}
			
			if (num1 >= 0 && num2 >= 0) return num1 - num2;
			if (num1 < 0 && num2 >= 0) return 1;
			if (num1 >= 0 && num2 < 0) return -1;
			return this.UDP_name.compareTo(o.UDP_name);
		}
		
		//Check sort mode
		if(sortClass == SORTCLASS_CONTIG && SORTMODE == SORTMODE_NORMAL_CTGLEN)
		{
			return (int)(this.length - o.length);
		}
		
		//Else, normal
		if (this.UDP_name == null)
		{
			System.err.println("Contig.compareTo || [this] UDP name field is null!!");
			System.err.println("Contig.compareTo || Contig info:");
			printInfo();
			throw new NullPointerException();
		}
		if (o.UDP_name == null)
		{
			System.err.println("Contig.compareTo || [o] UDP name field is null!!");
			System.err.println("Contig.compareTo || Contig info:");
			o.printInfo();
			throw new NullPointerException();
		}
		return this.UDP_name.compareTo(o.UDP_name); //Just keep it alphabetical
	}
	
	public String toString()
	{
		return this.UDP_name;
	}
	
	public String printInfo()
	{
		String s = "";
		s += this.getUDPName() + "\t";
		s += this.getUCSCName() + "\t";
		switch (this.sortClass)
		{
		case SORTCLASS_AUTOSOME: s += "AUTOSOME\t"; break;
		case SORTCLASS_SEXCHROM: s += "SEXCHROM\t"; break;
		case SORTCLASS_MITO: s += "MITO\t"; break;
		case SORTCLASS_CONTIG: s += "CONTIG\t"; break;
		case SORTCLASS_COMPOUND: s += "COMPOUND\t"; break;
		case SORTCLASS_UNKNOWN: s += "UNK\t"; break;
		}
		s += this.getLength() + "\t";
		s += String.format("0x%08x", guid) + "\t";
		s += "[";
		int nNames = names.size();
		int i = 0;
		for (String n : names)
		{
			s += n;
			if (i < nNames-1) s += ",";
			i++;
		}
		s += "]";
		
		return s;
	}
	
	public static void setSortMode(int sortModeEnum)
	{
		SORTMODE = sortModeEnum;
	}
	
	public boolean hasName(String name)
	{
		if(name == null || name.isEmpty()) return false;
		if (name.equals(UDP_name)) return true;
		if (name.equals(UCSC_name)) return true;
		for(String s : names)
		{
			if(name.equals(s)) return true;
		}
		return false;
	}
	
}
