package waffleoRai_Utils;

import java.util.*;

/*
 * UPDATES
 * 
 * 2017.09.21
 * 	1.0.0 -> 1.1.0 | Fixed for compatibility with FileBuffer changes.
 */

/**
 * Huffman encoding frequency table. To use with Huffman class.
 * @author Blythe Hospelhorn
 * @version 1.1.0
 * @since September 21, 2017
 *
 */
public class HuffTable 
{
	private class tableEntry
	{
		private long frequency;
		private String HuffCode;
		
		public tableEntry()
		{
			this.frequency = 0;
			this.HuffCode = "";
		}

		public long getFreq()
		{
			return this.frequency;
		}
		public void setFreq(long freq)
		{
			this.frequency = freq;
		}
		public void incrementFreq()
		{
			this.frequency++;
		}
		public String getCode()
		{
			return this.HuffCode;
		}
		public void setCode(String newCode)
		{
			/*Check validity of code!*/
			this.HuffCode = newCode;
		}
		
	}
	
	public class HuffPoint implements Comparator<HuffPoint>, Comparable<HuffPoint>
	{
		private long symbol;
		private long frequency;
		private String code;
		
		public HuffPoint(long s, long f, String c)
		{
			this.symbol = s;
			this.frequency = f;
			this.code = c;
		}
		
		public long getSymbol()
		{
			return this.symbol;
		}
		public long getFreq()
		{
			return this.frequency;
		}
		public String getHuffCode()
		{
			return this.code;
		}
		
		public void setSymbol(long s)
		{
			this.symbol = s;
		}
		public void setFreq(long f)
		{
			this.frequency = f;
		}
		public void setCode(String c)
		{
			this.code = c;
		}
	
		public int compareTo(HuffPoint other)
		{
			/*Sorts by symbol*/
			return (int)(this.symbol - other.getSymbol());
		}
		
		public int compare(HuffPoint p1, HuffPoint p2)
		{
			/*Comparator - sorts by frequency*/
			return (int)(p1.getFreq() - p2.getFreq());
		}
	
	}
	
	private Map<Long, tableEntry> contents;
	private int bitDepth;
	private long fileSize;
	
	public HuffTable(int bits)
	{
		if (bits > 64) bits = 64;
		if (bits < 2) bits = 2;
		this.contents = new HashMap<Long, tableEntry>();
		this.bitDepth = bits;
		this.fileSize = 0;
	}
	
	public int getBitDepth()
	{
		return this.bitDepth;
	}
	
	public long getFileSize()
	{
		return this.fileSize;
	}
	
	public void setFileSize(long newSize)
	{
		this.fileSize = newSize;
	}
	
	public long topValue()
	{
		long max = 1L;
		return max << this.bitDepth;
	}
	
	public boolean entryExists(long key)
	{
		if (this.contents.containsKey(key)) return true;
		return false;
	}
	
	public int numberValidEntries()
	{
		int tot = 0;
		for (Map.Entry<Long, tableEntry> e:contents.entrySet())
		{
			if (e.getValue().getFreq() > 0) tot++;
		}
		return tot;
	}
	
 	public long getFrequency(long key)
	{
		if (!entryExists(key)) return 0;
		
		return this.contents.get(key).getFreq();
	}
 	
	public void setFrequency(long key, long newFreq)
	{
		if (!entryExists(key))
		{
			if (key > this.topValue()) return;
			else
			{
				tableEntry te = new tableEntry();
				te.setFreq(newFreq);
				this.contents.put(key, te);
				return;
			}
		}
		
		this.contents.get(key).setFreq(newFreq);
	}
	
	public void incrementFrequency(long key)
	{
		if (!entryExists(key))
		{
			if (key > this.topValue()) return;
			else
			{
				tableEntry te = new tableEntry();
				te.incrementFreq();
				this.contents.put(key, te);
				return;
			}
		}
		
		this.contents.get(key).incrementFreq();
	}
	
	public String getHuffCode(long key)
	{
		if (!entryExists(key)) return "";
		
		return this.contents.get(key).getCode();	
	}
	
	public void setHuffCode(long key, String newCode)
	{
		if (!entryExists(key)) return;
		
		this.contents.get(key).setCode(newCode);
	}
	
	public boolean entryEmpty(long key)
	{
		/*Returns whether entry's frequency is above 0.*/
		if (!entryExists(key)) return true;
		
		if (this.contents.get(key).getFreq() <= 0) return true;
		return false;
	}
	
	public List<HuffPoint> contentsToList(boolean sortByFreq)
	{
		List<HuffPoint> myList = new ArrayList<HuffPoint>();
		
		HuffPoint aPoint = null;
		long s = 0;
		long f = 0;
		String c = null;
		
		for (Map.Entry<Long, tableEntry> entry:this.contents.entrySet())
		{
			s = entry.getKey();
			f = entry.getValue().getFreq();
			c = entry.getValue().getCode();
			
			aPoint = new HuffPoint(s, f, c);
			
			myList.add(aPoint);
		}
		
		if (sortByFreq)
		{
			Collections.sort(myList, new HuffPoint(0, 0, null));
		}
		else
		{
			Collections.sort(myList);	
		}
		
		return myList;
	}

	public String toString()
	{
		String s = "";
		List<HuffPoint> elist = this.contentsToList(false);
		
		for (HuffPoint p:elist)
		{
			s += Long.toHexString(p.getSymbol()) + "\t";
			s += "[" + Long.toHexString(p.getFreq()) + "]\t";
			s += "(" + p.getHuffCode() + ")\n";
		}
		
		return s;
	}

}
