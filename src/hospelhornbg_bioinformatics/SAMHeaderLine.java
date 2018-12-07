package hospelhornbg_bioinformatics;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;
import waffleoRai_Utils.TallyMap;

public class SAMHeaderLine {
	
	private String key;
	private Map<String, String> linefields;
	private TallyMap fieldCounts;
	private List<String> fieldOrder;
	
	public SAMHeaderLine()
	{
		key = null;
		linefields = new HashMap<String,String>();
		fieldCounts = new TallyMap();
		fieldOrder = new LinkedList<String>();
	}
	
	public SAMHeaderLine(String line) throws UnsupportedFileTypeException
	{
		if (line == null) throw new FileBuffer.UnsupportedFileTypeException();
		if (line.isEmpty()) throw new FileBuffer.UnsupportedFileTypeException();
		if (line.charAt(0) != '@') throw new FileBuffer.UnsupportedFileTypeException();
		
		String[] fields = line.split("\t");
		
		if (fields.length < 1) throw new FileBuffer.UnsupportedFileTypeException();
		linefields = new HashMap<String, String>();
		fieldCounts = new TallyMap();
		fieldOrder = new LinkedList<String>();
		
		key = fields[0].substring(1);
		int nfields = fields.length;
		for (int i = 1; i < nfields; i++)
		{
			String sf[] = fields[i].split(":");
			fieldCounts.increment(sf[0].hashCode());
			if (sf.length < 2) continue; //Eat
			linefields.put(sf[0], sf[1]);
			fieldOrder.add(sf[0]);
		}
	}
	
	public String getKey()
	{
		return key;
	}
	
	public String getValue(String fieldKey)
	{
		return linefields.get(fieldKey);
	}
	
	public int getFieldCount(String fieldKey)
	{
		return fieldCounts.getCount(fieldKey.hashCode());
	}

	public void setKey(String s)
	{
		key = s;
	}
	
	public void addField(String fkey, String fval)
	{
		linefields.put(fkey, fval);
		fieldCounts.increment(fkey.hashCode());
		fieldOrder.add(fkey);
	}

	public void setField(String fkey, String fval)
	{
		linefields.put(fkey, fval);
	}
	
	private void cleanFieldList()
	{
		List<String> nlist = new LinkedList<String>();
		for(String k : fieldOrder)
		{
			if(nlist.contains(k)) continue;
			nlist.add(k);
		}
		fieldOrder = nlist;
	}
	
	public String toString()
	{
		cleanFieldList();
		String s = "@" + key;
		for(String k : fieldOrder)
		{
			s += "\t" + k + ":" + linefields.get(k);
		}
		return s;
	}
	
}
