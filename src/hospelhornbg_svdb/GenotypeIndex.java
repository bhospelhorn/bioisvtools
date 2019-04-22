package hospelhornbg_svdb;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import hospelhornbg_svdb.GenotypeTable.IndexUnit;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class GenotypeIndex {
	
	//genox format...
	//Magic [8]
	//Version [4]
	//Reserved [4]
	
	//Records...
	//	VarUID [4]
	//	Offset [8]
	//	Length [4]
	
	public static final String MAGIC = "svdbGIDX";
	public static final int VERSION = 1;
	public static final String EXT = "genox";
	
	public static class IndexRecord
	{
		//private String tablePath;
		private long offset;
		private int recLen;
		
		public IndexRecord(long off, int len)
		{
			offset = off;
			recLen = len;
		}
		
		public long getOffset()
		{
			return offset;
		}
		
		public int getLength()
		{
			return recLen;
		}
	}
	
	//private String tablePath;
	private Map<Integer, IndexRecord> recordMap;
	
	private GenotypeIndex()
	{
		//tablePath = tblPath;
		recordMap = new TreeMap<Integer, IndexRecord>();
	}
	
	public static GenotypeIndex readIndex(String indexPath) throws IOException, UnsupportedFileTypeException
	{
		FileBuffer idxfile = FileBuffer.createBuffer(indexPath, true);
		long cpos = idxfile.findString(0, 0x10, MAGIC);
		if (cpos != 0) throw new FileBuffer.UnsupportedFileTypeException();
		
		GenotypeIndex idx = new GenotypeIndex();
		cpos = 16; //Skip rest of header for now
		long fsz = idxfile.getFileSize();
		while (cpos < fsz)
		{
			int varid = idxfile.intFromFile(cpos); cpos += 4;
			long off = idxfile.longFromFile(cpos); cpos += 8;
			int len = idxfile.intFromFile(cpos); cpos += 4;
			idx.recordMap.put(varid, new IndexRecord(off, len));
		}
		
		return idx;
	}
	
	public static GenotypeIndex generateIndex(GenotypeTable tbl) throws IOException
	{
		GenotypeIndex idx = new GenotypeIndex();
		
		tbl.openStream();
		IndexUnit iu = null;
		while((iu = tbl.getIndexingInfoForNextVariant()) != null)
		{
			idx.recordMap.put(iu.varID, new IndexRecord(iu.offset, iu.len));
		}
		tbl.closeStream();
		return idx;
	}
	
	public void writeIndex(String indexPath) throws IOException
	{
		int nrec = recordMap.size();
		int size = 16 + (nrec * 16);
		FileBuffer out = new FileBuffer(size, true);
		
		//Header
		out.printASCIIToFile(MAGIC);
		out.addToFile(VERSION);
		out.addToFile(0);
		
		//Order records...
		Set<Integer> keys = recordMap.keySet();
		List<Integer> okeys = new ArrayList<Integer>(nrec + 1);
		okeys.addAll(keys);
		Collections.sort(okeys);
		
		//Output records...
		for(Integer k : okeys)
		{
			IndexRecord r = recordMap.get(k);
			out.addToFile(k);
			out.addToFile(r.getOffset());
			out.addToFile(r.getLength());
		}
		
		out.writeFile(indexPath);
		
	}
	
	public IndexRecord getIndexRecord(int variantUID)
	{
		return recordMap.get(variantUID);
	}
	
	public List<Integer> getAllVariantIDs()
	{
		Set<Integer> keys = recordMap.keySet();
		List<Integer> okeys = new ArrayList<Integer>(recordMap.size() + 1);
		okeys.addAll(keys);
		Collections.sort(okeys);
		return okeys;
	}

}
