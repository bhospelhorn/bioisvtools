package hospelhornbg_svdb;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.ConcurrentSkipListMap;

import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.FileBuffer;

public class DBVariantTable {
	
	//Variant UID format (64-bit)
	//	Contig ID [2]
	//	Contig Section [2]
	//	Random [4]
	
	public static final int MAX_IDXREC_IN_MEM = 0xFFFF;
	public static final int MAX_VAR_IN_MEM = 2048;
	public static final String VINDEX_NAME = "vars.vidx";
	public static final String CINDEX_NAME = "vars.cidx";
	
	public static final String VAR_DIR = "var";
	
	public static final String VARTBL_EXT = "vdb";
	
	private static class LookupRecord
	{
		public SVType type;
		public long offset;
		public int recordSize;
	}
	
	private static class GenomeIndex
	{
		//For correlating variant UIDs to contig/pos!
		private ConcurrentMap<Contig, Integer> shortIDMap;
		
		public GenomeIndex(GenomeBuild gb)
		{
			//TODO
		}
		
		public long getVarIDPrefix(Contig c, int pos)
		{
			//TODO
			return 0;
		}
		
	}
	
	
	private static class VariantCache
	{
		//Index is sorted
		//This is so that in case the index gets too large to load
		//	fully into memory, it itself can easily be used like a tree
		//	for lookup.
		private GenomeBuild genome;
		
		private String srcDir;
		private boolean cachedIndex;
		
		private ConcurrentMap<Long, DBVariant> cache;
		private ConcurrentMap<Long, LookupRecord> index;
		
		private ConcurrentLinkedQueue<Long> cacheQueue;
		private ConcurrentLinkedQueue<Long> indexCacheQueue;
		
		public VariantCache(GenomeBuild gb, String dir, boolean cacheIndex) throws IOException
		{
			genome = gb;
			srcDir = dir;
			cache = new ConcurrentSkipListMap<Long, DBVariant>();
			index = new ConcurrentHashMap<Long, LookupRecord>();
			cacheQueue = new ConcurrentLinkedQueue<Long>();
			cachedIndex = cacheIndex;
			indexCacheQueue = new ConcurrentLinkedQueue<Long>();
			if(!cacheIndex) loadIndex();
		}
		
		private String getIndexPath()
		{
			return srcDir + File.separator + VINDEX_NAME;
		}
		
		private void loadIndex() throws IOException
		{
			//VIDX format (BE):
			//	Record [32]
			//		Var UID[8]
			//		FileOffset [8]
			//		Record Size [4]
			//		SV Type [4]
			//		Padding [8]
			
			String idxpath = getIndexPath();
			if(!FileBuffer.fileExists(idxpath)) return;
			
			long fsz = FileBuffer.fileSize(idxpath);
			long cpos = 0;
			FileBuffer idx = FileBuffer.createBuffer(idxpath, true);
			while(cpos < fsz)
			{
				long vid = idx.longFromFile(cpos); cpos += 8;
				long off = idx.longFromFile(cpos); cpos += 8;
				int sz = idx.intFromFile(cpos); cpos += 4;
				int etype = idx.intFromFile(cpos); cpos+=4;
				cpos += 8;
				LookupRecord lr = new LookupRecord();
				lr.offset = off;
				lr.recordSize = sz;
				lr.type = SVType.getTypeByID(etype);
				index.put(vid, lr);
			}
		}
		
		private LookupRecord searchForRecord(long varUID, long minOff, long maxOff, String idxpath)
		{
			long ssz = maxOff - minOff;
			long stoff = ssz / 2L;
			stoff &= ~(0x1FL); //To set to nearest 32
			stoff += minOff;
			
			//Read variant UID at that offset
			int idxid = -1;
			FileBuffer rec = null;
			try 
			{
				rec = new FileBuffer(idxpath, stoff, 32, true);
			} 
			catch (IOException e) 
			{
				e.printStackTrace();
				return null;
			}
			idxid = rec.intFromFile(0);
			
			//Compare retrieved id to target
			if(idxid == varUID)
			{
				//Parse and return this record
				long off = rec.longFromFile(8);
				int sz = rec.intFromFile(16);
				int etype = rec.intFromFile(20);
				LookupRecord lr = new LookupRecord();
				lr.offset = off;
				lr.recordSize = sz;
				lr.type = SVType.getTypeByID(etype);
				return lr;
			}
			else if (idxid < varUID)
			{
				//Must be after
				return searchForRecord(varUID, stoff+32, maxOff, idxpath);
			}
			else if (idxid > varUID)
			{
				//Must be before
				return searchForRecord(varUID, minOff, stoff, idxpath);
			}
			
			return null;
		}
		
		private LookupRecord getIndexRecord(long varUID)
		{
			//This is for disc only index lookup!
			String idxpath = getIndexPath();
			if (!FileBuffer.fileExists(idxpath)) return null;
			long fsz = FileBuffer.fileSize(idxpath);
			
			//Start at middle record
			//All records are 16 bytes
			return searchForRecord(varUID, 0, fsz, idxpath);
		}
		
		private String getFilename(SVType type)
		{
			String typeroot = type.name();
			return srcDir + File.separator + typeroot + "." + VARTBL_EXT;
		}
		
		public DBVariant getVariant(long varUID)
		{
			//TODO
			return null;
		}
		
		public List<DBVariant> getVariantsInRegion(Contig c, int start, int end)
		{
			//TODO
			return null;
		}
		
		public void writeMasterIndex()
		{
			//TODO
			//Write the vidx
		}
		
		public void writeRegionIndex()
		{
			//TODO
			//Write the cidx
		}
		
	}

}
