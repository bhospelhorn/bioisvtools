package hospelhornbg_svdb;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Paths;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.ConcurrentSkipListMap;

import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.VCFReadStreamer;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.Gene;
import hospelhornbg_genomeBuild.GeneFunc;
import hospelhornbg_genomeBuild.GeneSet;
import hospelhornbg_genomeBuild.GenomeBuild;
import hospelhornbg_segregation.Family;
import hospelhornbg_segregation.FamilyMember;
import hospelhornbg_segregation.Population;
import hospelhornbg_svdb.DBVariant.ParsedVariant;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;
import waffleoRai_Utils.StreamBuffer;

public class DBVariantTable implements VariantTable{
	
	//Variant UID format (64-bit)
	//	Contig ID [2]
	//	Contig Section [2]
	//	Random [4]
	
	/*----- Constants -----*/
	
	public static boolean breaker = true;
	
	public static final String GENOTABLE_MAGIC = "GENO_TBL";
	public static final int GENOTABLE_VERSION = 2;
	
	public static final int MAX_IDXREC_IN_MEM = 0xFFFF;
	public static final int MAX_VAR_IN_MEM = 2048;
	public static final int MAX_GENO_IN_MEM = 2048;
	public static final String VINDEX_NAME = "vars.vidx";
	public static final String CINDEX_NAME = "vars.cidx";
	
	public static final String VAR_DIR = "var";
	
	public static final String VARTBL_EXT = "vdb";
	public static final int VARTBL_VERSION = 1;
	public static final String VARTBL_MAGIC = "VDBt";
	
	public static final String GENOTBL_EXT = "genot";
	public static final String GENOTBL_FILE = "genotypes";
	
	//public static final int NEWVAR_RESULT_FAIL = -1;
	//public static final int NEWVAR_RESULT_MERGED = 1;
	//public static final int NEWVAR_RESULT_NEW = 0;
	
	/*----- Inner Classes (Minor) -----*/
	
	private static class LookupRecord
	{
		public SVType type;
		public long offset;
		public int recordSize;
		
		public long genoOffset;
		public int genoRecSize;
	}
	
	private static class RegionIndex
	{
		//Right now, load fully into mem
		
		//cidx format (BE)
		//Record [16]
		//	Prefix [4]
		//	Start Record # [8]
		//	Record Count [4]
		
		public static class RegionRecord
		{
			public long startIndex;
			public int recordCount;
		}
		
		private String cidx_path;
		private ConcurrentMap<Integer, RegionRecord> index;
		
		public RegionIndex(String cidxPath) throws IOException
		{
			cidx_path = cidxPath;
			index = new ConcurrentSkipListMap<Integer, RegionRecord>();
			loadIndex();
		}
		
		private void loadIndex() throws IOException
		{
			if(!FileBuffer.fileExists(cidx_path)) return;
			long fsz = FileBuffer.fileSize(cidx_path);
			FileBuffer idx = FileBuffer.createBuffer(cidx_path, true);
			long cpos = 0;
			while(cpos < fsz)
			{
				int prefix = idx.intFromFile(cpos); cpos += 4;
				RegionRecord r = new RegionRecord();
				r.startIndex = idx.longFromFile(cpos); cpos += 8;
				r.recordCount = idx.intFromFile(cpos); cpos += 4;
				index.put(prefix, r);
			}
		}
		
		public RegionRecord getRecordInterval(int chrPrefix, int stPrefix, int edPrefix)
		{
			int stpre = (chrPrefix << 16) | (stPrefix & 0xFFFF);
			int edpre = (chrPrefix << 16) | (edPrefix & 0xFFFF);
			
			RegionRecord str = index.get(stpre);
			if(str == null) return null;
			
			RegionRecord combo = new RegionRecord();
			combo.startIndex = str.startIndex;
			combo.recordCount = str.recordCount;
			
			//Now, count records...
			int pre = stpre;
			while(pre <= edpre)
			{
				RegionRecord r = index.get(pre);
				if(r != null)
				{
					combo.recordCount += r.recordCount;
				}
				pre++;
			}
			
			
			return combo;
		}
		
		public void reIndex(String vidxPath) throws IOException
		{
			//Clear current map to open up memory since it is now obsolete
			index.clear();
			
			//VIDX format (BE):
			//	Record [32]
			//		Var UID[8]
			//		FileOffset [8]
			//		Record Size [4]
			//		SV Type [4]
			//		Padding [8]
			
			//We are assuming everything is sorted correctly.
			if(!FileBuffer.fileExists(vidxPath)) return;
			if(FileBuffer.fileSize(vidxPath) < 1L) return;
			
			long fsz = FileBuffer.fileSize(vidxPath);
			long cpos = 0;
			FileBuffer idx = FileBuffer.createBuffer(vidxPath, true);
			int prefix = 0;
			RegionRecord now = null;
			while(cpos < fsz)
			{
				//Only need the var IDs
				long varid = idx.longFromFile(cpos); 
				
				//Decide what to do with it
				int varpre = (int)(varid >>> 32);
				if(varpre != prefix)
				{
					if(now != null) index.put(prefix, now);
					now = new RegionRecord();
					now.startIndex = cpos >>> 5; //Divide by 32
					now.recordCount++;
					prefix = varpre;
				}
				else
				{
					//In same block as previous variant
					if(now != null) now.recordCount++;
				}
				
				cpos += 32;
			}
			
			writeIndex();
		}
		
		private void writeIndex() throws IOException
		{
			List<Integer> keys = new ArrayList<Integer>(index.size() + 1);
			keys.addAll(index.keySet());
			Collections.sort(keys);
			BufferedOutputStream bw = new BufferedOutputStream(new FileOutputStream(cidx_path));
			for(Integer k : keys)
			{
				FileBuffer rec = new FileBuffer(16, true);
				rec.addToFile(k);
				RegionRecord rr = index.get(k);
				if (rr == null) continue;
				rec.addToFile(rr.startIndex);
				rec.addToFile(rr.recordCount);
				bw.write(rec.getBytes());
			}
			bw.close();
		}
		
	}
	
	public static class GeneHitCounter
	{
		private Set<Integer> total_hits_indiv;
		private Set<Integer> exon_hits_indiv;
		private int total_hits_var;
		private int exon_hits_var;
		
		public GeneHitCounter()
		{
			total_hits_indiv = new TreeSet<Integer>();
			exon_hits_indiv = new TreeSet<Integer>();
		}
	
		public void setTotalHits(int i){total_hits_var = i;}
		public void setExonHits(int i){exon_hits_var = i;}
		public void addIndivTotal(Integer uid) {total_hits_indiv.add(uid);}
		public void addIndivExon(Integer uid) {exon_hits_indiv.add(uid);}
		
		public synchronized void addIndivTotal_sync(Integer uid) {total_hits_indiv.add(uid);}
		public synchronized void addIndivExon_sync(Integer uid) {exon_hits_indiv.add(uid);}
		
		public int countIndivTotal() {return total_hits_indiv.size();}
		public int countIndivExon() {return exon_hits_indiv.size();}
		
		public int getTotalHitsVar() {return total_hits_var;}
		public int getExonHitsVar() {return exon_hits_var;}
		public Set<Integer> getTotalHitsIndiv_setref(){return total_hits_indiv;}
		public Set<Integer> getExonHitsIndiv_setref(){return exon_hits_indiv;}
		
		public void addIndividuals_Total(Collection<Integer> sids) {total_hits_indiv.addAll(sids);}
		public void addIndividuals_Exon(Collection<Integer> sids) {exon_hits_indiv.addAll(sids);}
		public synchronized void addIndividuals_Total_sync(Collection<Integer> sids) {total_hits_indiv.addAll(sids);}
		public synchronized void addIndividuals_Exon_sync(Collection<Integer> sids) {exon_hits_indiv.addAll(sids);}
		
		public boolean removeIndividual_Total(int sid) {return this.total_hits_indiv.remove(sid);}
		public boolean removeIndividual_Exon(int sid) {return this.exon_hits_indiv.remove(sid);}
		public synchronized boolean removeIndividual_Total_sync(int sid) {return this.total_hits_indiv.remove(sid);}
		public synchronized boolean removeIndividual_Exon_sync(int sid) {return this.exon_hits_indiv.remove(sid);}
		
		public void removeIndividuals_Total(Collection<Integer> sids) {total_hits_indiv.removeAll(sids);}
		public void removeIndividuals_Exon(Collection<Integer> sids) {exon_hits_indiv.removeAll(sids);}
		public synchronized void removeIndividuals_Total_sync(Collection<Integer> sids) {total_hits_indiv.removeAll(sids);}
		public synchronized void removeIndividuals_Exon_sync(Collection<Integer> sids) {exon_hits_indiv.removeAll(sids);}
		
		public void incrementTotalHits() {total_hits_var++;}
		public void decrementTotalHits() {total_hits_var--;}
		public void incrementExonHits() {exon_hits_var++;}
		public void decrementExonHits() {exon_hits_var--;}
		
		public synchronized void incrementTotalHits_sync() {total_hits_var++;}
		public synchronized void decrementTotalHits_sync() {total_hits_var--;}
		public synchronized void incrementExonHits_sync() {exon_hits_var++;}
		public synchronized void decrementExonHits_sync() {exon_hits_var--;}
	}
	
	/*----- Inner Classes (Major) -----*/
	
	private static class ListGroup
	{
		public List<Long> deleteVars;
		public List<Long> hetVars;
		public List<Long> homVars;
		
		public ListGroup(List<Long> l_del, List<Long> l_het, List<Long> l_hom)
		{
			deleteVars = l_del;
			hetVars = l_het;
			homVars = l_hom;
		}
	}
	
	private static class GenotypeCache implements Iterable<VariantGenotype>
	{
		
		public class GenotypeTableIterator implements Iterator<VariantGenotype>
		{

			private long cpos;
			private long fsz;
			
			public GenotypeTableIterator()
			{
				cpos = 16;
				fsz = FileBuffer.fileSize(gtbl_path);
			}
			
			@Override
			public boolean hasNext() 
			{
				return (cpos < fsz);
			}

			@Override
			public VariantGenotype next() 
			{
				VariantGenotype vg = VariantGenotype.readFromGENOT(openFile, cpos);
				if (vg == null) return vg;
				cpos += vg.calculateSerializedSize(true);
				return vg;
			}
			
		}
		
		private ConcurrentMap<Long, VariantGenotype> genoCache;
		private ConcurrentLinkedQueue<Long> cacheQueue;
		
		private String gtbl_path;
		private StreamBuffer openFile;
		
		private ConcurrentMap<Long, VariantGenotype> dirtyQueue;
		//private ConcurrentLinkedQueue<Long> removeQueue;
		
		private ConcurrentMap<Long, LookupRecord> tempIndex; //Usually null
		
		public GenotypeCache(String genoTablePath)
		{
			genoCache = new ConcurrentSkipListMap<Long, VariantGenotype>();
			dirtyQueue = new ConcurrentSkipListMap<Long, VariantGenotype>();
			cacheQueue = new ConcurrentLinkedQueue<Long>();
			//removeQueue = new ConcurrentLinkedQueue<Long>();
			gtbl_path = genoTablePath;
			try 
			{
				openFile = new StreamBuffer(gtbl_path, true);
			} 
			catch (IOException e) 
			{
				//e.printStackTrace();
				System.err.println("Genotype file does not currently exist! Nothing to read!");
				openFile = null;
			}
		}
		
		protected void buildTemporaryIndex()
		{
			//WILL EAT MEMORY!
			tempIndex = new ConcurrentHashMap<Long, LookupRecord>();
			
			//Scan 
			long fsz = FileBuffer.fileSize(gtbl_path);
			long cpos = 16; //GENOT header
			while(cpos < fsz)
			{
				VariantGenotype vg = VariantGenotype.readFromGENOT(openFile, cpos);
				int sz = vg.calculateSerializedSize(true);
				LookupRecord lr = new LookupRecord();
				lr.genoOffset = cpos;
				lr.genoRecSize = sz;
				tempIndex.put(vg.getVariantUID(), lr);
				cpos += sz;
			}
		}
		
		protected void updateWithGenoIndex(long vid, LookupRecord lr)
		{
			if(tempIndex == null) buildTemporaryIndex();
			LookupRecord gr = tempIndex.get(vid);
			if (gr == null) return;
			lr.genoOffset = gr.genoOffset;
			lr.genoRecSize = gr.genoRecSize;
		}
		
		protected void freeTemporaryIndex()
		{
			tempIndex = null;
		}
		
		private void freeOldestGT()
		{
			long targ = cacheQueue.poll();
			genoCache.remove(targ);
		}
		
		public VariantGenotype getVariantGenotypes(long vid, LookupRecord lr)
		{
			//Check cache
			VariantGenotype vg = genoCache.get(vid);
			if (vg != null) return vg;
			
			//Cache miss
			if(genoCache.size() >= MAX_GENO_IN_MEM)
			{
				//Free last
				freeOldestGT();
			}
			
			//Load into cache and return
			vg = VariantGenotype.readFromGENOT(openFile, lr.genoOffset);
			if(vg != null) {
				cacheQueue.add(vg.getVariantUID());
				genoCache.put(vg.getVariantUID(), vg);
			}
			
			return vg;
		}
		
		public boolean queueForWriting(VariantGenotype newGT)
		{
			if(dirtyQueue.size() >= MAX_GENO_IN_MEM) return false;
			dirtyQueue.put(newGT.getVariantUID(), newGT);
			return true;
		}
		
		public boolean writeGenotypeTable() throws IOException
		{
			if(dirtyQueue.isEmpty()) return true; //Nothing new to write
			
			//Clear the cache! It is now obsolete. Also saves memory.
			genoCache.clear();
			cacheQueue.clear();
			tempIndex = null;
			
			//Temp file
			String temp = FileBuffer.getTempDir() + File.separator + "svdb_genot.tmp";
			
			//Generate header
			FileBuffer genotHeader = new FileBuffer(16, true);
			genotHeader.printASCIIToFile(GENOTABLE_MAGIC);
			genotHeader.addToFile(GENOTABLE_VERSION);
			genotHeader.addToFile(0);
			
			//Stream out!
			BufferedOutputStream bw = new BufferedOutputStream(new FileOutputStream(temp));
			bw.write(genotHeader.getBytes());
			
			//Write clean variants from input file...
			long cpos = 16;
			long fsz = FileBuffer.fileSize(gtbl_path);
			while(cpos < fsz)
			{
				VariantGenotype vg = VariantGenotype.readFromGENOT(openFile, cpos);
				cpos += vg.calculateSerializedSize(true);
				//See if clean...
				long vid = vg.getVariantUID();
				if(!dirtyQueue.containsKey(vid))
				{
					//Just copy back
					bw.write(vg.serializeForGENOT(true).getBytes());
				}
				//Otherwise, eat.
			}
			
			//Write dirty variant info...
			List<Long> dirtyList = new ArrayList<Long>(dirtyQueue.size() + 1);
			dirtyList.addAll(dirtyQueue.keySet());
			for(Long vid : dirtyList)
			{
				VariantGenotype vg = dirtyQueue.get(vid);
				bw.write(vg.serializeForGENOT(true).getBytes());
			}
			
			bw.close();
			dirtyQueue.clear();
			
			//Close old file
			openFile = null;
			
			//Replace old file!
			Files.move(Paths.get(temp), Paths.get(gtbl_path));
			openFile = new StreamBuffer(gtbl_path, true);
			
			return true;
		}

		public ListGroup removeSample(int sampleUID) throws IOException
		{
			if(!dirtyQueue.isEmpty()) writeGenotypeTable();
			List<Long> deleteVars = new LinkedList<Long>();
			List<Long> homVars = new LinkedList<Long>();
			List<Long> hetVars = new LinkedList<Long>();
			String temp = FileBuffer.getTempDir() + File.separator + "genot_remSample.tmp";
			
			genoCache.clear();
			cacheQueue.clear();
			tempIndex = null;
			
			FileBuffer genotHeader = new FileBuffer(16, true);
			genotHeader.printASCIIToFile(GENOTABLE_MAGIC);
			genotHeader.addToFile(GENOTABLE_VERSION);
			genotHeader.addToFile(0);
			
			//Stream out!
			BufferedOutputStream bw = new BufferedOutputStream(new FileOutputStream(temp));
			bw.write(genotHeader.getBytes());
			
			//Write clean variants from input file...
			long cpos = 16;
			long fsz = FileBuffer.fileSize(gtbl_path);
			while(cpos < fsz)
			{
				VariantGenotype vg = VariantGenotype.readFromGENOT(openFile, cpos);
				cpos += vg.calculateSerializedSize(true);
				
				//Check to see if we need to delete a genotype record
				//	or the whole variant record...
				SVDBGenotype g = vg.removeGenotype(sampleUID);
				if(vg.genotypeRecords() > 0)
				{
					//Write back out
					bw.write(vg.serializeForGENOT(true).getBytes());
				}
				else deleteVars.add(vg.getVariantUID());
				
				if(g != null)
				{
					if(g.isHomozygous()) homVars.add(vg.getVariantUID());
					else hetVars.add(vg.getVariantUID());
				}
			}
			
			bw.close();
			openFile = null;
			
			//Replace old file!
			Files.move(Paths.get(temp), Paths.get(gtbl_path));
			openFile = new StreamBuffer(gtbl_path, true);

			return new ListGroup(deleteVars, hetVars, homVars);
		}
		
		@Override
		public Iterator<VariantGenotype> iterator() 
		{
			return new GenotypeTableIterator();
		}
		
		public ListGroup getVariantIDsForSample(int sampleUID)
		{
			List<Long> homVars = new LinkedList<Long>();
			List<Long> hetVars = new LinkedList<Long>();
			
			for(VariantGenotype vg : this)
			{
				if(vg.hasGenotypeFor(sampleUID))
				{
					SVDBGenotype g = vg.getGenotype(sampleUID);
					if(g.isHomozygous()) homVars.add(vg.getVariantUID());
					else hetVars.add(vg.getVariantUID());
				}
			}
			
			return new ListGroup(null, hetVars, homVars);
		}
		
	}
	
	private static class VariantCache
	{
		//Index is sorted
		//This is so that in case the index gets too large to load
		//	fully into memory, it itself can easily be used like a tree
		//	for lookup.
		
		//VIDX format (BE):
		//	Record [32]
		//		Var UID[8]
		//		FileOffset [8]
		//		Record Size [4]
		//		SV Type [1]
		//		Genotype Data Size [3]
		//		Genotype Table Offset [8]
		
		private GenomeBuild genome;
		private GeneSet genes;
		private GenomeIndex gindex;
		
		private RegionIndex cindex;
		
		private String srcDir;
		private boolean cachedIndex;
		
		private ConcurrentMap<Long, DBVariant> cache;
		private ConcurrentMap<Long, LookupRecord> index;
		
		private ConcurrentLinkedQueue<Long> cacheQueue;
		private ConcurrentLinkedQueue<Long> indexCacheQueue;
		
		private ConcurrentMap<Long, DBVariant> dirtyQueue;
		private ConcurrentLinkedQueue<Long> removeQueue;
		
		private GenotypeCache genoCache;
		
		public VariantCache(GenomeBuild gb, GeneSet gs, String dir, boolean cacheIndex) throws IOException
		{
			genome = gb;
			genes = gs;
			gindex = new GenomeIndex(gb);
			srcDir = dir;
			cache = new ConcurrentSkipListMap<Long, DBVariant>();
			index = new ConcurrentHashMap<Long, LookupRecord>();
			cacheQueue = new ConcurrentLinkedQueue<Long>();
			cachedIndex = cacheIndex;
			indexCacheQueue = new ConcurrentLinkedQueue<Long>();
			cindex = new RegionIndex(srcDir + File.separator + CINDEX_NAME);
			genoCache = new GenotypeCache(srcDir + File.separator + GENOTBL_FILE + "." + GENOTBL_EXT);
			dirtyQueue = new ConcurrentSkipListMap<Long, DBVariant>();
			removeQueue = new ConcurrentLinkedQueue<Long>();
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
			//		SV Type [1]
			//		Genotype Data Size [3]
			//		Genotype Table Offset [8]
			
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
				int etype = Byte.toUnsignedInt(idx.getByte(cpos)); cpos++;
				int gsz = idx.shortishFromFile(cpos); cpos += 3;
				long goff = idx.longFromFile(cpos); cpos += 8;
				LookupRecord lr = new LookupRecord();
				lr.offset = off;
				lr.recordSize = sz;
				lr.type = SVType.getTypeByID(etype);
				lr.genoOffset = goff;
				lr.genoRecSize = gsz;
				index.put(vid, lr);
			}
		}
		
		private void clearOldestIndexRecordFromCache()
		{
			long oldid = indexCacheQueue.poll();
			index.remove(oldid);
		}
		
		private void clearOldestVariantFromCache()
		{
			long vid = cacheQueue.poll();
			cache.remove(vid);
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
				int etype =  Byte.toUnsignedInt(rec.getByte(20));
				int gsz = rec.shortishFromFile(21);
				long goff = rec.longFromFile(24);
				LookupRecord lr = new LookupRecord();
				lr.offset = off;
				lr.recordSize = sz;
				lr.type = SVType.getTypeByID(etype);
				lr.genoOffset = goff;
				lr.genoRecSize = gsz;
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
		
		private LookupRecord getIndexRecordFromDisc(long varUID)
		{
			//This is for disc only index lookup!
			String idxpath = getIndexPath();
			if (!FileBuffer.fileExists(idxpath)) return null;
			long fsz = FileBuffer.fileSize(idxpath);
			if(fsz < 32L) return null;
			
			//Start at middle record
			//All records are 32 bytes
			return searchForRecord(varUID, 0, fsz, idxpath);
		}
		
		public LookupRecord getIndexRecord(long varUID)
		{
			if(!cachedIndex) return index.get(varUID);
			//See if there's a hit
			LookupRecord rec = index.get(varUID);
			if (rec != null) {
				indexCacheQueue.remove(varUID);
				indexCacheQueue.add(varUID);
				return rec;
			}
			//Otherwise, miss
			if(index.size() >= MAX_IDXREC_IN_MEM)
			{
				//Clear oldest record
				clearOldestIndexRecordFromCache();
			}
			rec = getIndexRecordFromDisc(varUID);
			if (rec != null)
			{
				indexCacheQueue.add(varUID);
				index.put(varUID, rec);
				return rec;
			}
			
			return null;
		}
		
		private String getFilename(SVType type)
		{
			String typeroot = type.name();
			return srcDir + File.separator + typeroot + "." + VARTBL_EXT;
		}
		
		public DBVariant getVariant(long varUID)
		{
			//Check cache...
			DBVariant hit = cache.get(varUID);
			if(hit != null)
			{
				cacheQueue.remove(varUID);
				cacheQueue.add(varUID);
				return hit;
			}
			
			//Cache miss
			if (cache.size() >= MAX_VAR_IN_MEM) clearOldestVariantFromCache();
			
			LookupRecord lookup = getIndexRecord(varUID);
			if (lookup == null) return null;
			if(lookup.type == null) return null;
			String tblpath = getFilename(lookup.type);
			long edoff = lookup.offset + Integer.toUnsignedLong(lookup.recordSize);
			try 
			{
				FileBuffer raw = FileBuffer.createBuffer(tblpath, lookup.offset, edoff, true);
				DBVariant var = DBVariant.getFromVDBRecord(raw, genome, genes, 0).getVariant();
				cache.put(varUID, var);
				cacheQueue.add(varUID);
				return var;
			} 
			catch (IOException e) 
			{
				e.printStackTrace();
				return null;
			}

		}
		
		public boolean hasVariant(long varUID)
		{
			return (getIndexRecord(varUID) != null);
		}
		
		public DBVariant getVariantNoCache(long varUID)
		{
			DBVariant hit = cache.get(varUID);
			if(hit != null) return hit;
			
			//Cache miss
			LookupRecord lookup = getIndexRecord(varUID);
			if (lookup == null) return null;
			if(lookup.type == null) return null;
			String tblpath = getFilename(lookup.type);
			long edoff = lookup.offset + Integer.toUnsignedLong(lookup.recordSize);
			try 
			{
				FileBuffer raw = FileBuffer.createBuffer(tblpath, lookup.offset, edoff, true);
				DBVariant var = DBVariant.getFromVDBRecord(raw, genome, genes, 0).getVariant();
				return var;
			} 
			catch (IOException e) 
			{
				e.printStackTrace();
				return null;
			}

		}
		
		public List<Long> getVariantIDsInApproximateRegion(Contig c, int start, int end) throws IOException
		{
			int cpre = gindex.getContigPrefix(c);
			int stpre = gindex.getPositionDivision(c, start);
			int edpre = gindex.getPositionDivision(c, end);
			List<Long> found = new LinkedList<Long>();
			if(cachedIndex)
			{
				RegionIndex.RegionRecord r = cindex.getRecordInterval(cpre, stpre, edpre);
				if(r == null) return found;
				
				//Scan the vidx, loading all lookup records in range into index cache
				long vidx_off = r.startIndex << 5;
				long vidx_end_off = vidx_off + ((long)r.recordCount * 32L);
				
				String vidx_path = getIndexPath();
				if(!FileBuffer.fileExists(vidx_path)) return found;
				FileBuffer vidx = FileBuffer.createBuffer(vidx_path, vidx_off, vidx_end_off, true);
				long fsz = vidx.getFileSize();
				long cpos = 0;
				while(cpos < fsz)
				{
					long vid = vidx.longFromFile(cpos); cpos += 8;
					//See if already cached...
					LookupRecord lr = index.get(vid);
					if(lr != null)
					{
						indexCacheQueue.remove(vid);
						indexCacheQueue.add(vid);
						cpos += 24;
					}
					else
					{
						lr = new LookupRecord();
						lr.offset = vidx.longFromFile(cpos); cpos += 8;
						lr.recordSize = vidx.intFromFile(cpos); cpos += 4;
						int etype = vidx.intFromFile(cpos); cpos+=4;
						cpos += 8;
						lr.type = SVType.getTypeByID(etype);
						//Cache
						if(index.size() >= MAX_IDXREC_IN_MEM) clearOldestIndexRecordFromCache();
						index.put(vid, lr);
						indexCacheQueue.add(vid);
					}
					found.add(vid);
				}
				
			}
			else
			{
				//Already in memory - no need to check the cidx
				int stp = (cpre << 16) | (stpre & 0xFFFF);
				int edp = (cpre << 16) | (edpre & 0xFFFF);
				
				long min = Integer.toUnsignedLong(stp) << 32;
				long max = (Integer.toUnsignedLong(edp) << 32) | 0xFFFFFFFF;
				
				List<Long> allids = new ArrayList<Long>(index.size() + 1);
				allids.addAll(index.keySet());
				Collections.sort(allids);
				for(Long l : allids)
				{
					if (l > max) break;
					if (l >= min) found.add(l);
				}
			}
			
			return found;
		}
		
		private List<DBVariant> scanForTRAs(Contig c, int start, int end) throws IOException, UnsupportedFileTypeException
		{
			String tratblpath = getFilename(SVType.TRA);
			List<DBVariant> vlist = new LinkedList<DBVariant>();
			if(!FileBuffer.fileExists(tratblpath)) return vlist;
			VDBFileScanner scanner = new VDBFileScanner(genome, genes, tratblpath);
			scanner.open();
			for(DBVariant v : scanner)
			{
				//Check the variant END to see if in range
				//Also check the end of the start 
				if(v.getChrom().equals(c))
				{
					int st2 = v.getStartPosition().getEnd();
					if(st2 >= start && st2 < end)
					{
						vlist.add(v);
						continue;
					}
				}
				Contig c2 = v.getEndChrom();
				if(c2 != null)
				{
					if(c2.equals(c))
					{
						int ed1 = v.getEndPosition().getStart();
						int ed2 = v.getEndPosition().getEnd();
						if(ed2 < start) continue;
						if(ed1 >= end) continue;
						vlist.add(v);
					}
				}
			}
			scanner.close();
			return vlist;
		}
		
		public Collection<DBVariant> getVariantsInRegion(Contig c, int start, int end) throws IOException
		{
			ConcurrentLinkedQueue<DBVariant> varlist = new ConcurrentLinkedQueue<DBVariant>();
			//Spawn sub threads for scanning all TRAs and checking earlier contigs
			Runnable trascanner = new Runnable() {

				@Override
				public void run() 
				{
					List<DBVariant> vlist;
					try 
					{
						vlist = scanForTRAs(c, start, end);
						varlist.addAll(vlist);
					} 
					catch (IOException | UnsupportedFileTypeException e) 
					{
						e.printStackTrace();
					}
			
				}
				
			};
			Thread traThread = new Thread(trascanner);
			traThread.start();
			
			Runnable endscanner = new Runnable() {
				@Override
				public void run()
				{
					//Figure out max pos
					int secsz = (int)(c.getLength()/0xFFFFL);
					int sec = start/secsz;
					int secst = sec * secsz;
					int max = secst - 1;
					try 
					{
						List<Long> idlist = getVariantIDsInApproximateRegion(c, 0, max);
						for(Long vid : idlist)
						{
							DBVariant var = getVariantNoCache(vid);
							//See if variant is actually in requested region
							if(var.getEndPosition().getEnd() < start) continue;
							//if(var.getStartPosition().getStart() > end) continue;
							varlist.add(var);
						}
					} 
					catch (IOException e) 
					{
						e.printStackTrace();
					}
				}
			};
			Thread endThread = new Thread(endscanner);
			endThread.start();
			
			List<Long> checkids = getVariantIDsInApproximateRegion(c, start, end);
			for(Long vid : checkids)
			{
				DBVariant var = getVariant(vid);
				//See if variant is actually in requested region
				if(var.getEndPosition().getEnd() < start) continue;
				if(var.getStartPosition().getStart() > end) continue;
				varlist.add(var);
			}
			
			//Wait for the other threads to finish
			while(traThread.isAlive() || endThread.isAlive())
			{
				try 
				{
					Thread.sleep(10);
				} 
				catch (InterruptedException e) 
				{
					// Eat, but print stack trace as a warning
					e.printStackTrace();
				}
			}
			
			return varlist;
		}
		
		public void generateMasterIndex() throws IOException
		{
			//Write the vidx
			
			//Wipe the index cache
			indexCacheQueue.clear();
			index.clear();
			
			//If we are set to cache the index, we may need some disc space...
			int tmp_idx = 0;
			String tmp_path_stem = FileBuffer.getTempDir() + File.separator + "svdb_varindextmp_";
			
			SVType[] types = SVType.values();
			for(SVType t : types)
			{
				String vdbpath = getFilename(t);
				if(!FileBuffer.fileExists(vdbpath)) continue;
				long fsz = FileBuffer.fileSize(vdbpath);
				FileBuffer file = new StreamBuffer(vdbpath, 8, fsz, true);
				
				//Scan the file
				long cpos = 0;
				while(cpos < fsz)
				{
					ParsedVariant pv = DBVariant.getFromVDBRecord(file, genome, genes, cpos);
					LookupRecord lr = new LookupRecord();
					lr.offset = cpos;
					lr.recordSize = (int)pv.getSize();
					lr.type = t;
					cpos += pv.getSize();
					DBVariant v = pv.getVariant();
					if (v != null)
					{
						long vid = v.getLongID();
						genoCache.updateWithGenoIndex(vid, lr);
						if(!cachedIndex) index.put(vid, lr);
						else
						{
							//We have to see if there is space in the cache
							if(index.size() >= MAX_IDXREC_IN_MEM)
							{
								//Dump to tmp file
								String tmppath = tmp_path_stem + tmp_idx + ".tmp";
								dumpIndexFromMemory(tmppath);
								index.clear();
								tmp_idx++;
							}
							index.put(vid, lr);
						}
					}
				}
				
			}
			
			genoCache.freeTemporaryIndex();
			
			//Scan done
			if(!cachedIndex) dumpIndexFromMemory(getIndexPath());
			else
			{
				if(tmp_idx < 1) dumpIndexFromMemory(getIndexPath()); //Never filled, don't have to combine temps
				else combineTempIndexDumps(tmp_path_stem, tmp_idx);
			}
		}
		
		private void combineTempIndexDumps(String stem, int filecount) throws IOException
		{
			//Also make sure to dump what remains in memory!
			
			List<Long> memids = new ArrayList<Long>(index.size() + 1);
			memids.addAll(index.keySet());
			Collections.sort(memids);
			int memind = 0;
			
			FileBuffer[] temps = new FileBuffer[filecount];
			long[] cposs = new long[filecount];
			long[] fszs = new long[filecount];
			for(int i = 0; i < filecount; i++)
			{
				String tmppath = stem + i + ".tmp";
				fszs[i] = FileBuffer.fileSize(tmppath);
				temps[i] = new StreamBuffer(tmppath, StreamBuffer.DEFO_SUBBUF_SIZE, 16);
				temps[i].setEndian(true);
				cposs[i] = 0;
			}
			
			String outpath = this.getIndexPath();
			BufferedOutputStream bw = new BufferedOutputStream(new FileOutputStream(outpath));
			
			//Scan all temps and the mem set to find the lowest var UID
			//Grab that and write record
			boolean done = false;
			long lowest = Long.MAX_VALUE;
			int lind = -1;
			while(!done)
			{
				lowest = Long.MAX_VALUE;
				lind = -1;
				//Scan the temps
				for(int i = 0; i < filecount; i++)
				{
					if (cposs[i] >= fszs[i]) continue;
					long vid = temps[i].longFromFile(cposs[i]);
					if (vid <= lowest)
					{
						lowest = vid;
						lind = i;
					}
				}
				//Compare to what's in memory
				if(memind < memids.size() && memids.get(memind) < lowest)
				{
					//Serialize and write the record in memory
					FileBuffer srec = new FileBuffer(32, true);
					LookupRecord lr = index.get(memids.get(memind));
					if (lr == null) continue;
					srec.addToFile(memids.get(memind));
					srec.addToFile(lr.offset);
					srec.addToFile(lr.recordSize);
					srec.addToFile((byte)lr.type.getID());
					srec.add24ToFile(lr.genoRecSize);
					srec.addToFile(lr.genoOffset);
					
					bw.write(srec.getBytes());
					memind++;
				}
				else
				{
					//Copy the data from the specified file and advance cpos
					//Just dump the next 32 bytes to bw
					bw.write(temps[lind].getBytes(cposs[lind], cposs[lind]+32));
					cposs[lind] += 32;
				}
				//Check for done condition
				done = true;
				if(memind < memids.size()) done = false;
				if(done)
				{
					for(int i = 0; i < filecount; i++)
					{
						if (cposs[i] < fszs[i])
						{
							done = false;
							break;
						}
					}
				}
			}
			
			bw.close();
		}
		
		private void dumpIndexFromMemory(String path) throws IOException
		{
			//Sort variant UIDs
			List<Long> allids = new ArrayList<Long>(index.size() + 1);
			allids.addAll(index.keySet());
			Collections.sort(allids);
			
			BufferedOutputStream bw = new BufferedOutputStream(new FileOutputStream(path));
			
			for(Long vid : allids)
			{
				FileBuffer srec = new FileBuffer(32, true);
				LookupRecord lr = index.get(vid);
				if (lr == null) continue;
				srec.addToFile(vid);
				srec.addToFile(lr.offset);
				srec.addToFile(lr.recordSize);
				srec.addToFile((byte)lr.type.getID());
				srec.add24ToFile(lr.genoRecSize);
				srec.addToFile(lr.genoOffset);
				
				bw.write(srec.getBytes());
			}
			
			bw.close();
		}
		
		public void generateRegionIndex() throws IOException
		{
			/*//Write the cidx
			//Record numbers in vidx for start end of each region
			String vidx_path = this.getIndexPath();
			String cidx_path = cindex.cidx_path;
			if(!FileBuffer.fileExists(vidx_path)) return;
			BufferedInputStream br = new BufferedInputStream(new FileInputStream(vidx_path));
			BufferedOutputStream bw = new BufferedOutputStream(new FileOutputStream(cidx_path));
			long read = 0;
			long fsz = FileBuffer.fileSize(vidx_path);
			long rindex = 0;
			long lastindex = 0;
			int rcount = 0;
			int lastPrefix = -1;
			
			while(read < fsz)
			{
				byte[] vrec = new byte[32];
				br.read(vrec);
				read+=32;
				//Read the first four bytes as an int
				int prefix = 0;
				prefix |= (Byte.toUnsignedInt(vrec[0]) << 24);
				prefix |= (Byte.toUnsignedInt(vrec[1]) << 16);
				prefix |= (Byte.toUnsignedInt(vrec[2]) << 8);
				prefix |= Byte.toUnsignedInt(vrec[3]);
				
				if(prefix != lastPrefix)
				{
					//Write record to cidx
					if(lastPrefix != -1)
					{
						FileBuffer crec = new FileBuffer(16, true);
						crec.addToFile(lastPrefix);
						crec.addToFile(lastindex);
						crec.addToFile(rcount);
						bw.write(crec.getBytes());	
					}
					
					lastindex = rindex;
					rcount = 0;
					lastPrefix = prefix;
				}
				rcount++;
			}
			
			br.close();
			bw.close();*/
			cindex.reIndex(getIndexPath());
		}

		public int calculateLeewayBP(StructuralVariant sv, double leewayPerc)
		{
			//double ldbl = (double)leewayRaw/1000.0;
			int vlen = 0;
			
			if(sv.getType() == SVType.TRA || sv.getType() == SVType.BND)
			{
				//Use the CIPOS90
				int st = sv.getCIPosition(false, false, false);
				int ed = sv.getCIPosition(false, false, true);
				vlen = ed-st;
			}
			else
			{
				vlen = sv.getEndPosition() - sv.getPosition();
			}
			
			return (int)Math.round((double)vlen * leewayPerc);
		}
		
		public long generateUID(StructuralVariant sv)
		{
			Contig c1 = sv.getChromosome();
			int pos = sv.getCIPosition(false, false, false);
			long id = gindex.getVarIDPrefix(c1, pos);
			
			Random r = new Random();
			long suffix = r.nextLong() & 0xFFFFFFFF;
			id |= suffix;
			while(dirtyQueue.containsKey(id) || hasVariant(id))
			{
				//Make sure it's unique
				id &= ~0xFFFFFFFF;
				suffix = r.nextLong() & 0xFFFFFFFF;
				id |= suffix;
			}
			
			return id;
		}
		
		public long addVariant(StructuralVariant sv, int mergeFactor, Map<String, FamilyMember> sampleMap)
		{
			if(sv == null) return -1;
			//Get some basic info on the variant
			double ldbl = (double)mergeFactor/1000.0;
			int lbp = calculateLeewayBP(sv, ldbl);

			//See if it merges to existing variant
			int stEnd = sv.getCIPosition(false, false, true);
			Collection<Long> cand;
			try {cand = this.getVariantIDsInApproximateRegion(sv.getChromosome(), 0, stEnd + lbp);} 
			catch (IOException e) 
			{
				e.printStackTrace();
				return -1;
			}
			
			if(cand != null && !cand.isEmpty())
			{
				for(Long vid : cand)
				{
					//Grab lookup record
					LookupRecord lr = getIndexRecord(vid);
					//Toss if not the same type
					if(lr.type != sv.getType()) continue;
					//Read variant if same type
					DBVariant dbv = getVariant(vid);
					if(dbv.svIsEquivalent(sv, ldbl))
					{
						if(mergeInVariant(sv, dbv, sampleMap)) return vid;
						else return -1;	
					}
				}
			}
			

			//If not, add anew
			//If it gets here, did not find a merge partner
			//Generate new UID
			long nuid = addNewVariant(sv, sampleMap);
			return nuid;
		}
		
		private long addNewVariant(StructuralVariant sv, Map<String, FamilyMember> sampleMap)
		{
			long nuid = generateUID(sv);
			String vname = sv.getType().getString() + "_" + Long.toHexString(nuid);
			DBVariant dbv = DBVariant.getFromVariant(sv, vname);
			dbv.setLongUID(nuid);
			
			//Try adding genotypes...
			VariantGenotype geno = null;
			if(sampleMap != null && !sampleMap.isEmpty())
			{
				Collection<FamilyMember> samples = sampleMap.values();
				geno = new VariantGenotype(nuid);
				for(FamilyMember m : samples)
				{
					SVDBGenotype g = SVDBGenotype.generateGenotype(m, sv);
					if(g != null)
					{
						Collection<Population> ptags = m.getPopulationTags();
						//Update counts in variant
						dbv.incrementTotalCount();
						for(Population p : ptags) dbv.incrementTotalCount(p);
						//Check if hom
						if(g.isHomozygous())
						{
							dbv.incrementHomozygoteCount();
							for(Population p : ptags) dbv.incrementHomozygoteCount(p);
						}
						//Pull genotype and update as well
						geno.addGenotype(g);
					}
				}
			}
			
			if(!queueVariantForUpdate(dbv))
			{
				//Dump and restart
				try 
				{
					updateTable();
				} 
				catch (IOException | UnsupportedFileTypeException e) 
				{
					e.printStackTrace();
					return -1;
				}
				queueVariantForUpdate(dbv);
			}
			
			if(geno != null)
			{
				if(!queueGenotypeForUpdate(geno))
				{
					try 
					{
						updateTablesAndReindex();
					} 
					catch (IOException | UnsupportedFileTypeException e) 
					{
						e.printStackTrace();
						return -1;
					}
					if(queueGenotypeForUpdate(geno)) return nuid;
					return -1;
				}
			}
			
			return nuid;

		}
		
		private boolean mergeInVariant(StructuralVariant src, DBVariant target, Map<String, FamilyMember> sampleMap)
		{
			if(src == null) return false;
			
			//Bring in target...
			//LookupRecord lr = getIndexRecord(targetUID); //We may reuse this
			//DBVariant target = getVariant(targetUID);
			if(target == null) return false;

			int nPos = src.getPosition();
			int nEnd = src.getEndPosition();
			boolean change = false;
			if(!target.getStartPosition().contains(nPos))
			{
				change = true;
				if(nPos < target.getStartPosition().getStart()) target.getStartPosition().setStart(nPos);
				else if (nPos > target.getStartPosition().getEnd()) target.getStartPosition().setEnd(nPos);
			}
			
			if(!target.getEndPosition().contains(nEnd))
			{
				change = true;
				if(nEnd < target.getEndPosition().getStart()) target.getEndPosition().setStart(nEnd);
				else if (nEnd > target.getEndPosition().getEnd()) target.getEndPosition().setEnd(nEnd);
			}
		
			//Update allele counts (or check if need to)
			boolean genoChange = false;
			VariantGenotype geno = null;
			if(sampleMap != null && !sampleMap.isEmpty())
			{
				Collection<FamilyMember> samples = sampleMap.values();
				long vid = target.getLongID();
				geno = genoCache.getVariantGenotypes(vid, getIndexRecord(vid));
				if (geno == null) geno = new VariantGenotype(vid);
				for(FamilyMember m : samples)
				{
					SVDBGenotype g = SVDBGenotype.generateGenotype(m, src);
					if(g != null)
					{
						change = true;
						Collection<Population> ptags = m.getPopulationTags();
						//Update counts in variant
						target.incrementTotalCount();
						for(Population p : ptags) target.incrementTotalCount(p);
						//Check if hom
						if(g.isHomozygous())
						{
							target.incrementHomozygoteCount();
							for(Population p : ptags) target.incrementHomozygoteCount(p);
						}
						//Pull genotype and update as well
						geno.addGenotype(g);
						genoChange = true;
					}
				}
			}
			
			if(!change) return true; //Don't need to update record
			if(!queueVariantForUpdate(target))
			{
				//Dump and restart
				try 
				{
					updateTablesAndReindex();
				} 
				catch (IOException | UnsupportedFileTypeException e) 
				{
					e.printStackTrace();
					return false;
				}
				return queueVariantForUpdate(target);
			}
			
			if(genoChange && geno != null)
			{
				if(!genoCache.queueForWriting(geno))
				{
					try 
					{
						updateTablesAndReindex();
					} 
					catch (IOException | UnsupportedFileTypeException e) 
					{
						e.printStackTrace();
						return false;
					}
					return genoCache.queueForWriting(geno);
				}
			}
			
			return true;
		}
				
		private boolean addVariantOnDisk(DBVariant dbv)
		{
			String tpath = getFilename(dbv.getType());
			FileBuffer rec = dbv.toVDBRecord();
			try
			{
				rec.appendToFile(tpath);
			} 
			catch (NoSuchFileException e) 
			{
				try 
				{
					rec.writeFile(tpath);
				} 
				catch (IOException e1) 
				{
					e1.printStackTrace();
					return false;
				}
			} 
			catch (IOException e) 
			{
				e.printStackTrace();
				return false;
			}
			return true;
		}
		
		public boolean queueVariantForUpdate(DBVariant dbv)
		{
			if(dbv == null) return false;
			if(dirtyQueue.size() >= MAX_VAR_IN_MEM) return false;
			dirtyQueue.put(dbv.getLongID(), dbv);
			return true;
		}
		
		public void updateTable() throws IOException, UnsupportedFileTypeException
		{
			System.err.println("UpdateTable called!");
			if(dirtyQueue.isEmpty() && removeQueue.isEmpty()) return;
			
			//168.1 MB
			
			cache.clear();
			index.clear();
			cacheQueue.clear();
			indexCacheQueue.clear();
			
			//167.3 MB
			
			FileBuffer header = new FileBuffer(8, true);
			header.printASCIIToFile(VARTBL_MAGIC);
			header.addToFile(VARTBL_VERSION);
			
			SVType[] types = SVType.values();
			for(SVType t : types)
			{
				String tpath = getFilename(t);
				if(FileBuffer.fileExists(tpath))
				{
					String temp = FileBuffer.getTempDir() + File.separator + "vtableupdate.tmp";
					BufferedOutputStream bw = new BufferedOutputStream(new FileOutputStream(temp));
					bw.write(header.getBytes());
					//168.2
					VDBFileScanner streamer = new VDBFileScanner(genome, genes, tpath);
					streamer.open();
					//166.7
					for(DBVariant dbv : streamer)
					{
						System.err.println("DBVariantTable.VariantCache.updateTable || -DEBUG- Variant read from old file!");
						//See if dirty
						if(!dirtyQueue.containsKey(dbv.getLongID()) && !removeQueue.contains(dbv.getLongID()))
						{
							System.err.println("DBVariantTable.VariantCache.updateTable || -DEBUG- Variant good to copy!");
							//Copy!
							FileBuffer rec = dbv.toVDBRecord();
							bw.write(rec.getBytes(0, rec.getFileSize()));
							//System.err.println("Hitting while breaker");
							//while(breaker);
						}
					}
					bw.close();
					
					//Copy back
					//Files.deleteIfExists(Paths.get(tpath));
					Files.move(Paths.get(tpath), Paths.get(tpath + ".old"));
					Files.move(Paths.get(temp), Paths.get(tpath)); //Looks like it won't move if target exists >.>
				}
			}
			
			//Write the dirty records
			List<Long> dirtyids = new ArrayList<Long>(dirtyQueue.size() + 1);
			dirtyids.addAll(dirtyQueue.keySet());
			for(Long id : dirtyids)
			{
				if(removeQueue.contains(id)) continue;
				if(!this.addVariantOnDisk(dirtyQueue.get(id))) throw new IOException();
			}
			
			dirtyQueue.clear();
			removeQueue.clear();

		}
	
		public VariantGenotype getVariantGenotypes(long varUID)
		{
			//Lookup
			LookupRecord lr = getIndexRecord(varUID);
			if(lr == null) return null;
			return genoCache.getVariantGenotypes(varUID, lr);
		}
		
		public boolean queueGenotypeForUpdate(VariantGenotype vg)
		{
			return genoCache.queueForWriting(vg);
		}
		
		public void updateTables() throws IOException, UnsupportedFileTypeException
		{
			updateTable();
			genoCache.writeGenotypeTable();
		}
		
		public void updateTablesAndReindex() throws IOException, UnsupportedFileTypeException
		{
			updateTables();
			generateMasterIndex();
			generateRegionIndex();
		}
		
		public Iterator<VariantGenotype> getGenotypeIterator()
		{
			return genoCache.iterator();
		}

		public boolean removeSample(FamilyMember mem)
		{
			int sampleUID = mem.getUID();
			ListGroup vlist = null;
			try {vlist = genoCache.removeSample(sampleUID);}
			catch (IOException e) 
			{
				e.printStackTrace();
				return false;
			}
			
			//First fill the remove queue...
			for(Long r : vlist.deleteVars) removeQueue.add(r);
			
			//Write and index
			try {updateTablesAndReindex();} 
			catch (IOException | UnsupportedFileTypeException e) {
				e.printStackTrace();
				return false;
			}
			
			//Now, deal with the allele count updating...
			Collection<Population> ptags = mem.getPopulationTags();
			for(Long v : vlist.hetVars)
			{
				DBVariant target = getVariantNoCache(v);
				target.decrementTotalCount();
				for(Population p : ptags) target.decrementTotalCount(p);
				dirtyQueue.put(v, target);
			}
			
			for(Long v : vlist.homVars)
			{
				DBVariant target = getVariantNoCache(v);
				target.decrementTotalCount();
				target.decrementHomozygoteCount();
				for(Population p : ptags) {
					target.decrementTotalCount(p);
					target.decrementHomozygoteCount(p);
				}
				dirtyQueue.put(v, target);
			}
			
			try {updateTablesAndReindex();} 
			catch (IOException | UnsupportedFileTypeException e) {
				e.printStackTrace();
				return false;
			}
			
			return true;
		}
		
		public boolean updatePopFlags(int sampleUID, Collection<Population> added, Collection<Population> removed)
		{
			ListGroup vlist = genoCache.getVariantIDsForSample(sampleUID);
			
			for(Long vid : vlist.hetVars)
			{
				DBVariant dbv = this.getVariantNoCache(vid);
				if(dbv == null) return false;
				for(Population p : added) dbv.incrementTotalCount(p);
				for(Population p : removed) dbv.decrementTotalCount(p);
			}
			
			try {updateTablesAndReindex();} 
			catch (IOException | UnsupportedFileTypeException e) {
				e.printStackTrace();
				return false;
			}
			
			for(Long vid : vlist.homVars)
			{
				DBVariant dbv = this.getVariantNoCache(vid);
				if(dbv == null) return false;
				for(Population p : added) {
					dbv.incrementTotalCount(p);
					dbv.incrementHomozygoteCount(p);
				}
				for(Population p : removed) {
					dbv.decrementTotalCount(p);
					dbv.decrementHomozygoteCount(p);
				}
			}
			
			try {updateTablesAndReindex();} 
			catch (IOException | UnsupportedFileTypeException e) {
				e.printStackTrace();
				return false;
			}
			
			return true;
		}
		
		public SVType getVariantType(long varUID)
		{
			LookupRecord lr = this.getIndexRecord(varUID);
			if(lr == null) return null;
			return lr.type;
		}
		
		public Map<String, GeneHitCounter> generateGeneHitMap()
		{
			Map<String, GeneHitCounter> map = new HashMap<String, GeneHitCounter>();
			for(VariantGenotype vg : genoCache)
			{
				//See what genes the variant overlaps with
				long vid = vg.getVariantUID();
				DBVariant var = getVariantNoCache(vid);
				List<Gene> glist = var.getGeneListReference();
				if(glist == null || glist.isEmpty()) continue; //Nothing to count!
				
				//Iterate through genes and update counts.
				//Remember to check if exonic
				for(Gene g : glist)
				{
					//See if there is already a record
					GeneHitCounter ghc = map.get(g.getID());
					if(ghc == null)
					{
						ghc = new GeneHitCounter();
						map.put(g.getID(), ghc);
					}
					//Total var always gets incremented
					ghc.total_hits_var++;
					//See if exonic
					boolean isExonic = g.getRelativeRegionLocationEffect(var.getStartPosition().getStart(), var.getEndPosition().getEnd()) == GeneFunc.EXONIC;
					if(isExonic) ghc.exon_hits_var++;
					Collection<SVDBGenotype> gtlist = vg.getGenotypes();
					for(SVDBGenotype gt : gtlist)
					{
						ghc.total_hits_indiv.add(gt.getIndividualUID());
						if(isExonic) ghc.exon_hits_indiv.add(gt.getIndividualUID());
					}
				}
			}
			return map;
		}
		
	}

	/*----- Instance Variables -----*/
	
	private VariantCache varCache;
	
	/*----- Construction -----*/
	
	public DBVariantTable(GenomeBuild gb, GeneSet gs, String srcDir, boolean limitVarIndexRecords) throws IOException
	{
		String vardirpath = srcDir + File.separator + VAR_DIR;
		if(!FileBuffer.directoryExists(vardirpath)) Files.createDirectories(Paths.get(vardirpath));
		varCache = new VariantCache(gb, gs, vardirpath, limitVarIndexRecords);
	}
	
	/*----- Read -----*/
	
	public DBVariant getVariant(long varUID)
	{
		return varCache.getVariant(varUID);
	}
	
	public VariantGenotype getGenotype(long varUID)
	{
		return varCache.getVariantGenotypes(varUID);
	}
	
	public List<Long> getVariantIDsForSample(int sampleUID)
	{
		List<Long> idlist = new LinkedList<Long>();
		Iterator<VariantGenotype> genoIterator = varCache.getGenotypeIterator();
		while(genoIterator.hasNext())
		{
			VariantGenotype vg = genoIterator.next();
			//See if target sample is in list.
			if (vg.hasGenotypeFor(sampleUID)) idlist.add(vg.getVariantUID());
		}
		
		return idlist;
	}
	
	public List<Long> getVariantIDsForFamily(Family fam)
	{
		List<Long> idlist = new LinkedList<Long>();
		if(fam == null) return idlist;
		List<FamilyMember> members = fam.getAllFamilyMembers();
		Iterator<VariantGenotype> genoIterator = varCache.getGenotypeIterator();
		while(genoIterator.hasNext())
		{
			VariantGenotype vg = genoIterator.next();
			//See if target sample is in list.
			for(FamilyMember m : members)
			{
				if (vg.hasGenotypeFor(m.getUID())) idlist.add(vg.getVariantUID());
				break;
			}
		}
		
		return idlist;
	}
	
	public List<Long> getVariantIDsForSampleOfType(int sampleUID, SVType type)
	{
		List<Long> list = new LinkedList<Long>();
		List<Long> allvar = getVariantIDsForSample(sampleUID);
		if(allvar == null) return list;
		for(Long vid : allvar)
		{
			if(varCache.getVariantType(vid) == type) list.add(vid);
		}
		return list;
	}
	
	public Collection<DBVariant> getVariantsInRegion(Contig c, int start, int end)
	{
		try 
		{
			return varCache.getVariantsInRegion(c, start, end);
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
			return null;
		}
	}
	
	public Map<String, GeneHitCounter> generateGeneHitMap()
	{
		return varCache.generateGeneHitMap();
	}
	
	/*----- Write -----*/
	
	public boolean addVCF(String vcfpath, Map<String, FamilyMember> sampleMap, int mergeFactor, boolean ignoreTRA, int threads) throws IOException
	{
		VCFReadStreamer vcfStream = new VCFReadStreamer(vcfpath, varCache.genome);
		
		vcfStream.open();
		Iterator<StructuralVariant> iterator = vcfStream.getSVIterator();
		
		StructuralVariant sv = null;
		while((sv = iterator.next()) != null)
		{
			long added = varCache.addVariant(sv, mergeFactor, sampleMap);
			if(added == -1)
			{
				System.err.println("DBVariantTable.addVCF || There was an error adding a variant! Returning...");
				vcfStream.close();
				return false;
			}
		}
		
		vcfStream.close();
		return true;
	}
	
	public boolean removeSample(FamilyMember sample)
	{
		if(sample == null) return false;
		return varCache.removeSample(sample);
	}
	
	public boolean removeFamily(Family fam)
	{
		if(fam == null) return false;
		List<FamilyMember> members = fam.getAllFamilyMembers();
		boolean b = true;
		for(FamilyMember m : members)
		{
			b = b && removeSample(m);
		}
		return b;
	}
	
	public boolean updateSampleCounts(DBSampleTable sampleTable)
	{
		//TODO
		//Updates the counts of samples with a variant and homozygotes
		//	stored in the data for each variant.
		//Sample table is needed to check population tags for each sample
		//Population frequencies are not stored in these variant records
		return false;
	}
	
	public boolean updateSamplePopulationFlags(FamilyMember sample, Collection<Population> oldFlags)
	{
		Collection<Population> nflags = sample.getPopulationTags();
		
		List<Population> addedFlags = new LinkedList<Population>();
		List<Population> removedFlags = new LinkedList<Population>();
		
		if (oldFlags != null)
		{
			for(Population p : oldFlags)
			{
				if(!nflags.contains(p)) removedFlags.add(p);
			}	
			for(Population p : nflags)
			{
				if(!oldFlags.contains(p)) addedFlags.add(p);
			}
		}
		else addedFlags.addAll(nflags);
		
		return varCache.updatePopFlags(sample.getUID(), addedFlags, removedFlags);
	}
	
	public void save() throws IOException
	{
		try 
		{
			varCache.updateTablesAndReindex();
		} 
		catch (UnsupportedFileTypeException e) 
		{
			e.printStackTrace();
		}
	}
	
	public void clearVariantTable()
	{
		//TODO
	}
	
	public void dumpTable(String directory)
	{
		//TODO
	}
	
	public void close()
	{
		//TODO
	}
	
	public void indexByRegion() throws IOException
	{
		//TODO
	}
	
	public boolean updateSampleGenotypeTable()
	{
		return false;
	}
	
}
