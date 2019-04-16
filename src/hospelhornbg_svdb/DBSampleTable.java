package hospelhornbg_svdb;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.ConcurrentSkipListMap;

import hospelhornbg_segregation.Family;
import hospelhornbg_segregation.FamilyMember;
import hospelhornbg_segregation.Population;
import waffleoRai_Utils.CompositeBuffer;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class DBSampleTable implements Iterable<Family>{
	
	public static final int MAX_FAM_IN_MEM = 256;
	public static final String FINDEX_NAME = "famid.fidx";
	public static final String SINDEX_NAME = "smplid.sidx";
	
	public static final String FAMI_DIR = "fam";
	
	public static final int FAMILY_NAME_MAX_LEN = 28;
	public static final int SAMPLE_NAME_MAX_LEN = 24;
	
	private static class FamilyCache implements Iterable<Family>
	{	
		private String srcDir;
		
		private ConcurrentLinkedQueue<Integer> usedQueue;
		
		//private ConcurrentMap<String, Family> nameMap;
		private ConcurrentMap<Integer, Family> idMap;
		
		private ConcurrentMap<String, Integer> nameIndex;
		
		public FamilyCache(String dirpath) throws IOException
		{
			srcDir = dirpath;
			usedQueue = new ConcurrentLinkedQueue<Integer>();
			//nameMap = new ConcurrentHashMap<String, Family>();
			idMap = new ConcurrentSkipListMap<Integer, Family>();
			nameIndex = new ConcurrentHashMap<String, Integer>();
			readIndex();
		}

		private void readIndex() throws IOException
		{
			String idxPath = srcDir + File.separator + FINDEX_NAME;
			//FamName [28]
			//Fam UID [4]
			
			FileBuffer idxFile = FileBuffer.createBuffer(idxPath, true);
			long cpos = 0;
			long fsz = idxFile.getFileSize();
			while(cpos < fsz)
			{
				String fname = idxFile.getASCII_string(cpos, FAMILY_NAME_MAX_LEN); cpos += FAMILY_NAME_MAX_LEN;
				int uid = idxFile.intFromFile(cpos); cpos += 4;
				nameIndex.put(fname, uid);
			}
		}
		
		private Family loadFamily(int uid)
		{
			String fpath = srcDir + File.separator + Integer.toHexString(uid) + ".fam";
			if (!FileBuffer.fileExists(fpath)) return null;
			if (idMap.size() >= MAX_FAM_IN_MEM)
			{
				int freeid = usedQueue.poll();
				idMap.remove(freeid);
				//if (freefam != null) nameMap.remove(freefam.getFamilyName());
			}
			
			try 
			{
				Family f = Family.readFromFAMI(fpath);
				//nameMap.put(f.getFamilyName(), f);
				idMap.put(uid, f);
				usedQueue.add(uid);
				return f;
			} 
			catch (IOException e) 
			{
				e.printStackTrace();
				return null;
			} 
			catch (UnsupportedFileTypeException e) 
			{
				e.printStackTrace();
				return null;
			}
		}
		
		public Family getFamily(int uid)
		{
			Family f = idMap.get(uid);
			if (f != null)
			{
				usedQueue.remove(uid);
				usedQueue.add(uid);
				return f;
			}
			//Try to load into memory
			return loadFamily(uid);
		}
		
		public Family getFamily(String famName)
		{
			//Family f = nameMap.get(famName);
			int uid = nameIndex.get(famName);
			return getFamily(uid);
		}
			
		public void writeIndex() throws IOException
		{
			String idxPath = srcDir + File.separator + FINDEX_NAME;
			
			int nfam = nameIndex.size();
			FileBuffer out = new CompositeBuffer(nfam);
			
			List<String> fnames = new ArrayList<String>(nfam+1);
			fnames.addAll(nameIndex.keySet());
			Collections.sort(fnames);
			
			for(String s : fnames)
			{
				FileBuffer rec = new FileBuffer(32, true);
				String fname = s;
				int nlen = s.length();
				if (nlen > FAMILY_NAME_MAX_LEN) {
					fname = s.substring(0, FAMILY_NAME_MAX_LEN);
					nlen = FAMILY_NAME_MAX_LEN;
				}
				rec.printASCIIToFile(fname);
				while(nlen < FAMILY_NAME_MAX_LEN)
				{
					rec.addToFile((byte)0x00);
					nlen++;
				}
				int fid = nameIndex.get(s);
				rec.addToFile(fid);
				
				out.addToFile(rec);
			}
			
			out.writeFile(idxPath);
		}
	
		private int generateFamilyUID(String fname)
		{
			Collection<Integer> used = nameIndex.values();
			Random r = new Random();
			int id = fname.hashCode() ^ r.nextInt();
			while(used.contains(id) || id == -1) id = fname.hashCode() ^ r.nextInt();
			return id;
		}
		
		public int addOrReplaceFamily(Family f) throws IOException
		{
			//See if we're adding or replacing...
			if (f == null) return -1;
			String fname = f.getFamilyName();
			if(nameIndex.containsKey(fname))
			{
				//Replace
				int fuid = nameIndex.get(fname);
				//Clear if loaded
				usedQueue.remove(fuid);
				idMap.remove(fuid);
				//Overwrite file
				String fampath = srcDir + File.separator + Integer.toHexString(fuid) + ".fam";
				Family.writeToFAMI(f, fampath, true);
				//Reload
				idMap.put(fuid, f);
				usedQueue.add(fuid);
				return fuid;
			}
			else
			{
				//Add - Pretty straightforward
				int fuid = generateFamilyUID(fname);
				//Map Values
				nameIndex.put(f.getFamilyName(), fuid);
				if(idMap.size() >= MAX_FAM_IN_MEM)
				{
					//Unmap the family loaded the longest...
					int oldest = usedQueue.poll();
					idMap.remove(oldest);
				}
				idMap.put(fuid, f);
				usedQueue.add(fuid);
				//Write to file
				String fampath = srcDir + File.separator + Integer.toHexString(fuid) + ".fam";
				Family.writeToFAMI(f, fampath, true);
				return fuid;
			}
		}
		
		public int getFamilyID(String familyName)
		{
			Integer id = nameIndex.get(familyName);
			if (id == null) return -1;
			return id;
		}
		
		public int countFamilies()
		{
			return nameIndex.size();
		}

		
		public class FCIterator implements Iterator<Family>
		{

			private int nowIndex;
			private int size;
			private List<Integer> idList;
			
			public FCIterator()
			{
				nowIndex = 0;
				size = nameIndex.size();
				idList = new ArrayList<Integer>(size+1);
				idList.addAll(nameIndex.values());
				Collections.sort(idList);
			}
			
			@Override
			public boolean hasNext() 
			{
				return (nowIndex < size);
			}

			@Override
			public Family next() 
			{
				int fuid = idList.get(nowIndex);
				nowIndex++;
				Family f = getFamily(fuid);
				return f;
			}
			
		}
		
		public Iterator<Family> iterator() 
		{
			return new FCIterator();
		}
		
	}
	
	private FamilyCache fCache;
	
	private String sidx_path;
	private ConcurrentMap<String, Integer> nameMap;
	private ConcurrentMap<Integer, Integer> famIdMap;
	
	public DBSampleTable(String projectDir) throws IOException
	{
		String famdir = projectDir + File.separator + FAMI_DIR;
		fCache = new FamilyCache(famdir);
		sidx_path = famdir + File.pathSeparator + SINDEX_NAME;
		nameMap = new ConcurrentHashMap<String, Integer>();
		famIdMap = new ConcurrentSkipListMap<Integer, Integer>();
		readSampleIndex();
	}
	
	private void readSampleIndex() throws IOException
	{
		//SIDX Format (BE):
		//	Record [32]
		//		Sample UID [4]
		//		Family UID [4]
		//		Sample Name [24]
		
		nameMap.clear();
		famIdMap.clear();
		
		if (!FileBuffer.fileExists(sidx_path)) return;
		FileBuffer sidx = FileBuffer.createBuffer(sidx_path, true);
		long sz = sidx.getFileSize();
		long cpos = 0;
		while (cpos < sz)
		{
			int suid = sidx.intFromFile(cpos); cpos += 4;
			int fuid = sidx.intFromFile(cpos); cpos += 4;
			String sname = sidx.getASCII_string(cpos, SAMPLE_NAME_MAX_LEN); cpos += SAMPLE_NAME_MAX_LEN;
			nameMap.put(sname, suid);
			famIdMap.put(suid, fuid);
		}
		
	}
	
	private void writeSampleIndex() throws IOException
	{
		int scount = nameMap.size();
		FileBuffer out = new CompositeBuffer(scount);
		
		List<String> names = new ArrayList<String>(scount);
		names.addAll(nameMap.keySet());
		Collections.sort(names);
		
		for(String s : names)
		{
			int suid = nameMap.get(s);
			int fuid = famIdMap.get(suid);
			String sname = s;
			int nlen = sname.length();
			if (nlen > SAMPLE_NAME_MAX_LEN)			
			{
				sname = s.substring(0, SAMPLE_NAME_MAX_LEN);
				nlen = SAMPLE_NAME_MAX_LEN;
			}
			
			FileBuffer rec = new FileBuffer(32, true);
			rec.addToFile(suid);
			rec.addToFile(fuid);
			rec.printASCIIToFile(sname);
			while(nlen < SAMPLE_NAME_MAX_LEN)
			{
				rec.addToFile((byte)0x00);
				nlen++;
			}
			
			out.addToFile(rec);
		}
		
		out.writeFile(sidx_path);
	}
	
	public void saveTable() throws IOException
	{
		//Writes sample and family index files
		//Fami files should be written as they are modified...
		fCache.writeIndex();
		writeSampleIndex();
	}
	
	public Family getFamily(int uid)
	{
		return fCache.getFamily(uid);
	}
	
	public Family getFamily(String famName)
	{
		return fCache.getFamily(famName);
	}
	
	public FamilyMember getSample(int uid)
	{
		if(!famIdMap.containsKey(uid)) return null;
		int fuid = famIdMap.get(uid);
		Family f = fCache.getFamily(fuid);
		if (f == null) return null;
		
		return f.getMember(uid);
	}
	
	public FamilyMember getSample(String sampleName)
	{
		Integer suid = nameMap.get(sampleName);
		if(suid == null) return null;
		return getSample(suid);
	}
	
	private int generateNewSampleID(String sampleName)
	{
		Random r = new Random();
		Set<Integer> used = new TreeSet<Integer>();
		used.addAll(nameMap.values());
		int id = sampleName.hashCode() ^ r.nextInt();
		while(id == -1 || used.contains(id)) id = sampleName.hashCode() ^ r.nextInt();
		return id;
	}
	
	public boolean addOrReplaceFamily(Family f)
	{
		if(f == null) return false;
		//Before sending it to the cache method,
		//	make sure all samples have UNIQUE ids...
		//Also map to sample index!
		String fname = f.getFamilyName();
		int oldid = fCache.getFamilyID(fname);
		boolean newfam = (oldid != -1);
		List<FamilyMember> members = f.getAllFamilyMembers();
		for(FamilyMember mem : members)
		{
			//See if there is a sampleID for this person
			String mname = mem.getName();
			//The UID should be unique already if in table...
			if(!newfam)
			{
				if(!nameMap.containsKey(mname))
				{
					//If no hit, try looking up the sample ID directly...
					//	In case sample name was changed...
					if(famIdMap.containsKey(mem.getUID()))
					{
						//Either new family member with conflicting UID, or old one renamed
						int mappedid = famIdMap.get(mem.getUID());
						if (mappedid == oldid)
						{
							//Assume it's the same person renamed
							//Go through all keys until find one mapped to this UID
							//	and remove
							List<String> allsnames = new LinkedList<String>();
							allsnames.addAll(nameMap.keySet());
							for(String k : allsnames)
							{
								int v = nameMap.get(k);
								if (v == mem.getUID())
								{
									nameMap.remove(k);
									break;
								}
							}
							nameMap.put(mname, mem.getUID());
						}
						else
						{
							//Generate new UID until truly unique and add to maps
							//(We're assuming it's a new fam member)
							int id = generateNewSampleID(mname);
							mem.setUID(id);
							nameMap.put(mname, id);
							famIdMap.put(id, oldid);
						}
					}
					else
					{
						//Assume new family member
						//We already checked to see if map had ID, so we know
						//	that UID isn't currently in use
						//Add to maps
						nameMap.put(mname, mem.getUID());
						famIdMap.put(mem.getUID(), oldid);
					}
				}	
			}
			else
			{
				//Entire family is new.
				//Make sure ids are unique, then put in map
				if(famIdMap.containsKey(mem.getUID()))
				{
					int id = generateNewSampleID(mname);
					mem.setUID(id);
				}
				nameMap.put(mname, mem.getUID());
			}
		}
		
		int fuid = -1;
		try 
		{
			fuid = fCache.addOrReplaceFamily(f);
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
			return false;
		}
		
		if(fuid == -1) return false;
		
		if(newfam)
		{
			//Now we can add family IDs to sample index if new fam
			//(Added before for existing fam)
			for(FamilyMember mem : members)
			{
				famIdMap.put(mem.getUID(), fuid);
			}
		}
		
		return true;
		
	}

	public int countFamilies()
	{
		return fCache.countFamilies();
	}
	
	public int countSamples()
	{
		return famIdMap.size();
	}
	
	public int countSamplesInPopulation(Population p)
	{
		int tot = 0;
		for(Family f : fCache)
		{
			if (f != null)
			{
				List<FamilyMember> members = f.getAllFamilyMembers();
				for(FamilyMember mem : members)
				{
					if (mem.hasPopulationTag(p)) tot++;
				}
			}
		}
		
		return tot;
	}

	@Override
	public Iterator<Family> iterator() 
	{
		return fCache.iterator();
	}
	
}
