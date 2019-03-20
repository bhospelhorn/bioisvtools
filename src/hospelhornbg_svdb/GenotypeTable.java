package hospelhornbg_svdb;

import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collection;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import hospelhornbg_svdb.GenotypeIndex.IndexRecord;
import hospelhornbg_svdb.SVDBGenotype.SVDBAllele;
import waffleoRai_Utils.CompositeBuffer;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.StreamBuffer;


public class GenotypeTable {
	
	//Magic [8]
	//Version [4]
	//Reserved [4]
	//Records...
	
	//Format for non-TRA record:
	//VarID [4]
	//Number Indivs [4]
	
	//Genotype:
	//	Sample ID [4]
	//	# Unique Alleles [2]
	
	//	Alleles:
	//		Allele Copies [2]
	//		Start [4]
	//		End [4]
	
	public static final String MAGIC = "svdbGENO";
	public static final int VERSION = 1;
	public static final String EXT = "genot";
	
	private String tablePath;
	private GenotypeIndex index;
	
	private StreamBuffer openStream;
	private long streamPosition;
	//private long fileSize;
	
	public GenotypeTable(String filePath)
	{
		tablePath = filePath;
		streamPosition = -1;
		//First, try to read index...
		try 
		{
			String ipath = this.generateIndexPath();
			if (FileBuffer.fileExists(ipath)) index = GenotypeIndex.readIndex(ipath);
			else index = GenotypeIndex.generateIndex(this);
		}
		catch(Exception e)
		{
			System.err.println("GenotypeTable.<init> || WARNING: Table could not be indexed! Index will remain null!");
			e.printStackTrace();
		}
	}
	
	protected class IndexUnit
	{
		public int varID;
		public long offset;
		public int len;
	}
	
	protected void openStream() throws IOException
	{
		openStream = new StreamBuffer(tablePath, true);
		streamPosition = 16;
	}
	
	protected void closeStream()
	{
		openStream = null;
		streamPosition = -1;
	}
	
	protected IndexUnit getIndexingInfoForNextVariant()
	{
		if (openStream == null) return null;
		long fsz = FileBuffer.fileSize(tablePath);
		if (streamPosition >= fsz) return null;
		
		IndexUnit iu = new IndexUnit();
		iu.offset = streamPosition;
		iu.len = 0;
		
		iu.varID = openStream.intFromFile(streamPosition); 
		streamPosition += 4;
		
		int gcount = openStream.intFromFile(streamPosition);
		for (int i = 0; i < gcount; i++)
		{
			//Skip sampleID
			streamPosition += 4;
			short nall = openStream.shortFromFile(streamPosition); streamPosition += 2;
			//Skip nall * 10
			streamPosition += (nall * 10);
		}
		iu.len = (int)(streamPosition - iu.offset);
		
		return iu;
	}
	
	public List<SVDBGenotype> getGenotypesForVariant(int vuid) throws IOException
	{
		List<SVDBGenotype> list = new LinkedList<SVDBGenotype>();
		if (index == null) return list;
		IndexRecord r = index.getIndexRecord(vuid);
		if (r == null) return list;
		
		long stpos = r.getOffset();
		long edpos = r.getOffset() + r.getLength();
		FileBuffer rec = new FileBuffer(tablePath, stpos, edpos, true);
		long cpos = 4;
		int gcount = rec.intFromFile(cpos); cpos += 4;
		
		for (int i = 0; i < gcount; i++)
		{
			int sid = rec.intFromFile(cpos); cpos += 4;
			int acount = (int)rec.shortFromFile(cpos); cpos += 2;
			SVDBGenotype g = new SVDBGenotype(sid, acount);
			for (int j = 0; j < acount; j++)
			{
				int count = (int)rec.shortFromFile(cpos); cpos += 2;
				int st = rec.intFromFile(cpos); cpos += 4;
				int ed = rec.intFromFile(cpos); cpos += 4;
				g.addAllele(count, st, ed);
			}
			list.add(g);
		}
		
		return list;
	}
	
	public Map<Integer, SVDBGenotype> getGenotypesForSample(int sampleUID)
	{
		if (index == null) return null;
		try {openStream();}
		catch(IOException e)
		{
			e.printStackTrace();
			return null;
		}
		
		Map<Integer, SVDBGenotype> gmap = new TreeMap<Integer, SVDBGenotype>();
		//Scan entire file...
		long fsz = FileBuffer.fileSize(tablePath);
		while(streamPosition < fsz)
		{
			int varid = openStream.intFromFile(streamPosition); streamPosition += 4;
			int icount = openStream.intFromFile(streamPosition); streamPosition += 4;
			for (int i = 0; i < icount; i++)
			{
				int sid = openStream.intFromFile(streamPosition); streamPosition += 4;
				short acount = openStream.shortFromFile(streamPosition); streamPosition += 2;
				if (sid == sampleUID)
				{
					//Read genotype data
					SVDBGenotype g = new SVDBGenotype(sampleUID, acount);
					for(short j = 0; j < acount; j++)
					{
						int c = (int)openStream.shortFromFile(streamPosition); streamPosition += 2;
						int st = openStream.intFromFile(streamPosition); streamPosition += 4;
						int ed = openStream.intFromFile(streamPosition); streamPosition += 4;
						g.addAllele(c, st, ed);
					}
					gmap.put(varid, g);
				}
				else streamPosition += (10 * acount);
			}
		}
		
		
		closeStream();
		return gmap;
	}
	
	private byte[] serializeVariantRecord(int varID, List<SVDBGenotype> genoList)
	{
		int gcount = genoList.size();
		CompositeBuffer rec = new CompositeBuffer(1 + gcount);
		
		//VID
		FileBuffer idinfo = new FileBuffer(8, true);
		idinfo.addToFile(varID);
		idinfo.addToFile(gcount);
		rec.addToFile(idinfo);
		
		//Genotypes
		for(SVDBGenotype g : genoList)
		{
			int nall = g.getUniqueAlleleCount();
			FileBuffer sgeno = new FileBuffer(6 + (nall * 10), true);
			sgeno.addToFile(g.getIndividualUID());
			sgeno.addToFile((short)nall);
			Collection<SVDBAllele> all = g.getAlleles();
			for (SVDBAllele a : all)
			{
				sgeno.addToFile((short)a.getAlleleCount());
				sgeno.addToFile(a.getAllele().getStart());
				sgeno.addToFile(a.getAllele().getEnd());
			}
			rec.addToFile(sgeno);
		}
		
		//Convert to byte array
		int size = (int)rec.getFileSize();
		byte[] arr = new byte[size];
		
		for (int i = 0; i < size; i++)
		{
			arr[i] = rec.getByte(i);
		}
		
		return arr;
	}
	
	private static byte[] serializeHeader()
	{
		FileBuffer header = new FileBuffer(16, true);
		header.printASCIIToFile(MAGIC);
		header.addToFile(VERSION);
		header.addToFile(0);
		byte[] arr = new byte[16];
		for (int i = 0; i < 16; i++)
		{
			arr[i] = header.getByte(i);
		}
		
		return arr;
	}
	
	public List<Integer> removeSampleFromGenotypeTable(int sampleUID) throws IOException
	{
		if (index == null) return null;
		//returns the list of variant IDs where requested sample was only sample
		
		String temppath = FileBuffer.generateTemporaryPath("GenotypeTableTemp");
		
		//Copies the file to the temp while removing requested data.
		//Moves temp file back to original path when done and re-indexes.
		
		//Uses index to read variant by variant :)
		List<Integer> outids = new LinkedList<Integer>();
		List<Integer> allids = index.getAllVariantIDs();
		
		FileOutputStream myStream = new FileOutputStream(temppath);
		myStream.write(serializeHeader());
		
		for(Integer vid : allids)
		{
			List<SVDBGenotype> genos = this.getGenotypesForVariant(vid);
			List<SVDBGenotype> ngenos = new LinkedList<SVDBGenotype>();
			//See if this sample has a geno for this variant...
			for(SVDBGenotype g : genos)
			{
				if (g.getIndividualUID() != sampleUID)
				{
					ngenos.add(g);	
				}
			}
			if(ngenos.isEmpty())
			{
				outids.add(vid);
			}
			else
			{
				myStream.write(serializeVariantRecord(vid, ngenos));
			}
		}
		
		myStream.close();
		
		//Move file back...
		Files.move(Paths.get(temppath), Paths.get(tablePath));
		
		//Reindex
		index = GenotypeIndex.generateIndex(this);
		index.writeIndex(generateIndexPath());
		
		return outids;
	}
	
	public void addGenotypes(Map<Integer, Collection<SVDBGenotype>> additions) throws IOException
	{
		String temppath = FileBuffer.generateTemporaryPath("GenotypeTableTemp");
		Set<Integer> set = new TreeSet<Integer>();
		set.addAll(index.getAllVariantIDs());
		set.addAll(additions.keySet());
		
		List<Integer> allids = new LinkedList<Integer>();
		allids.addAll(set);
		Collections.sort(allids);
		
		FileOutputStream myStream = new FileOutputStream(temppath);
		myStream.write(serializeHeader());
		
		for(Integer vid : allids)
		{
			List<SVDBGenotype> genos = this.getGenotypesForVariant(vid);
			if (additions.containsKey(vid))
			{
				//Add or replace with new genotypes!
				List<SVDBGenotype> ngenos = new LinkedList<SVDBGenotype>();
				ngenos.addAll(genos);
				Collection<SVDBGenotype> addlist = additions.get(vid);
				for(SVDBGenotype g : addlist)
				{
					//Make sure sample isn't in current list...
					int i = -1;
					int ind = 0;
					for (SVDBGenotype og : ngenos)
					{
						if (g.getIndividualUID() == og.getIndividualUID())
						{
							i = ind;
							break;
						}
						ind++;
					}
					if(i >= 0) ngenos.remove(i);
					ngenos.add(g);
				}
				myStream.write(serializeVariantRecord(vid, ngenos));
			}
			else
			{
				//Just rewrite!
				myStream.write(serializeVariantRecord(vid, genos));
			}
		}
		
		myStream.close();

		//Move file back...
		Files.move(Paths.get(temppath), Paths.get(tablePath));
		
		//Reindex
		index = GenotypeIndex.generateIndex(this);
		index.writeIndex(generateIndexPath());
		
	}

	public String getFilePath()
	{
		return tablePath;
	}
	
	public String generateIndexPath()
	{
		int dot = tablePath.lastIndexOf('.');
		if (dot >= 0) return tablePath.substring(0, dot) + "." + GenotypeIndex.EXT;
		return tablePath + "." + GenotypeIndex.EXT;
	}

	public void generateEmptyTable() throws IOException
	{
		FileBuffer t = new FileBuffer(16, true);
		t.printASCIIToFile(MAGIC);
		t.addToFile(VERSION);
		t.addToFile(0);
		t.writeFile(tablePath);
		
		index.writeIndex(generateIndexPath());
	}
	
}
