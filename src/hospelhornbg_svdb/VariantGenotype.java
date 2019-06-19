package hospelhornbg_svdb;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;

import hospelhornbg_svdb.SVDBGenotype.SVDBAllele;
import waffleoRai_Utils.FileBuffer;

public class VariantGenotype {
	
	//GENOT Format
	
	//Magic [8]
	//Version [4]
	//Reserved [4]
	//Records...
		
	//Format for non-TRA record:
		//VarID [8]
		//Number Indivs [4]
		
	//Genotype:
		//	Sample ID [4]
		//	# Unique Alleles [2]
		
	//	Alleles:
		//		Allele Copies [2]
		//		Start [4]
		//		End [4]
	
	private long varUID;
	private ConcurrentMap<Integer, SVDBGenotype> gMap;
	
	public VariantGenotype(long uid)
	{
		varUID = uid;
		gMap = new ConcurrentHashMap<Integer, SVDBGenotype>();
	}
	
	public static VariantGenotype readFromGENOT(FileBuffer src, long stoff)
	{
		if(src == null) return null;
		
		long cpos = stoff;
		long vid = src.longFromFile(cpos); cpos += 8;
		VariantGenotype vg = new VariantGenotype(vid);
		int icount = src.intFromFile(cpos); cpos += 4;
		for(int i = 0; i < icount; i++)
		{
			SVDBGenotype gt = readSingleGenotypeFromGENOT(src, cpos);
			cpos += 6 + (10 * gt.getUniqueAlleleCount());
			vg.gMap.put(gt.getIndividualUID(), gt);
		}
		
		return vg;
	}
	
	public static SVDBGenotype readSingleGenotypeFromGENOT(FileBuffer src, long stoff)
	{
		long cpos = stoff;
		int sid = src.intFromFile(cpos); cpos += 4;
		int uniqueAlleles = Short.toUnsignedInt(src.shortFromFile(cpos)); cpos += 2;
		
		SVDBGenotype gt = new SVDBGenotype(sid, uniqueAlleles);
		for(int i = 0; i < uniqueAlleles; i++)
		{
			int c = Short.toUnsignedInt(src.shortFromFile(cpos)); cpos += 2;
			int st = src.intFromFile(cpos); cpos += 4;
			int ed = src.intFromFile(cpos); cpos += 4;
			gt.addAllele(c, st, ed);
		}
		
		return gt;
	}
	
	public void readDataFromBLOB(byte[] data)
	{
		if(data == null) return;
		//Wrap into FileBuffer
		FileBuffer loader = new FileBuffer(data.length, true);
		for(byte b : data) loader.addToFile(b);
		readDataFromBLOB(loader);
	}
	
	public void readDataFromBLOB(FileBuffer data)
	{
		if(data == null) return;
		long cpos = 0;
		long fsz = data.getFileSize();
		while(cpos < fsz)
		{
			int sid = data.intFromFile(cpos); cpos += 4;
			int acount = Short.toUnsignedInt(data.shortFromFile(cpos)); cpos += 2;
			
			SVDBGenotype sgeno = new SVDBGenotype(sid, acount);
			for(int i = 0; i < acount; i++)
			{
				int acopies = Short.toUnsignedInt(data.shortFromFile(cpos)); cpos += 2;
				int astart = data.intFromFile(cpos); cpos += 4;
				int aend = data.intFromFile(cpos); cpos += 4;
				sgeno.addAllele(acopies, astart, aend);
			}
			gMap.put(sid, sgeno);
		}
	}
	
	public int calculateSerializedSize()
	{
		int sz = 12;
		for(SVDBGenotype gt : gMap.values())
		{
			sz += 6 + (10 * gt.getUniqueAlleleCount());
		}
		
		return sz;
	}
	
	public FileBuffer serializeForGENOT()
	{
		int sz = calculateSerializedSize();
		FileBuffer outbuff = new FileBuffer(sz, true);
		
		//Header
		int icount = gMap.size();
		outbuff.addToFile(varUID);
		outbuff.addToFile(icount);
		
		//Genotypes
		List<Integer> sids = new ArrayList<Integer>(icount+1);
		sids.addAll(gMap.keySet());
		Collections.sort(sids);
		
		for(Integer sid : sids)
		{
			SVDBGenotype gt = gMap.get(sid);
			int acount = gt.getUniqueAlleleCount();
			outbuff.addToFile(sid);
			outbuff.addToFile((short)acount);
			Collection<SVDBAllele> alist = gt.getAlleles();
			for(SVDBAllele a : alist)
			{
				outbuff.addToFile((short)a.getAlleleCount());
				outbuff.addToFile(a.getAllele().getStart());
				outbuff.addToFile(a.getAllele().getEnd());
			}
		}
		
		return outbuff;
	}
	
	public byte[] getGenotypesAsBLOB()
	{
		if(gMap.isEmpty()) {byte[] barr = {-1}; return barr;}
		FileBuffer me = serializeForGENOT();
		//Knock off the first 12 bytes...
		return me.getBytes(12, me.getFileSize());
	}
	
	public long getVariantUID()
	{
		return varUID;
	}
	
	public SVDBGenotype getGenotype(int sampleUID)
	{
		return gMap.get(sampleUID);
	}

	public boolean hasGenotypeFor(int sampleUID)
	{
		return (gMap.containsKey(sampleUID));
	}
	
	public void addGenotype(SVDBGenotype gt)
	{
		if(gt == null) return;
		gMap.put(gt.getIndividualUID(), gt);
	}

	public int genotypeRecords()
	{
		return gMap.size();
	}
	
	public SVDBGenotype removeGenotype(int sampleUID)
	{
		return gMap.remove(sampleUID);
	}
	
	public Collection<SVDBGenotype> getGenotypes()
	{
		List<SVDBGenotype> list = new ArrayList<SVDBGenotype>(gMap.size() + 1);
		list.addAll(gMap.values());
		return list;
	}

	public Set<Integer> getAllIndividuals()
	{
		Set<Integer> set = new HashSet<Integer>();
		set.addAll(gMap.keySet());
		return set;
	}
	
}
