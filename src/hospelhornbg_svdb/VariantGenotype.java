package hospelhornbg_svdb;

import java.io.IOException;
import java.io.InputStream;
import java.sql.Blob;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;

import hospelhornbg_bioinformatics.Genotype;
import hospelhornbg_bioinformatics.SVType;
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
	
	private boolean corrupted_flag;
	
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
	
	public void readDataFromBLOB(Blob blob) throws SQLException, IOException
	{
		if(blob == null) return;
		//Wrap into FileBuffer
		FileBuffer loader = new FileBuffer((int)blob.length(), true);
		//for(byte b : data) loader.addToFile(b);
		InputStream is = blob.getBinaryStream();
		int b = -1;
		while((b = is.read()) != -1) loader.addToFile((byte)b);
		is.close();
		readDataFromBLOB(loader);
	}
	
	public void readDataFromBLOB(FileBuffer data)
	{
		if(data == null) return;
		long cpos = 0;
		long fsz = data.getFileSize();
		while(cpos < fsz)
		{
			int sid = -1;
			try{
				sid = data.intFromFile(cpos); cpos += 4;
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
			catch(IndexOutOfBoundsException x)
			{
				System.err.println("Corrupted genotype record: cpos = 0x" + Long.toHexString(cpos));
				System.err.println("sampleID = 0x" + Integer.toHexString(sid));
				//x.printStackTrace();
				corrupted_flag = true;
			}
		}
	}
	
	public int calculateSerializedSize(boolean include_vid)
	{
		int sz = 4;
		if(include_vid) sz += 8;
		for(SVDBGenotype gt : gMap.values())
		{
			sz += 6 + (10 * gt.getUniqueAlleleCount());
		}
		
		return sz;
	}
	
	public FileBuffer serializeForGENOT(boolean include_vid)
	{
		int sz = calculateSerializedSize(include_vid);
		FileBuffer outbuff = new FileBuffer(sz, true);
		
		//Header
		int icount = gMap.size();
		if(include_vid) outbuff.addToFile(varUID);
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
	
	public byte[] getGenotypesAsBLOBBytes()
	{
		if(gMap.isEmpty()) {byte[] barr = {-1}; return barr;}
		FileBuffer me = serializeForGENOT(true); //...For now. Just for backward compatibility
		//Knock off the first 12 bytes...
		return me.getBytes(12, me.getFileSize());
	}
	
	public FileBuffer getGenotypesAsBLOB()
	{
		FileBuffer me = serializeForGENOT(true);
		try {return me.createReadOnlyCopy(12, me.getFileSize());} 
		catch (IOException e) {e.printStackTrace();}
		return null;
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

	public boolean isEmpty()
	{
		return gMap.isEmpty();
	}
	
	public Genotype getAsGenotypeObject(int sampleUID, SVType svtype, boolean hasData)
	{
		SVDBGenotype localg = gMap.get(sampleUID);
		Genotype g = new Genotype();
		if(localg == null)
		{
			//Ref/Ref or ./.
			if(hasData) g.setAlleles("0/0");
			else g.setAlleles("./.");	
			g.setCopyNumber(2);
		}
		else
		{
			Collection<SVDBAllele> alist = localg.getAlleles();
			//TODO maybe something more elegant some day...
			//Right now, just make all alt counts a "1" allele?
			int altcount = 0;
			for(SVDBAllele a : alist)
			{
				altcount += a.getAlleleCount();
			}
			
		
			if (altcount == 1)
			{
				g.setAlleles("0/1");
				if(svtype == SVType.DEL || svtype == SVType.DELME) g.setCopyNumber(1);
				else if (svtype == SVType.DUP || svtype == SVType.TANDEM) g.setCopyNumber(3);
				else g.setCopyNumber(2);
			}
			else if(altcount == 2)
			{
				g.setAlleles("1/1");
				if(svtype == SVType.DEL || svtype == SVType.DELME) g.setCopyNumber(0);
				else if (svtype == SVType.DUP || svtype == SVType.TANDEM) g.setCopyNumber(4);
				else g.setCopyNumber(2);
			}
			else if (altcount > 2)
			{
				//Weird. Should probably refine if this ever comes up.
				int[] aarr = new int[altcount];
				for(int i = 0; i < altcount; i++) aarr[i] = 1;
				g.setAlleles(aarr);
				g.setCopyNumber(altcount);
			}
			else if(altcount == 0)
			{
				//Ref/ref
				//Shouldn't come up if we get to this block, but just in case.
				g.setAlleles("0/0");
				g.setCopyNumber(2);
			}
		}
		
		return g;
	}
	
	public boolean isCorrupted(){return this.corrupted_flag;}
	
}
