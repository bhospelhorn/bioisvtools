package hospelhornbg_svdb;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;

import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.FileBuffer;

public class RegionIndex {
	
	//Regions are 8192 bp
	//8192 = 2^13 (so shift right 13 from position to get region)
	
	/*
	 * Disk bin format (bin for each type)
	 * #Contigs [4]
	 * Contig Entry
	 * 		Contig ID [4]
	 * 		#Regions [2]
	 * 		Region Entry
	 * 			#Variants [4]
	 * 			Variant UIDs [8] * n	
	 * 
	 */
	
	private TypeMap map;
	
	public static class TypeMap
	{
		private boolean concurrent;
		private Map<SVType, ContigMap> map;
		
		public TypeMap(boolean threadsafe)
		{
			if(threadsafe) map = new ConcurrentHashMap<SVType, ContigMap>();
			else map = new HashMap<SVType, ContigMap>();
			concurrent = threadsafe;
		}
		
		public Collection<Long> getVariants(SVType type, Contig c, int block)
		{
			ContigMap cm = map.get(type);
			if(cm == null) return null;
			return cm.getVariants(c, block);
		}
		
		public Collection<Long> getBlock(SVType type, Contig c, int block)
		{
			ContigMap cm = map.get(type);
			if(cm == null) return null;
			return cm.getVariantsDirect(c, block);
		}
		
		public void addVariant(DBVariant var)
		{
			if(var == null) return;
			SVType t = var.getType();
			if(t == null) return;
			ContigMap cm = map.get(t);
			if(cm == null)
			{
				cm = new ContigMap(concurrent);
				map.put(t, cm);
			}
			cm.addVariant(var);
		}
		
		public static TypeMap readIndex(String pathStem, GenomeBuild gb, boolean threadsafe) throws IOException
		{
			TypeMap tm = new TypeMap(threadsafe);
			for(SVType t : SVType.values())
			{
				String tpath = pathStem + t.getString() + ".bin";
				if(FileBuffer.fileExists(tpath))
				{
					FileBuffer in = FileBuffer.createBuffer(tpath, true);
					ContigMap cm = ContigMap.readIndex(in, 0, threadsafe, gb);
					tm.map.put(t, cm);
				}
			}
			
			return tm;
		}
		
		public void saveIndex(String pathStem) throws IOException
		{
			for(SVType t : map.keySet())
			{
				String tpath = pathStem + t.getString() + ".bin";
				ContigMap cm = map.get(t);
				BufferedOutputStream bos = new BufferedOutputStream(new FileOutputStream(tpath));
				cm.serializeMe(bos);
				bos.close();
			}
		}
		
	}
	
	public static class ContigMap
	{
		private boolean concurrent;
		private Map<Contig, RegionMap> map;
		
		public ContigMap(boolean threadsafe)
		{
			if(threadsafe) map = new ConcurrentHashMap<Contig, RegionMap>();
			else map = new HashMap<Contig, RegionMap>();
			concurrent = threadsafe;
		}
		
		public Collection<Long> getVariants(Contig c, int block)
		{
			if(c == null) return null;
			RegionMap rm = map.get(c);
			if(rm == null) return null;
			return rm.getVariantsInBlock(block);
		}
		
		public Collection<Long> getVariantsDirect(Contig c, int block)
		{
			if(c == null) return null;
			RegionMap rm = map.get(c);
			if(rm == null) return null;
			return rm.getBlock(block);
		}

		private void addVUID(long vuid, Contig c, int start, int end)
		{
			int stblock = start >>> 13;
			int edblock = end >>> 13;
			RegionMap rm = map.get(c);
			if(rm == null)
			{
				rm = new RegionMap(c.getLength(), concurrent);
				map.put(c, rm);
			}
			for(int b = stblock; b <= edblock; b++) rm.addVariant(b, vuid);
		}
		
		public void addVariant(DBVariant var)
		{
			if(var.getType() == SVType.BND || var.getType() == SVType.TRA)
			{
				Contig c = var.getChrom();
				addVUID(var.getLongID(), c, var.getStartPosition().getStart(), var.getStartPosition().getEnd());
				addVUID(var.getLongID(), var.getEndChrom(), var.getEndPosition().getStart(), var.getEndPosition().getEnd());
			}
			else
			{
				Contig c = var.getChrom();
				int st = var.getStartPosition().getStart();
				int ed = var.getEndPosition().getEnd();
				addVUID(var.getLongID(), c, st, ed);
			}
		}
		
		public void serializeMe(OutputStream output) throws IOException
		{
			if(output == null) return;
			FileBuffer ibuff = new FileBuffer(4, true);
			ibuff.addToFile(map.size());
			
			output.write(ibuff.getBytes());
			List<Contig> clist = new ArrayList<Contig>(map.size()+1);
			clist.addAll(map.keySet());
			Collections.sort(clist);
			
			for(Contig c : clist)
			{
				ibuff.replaceInt(c.getUID(), 0);
				output.write(ibuff.getBytes());
				RegionMap rm = map.get(c);
				output.write(rm.serializeMe().getBytes());
			}
			
		}
		
		public static ContigMap readIndex(FileBuffer src, long stpos, boolean threadsafe, GenomeBuild gb)
		{
			if(src == null) return null;
			long cpos = stpos;
			ContigMap cm = new ContigMap(threadsafe);
			
			int ccount = src.intFromFile(cpos); cpos += 4;
			for(int i = 0; i < ccount; i++)
			{
				int cuid = src.intFromFile(cpos); cpos += 4;
				Contig c = gb.getContigByUID(cuid);
				if(c == null) return null;
				
				RegionMap rm = RegionMap.readIndex(src, cpos, threadsafe);
				cpos += rm.getSerializedSize();
				cm.map.put(c, rm);
			}
			
			return cm;
		}
		
	}
	
	public static class RegionMap
	{
		private List<Collection<Long>> var_array;
		
		public RegionMap(int reg_count, boolean threadsafe)
		{
			var_array = new ArrayList<Collection<Long>>(reg_count);
			if(threadsafe) for(int i = 0; i < reg_count; i++){
				ConcurrentHashMap<Long, Boolean> dummy = new ConcurrentHashMap<Long,Boolean>();
				var_array.add(dummy.keySet(true));
			}
			else for(int i = 0; i < reg_count; i++) var_array.add(new TreeSet<Long>());
		}
		
		public RegionMap(long contig_size, boolean threadsafe)
		{
			long sections = (contig_size >>> 13) + 1;
			int sec = (int)sections;
			var_array = new ArrayList<Collection<Long>>(sec);
			if(threadsafe) for(int i = 0; i < sec; i++){
				ConcurrentHashMap<Long, Boolean> dummy = new ConcurrentHashMap<Long,Boolean>();
				var_array.add(dummy.keySet(true));
			}
			else for(int i = 0; i < sec; i++) var_array.add(new TreeSet<Long>());
		}
		
		public Collection<Long> getVariantsInBlock(int block)
		{
			if(block < 0 || block >= var_array.size()) return null;
			Collection<Long> vars = new TreeSet<Long>();
			vars.addAll(var_array.get(block));
			return vars;
		}
		
		public Collection<Long> getBlock(int block)
		{
			if(block < 0 || block >= var_array.size()) return null;
			return var_array.get(block);
		}
		
		public void addVariant(int block, long varUID)
		{
			if(block < 0 || block >= var_array.size()) return;
			var_array.get(block).add(varUID);
		}
		
		public long getSerializedSize()
		{
			int bcount = var_array.size();
			int vcount = 0;
			for(Collection<Long> c : var_array) vcount += c.size();
			long size = 4; //#blocks
			size += (bcount << 2); //4 * #blocks (#var field for each block)
			size += (vcount << 3); //8 * total vars
			
			return size;
		}
		
		public FileBuffer serializeMe()
		{
			int bcount = var_array.size();
			//Count
			int vcount = 0;
			for(Collection<Long> c : var_array) vcount += c.size();
			
			//Size
			int size = 4; //#blocks
			size += (bcount << 2); //4 * #blocks (#var field for each block)
			size += (vcount << 3); //8 * total vars
			
			//Allocate
			FileBuffer out = new FileBuffer(size, true);
			
			//Add
			out.addToFile((short)bcount);
			for(Collection<Long> col : var_array)
			{
				out.addToFile(col.size());
				for(Long vuid : col) out.addToFile(vuid);
			}
				
			return out;
		}

		public static RegionMap readIndex(FileBuffer src, long stpos, boolean threadsafe)
		{
			long cpos = stpos;
			
			int sec = Short.toUnsignedInt(src.shortFromFile(cpos)); cpos+=2;
			RegionMap rm = new RegionMap(sec, threadsafe);
			for(int i = 0; i < sec; i++)
			{
				Collection<Long> reg = rm.var_array.get(i);
				int vcount = src.intFromFile(cpos); cpos += 4;
				for(int j = 0; j < vcount; j++)
				{
					reg.add(src.longFromFile(cpos)); cpos += 8;
				}
			}
			
			return rm;
		}
		
	}

	public RegionIndex(boolean threadsafe)
	{
		map = new TypeMap(threadsafe);
	}
	
	public static RegionIndex loadIndex(String pathstem, GenomeBuild genome, boolean threadsafe) throws IOException
	{
		RegionIndex index = new RegionIndex(threadsafe);
		index.map = TypeMap.readIndex(pathstem, genome, threadsafe);
		
		return index;
	}
	
	public void saveIndex(String pathstem) throws IOException
	{
		map.saveIndex(pathstem);
	}
	
	public Collection<Long> getVariantIDsInApproximateRegion(SVType type, Contig c, int start, int end)
	{
		Set<Long> set = new TreeSet<Long>();
		
		//Determine blocks...
		int st_block = start >>> 13;
		int ed_block = end >>> 13;
		for(int b = st_block; b <= ed_block; b++)
		{
			Collection<Long> col = map.getBlock(type, c, b);
			if(col != null) set.addAll(col);
		}
		
		return set;
	}
	
	public Collection<Long> getVariantIDsInApproximateRegion(Contig c, int start, int end, boolean ignoreTRA)
	{
		Set<Long> set = new TreeSet<Long>();
		
		//Determine blocks...
		int st_block = start >>> 13;
		int ed_block = end >>> 13;
		for (SVType type : SVType.values())
		{
			if(ignoreTRA && (type == SVType.TRA)) continue;
			if(ignoreTRA && (type == SVType.BND)) continue;
			for(int b = st_block; b <= ed_block; b++)
			{
				Collection<Long> col = map.getBlock(type, c, b);
				if(col != null) set.addAll(col);
			}
		}
		
		return set;
	}
	
	
	public void indexVariant(DBVariant var)
	{
		map.addVariant(var);
	}
	
}
