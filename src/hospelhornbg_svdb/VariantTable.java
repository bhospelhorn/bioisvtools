package hospelhornbg_svdb;

import java.io.IOException;
import java.sql.SQLException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;

import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.GenomeBuild;
import hospelhornbg_segregation.Family;
import hospelhornbg_segregation.FamilyMember;
import hospelhornbg_segregation.Population;
import hospelhornbg_svdb.DBVariantTable.GeneHitCounter;

public interface VariantTable {
	
	public static class GenomeIndex
	{
		//For correlating variant UIDs to contig/pos!
		private ConcurrentMap<Contig, Integer> shortIDMap;
		
		public GenomeIndex(GenomeBuild gb)
		{
			shortIDMap = new ConcurrentHashMap<Contig, Integer>();
			if(gb == null) return;
			List<Contig> clist = gb.getChromosomes();
			Set<Integer> usedids = new TreeSet<Integer>();
			for(Contig c : clist)
			{
				int id = (int)(c.getLength() & 0x7FFFL);
				while(usedids.contains(id)) {
					id += 1;
					if(id > 0x7FFF) id = 0;
				}
				usedids.add(id);
				shortIDMap.put(c, id);
			}
		}
		
		public int getPositionDivision(Contig c, int pos)
		{
			int divsz = (int)c.getLength()/0xFFFF;
			if(divsz == 0) return 0;
			return pos/divsz;
		}
		
		public long getVarIDPrefix(Contig c, int pos)
		{
			if (c == null) return 0;
			if (pos < 0) return 0;
			long prefix = 0;
			int cid = shortIDMap.get(c);
			prefix = Integer.toUnsignedLong(cid) << 48;
			prefix &= 0xFFFF000000000000L;
			long pdiv = Integer.toUnsignedLong(getPositionDivision(c, pos) & 0xFFFF);
			pdiv = pdiv << 32;
			prefix |= pdiv;
			return prefix;
		}
		
		public int getContigPrefix(Contig c)
		{
			return shortIDMap.get(c);
		}
		
	}
	
	public DBVariant getVariant(long varUID);
	public VariantGenotype getGenotype(long varUID);
	public List<Long> getVariantIDsForSample(int sampleUID);
	public List<Long> getVariantIDsForFamily(Family fam);
	public List<Long> getVariantIDsForSampleOfType(int sampleUID, SVType type);
	public Collection<DBVariant> getVariantsInRegion(Contig c, int start, int end);
	public Map<String, GeneHitCounter> generateGeneHitMap();
	
	public boolean addVCF(String vcfpath, Map<String, FamilyMember> sampleMap, int mergeFactor, boolean ignoreTRA, int threads) throws IOException;
	public boolean removeSample(FamilyMember sample);
	public boolean removeFamily(Family fam);
	public boolean updateSampleCounts(DBSampleTable sampleTable);
	public boolean updateSamplePopulationFlags(FamilyMember sample, Collection<Population> oldFlags);
	public void save() throws IOException;
	public void close() throws SQLException;
	
	public void dumpTable(String directory);
	public void clearVariantTable();
	
	
}
