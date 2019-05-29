package hospelhornbg_svdb;

import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_segregation.Family;
import hospelhornbg_segregation.FamilyMember;
import hospelhornbg_segregation.Population;
import hospelhornbg_svdb.DBVariantTable.GeneHitCounter;

public interface VariantTable {
	
	public DBVariant getVariant(long varUID);
	public VariantGenotype getGenotype(long varUID);
	public List<Long> getVariantIDsForSample(int sampleUID);
	public List<Long> getVariantIDsForFamily(Family fam);
	public List<Long> getVariantIDsForSampleOfType(int sampleUID, SVType type);
	public Collection<DBVariant> getVariantsInRegion(Contig c, int start, int end);
	public Map<String, GeneHitCounter> generateGeneHitMap();
	
	public boolean addVCF(String vcfpath, Map<String, FamilyMember> sampleMap, int mergeFactor) throws IOException;
	public boolean removeSample(FamilyMember sample);
	public boolean removeFamily(Family fam);
	public boolean updateSampleCounts(DBSampleTable sampleTable);
	public boolean updateSamplePopulationFlags(FamilyMember sample, Collection<Population> oldFlags);
	public void save() throws IOException;

}
