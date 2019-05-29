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

public class SQLVariantTable implements VariantTable{
	
	private String dbPath;
	
	public SQLVariantTable(String path)
	{
		//TODO
	}

	@Override
	public DBVariant getVariant(long varUID) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public VariantGenotype getGenotype(long varUID) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<Long> getVariantIDsForSample(int sampleUID) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<Long> getVariantIDsForFamily(Family fam) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<Long> getVariantIDsForSampleOfType(int sampleUID, SVType type) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Collection<DBVariant> getVariantsInRegion(Contig c, int start, int end) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Map<String, GeneHitCounter> generateGeneHitMap() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean addVCF(String vcfpath, Map<String, FamilyMember> sampleMap, int mergeFactor) throws IOException {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean removeSample(FamilyMember sample) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean removeFamily(Family fam) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean updateSampleCounts(DBSampleTable sampleTable) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public boolean updateSamplePopulationFlags(FamilyMember sample, Collection<Population> oldFlags) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public void save() throws IOException {
		// TODO Auto-generated method stub
		
	}

}
