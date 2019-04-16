package hospelhornbg_svdb;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import hospelhornbg_bioinformatics.Interval;

public class SVDBGenotype {

	public static class SVDBAllele
	{
		private int count;
		private Interval range;
		
		public SVDBAllele(int aCount, int stPos, int edPos)
		{
			count = aCount;
			range = new Interval(stPos, edPos);
		}
		
		public int getAlleleCount()
		{
			return count;
		}
		
		public Interval getAllele()
		{
			return range;
		}
	}

	//private int variantUID;
	private int indivUID;
	private List<SVDBAllele> alleles;
	
	public SVDBGenotype(int individualUID, int uniqueAlleles)
	{
		indivUID = individualUID;
		alleles = new ArrayList<SVDBAllele>(uniqueAlleles+1);
	}
	
	public int getIndividualUID()
	{
		return indivUID;
	}
	
	public void addAllele(int count, int stPos, int edPos)
	{
		alleles.add(new SVDBAllele(count, stPos, edPos));
	}
	
	public Collection<SVDBAllele> getAlleles()
	{
		List<SVDBAllele> alist = new ArrayList<SVDBAllele>(alleles.size() + 1);
		alist.addAll(alleles);
		return alist;
	}
	
	public int getUniqueAlleleCount()
	{
		return alleles.size();
	}
	
}
