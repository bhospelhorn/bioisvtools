package hospelhornbg_svdb;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import hospelhornbg_bioinformatics.Genotype;
import hospelhornbg_bioinformatics.Interval;
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_segregation.FamilyMember;

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
	
	public boolean isHomozygous()
	{
		for(SVDBAllele a : alleles)
		{
			if (a.getAlleleCount() > 1) return true;
		}
		return false;
	}
	
	public static SVDBGenotype generateGenotype(FamilyMember sample, StructuralVariant sv)
	{
		String name = sample.getName();
		Genotype g = sv.getSampleGenotype(name);
		if(g != null)
		{
			SVDBGenotype geno = new SVDBGenotype(sample.getUID(), 1);
			//Get GT field...
			int[] calls = g.getAlleles();
			//We're really only looking at 1's and 0's
			int oneCount = 0;
			for(int a : calls) if (a == 1) oneCount++;
			geno.addAllele(oneCount, sv.getPosition(), sv.getEndPosition());
			return geno;
		}
		return null;
	}
	
}
