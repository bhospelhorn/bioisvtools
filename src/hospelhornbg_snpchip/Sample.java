package hospelhornbg_snpchip;


import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import hospelhornbg_bioinformatics.AffectedStatus;
import hospelhornbg_bioinformatics.Sex;

public class Sample {
	
	/* ----- Instance Variables ----- */
	
	private String name;
	
	private Sex sex;
	private AffectedStatus affected;
	
	private Sample parent1;
	private Sample parent2;
	
	private Set<Sample> children;
	
	private Map<SNP, SNPGeno> snps;
	
	/* ----- Construction ----- */
	
	public Sample(String sampleName, Sex s)
	{
		name = sampleName;
		sex = s;
		snps = new HashMap<SNP, SNPGeno>();
		parent1 = null;
		parent2 = null;
		affected = AffectedStatus.UNAFFECTED;
		children = new HashSet<Sample>();
	}
	
	/* ----- Getters ----- */
	
	public List<SNP> getSNPList()
	{
		int nsnps = snps.size();
		if (nsnps < 1) return null;
		List<SNP> mylist = new ArrayList<SNP>(snps.size());
		mylist.addAll(snps.keySet());
		Collections.sort(mylist);
		return mylist;
	}
	
	public SNPGeno getGeno(SNP snp)
	{
		return snps.get(snp);
	}
	
	/* ----- Setters ----- */
	
	public void setName(String sampleName)
	{
		name = sampleName;
	}
	
	public void setSex(Sex s)
	{
		sex = s;
	}
	
	public void setAffected(AffectedStatus as)
	{
		affected = as;
	}
	
	public void setParent1(Sample p1)
	{
		parent1.removeChild(this);
		parent1 = p1;
		parent1.addChild(this);
	}
	
	public void setParent2(Sample p2)
	{
		parent2.removeChild(this);
		parent2 = p2;
		parent2.addChild(this);
	}
	
	public void addSNP(SNP snp, SNPGeno geno)
	{
		snps.put(snp, geno);
	}
	
	public void clearSNPGenotypes()
	{
		snps.clear();
	}
	
	private void addChild(Sample c)
	{
		children.add(c);
	}
	
	private void removeChild(Sample c)
	{
		children.remove(c);
	}
	
	/* ----- Phasing ----- */
	
	public void phaseWithParents(boolean overwrite)
	{
		if (parent1 == null)
		{
			if (parent2 == null) return;
			parent1 = parent2;
			parent2 = null;
			phaseWithOneParent(overwrite);
			return;
		}
		//Only reach this point if parent1 != null
		if (parent2 == null)
		{
			phaseWithOneParent(overwrite);
			return;
		}
		phaseWithBothParents(overwrite);
	}
	
	private void phaseWithOneParent(boolean overwrite)
	{
		List<SNP> snplist = this.getSNPList();
		if (snplist == null) return;
		if (snplist.isEmpty()) return;
		
		//This method is private. If it's called, the only parent should be set as parent 1
		//It also doesn't present terribly good results.
		//All it can do is mark SNPs that could not have been inherited from known parent!
		for (SNP snp : snplist)
		{
			//Pull SNP geno from parent
			SNPGeno pgeno = parent1.getGeno(snp);
			//Pull SNP geno from child
			SNPGeno mygeno = getGeno(snp);
			//Erase any phasing marks currently present on child geno (if overwrite)
			if (overwrite) mygeno.clearPhasingMarks();
			//Note: Needs to work differently for MT and X and Y
		}
	}
	
	private void phaseWithBothParents(boolean overwrite)
	{
		//NOTE!! Write a method that looks for UPD? Or include here?
	}
	
	public void phaseWithChildren()
	{
		//Must have at least two children!
	}

}
