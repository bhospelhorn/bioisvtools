package hospelhornbg_segregation;

import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import hospelhornbg_bioinformatics.Genotype;

//TODO: What about pseudo autosomal regions and regions that pair up between X and Y?
public class SexChromInheritor {
	
	public void checkMammalXCandidate(Map<Individual, Genotype> genomap, Individual pb, Candidate c)
	{
		Genotype pbg = genomap.get(pb);
		if (pbg == null) c.setInheritancePattern(pb, Inheritance.UNRESOLVED);
		
		boolean het = false;
		if (pb.getExpectedXCount() > 1)
		{
			het = !pbg.isHomozygous();
		}
		
		List<Individual> unaff = new LinkedList<Individual>();
		List<Individual> aff = new LinkedList<Individual>();
		for(Individual i : genomap.keySet())
		{
			if(i.isAffected()) aff.add(i);
			else unaff.add(i);
		}
		
		int allele = c.getAllele();
		//Do any unaffected have this allele?
		boolean uhave = false;
		boolean uhom = false;
		boolean uhemi = false;
		for(Individual i : unaff)
		{
			Genotype g = genomap.get(i);
			if (g == null) continue;
			if (g.hasAllele(allele))
			{
				uhave = true;
				if (i.getExpectedXCount() == 1) uhemi = true;
				else if (g.isHomozygous()) uhom = true;
			}
		}
		if (het)
		{	
			if(uhave)
			{
				//At least one unaff has this allele
				if(uhom)
				{
					//There is at least one unaff XX homozygote for this allele
					c.setInheritancePattern(pb, Inheritance.UNRESOLVED);
				}
				else if (!uhom && uhemi)
				{
					//There is at least one X(N) hemizygote for this allele
					if (patXRescue(genomap, aff, c)) {
						if (patXRescueDom(genomap, aff, c)) c.setInheritancePattern(pb, Inheritance.X_PATIMPRINT_DOM);
						else c.setInheritancePattern(pb, Inheritance.X_PATIMPRINT_HH);
					}
					else c.setInheritancePattern(pb, Inheritance.UNRESOLVED);
				}
				else if (!uhom && !uhemi)
				{
					//Unaffs are all hets
					//Check affs.
					boolean allhave = true;
					for (Individual i : aff)
					{
						Genotype g = genomap.get(i);
						if (g == null) continue;
						if (!g.hasAllele(allele))
						{
							allhave = false;
							break;
						}
					}
					//Could allele have come from parents? (Check De Novo)
					Individual mom = pb.getMother(); //No hemis, so assume didn't come from dad.
					if (mom != null)
					{
						Genotype mgeno = genomap.get(mom);
						if (mgeno != null)
						{
							if (mgeno.hasAllele(allele)) {
								if(allhave)c.setInheritancePattern(pb, Inheritance.X_LINKED_DOM);
								else c.setInheritancePattern(pb, Inheritance.HALF_HET);
							}
							else {
								if(allhave)c.setInheritancePattern(pb, Inheritance.X_LINKED_DN);
								else c.setInheritancePattern(pb, Inheritance.DENOVO_HET);
							}
						}
						else {
							if(allhave)c.setInheritancePattern(pb, Inheritance.X_LINKED_DN);
							else c.setInheritancePattern(pb, Inheritance.DENOVO_HET);
						}
					}
					else {
						if(allhave)c.setInheritancePattern(pb, Inheritance.X_LINKED_DOM);
						else c.setInheritancePattern(pb, Inheritance.HALF_HET);
					}
				}
			}
			else
			{
				//No unaffs have this allele
				//Do any affs NOT have this allele?
				boolean allhave = true;
				for (Individual i : aff)
				{
					Genotype g = genomap.get(i);
					if (g == null) continue;
					if (!g.hasAllele(allele))
					{
						allhave = false;
						break;
					}
				}
				//Check for DN...
				boolean isDN = false;
				Individual mom = pb.getMother(); //No hemis, so assume didn't come from dad.
				if (mom != null)
				{
					Genotype mgeno = genomap.get(mom);
					if (mgeno != null)
					{
						isDN = !mgeno.hasAllele(allele);
					}
				}
				if(isDN)
				{
					if(allhave) c.setInheritancePattern(pb, Inheritance.X_LINKED_DN);
					else c.setInheritancePattern(pb, Inheritance.DENOVO_HET);
				}
				else
				{
					if(allhave) c.setInheritancePattern(pb, Inheritance.X_LINKED_DOM);
					else c.setInheritancePattern(pb, Inheritance.HALF_HET);	
				}
			}
			
		}
		else
		{
			//PB is homozygous or hemizygous for this allele
			if (uhom)
			{
				c.setInheritancePattern(pb, Inheritance.UNRESOLVED);
			}
			else if (!uhom && uhemi)
			{
				//If mother is het and mother's X is maternal or PB has MV, could be imprint rescued
				if(pb.getExpectedXCount() > 1)
				{
					if(patXRescue(genomap, aff, c))
					{
						if (patXRescueDom(genomap, aff, c)) c.setInheritancePattern(pb, Inheritance.X_PATIMPRINT_DOM);
						else c.setInheritancePattern(pb, Inheritance.X_PATIMPRINT_HH);
					}
					else c.setInheritancePattern(pb, Inheritance.UNRESOLVED);
				}
				else c.setInheritancePattern(pb, Inheritance.UNRESOLVED);
			}
			else
			{
				//There are no unaffected homozygotes or hemizygotes
				boolean allhave = true;
				boolean allhom = true;
				for (Individual i : aff)
				{
					Genotype g = genomap.get(i);
					if (g == null) continue;
					if (!g.hasAllele(allele))
					{
						allhave = false;
						allhom = false;
						break;
					}
					else
					{
						if (i.getExpectedXCount() > 1 && !g.isHomozygous()) allhom = false;
					}
				}
				
				if (uhave)
				{
					//At least one unaff has this allele (is het, not hemi)
					if(!allhom)
					{
						//Halfhet (check DN)
						if (checkXDeNovo(genomap, pb, allele)) c.setInheritancePattern(pb, Inheritance.DENOVO_HET);
						else c.setInheritancePattern(pb, Inheritance.HALF_HET);
					}
					else
					{
						//X-rec (check DN)
						if (checkXDeNovo(genomap, pb, allele)) c.setInheritancePattern(pb, Inheritance.X_LINKED_MV);
						else c.setInheritancePattern(pb, Inheritance.X_LINKED_REC);
					}
				}
				else
				{
					//No unaff have this allele
				}
			}
		}
		
	}
	
	private boolean patXRescue(Map<Individual, Genotype> genomap, List<Individual> aff, Candidate c)
	{
		for(Individual a : aff)
		{
			if (a.getExpectedXCount() == 1)
			{
				Genotype ag = genomap.get(a);
				if (ag != null)
				{
					if (ag.hasAllele(c.getAllele())) return false;
				}
			}
		}
		
		return true;
	}
	
	private boolean patXRescueDom(Map<Individual, Genotype> genomap, List<Individual> aff, Candidate c)
	{
		for(Individual a : aff)
		{
			if (a.getExpectedXCount() > 1)
			{
				Genotype ag = genomap.get(a);
				if (ag != null)
				{
					if (!ag.hasAllele(c.getAllele())) return false;
				}
			}
		}
		
		return true;
	}
	
	private boolean checkXDeNovo(Map<Individual, Genotype> genomap, Individual pb, int allele)
	{
		//TODO: Write
		return false;
	}
	
	public void checkMammalYCandidate(Map<Individual, Genotype> genomap, Individual pb, Candidate c)
	{
		//TODO: Write
	}

}
