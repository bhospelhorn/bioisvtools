package hospelhornbg_segregation;

import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import hospelhornbg_bioinformatics.Genotype;

public class SexChromInheritor {
	
	public static void checkMammalXCandidate(Map<Individual, Genotype> genomap, Individual pb, Candidate c)
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
					if (!allhave)
					{
						c.setInheritancePattern(pb, Inheritance.UNRESOLVED);
					}
					else
					{
						if(allhom)
						{
							//Call homrec
							if (checkXDeNovo(genomap, pb, allele)) c.setInheritancePattern(pb, Inheritance.X_LINKED_MV);
							else c.setInheritancePattern(pb, Inheritance.X_LINKED_REC);
						}
						else
						{
							//Call dom
							if (checkXDeNovo(genomap, pb, allele)) c.setInheritancePattern(pb, Inheritance.X_LINKED_DN);
							else c.setInheritancePattern(pb, Inheritance.X_LINKED_DOM);
						}
					}
				}
			}
		}
		
	}
	
	private static boolean patXRescue(Map<Individual, Genotype> genomap, List<Individual> aff, Candidate c)
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
	
	private static boolean patXRescueDom(Map<Individual, Genotype> genomap, List<Individual> aff, Candidate c)
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
	
	private static boolean checkXDeNovo(Map<Individual, Genotype> genomap, Individual pb, int allele)
	{
		//Get parents and genotypes
		Individual mother = pb.getMother();
		Individual father = pb.getFather();
		
		Genotype pbgeno = genomap.get(pb);
		Genotype mgeno = null;
		Genotype pgeno = null;
		
		if (mother != null) mgeno = genomap.get(mother);
		if (father != null) pgeno = genomap.get(father);
		
		int pbcn = pb.getExpectedXCount();
		switch(pbcn)
		{
		case 0:
			//Should not happen
			if (mgeno == null) return false;
			return (mgeno.getAlleleCallCount() > 0);
		case 1:
			if (mgeno == null) return false;
			return (!mgeno.hasAllele(allele));
		case 2:
			int acount = pbgeno.countAlleleOccurrences(allele);
			if (acount < 2)
			{
				//Only need one parent
				if (mgeno == null || pgeno == null) return false;
				if (mgeno.hasAllele(allele)) return false;
				if (pgeno.hasAllele(allele)) return false;
				return true;
			}
			else
			{
				//Both parents must have
				boolean mhas = true;
				boolean fhas = true;
				if (mgeno != null) mhas = mgeno.hasAllele(allele);
				if (pgeno != null) fhas = pgeno.hasAllele(allele);
				return !(mhas && fhas);
			}
		default:
			break;
		}
		
		return false;
	}
	
	public static void checkMammalYCandidate(Map<Individual, Genotype> genomap, Individual pb, Candidate c)
	{
		//Get expected Y count
		int cn = pb.getExpectedYCount();
		if (cn == 0)
		{
			c.setInheritancePattern(pb, Inheritance.UNRESOLVED);
			return;
		}
		
		int allele = c.getAllele();
		List<Individual> unaff = new LinkedList<Individual>();
		List<Individual> aff = new LinkedList<Individual>();
		for(Individual i : genomap.keySet())
		{
			if(i.isAffected()) aff.add(i);
			else unaff.add(i);
		}
		if (cn == 1)
		{
			//Do any unaff have at all?
			boolean anyunaff = false;
			for(Individual u : unaff)
			{
				if (u.getExpectedYCount() > 0)
				{
					Genotype g = genomap.get(u);
					if (g.hasAllele(allele))
					{
						anyunaff = true;
						break;
					}
				}
			}
			
			if(anyunaff)
			{
				//At least one unaff has this allele
				//Are any unaff CN2+?
				boolean u2plus = false;
				for(Individual u : unaff)
				{
					if (u.getExpectedYCount() >= 2)
					{
						u2plus = true;
						break;
					}
				}
				
				if (u2plus)
				{
					//At least one unaff has 2 or more Y chroms
					//Are any of these 2+ Y unaffs homozygous for this allele?
					for(Individual u : unaff)
					{
						if (u.getExpectedYCount() >= 2)
						{
							Genotype g = genomap.get(u);
							if(g.hasAllele(allele) && g.isHomozygous())
							{
								c.setInheritancePattern(pb, Inheritance.UNRESOLVED);
								return;
							}
						}
					}
					c.setInheritancePattern(pb, Inheritance.Y_LINKED_SPECIAL);
					return;
				}
				else
				{
					//DNS
					c.setInheritancePattern(pb, Inheritance.UNRESOLVED);
					return;
				}
				
			}
			else
			{
				//No unaff have this allele
				//Do all affected have this allele?
				for(Individual a : aff)
				{
					if (a.getExpectedYCount() < 1)
					{
						c.setInheritancePattern(pb, Inheritance.UNRESOLVED);
						return;
					}
					Genotype g = genomap.get(a);
					if(!g.hasAllele(allele))
					{
						c.setInheritancePattern(pb, Inheritance.UNRESOLVED);
						return;
					}
				}
				//Check if de novo
				if(checkYDeNovo(genomap, pb, allele))
				{
					c.setInheritancePattern(pb, Inheritance.Y_LINKED_DN);
					return;
				}
				else
				{
					c.setInheritancePattern(pb, Inheritance.Y_LINKED_DOM);
					return;
				}
			}
			
		}
		else
		{
			//cn >= 2
			c.setInheritancePattern(pb, Inheritance.Y_LINKED_SPECIAL);
			return;
		}
	}
	
	private static boolean checkYDeNovo(Map<Individual, Genotype> genomap, Individual pb, int allele)
	{
		if (pb.getExpectedYCount() < 1) return false;
		Individual father = pb.getFather();
		
		//Genotype pbgeno = genomap.get(pb);
		Genotype pgeno = null;
		if (father != null) pgeno = genomap.get(father);
		else return false;
		
		//Does dad have allele?
		return !pgeno.hasAllele(allele);
	}

}
