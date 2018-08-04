package hospelhornbg_segregation;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import hospelhornbg_bioinformatics.Genotype;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_genomeBuild.Gene;
import hospelhornbg_genomeBuild.GeneSet;

public class Inheritor {
	
	/* --- Check Parental Relationships --- */
	
	public static final int HERITABILITY_NONE = 0;
	public static final int HERITABILITY_STANDARD = 1;
	public static final int HERITABILITY_DELETION_ONLY = 2;
	public static final int HERITABILITY_DUPLICATION_ONLY = 3;
	public static final int HERITABILITY_UNKNOWN = -1;
	
	private static boolean arrayContains(int[] array, int val)
	{
		if (array == null || array.length < 1) return false;
		for (int i : array) if (i == val) return true;
		return false;
	}
	
	public static int checkPC(Genotype parent, Genotype child)
	{
		if (parent == null || child == null) return HERITABILITY_UNKNOWN;
		if (parent.isGenotypeUnknown() || child.isGenotypeUnknown()) return HERITABILITY_UNKNOWN;
		//Do parent and child share at least one allele?
		boolean shared = false;
		int[] c_all = child.getAlleles();
		int[] p_all = parent.getAlleles();
		for (int a : c_all)
		{
			for (int pa : p_all)
			{
				if (pa == a) shared = true;
				break;
			}
			if (shared) return HERITABILITY_STANDARD;
		}
		
		//Check for CNVs
		int ccnv = child.getCopyNumber();
		int pcnv = parent.getCopyNumber();
		if (ccnv < 0 || ccnv >= 2) return HERITABILITY_NONE; //PC Error
		if (ccnv <= 1)
		{
			if (pcnv < 2) return HERITABILITY_DELETION_ONLY;
		}
		
		return HERITABILITY_NONE;
	}
	
	public static boolean checkPPC(Genotype parent1, Genotype parent2, Genotype child)
	{
		if (child == null || parent1 == null || parent2 == null) return false;
		int p1_pc = checkPC(parent1, child);
		int p2_pc = checkPC(parent2, child);
		if (p1_pc == 0 || p2_pc == 0) return false;
		if (p1_pc < 0 || p2_pc < 0) return true; //One or both unknown...
		
		if (p1_pc == HERITABILITY_STANDARD && p2_pc == HERITABILITY_STANDARD)
		{
			int[] c_all = child.getAlleles();
			int[] p1_all = parent1.getAlleles();
			int[] p2_all = parent2.getAlleles();
			
			//Are there any alleles neither parent has?
			for (int a : c_all)
			{
				boolean p1has = arrayContains(p1_all, a);
				boolean p2has = arrayContains(p2_all, a);
				if (!p1has && !p2has) return false;
			}
			
			//Can each parent have contributed an allele?
			Set<Integer> p1set = new HashSet<Integer>();
			for (int a : p1_all) p1set.add(a);
			for (int a : p1set)
			{
				//Find a match in c_all and make a version of c_all without it
				int[] ctemp = new int[c_all.length - 1];
				int ind = 0;
				for (int ca : c_all)
				{
					if (ca == a) break;
					ind++;
				}
				int j = 0;
				for (int i = 0; i < ctemp.length; i++)
				{
					if (j != ind) ctemp[i] = c_all[j];
					else
					{
						j++;
						ctemp[i] = c_all[j];
					}
					j++;
				}
				
				//See if there is any overlap between the p2 alleles and ctemp
				for (int pa : p2_all)
				{
					if (arrayContains(ctemp, pa)) return true;
				}
				
			}
			
		}
		else if (p1_pc == HERITABILITY_DELETION_ONLY && p2_pc == HERITABILITY_STANDARD)
		{
			//Make sure there's no alleles that came from neither parent...
			int[] c_all = child.getAlleles();
			int[] p1_all = parent1.getAlleles();
			int[] p2_all = parent2.getAlleles();
			
			for (int a : c_all)
			{
				boolean p1has = arrayContains(p1_all, a);
				boolean p2has = arrayContains(p2_all, a);
				if (!p1has && !p2has) return false;
			}
		}
		else if (p2_pc == HERITABILITY_DELETION_ONLY && p1_pc == HERITABILITY_STANDARD)
		{
			//Make sure there's no alleles that came from neither parent...
			int[] c_all = child.getAlleles();
			int[] p1_all = parent1.getAlleles();
			int[] p2_all = parent2.getAlleles();
			
			for (int a : c_all)
			{
				boolean p1has = arrayContains(p1_all, a);
				boolean p2has = arrayContains(p2_all, a);
				if (!p1has && !p2has) return false;
			}
		}
		else if (p1_pc == HERITABILITY_DELETION_ONLY && p2_pc == HERITABILITY_DELETION_ONLY)
		{
			//Make sure there's no alleles!
			return (child.getCopyNumber() == 0 || child.isGenotypeUnknown());
		}
		
		return false;
	}
	
	/* --- General Inheritance Checks --- */
	
	public static boolean checkCompHetCandidacy(Genotype parent1, Genotype parent2, Genotype child, int allele)
	{
		if (child == null) return false;
		if (parent1 == null && parent2 == null) return false;
		if (parent1 == null || parent2 == null) return true;
		
		boolean cHas = arrayContains(child.getAlleles(), allele);
		if (!cHas) return false;
		boolean p1Has = arrayContains(parent1.getAlleles(), allele);
		boolean p2Has = arrayContains(parent2.getAlleles(), allele);
		if (p1Has && p2Has) return false;
		
		return (p1Has || p2Has);
	}

	public static Inheritance checkHomozygousCandidate(Individual target, Genotype tgeno, Individual other, Genotype ogeno, int allele)
	{
		//NULL means inconclusive
		if (target == null) return null;
		if (other == null) return null;
		if (!target.isAffected()) return null;
		if (tgeno.isGenotypeUnknown()) return null;
		if (ogeno.isGenotypeUnknown()) return null;
		
		int[] t_all = tgeno.getAlleles();
		if (t_all == null || t_all.length < 1) return null;
		//int tallele = tgeno.getAlleles()[0];
		if (!arrayContains(t_all, allele)) return null;
		
		int[] o_all = ogeno.getAlleles();
		if (o_all == null || o_all.length < 1) return null;
		//if (!arrayContains(o_all, tallele)) return null;
		if (!arrayContains(o_all, allele)) return null;
		
		//Now for the actual analysis
		if (other.isAffected()) return Inheritance.DOMINANT;
		else
		{
			if(ogeno.isHomozygous())
			{
				//Check for unusual CN
				int tcnv = t_all.length;
				int ocnv = o_all.length;
				if (tcnv != 2 && tcnv != ocnv) return Inheritance.HALF_HET_SV;
			}
			else return Inheritance.HOM_REC;
		}
		
		return Inheritance.UNRESOLVED; //DNS
	}

	public static void checkHeterozygousCandidate(Map<Individual, Genotype> genomap, Candidate c)
	{
		if (genomap == null || genomap.isEmpty()) return;
		Set<Individual> allindivs = genomap.keySet();
		List<Individual> unaffected = new LinkedList<Individual>();
		List<Individual> affected = new LinkedList<Individual>();
		Set<Integer> allelePool = new HashSet<Integer>();
		for (Individual indiv : allindivs)
		{
			if (indiv.isAffected()) affected.add(indiv);
			else unaffected.add(indiv);
			int[] alleles = genomap.get(indiv).getAlleles();
			if (alleles != null && alleles.length > 0)
			{
				for (int a : alleles) allelePool.add(a);
			}
		}
		
		//Eliminate all alleles that are seen homozygous in unaffected indivs
		for (Individual indiv : unaffected)
		{
			Genotype g = genomap.get(indiv);
			if (g == null) continue;
			if (g.isGenotypeUnknown()) continue;
			if (g.isHomozygous())
			{
				//Get the allele
				int[] i_all = g.getAlleles();
				if (i_all != null && i_all.length > 0)
				{
					allelePool.remove(i_all[0]);
				}
			}
		}
		if(allelePool.isEmpty()){
			c.setInheritancePattern(Inheritance.UNRESOLVED);
			return;
		}
		
		//More than 1 affected?
		if (affected.size() < 2){
			c.setInheritancePattern(Inheritance.HALF_HET);
			//Get the allele(s) in question...
			for (int a : allelePool) c.addAllele(a);
			return;
		}
		
		//Are there any alleles that are only in unaffected individuals?
		//If so, remove
		Set<Integer> removelist = new HashSet<Integer>();
		for (int a : allelePool)
		{
			for (Individual indiv : unaffected)
			{
				Genotype g = genomap.get(indiv);
				if (g == null) continue;
				if (g.isGenotypeUnknown()) continue;
				int[] i_all = g.getAlleles();
				if (arrayContains(i_all, a)) removelist.add(a);
			}
		}
		allelePool.removeAll(removelist);
		
		//Any left?
		if (allelePool.isEmpty()){
			c.setInheritancePattern(Inheritance.HALF_HET);
			//Get the allele(s) in question...
			for (int a : allelePool) c.addAllele(a);
			return;
		}
		
		//Is there a remaining allele that's in all affecteds?
		removelist.clear();
		for (int a : allelePool)
		{
			for (Individual indiv : affected)
			{
				Genotype g = genomap.get(indiv);
				if (g == null) continue;
				if (g.isGenotypeUnknown()) continue;
				int[] i_all = g.getAlleles();
				if (!arrayContains(i_all, a)) removelist.add(a);
			}
		}
		allelePool.removeAll(removelist);
		
		if (allelePool.isEmpty()){
			//Nothing left in allele pool, so need to retrieve potential halfhets
			c.setInheritancePattern(Inheritance.HALF_HET);
			for (int a : removelist) c.addAllele(a);
			return;
		}
		
		c.setInheritancePattern(Inheritance.DOMINANT);
		//Get the allele(s) in question...
		for (int a : allelePool) c.addAllele(a);
	}

	public static void checkHomozygousCandidate(Map<Individual, Genotype> genomap, Candidate c)
	{
		if (genomap == null || genomap.isEmpty()) return;
		if (c == null) return;
		List<Integer> alist = c.getAlleles();
		if (alist.isEmpty()) return;
		int allele = alist.get(0);
		Set<Individual> allindivs = genomap.keySet();
		List<Individual> unaffected = new LinkedList<Individual>();
		List<Individual> affected = new LinkedList<Individual>();
		
		for (Individual indiv : allindivs)
		{
			if (indiv.isAffected()) affected.add(indiv);
			else unaffected.add(indiv);
		}
		
		int a_ind = 0;
		Inheritance ih = null;
		int naff = affected.size();
		for (Individual aff : affected)
		{
			a_ind++;
			Genotype a_geno = genomap.get(aff);
			Inheritance ih_aff = null;
			for (Individual unaff : unaffected)
			{
				Genotype u_geno = genomap.get(unaff);
				Inheritance ih_comp = checkHomozygousCandidate(aff, a_geno, unaff, u_geno, allele);
				if (ih_comp == Inheritance.UNRESOLVED){
					c.setInheritancePattern(Inheritance.UNRESOLVED);
					return;
				}
				if (ih_comp == null) continue;
				//comphetSV > dom > rec
				if (ih_comp == Inheritance.HALF_HET_SV) ih_aff =  Inheritance.HALF_HET_SV;
				else if (ih_comp == Inheritance.DOMINANT)
				{
					if (ih_aff == Inheritance.HOM_REC || ih_aff == null) ih_aff = Inheritance.DOMINANT;
				}
				else if (ih_comp == Inheritance.HOM_REC)
				{
					if (ih_aff == null) ih_aff = Inheritance.HOM_REC;
				}
			}
			if (a_ind < naff)
			{
				for (int i = a_ind; i < naff; i++)
				{
					Individual oaff = affected.get(a_ind);
					Genotype o_geno = genomap.get(oaff);
					Inheritance ih_comp = checkHomozygousCandidate(aff, a_geno, oaff, o_geno, allele);
					if (ih_comp == Inheritance.UNRESOLVED){
						c.setInheritancePattern(Inheritance.UNRESOLVED);
						return;
					}
					if (ih_comp == null) continue;
					//comphetSV > dom > rec
					if (ih_comp == Inheritance.HALF_HET_SV) ih_aff =  Inheritance.HALF_HET_SV;
					else if (ih_comp == Inheritance.DOMINANT)
					{
						if (ih_aff == Inheritance.HOM_REC || ih_aff == null) ih_aff = Inheritance.DOMINANT;
					}
					else if (ih_comp == Inheritance.HOM_REC)
					{
						if (ih_aff == null) ih_aff = Inheritance.HOM_REC;
					}
				}
			}
			if (ih_aff == Inheritance.UNRESOLVED){
				c.setInheritancePattern(Inheritance.UNRESOLVED);
				return;
			}
			if (ih_aff == null) continue;
			//comphetSV > dom > rec
			if (ih_aff == Inheritance.HALF_HET_SV) ih =  Inheritance.HALF_HET_SV;
			else if (ih_aff == Inheritance.DOMINANT)
			{
				if (ih == Inheritance.HOM_REC || ih == null) ih = Inheritance.DOMINANT;
			}
			else if (ih_aff == Inheritance.HOM_REC)
			{
				if (ih == null) ih = Inheritance.HOM_REC;
			}
			
		}
		
		if (ih == Inheritance.HALF_HET_SV){
			//Get CNs of all affecteds, compare to CNs of all unaffecteds...
			//Note suspicious CNs
			//Check for deletions...
			boolean delpass = true;
			for (Individual a : affected)
			{
				Genotype ag = genomap.get(a);
				if (ag == null) continue;
				if (ag.getCopyNumber() >= 2){
					delpass = false;
					break;
				}
			}
			if (delpass)
			{
				for (Individual u : unaffected)
				{
					Genotype ug = genomap.get(u);
					if (ug == null) continue;
					if (ug.getCopyNumber() < 2)
					{
						delpass = false;
						break;
					}
				}
			}
			if (delpass) c.flagDEL_candidate(true);
			//Check for duplications...
			boolean duppass = true;
			for (Individual a : affected)
			{
				Genotype ag = genomap.get(a);
				if (ag == null) continue;
				if (ag.getCopyNumber() <= 2){
					duppass = false;
					break;
				}
			}
			if (duppass)
			{
				for (Individual u : unaffected)
				{
					Genotype ug = genomap.get(u);
					if (ug == null) continue;
					if (ug.getCopyNumber() > 2)
					{
						duppass = false;
						break;
					}
				}
			}
			if (duppass) c.flagDUP_candidate(true);
		}
		
		c.setInheritancePattern(ih);
	}
	
	/* --- Segregation --- */
	
	public static List<Gene> getVariantGenes(Variant v, GeneSet genes)
	{
		List<Gene> glist = genes.annotateVariant(v);
		return glist;
	}
	
	private static boolean isInheritancePairable(Inheritance ip)
	{
		switch(ip)
		{
		case COMP_HET:
			break;
		case COMP_HET_DEL:
			break;
		case COMP_HET_SV:
			break;
		case DENOVO_DOM:
			break;
		case DENOVO_HET:
			return true;
		case DENOVO_HET_SV:
			return true;
		case DENOVO_REC:
			break;
		case DOMINANT:
			break;
		case HALF_HET:
			return true;
		case HALF_HET_SV:
			return true;
		case HOM_REC:
			break;
		case MVIOL:
			break;
		case UNRESOLVED:
			break;
		default:
			break;
		
		}
		return false;
	}
	
	public static List<Candidate> getCandidates(VariantPool varpool, Pedigree family, GeneSet genes)
	{
		//Checks
		if (varpool == null) return null;
		if (family == null) return null;
		if (genes == null) return null;
		
		//Prep
		List<Variant> variants = varpool.getVariants();
		Map<Gene, List<Candidate>> genemap = new HashMap<Gene, List<Candidate>>();
		List<Candidate> intergenics = new LinkedList<Candidate>();
		
		//Mark initial inheritance, gene annotate, generate initial candidate set
		for (Variant v : variants)
		{
			Candidate cinit = family.getCandidate(v);
			//Match to genes...
			List<Gene> glist = getVariantGenes(v, genes);
			List<Candidate> split_cand = cinit.addGenes(glist);
			if (split_cand == null || split_cand.isEmpty()) intergenics.add(cinit);
			else
			{
				for (Candidate c : split_cand)
				{
					Gene cgene = c.getGene();
					List<Candidate> clist = genemap.get(cgene);
					if (clist == null)
					{
						clist = new LinkedList<Candidate>();
						clist.add(c);
						genemap.put(cgene, clist);
					}
					else clist.add(c);
				}
			}
		}
		
		//Attempt to pair half-hets
		Set<Gene> allgenes = genemap.keySet();
		for (Gene g : allgenes)
		{
			List<Candidate> clist = genemap.get(g);
			if (clist == null) continue;
			if (clist.isEmpty()) continue;
			int num = clist.size();
			if (num < 2) continue;
			//Scan for comphets and unroll...
			List<Candidate> newlist = new ArrayList<Candidate>(num*2);
			for (Candidate c : clist)
			{
				if (isInheritancePairable(c.getInheritancePattern()))
				{
					List<Candidate> unrolled = c.unroll();
					if (unrolled != null) newlist.addAll(unrolled);
					else newlist.add(c);
				}
				else newlist.add(c);
			}
			num = newlist.size();
			genemap.put(g, newlist);
			//List<Candidate> finallist = new LinkedList<Candidate>();
			
			for (int i = 0; i < num - 1; i++)
			{
				Candidate c1 = newlist.get(i);
				if (isInheritancePairable(c1.getInheritancePattern()))
				{
					for (int j = i+1; j < num; j++)
					{
						Candidate c2 = newlist.get(j);
						if (isInheritancePairable(c2.getInheritancePattern()))
						{
							//Try to pair!
							boolean pairs = family.passCompHetPair(c1, c2);
							if (pairs)
							{
								c1.addPartner(c2, c2.getAlleles());
								c2.addPartner(c1, c1.getAlleles());
							}
						}
					}	
				}
			}
		}
		
		//Compile final massive list...
		List<Candidate> biglist = new LinkedList<Candidate>();
		biglist.addAll(intergenics);
		for (Gene g : allgenes)
		{
			List<Candidate> clist = genemap.get(g);
			biglist.addAll(clist);
		}
		Collections.sort(biglist);
		
		return biglist;
	}

}
