package hospelhornbg_segregation;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import hospelhornbg_bioinformatics.AffectedStatus;
import hospelhornbg_bioinformatics.Genotype;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_genomeBuild.Gene;
import hospelhornbg_genomeBuild.GeneSet;

//TODO: Need 2 fixes:
	// Check comphet candidate pairs against other affecteds
	// XY chrom should behave differently!
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

	public static void checkHeterozygousCandidate(Map<Individual, Genotype> genomap, Individual proband, Candidate c)
	{
		if (genomap == null || genomap.isEmpty()) return;
		if (c == null) return;
		int allele = c.getAllele();
		Set<Individual> allindivs = genomap.keySet();
		List<Individual> unaffected = new LinkedList<Individual>();
		//List<Individual> affected = new LinkedList<Individual>();
		
		//Load lists
		for (Individual indiv : allindivs)
		{
			if (indiv.getAffectedStatus() == AffectedStatus.UNAFFECTED) unaffected.add(indiv);
			//else if (indiv.getAffectedStatus() == AffectedStatus.AFFECTED) affected.add(indiv);
		}
		
		//1. Are any unaffected hom for this allele?
		for (Individual indiv : unaffected)
		{
			Genotype g = genomap.get(indiv);
			if (g == null) continue;
			if (g.isHomozygous())
			{
				if (g.hasAllele(allele))
				{
					if (cnvhh_rescue(genomap, proband)) c.setInheritancePattern(proband, Inheritance.HALF_HET_SV);
					else c.setInheritancePattern(proband, Inheritance.UNRESOLVED);
					return;
				}
			}
		}

		//2. Do any unaffected have that allele at all?
		for (Individual indiv : unaffected)
		{
			Genotype g = genomap.get(indiv);
			if (g == null) continue;
			if (g.hasAllele(allele))
			{
				if (halfHetCheck(c, proband))
				{
					c.setInheritancePattern(proband, Inheritance.HALF_HET);
					return;
				}
				if (cnvhh_rescue(genomap, proband)) c.setInheritancePattern(proband, Inheritance.HALF_HET_SV);
				else c.setInheritancePattern(proband, Inheritance.UNRESOLVED);
				return;
			}
		}
		

		//3. Any other affected have allele? Any het?
		c.setInheritancePattern(proband, Inheritance.DOMINANT);
	}

	public static void checkHomozygousCandidate(Map<Individual, Genotype> genomap, Individual proband, Candidate c)
	{
		if (genomap == null || genomap.isEmpty()) return;
		if (c == null) return;
		int allele = c.getAllele();
		Set<Individual> allindivs = genomap.keySet();
		List<Individual> unaffected = new LinkedList<Individual>();
		//List<Individual> affected = new LinkedList<Individual>();
		
		//Load lists
		for (Individual indiv : allindivs)
		{
			if (indiv.getAffectedStatus() == AffectedStatus.UNAFFECTED) unaffected.add(indiv);
			//else if (indiv.getAffectedStatus() == AffectedStatus.AFFECTED) affected.add(indiv);
		}
		
		//1. Are any unaffected hom for this allele?
		for (Individual indiv : unaffected)
		{
			Genotype g = genomap.get(indiv);
			if (g == null) continue;
			if (g.isHomozygous())
			{
				if (g.hasAllele(allele))
				{
					if (cnvhh_rescue(genomap, proband)) c.setInheritancePattern(proband, Inheritance.HALF_HET_SV);
					else c.setInheritancePattern(proband, Inheritance.UNRESOLVED);
					return;
				}
			}
		}
		
		//2. Do any unaffected have that allele at all?
		for (Individual indiv : unaffected)
		{
			Genotype g = genomap.get(indiv);
			if (g == null) continue;
			if (g.hasAllele(allele))
			{
				c.setInheritancePattern(proband, Inheritance.HOM_REC);
				return;
			}
		}
		
		//3. Any other affected have allele? Any het?
		c.setInheritancePattern(proband, Inheritance.DOMINANT);
		
	}
	
	private static boolean cnvhh_rescue(Map<Individual, Genotype> genomap, Individual proband)
	{
		Set<Individual> allindivs = genomap.keySet();
		List<Individual> unaffected = new LinkedList<Individual>();
		
		for (Individual indiv : allindivs)
		{
			if (indiv.getAffectedStatus() == AffectedStatus.UNAFFECTED) unaffected.add(indiv);
			//else if (indiv.getAffectedStatus() == AffectedStatus.AFFECTED) affected.add(indiv);
		}
		
		int pbcnv = genomap.get(proband).getCopyNumber();
		
		for (Individual indiv : unaffected)
		{
			Genotype g = genomap.get(indiv);
			int cnv = g.getCopyNumber();
			if (pbcnv == cnv) return false;
		}
		
		return true;
	}
	
	public static boolean isDenovo(Map<Individual, Genotype> genomap, Individual proband, int allele, boolean hom)
	{
		if (proband == null) return false;
		Individual p1 = proband.getParent1();
		Individual p2 = proband.getParent2();
		
		if (hom)
		{
			boolean p1has = true;
			boolean p2has = true;
			
			//Does parent1 have the allele?
			if (p1 != null)
			{
				Genotype g1 = genomap.get(p1);
				if (!g1.hasAllele(allele)) p1has = false;
			}

			//Does parent2 have the allele?
			if (p2 != null)
			{
				Genotype g2 = genomap.get(p2);
				if (!g2.hasAllele(allele)) p2has = false;
			}
		
			return !(p1has && p2has);
		}
		else
		{
			//Does parent1 have the allele?
			if (p1 != null)
			{
				Genotype g1 = genomap.get(p1);
				if (g1.hasAllele(allele)) return false;
			}

			//Does parent2 have the allele?
			if (p2 != null)
			{
				Genotype g2 = genomap.get(p2);
				if (g2.hasAllele(allele)) return false;
			}
		
			return true;
		}
		
	}
	
	public static void adjustInheritance(Map<Individual, Genotype> genomap, Individual proband, Candidate c)
	{
		if (genomap == null) return;
		if (proband == null) return;
		if (c == null) return;
		
		Inheritance ip = c.getInheritancePattern(proband);
		switch(ip)
		{
		case DOMINANT:
			boolean dn1 = isDenovo(genomap, proband, c.getAllele(), false);
			if(!dn1) return;
			c.setInheritancePattern(proband, Inheritance.DENOVO_DOM);
			break;
		case HALF_HET:
			boolean dn2 = isDenovo(genomap, proband, c.getAllele(), false);
			if(!dn2) return;
			c.setInheritancePattern(proband, Inheritance.DENOVO_HET);
			break;
		case HALF_HET_SV:
			boolean dn3 = isDenovo(genomap, proband, c.getAllele(), false);
			if(!dn3) return;
			c.setInheritancePattern(proband, Inheritance.DENOVO_HET_SV);
			break;
		case HOM_REC:
			boolean dn4 = isDenovo(genomap, proband, c.getAllele(), true);
			if(!dn4) return;
			c.setInheritancePattern(proband, Inheritance.DENOVO_REC);
			break;
		default:
			break;
		}
	}
	
	public static boolean halfHetCheck(Candidate c, Individual affected)
	{
		if (affected == null) return false;
		Individual p1 = affected.getParent1();
		Individual p2 = affected.getParent2();
		
		if (p1 == null) return true;
		if (p2 == null) return true;
		
		Variant v = c.getVariant();
		Genotype g1 = v.getSampleGenotype(p1.getName());
		Genotype g2 = v.getSampleGenotype(p2.getName());
		
		int a = c.getAllele();
		
		if (g1.hasAllele(a) && g2.hasAllele(a)) return false;
		
		return true;
	}
	
	public static boolean hetPairCheck(Candidate c1, Candidate c2, Collection<Individual> unaffected)
	{
		if (c1 == null) return true;
		if (c2 == null) return true;
		if (unaffected == null) return true;
		
		Variant v1 = c1.getVariant();
		Variant v2 = c2.getVariant();
		
		//Make sure no unaffected individuals have both
		for (Individual u : unaffected)
		{
			Genotype g1 = v1.getSampleGenotype(u.getName());
			Genotype g2 = v2.getSampleGenotype(u.getName());
			if (g1 == null || g2 == null) continue;
			boolean h1 = g1.hasAllele(c1.getAllele());
			boolean h2 = g2.hasAllele(c2.getAllele());
			if (h1 && h2) return false;
		}
		
		return true;
	}
	
	/* --- Segregation --- */
	
	public static List<Gene> getVariantGenes(Variant v, GeneSet genes)
	{
		List<Gene> glist = genes.annotateVariant(v);
		return glist;
	}
	
	private static boolean isInheritancePairable(Inheritance ip)
	{
		if (ip == null) return false;
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
			List<Candidate> cinitl = family.toCandidates(v);
			//Match to genes...
			List<Gene> glist = getVariantGenes(v, genes);
			for (Candidate cinit : cinitl)
			{
				List<Candidate> split_cand = cinit.addGenes(glist);
				if (split_cand == null || split_cand.isEmpty()) {
					intergenics.add(cinit);
					cinit.setIntergenic();
				}
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
		}
		
		//Get affected indiv list
		List<Individual> affected = family.getAllAffected();
		List<Individual> unaffected = family.getAllUnaffected();
		
		//Attempt to pair half-hets
		//Have to have both candidates in an individual
		Set<Gene> allgenes = genemap.keySet();
		for (Gene g : allgenes)
		{
			List<Candidate> llist = genemap.get(g);
			if (llist == null) continue;
			List<Candidate> clist = new ArrayList<Candidate>(llist.size() + 1);
			clist.addAll(llist);
			if (clist.isEmpty()) continue;
			int num = clist.size();
			if (num < 2) continue;
			
			//newlist is obsolete now - no unrolling needed
			for (int i = 0; i < num - 1; i++)
			{
				Candidate c1 = clist.get(i);
				for (Individual aff : affected)
				{
					Inheritance ip = c1.getInheritancePattern(aff);
					if (ip == null) continue;
					if (isInheritancePairable(ip))
					{
						for (int j = i+1; j < num; j++)
						{
							Candidate c2 = clist.get(j);
							if (isInheritancePairable(c2.getInheritancePattern(aff)))
							{
								//Try to pair!
								boolean pairs = family.checkPair(c1, c2);
								boolean stillpairs = Inheritor.hetPairCheck(c1, c2, unaffected);
								if (pairs && stillpairs)
								{
									c1.addPartner(aff, c2);
									c2.addPartner(aff, c1);
								}
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
