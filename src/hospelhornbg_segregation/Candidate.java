package hospelhornbg_segregation;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import hospelhornbg_bioinformatics.BreakendPair;
import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.Translocation;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.Gene;
import hospelhornbg_genomeBuild.GeneFunc;

public class Candidate implements Comparable<Candidate>{
	
	//private Inheritance eIPattern;
	private Map<Individual, Inheritance> iIPatternMap;

	//private Set<Integer> iAlleles;
	//private Map<Individual, Integer> iAlleleMap;
	private int iAllele;
	private boolean bDEL;
	private boolean bDUP;
	
	private Variant iVariant;
	//private Map<Candidate, Set<Integer>> iPartners;
	private Map<Individual, Set<Candidate>> iPartners;
	
	private Gene iGene;
	private GeneFunc iEffect;
	
	public Candidate(Variant variant, int allele)
	{
		iVariant = variant;
		//iVariants = null;
		//bMultiVar = false;
		//eIPattern = null;
		iIPatternMap = new HashMap<Individual, Inheritance>();
		iGene = null;
		//iHalfHetPhase = Inheritance.HALFHET_PHASE_NONE;
		iAllele = -1;
		//iAlleles = new LinkedList<Integer>();
		//iAlleles = new HashSet<Integer>();
		//iAlleleMap = new HashMap<Individual, Integer>();
		bDEL = false;
		bDUP = false;
		//iPartners = new HashMap<Candidate, Set<Integer>>();
		iPartners = new HashMap<Individual, Set<Candidate>>();
		iAllele = allele;
		iEffect = null;
	}

	public Inheritance getInheritancePattern(Individual indiv)
	{
		return iIPatternMap.get(indiv);
	}
	
	public int getAllele()
	{
		return iAllele;
	}
	
	public boolean isDEL_candidate()
	{
		return bDEL;
	}
	
	public boolean isDUP_candidate()
	{
		return bDUP;
	}
	
	public Gene getGene()
	{
		return iGene;
	}
	
	public GeneFunc getPositionEffect()
	{
		return iEffect;
	}
	
	public Variant getVariant()
	{
		return iVariant;
	}
	
	public List<Candidate> getAllPartners(Individual indiv)
	{
		Set<Candidate> pset = iPartners.get(indiv);
		int sz = 1;
		if (pset != null) sz += pset.size();
		List<Candidate> plist = new ArrayList<Candidate>(sz);
		if (pset != null) plist.addAll(pset);
		Collections.sort(plist);
		return plist;
	}
	
	public void setInheritancePattern(Individual indiv, Inheritance ip)
	{
		iIPatternMap.put(indiv, ip);
	}

	public void addAllele(int allele)
	{
		iAllele = allele;
	}
	
	public void flagDEL_candidate(boolean flag)
	{
		bDEL = flag;
	}
	
	public void flagDUP_candidate(boolean flag)
	{
		bDUP = flag;
	}
	
	public void setGene(Gene gene)
	{
		iGene = gene;
		//Determine position effect
		if (iGene == null)
		{
			iEffect = GeneFunc.INTERGENIC;
			return;
		}
		if (iVariant instanceof StructuralVariant)
		{
			if (iVariant instanceof BreakendPair)
			{
				BreakendPair bp = (BreakendPair)iVariant;
				StructuralVariant v1 = bp.getBNDVariant(false);
				StructuralVariant v2 = bp.getBNDVariant(true);
				Contig c1 = v1.getChromosome();
				Contig c2 = v2.getChromosome();
				Contig gc = iGene.getChromosome();
				GeneFunc e1 = GeneFunc.INTERGENIC;
				GeneFunc e2 = GeneFunc.INTERGENIC;
				if (c1.equals(gc))
				{
					e1 = iGene.getRelativeRegionLocationEffect(bp.getCIPosition(false, false, false), bp.getCIPosition(false, false, true));
				}
				if (c2.equals(gc))
				{
					e2 = iGene.getRelativeRegionLocationEffect(bp.getCIPosition(true, false, false), bp.getCIPosition(true, false, true));
				}
				if (e1.getPriority() > e2.getPriority()) iEffect = e2;
				else iEffect = e1;
			}
			else if (iVariant instanceof Translocation)
			{
				Translocation tra = (Translocation)iVariant;
				Contig c1 = tra.getChromosome1();
				Contig c2 = tra.getChromosome2();
				Contig gc = iGene.getChromosome();
				GeneFunc e1 = GeneFunc.INTERGENIC;
				GeneFunc e2 = GeneFunc.INTERGENIC;
				if (c1.equals(gc))
				{
					e1 = iGene.getRelativeRegionLocationEffect(tra.getCIPosition(false, false, false), tra.getCIPosition(false, false, true));
				}
				if (c2.equals(gc))
				{
					e2 = iGene.getRelativeRegionLocationEffect(tra.getCIPosition(true, false, false), tra.getCIPosition(true, false, true));
				}
				if (e1.getPriority() > e2.getPriority()) iEffect = e2;
				else iEffect = e1;
			}
			else
			{
				StructuralVariant sv = (StructuralVariant)iVariant;
				if(sv.getType() != SVType.INV) iEffect = iGene.getRelativeRegionLocationEffect(sv.getCIPosition(false, false, false), sv.getCIPosition(true, false, true));
				else
				{
					GeneFunc e1 = iGene.getRelativeRegionLocationEffect(sv.getCIPosition(false, false, false), sv.getCIPosition(false, false, true));
					GeneFunc e2 = iGene.getRelativeRegionLocationEffect(sv.getCIPosition(true, false, false), sv.getCIPosition(true, false, true));
					if (e1.getPriority() > e2.getPriority()) iEffect = e2;
					else iEffect = e1;
				}
			}
		}
		else iEffect = iGene.getRelativeLocationEffect(iVariant.getPosition());
	}
	
	public void setIntergenic()
	{
		iGene = null;
		iEffect = GeneFunc.INTERGENIC;
	}
	
	public void addPartner(Individual indiv, Candidate partner)
	{
		if (indiv == null) return;
		Set<Candidate> pset = iPartners.get(indiv);
		if (pset == null)
		{
			pset = new HashSet<Candidate>();
			pset.add(partner);
			iPartners.put(indiv, pset);
		}
		else pset.add(partner);
	}
	
	public void removePartner(Individual indiv, Candidate partner)
	{
		Set<Candidate> pset = iPartners.get(indiv);
		if (pset != null) pset.remove(partner);
	}
	
	public boolean hasPartners(Individual indiv)
	{
		Set<Candidate> pset = iPartners.get(indiv);
		if (pset != null) return pset.isEmpty();
		return false;
	}
	
	public boolean equals(Object o)
	{
		if (o == null) return false;
		if (o == this) return true;
		if (!(o instanceof Candidate)) return false;
		Candidate c = (Candidate)o;
		//if (this.getInheritancePattern() != c.getInheritancePattern()) return false;
		//if (this.isMultiVariant() != c.isMultiVariant()) return false;
		if (!(this.getGene().equals(c.getGene()))) return false;
		/*List<Variant> tvar = this.getVariants();
		List<Variant> cvar = c.getVariants();
		if (tvar == null && cvar == null) return true;
		if (tvar == null || cvar == null) return false; //Only one is null since already returned if both are
		int tsz = tvar.size();
		int csz = cvar.size();
		if (tsz != csz) return false;
		for (int i = 0; i < tsz; i++)
		{
			Variant tv = tvar.get(i);
			Variant cv = cvar.get(i);
			if (tv == null && cv == null) continue;
			if (tv == null || cv == null) return false;
			if (!tv.equals(cv)) return false;
		}*/
		if (this.iVariant == null && c.iVariant != null) return false;
		if (!(this.getVariant().equals(c.getVariant()))) return false;
		return true;
	}
	
	public int hashCode()
	{
		int vcode = 0x12345678;
		int gcode = 0x00ABCDEF;
		Variant v = this.getVariant();
		if (v != null) vcode = v.hashCode();
		Gene g = this.getGene();
		if (g != null) gcode = g.hashCode();
		return vcode ^ gcode;
	}

	@Override
	public int compareTo(Candidate o) 
	{
		if (o == null) return 1;
		if (o == this) return 0;
		//Simply [first] variant position.
		Variant tv = this.getVariant();
		Variant ov = o.getVariant();
		if (tv == null && ov == null) return 0;
		if (tv == null && ov != null) return -1;
		if (tv != null && ov == null) return 1;
		return tv.compareTo(ov);
	}
	
	private Map<Individual, Inheritance> copyIPMap()
	{
		Map<Individual, Inheritance> copy = new HashMap<Individual, Inheritance>();
		Set<Individual> keyset = iIPatternMap.keySet();
		for (Individual indiv : keyset)
		{
			copy.put(indiv, iIPatternMap.get(indiv));
		}
		return copy;
	}
	
	private Map<Individual, Set<Candidate>> copyPartnerMap()
	{
		Map<Individual, Set<Candidate>> copy = new HashMap<Individual, Set<Candidate>>();
		Set<Individual> keyset = iPartners.keySet();
		for (Individual indiv : keyset)
		{
			Set<Candidate> pset = iPartners.get(indiv);
			if (pset == null) continue;
			Set<Candidate> copyset = new HashSet<Candidate>();
			copyset.addAll(pset);
			copy.put(indiv, copyset);
		}
		return copy;
	}
	
	/**
	 * Get a list of new candidates otherwise identical to this candidate each
	 * linked to a provided gene.
	 * <br>!! IMPORTANT !! - This candidate's gene will NOT be changed!
	 * @param genes Genes to spawn new candidates linked to.
	 * @return Sorted list of new candidates, each one linked to a specified gene.
	 */
	public List<Candidate> addGenes(Collection<Gene> genes)
	{
		if (genes == null) return null;
		if (genes.isEmpty()) return null;
		List<Candidate> copies = new ArrayList<Candidate>(genes.size());
		for (Gene g : genes)
		{
			Candidate c = new Candidate(this.iVariant, this.iAllele);
			/*if (this.isMultiVariant()) c = new Candidate(this.iVariants);
			else c = new Candidate(this.iVariant);*/
			c.iIPatternMap = this.copyIPMap();
			c.setGene(g);
			c.flagDEL_candidate(bDEL);
			c.flagDUP_candidate(bDUP);
			//c.iAlleles.addAll(iAlleles);
			//c.iAllele = this.iAllele;
			c.iPartners = this.copyPartnerMap();
			copies.add(c);
		}
		Collections.sort(copies);
		return copies;
	}

	public boolean segregatesWithAny()
	{
		Collection<Inheritance> allip = iIPatternMap.values();
		for (Inheritance i : allip)
		{
			if (i != null && i != Inheritance.UNRESOLVED) return true;
		}
		return false;
	}
	
	public boolean isUnpairedHalfHet()
	{
		Set<Individual> aff = iIPatternMap.keySet();
		for (Individual i : aff)
		{
			Inheritance ip = iIPatternMap.get(i);
			if (i == null) continue;
			if (!Inheritance.isHalfhet(ip)) return false;
			//Check for partners
			Set<Candidate> hetpartners = iPartners.get(i);
			if (hetpartners == null) continue;
			if (hetpartners.isEmpty()) continue;
			return false;
		}
		return true;
	}

	public Collection<Variant> getAllPartnerVariants()
	{
		Set<Variant> vset = new HashSet<Variant>();
		Set<Individual> aff = iIPatternMap.keySet();
		for (Individual i : aff)
		{
			Set<Candidate> hetpartners = iPartners.get(i);
			if (hetpartners == null) continue;
			for (Candidate c : hetpartners)
			{
				vset.add(c.getVariant());
			}
		}
		return vset;
	}
	
	
}
