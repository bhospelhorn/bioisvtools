package hospelhornbg_segregation;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_genomeBuild.Gene;

public class Candidate implements Comparable<Candidate>{
	
	private Inheritance eIPattern;
	//private int iHalfHetPhase;
	//private int iAllele; //Index number for indexed, CN for CNVs, -1 for unknown
	private Set<Integer> iAlleles;
	private boolean bDEL;
	private boolean bDUP;
	
	//private boolean bMultiVar;
	private Variant iVariant;
	private Map<Candidate, Set<Integer>> iPartners;
	//private List<Variant> iVariants;
	
	private Gene iGene;
	
	public Candidate(Variant variant)
	{
		iVariant = variant;
		//iVariants = null;
		//bMultiVar = false;
		eIPattern = null;
		iGene = null;
		//iHalfHetPhase = Inheritance.HALFHET_PHASE_NONE;
		//iAllele = -1;
		//iAlleles = new LinkedList<Integer>();
		iAlleles = new HashSet<Integer>();
		bDEL = false;
		bDUP = false;
		iPartners = new HashMap<Candidate, Set<Integer>>();
	}

	public Inheritance getInheritancePattern()
	{
		return eIPattern;
	}
	
	public int getAllele()
	{
		if (iAlleles == null) return -1;
		if (iAlleles.isEmpty()) return -1;
		for (int i : iAlleles) return i;
		return -1;
	}
	
	public List<Integer> getAlleles()
	{
		//return iAllele;
		if (iAlleles.isEmpty()) return null;
		List<Integer> copy = new ArrayList<Integer>(iAlleles.size());
		copy.addAll(iAlleles);
		return copy;
	}
	
	public boolean isMonoallelic()
	{
		if (iAlleles == null) return false;
		if (iAlleles.size() == 1) return true;
		return false;
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
	
	public Variant getVariant()
	{
		return iVariant;
	}
	
	public List<Candidate> getAllPartners()
	{
		if (iPartners.isEmpty()) return null;
		List<Candidate> plist = new ArrayList<Candidate>(iPartners.size());
		plist.addAll(iPartners.keySet());
		Collections.sort(plist);
		return plist;
	}

	public Collection<Integer> getPartnerCompatibleAlleles(Candidate partner)
	{
		if (partner == null) return null;
		Set<Integer> pall = iPartners.get(partner);
		if (pall == null) return null;
		if (pall.isEmpty()) return null;
		Set<Integer> copy = new HashSet<Integer>();
		copy.addAll(pall);
		return copy;
	}
	
	public void setInheritancePattern(Inheritance ip)
	{
		eIPattern = ip;
	}

	public void addAllele(int allele)
	{
		iAlleles.add(allele);
	}

	public void clearAlleles()
	{
		iAlleles.clear();
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
	}
	
	public void addPartner(Candidate partner, Collection<Integer> goodalleles)
	{
		Set<Integer> alleles = new HashSet<Integer>();
		if (goodalleles != null) alleles.addAll(goodalleles);
		iPartners.put(partner, alleles);
	}
	
	public void addCompatibleAlleleToPartner(Candidate partner, int allele)
	{
		if (partner == null) return;
		Set<Integer> alleles = iPartners.get(partner);
		if (alleles == null)
		{
			alleles = new HashSet<Integer>();
			alleles.add(allele);
			iPartners.put(partner, alleles);
			return;
		}
		alleles.add(allele);
	}
	
	public void removePartner(Candidate partner)
	{
		iPartners.remove(partner);
	}
	
	public boolean equals(Object o)
	{
		if (o == null) return false;
		if (o == this) return true;
		if (!(o instanceof Candidate)) return false;
		Candidate c = (Candidate)o;
		if (this.getInheritancePattern() != c.getInheritancePattern()) return false;
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
			Candidate c = new Candidate(this.iVariant);
			/*if (this.isMultiVariant()) c = new Candidate(this.iVariants);
			else c = new Candidate(this.iVariant);*/
			c.setInheritancePattern(eIPattern);
			c.setGene(g);
			c.flagDEL_candidate(bDEL);
			c.flagDUP_candidate(bDUP);
			c.iAlleles.addAll(iAlleles);
			Set<Candidate> pset = iPartners.keySet();
			for (Candidate p : pset)
			{
				c.addPartner(p, iPartners.get(p));
			}
			copies.add(c);
		}
		Collections.sort(copies);
		return copies;
	}
	
	/**
	 * Split a multiallelic candidate into a series of monoallelic candidates.
	 * This is required for certain operations.
	 * <br>WARNING: Because unrolling should occur before comp het partnering, any partners
	 * this candidate has will NOT be carried over to its copies!
	 * @return A list of new monoallelic candidates derived from this candidate, if
	 * possible. Null if this candidate is already monoallelic (or somehow has no alleles...).
	 */
	public List<Candidate> unroll()
	{
		int nalleles = 0;
		if (iAlleles == null) return null;
		nalleles = iAlleles.size();
		if (nalleles == 1) return null;
		List<Candidate> clist = new ArrayList<Candidate>(nalleles);
		for (int i : iAlleles)
		{
			Candidate mac = new Candidate(this.getVariant());
			mac.setGene(iGene);
			mac.setInheritancePattern(eIPattern);
			mac.flagDEL_candidate(bDEL);
			mac.flagDUP_candidate(bDUP);
			mac.iAlleles.add(i);
			clist.add(mac);
		}
		return clist;
	}
	
}
