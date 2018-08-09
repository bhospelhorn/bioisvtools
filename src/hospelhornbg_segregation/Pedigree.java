package hospelhornbg_segregation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import hospelhornbg_bioinformatics.AffectedStatus;
import hospelhornbg_bioinformatics.Genotype;
import hospelhornbg_bioinformatics.Sex;
import hospelhornbg_bioinformatics.Variant;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class Pedigree {

	//Can read from PED file
	
	private String familyName;
	private Individual proband;
	
	private Map<String, Individual> indivMap;
	
	private List<Individual> iUnaffected;
	private List<Individual> iAffected;
	
 	public Pedigree(Individual PB)
	{
		proband = PB;
		indivMap = new HashMap<String, Individual>();
		familyName = proband.getName();
	}
	
	public Pedigree(String pedfile) throws IOException, UnsupportedFileTypeException
	{
		//Fam Indiv Father Mother Sex Affected
		indivMap = new HashMap<String, Individual>();
		Map<String, String[]> pmap = new HashMap<String, String[]>();
		
		FileReader fr = new FileReader(pedfile);
		BufferedReader br = new BufferedReader(fr);
		String line = null;
		while ((line = br.readLine()) != null)
		{
			String[] fields = line.split("\t");
			if (fields == null || fields.length < 6)
			{
				br.close();
				fr.close();
				throw new FileBuffer.UnsupportedFileTypeException();
			}
			String fam = fields[0];
			String name = fields[1];
			String father = fields[2];
			String mother = fields[3];
			String sexstr = fields[4];
			String affstr = fields[5];
			
			if (familyName == null || familyName.isEmpty()) familyName = fam;
			Individual i = new Individual(name);
			String[] parr = {father, mother};
			pmap.put(name, parr);
			if (sexstr.equals("1")) i.setSex(Sex.MALE);
			else if (sexstr.equals("2")) i.setSex(Sex.FEMALE);
			if (affstr.equals("1")) i.setAffectedStatus(AffectedStatus.UNAFFECTED);
			else if (affstr.equals("2")) i.setAffectedStatus(AffectedStatus.AFFECTED);
			indivMap.put(name, i);
			if (name.equals(familyName)) proband = i;
		}
		
		br.close();
		fr.close();
		
		//Now, link parents
		Collection<Individual> indivs = indivMap.values();
		for (Individual i : indivs)
		{
			String[] parents = pmap.get(i.getName());
			if (parents == null) continue;
			if (parents.length >= 1)
			{
				String father = parents[0];
				if (!father.equals("0"))
				{
					Individual p2 = indivMap.get(father);
					i.setParent2(p2);	
				}
			}
			if (parents.length >= 2)
			{
				String mother = parents[1];
				if (!mother.equals("0"))
				{
					Individual p1 = indivMap.get(mother);
					i.setParent1(p1);	
				}
			}
		}
		
	}

	/* --- Getters --- */
	
	public String getFamilyName()
	{
		return familyName;
	}
	
	public List<Individual> getAllMembers()
	{
		List<Individual> ilist = new LinkedList<Individual>();
		ilist.addAll(indivMap.values());
		return ilist;
	}
	
	/* --- Setters --- */
	
	public void setProband(String pbID)
	{
		Individual newpb = indivMap.get(pbID);
		if (newpb == null) return;
		if (!newpb.isAffected()) return;
		proband = newpb;
	}
	
	/* --- Lists --- */
	
	public void regenerateAffectedLists()
	{
		//Affected
		List<Individual> aff = new LinkedList<Individual>();
		Collection<Individual> all = indivMap.values();
		for (Individual i : all)
		{
			if (i.isAffected()) aff.add(i);
		}
		iAffected = aff;
		//Unaffected
		List<Individual> unaff = new LinkedList<Individual>();
		for (Individual i : all)
		{
			if (!i.isAffected()) unaff.add(i);
		}
		iUnaffected = unaff;
	}
	
	public List<Individual> getAllAffected()
	{
		if (iAffected != null) return iAffected;
		regenerateAffectedLists();
		return iAffected;
	}
	
	public List<Individual> getAllUnaffected()
	{
		if (iUnaffected != null) return iUnaffected;
		regenerateAffectedLists();
		return iUnaffected;
	}
	
	/* --- Inheritance Patterns --- */
	
	private void adjustInheritance(Individual i, Candidate c)
	{
		if (i == null) return;
		if (c == null) return;
		Individual p1 = i.getParent1();
		Individual p2 = i.getParent2();
		Inheritance ip = c.getInheritancePattern();
		Variant v = c.getVariant();
		switch(ip)
		{
		case DOMINANT:
			//One or the other parent must
				// A. Be affected
				// B. AND have the allele in question
			//If these conditions are not met, candidate is relabeled denovo_dom
			if (p1 != null)
			{
				if (p1.isAffected())
				{
					Genotype p1g = v.getSampleGenotype(p1.getName());
					if (p1g != null)
					{
						if (p1g.hasAllele(c.getAllele())) return;
					}
				}
			}
			if (p2 != null)
			{
				if (p2.isAffected())
				{
					Genotype p2g = v.getSampleGenotype(p2.getName());
					if (p2g != null)
					{
						if (p2g.hasAllele(c.getAllele())) return;
					}
				}
			}
			c.setInheritancePattern(Inheritance.DENOVO_DOM);
			break;
		case HALF_HET:
			//Exactly one parent must have this allele
			//If this condition is not met, candidate is relabeled denovo_het
			Genotype p1g = null;
			Genotype p2g = null;
			if (p1 != null) p1g = v.getSampleGenotype(p1.getName());
			if (p2 != null) p2g = v.getSampleGenotype(p2.getName());
			int a = c.getAllele();
			if (p1g != null && p2g != null)
			{
				boolean p1h = p1g.hasAllele(a);
				boolean p2h = p2g.hasAllele(a);
				if (p1h && !p2h) return;
				if (!p1h && p2h) return;
				if (p1h && p2h) c.setInheritancePattern(Inheritance.UNRESOLVED);
			}
			else if (p1g != null && p2g == null) return; //Inconclusive
			else if (p1g == null && p2g != null) return;  //Inconclusive
			c.setInheritancePattern(Inheritance.DENOVO_HET);
			return;
		case HALF_HET_SV:
			//Exactly one parent must have this CNV
			//If this condition is not met, candidate is relabeled denovo_het_sv
			Genotype p1g_hhsv = null;
			Genotype p2g_hhsv = null;
			if (p1 != null) p1g_hhsv = v.getSampleGenotype(p1.getName());
			if (p2 != null) p2g_hhsv = v.getSampleGenotype(p2.getName());
			boolean isdel = c.isDEL_candidate();
			boolean isdup = c.isDUP_candidate();
			if (p1g_hhsv != null && p2g_hhsv != null)
			{
				if (isdel)
				{
					boolean p1h = p1g_hhsv.getCopyNumber() < 2;
					boolean p2h = p2g_hhsv.getCopyNumber() < 2;
					if (p1h && !p2h) return;
					if (!p1h && p2h) return;
					if (p1h && p2h) c.setInheritancePattern(Inheritance.UNRESOLVED);
				}
				else if (isdup)
				{
					boolean p1h = p1g_hhsv.getCopyNumber() > 2;
					boolean p2h = p2g_hhsv.getCopyNumber() > 2;
					if (p1h && !p2h) return;
					if (!p1h && p2h) return;
					if (p1h && p2h) c.setInheritancePattern(Inheritance.UNRESOLVED);
				}
				else c.setInheritancePattern(Inheritance.UNRESOLVED);
			}
			else if (p1g_hhsv!= null && p2g_hhsv == null) return; //Inconclusive
			else if (p1g_hhsv == null && p2g_hhsv != null) return;  //Inconclusive
			c.setInheritancePattern(Inheritance.DENOVO_HET_SV);
			break;
		case HOM_REC:
			//Unaffected parents must be het for this allele
			//Affected parents must be hom for this allele
			//If these conditions are not met, candidate is relabeled denovo_rec
			int a_hr = c.getAllele();
			if (p1 != null)
			{
				Genotype p1g_hr = v.getSampleGenotype(p1.getName());
				if(p1g_hr != null)
				{
					if (!p1g_hr.hasAllele(a_hr))
					{
						c.setInheritancePattern(Inheritance.DENOVO_REC);
						return;
					}
					if (p1.isAffected())
					{
						if (!p1g_hr.isHomozygous())
						{
							c.setInheritancePattern(Inheritance.DENOVO_REC);
							return;
						}
					}
					else
					{
						if (p1g_hr.isHomozygous())
						{
							c.setInheritancePattern(Inheritance.DENOVO_REC);
							return;
						}
					}
				}
			}
			if (p2 != null)
			{
				Genotype p2g_hr = v.getSampleGenotype(p2.getName());
				if(p2g_hr != null)
				{
					if (!p2g_hr.hasAllele(a_hr))
					{
						c.setInheritancePattern(Inheritance.DENOVO_REC);
						return;
					}
					if (p2.isAffected())
					{
						if (!p2g_hr.isHomozygous())
						{
							c.setInheritancePattern(Inheritance.DENOVO_REC);
							return;
						}
					}
					else
					{
						if (p2g_hr.isHomozygous())
						{
							c.setInheritancePattern(Inheritance.DENOVO_REC);
							return;
						}
					}
				}
			}
			//c.setInheritancePattern(Inheritance.DENOVO_REC);
			break;
		case UNRESOLVED:
			//DNS - Ignored for now
			break;
		default:
			break;
		}
	}
	
	private void setInheritance_HOMOZYGOUS(Candidate c)
	{
		//Build genotype map
		Map<Individual, Genotype> genomap = new HashMap<Individual, Genotype>();
		if (indivMap == null) return;
		Collection<Individual> fam = indivMap.values();
		if (fam == null || fam.isEmpty()) return;
		Variant v = c.getVariant();
		if (v == null) return;
		for (Individual i : fam)
		{
			Genotype g = v.getSampleGenotype(i.getName());
			if (g == null) continue;
			genomap.put(i, g);
		}
		
		Inheritor.checkHeterozygousCandidate(genomap, c);
		//Now, check against local pedigree...
		adjustInheritance(proband, c);
		
	}
	
	private void setInheritance_HETEROZYGOUS(Candidate c)
	{
		Map<Individual, Genotype> genomap = new HashMap<Individual, Genotype>();
		if (indivMap == null) return;
		Collection<Individual> fam = indivMap.values();
		if (fam == null || fam.isEmpty()) return;
		Variant v = c.getVariant();
		if (v == null) return;
		for (Individual i : fam)
		{
			Genotype g = v.getSampleGenotype(i.getName());
			if (g == null) continue;
			genomap.put(i, g);
		}
		
		Inheritor.checkHomozygousCandidate(genomap, c);
		//Now, check against local pedigree...
		adjustInheritance(proband, c);
	}
	
	public Candidate getCandidate(Variant v)
	{
		Genotype pbgeno = v.getSampleGenotype(proband.getName());
		if (pbgeno == null) return null;
		
		Candidate c = new Candidate(v);
		
		boolean pbhomo = pbgeno.isHomozygous();
		if (pbhomo) setInheritance_HOMOZYGOUS(c);
		else setInheritance_HETEROZYGOUS(c);
		
		return c;
	}
	
	public boolean passCompHetPair(Candidate hh1, Candidate hh2)
	{
		if (!hh1.isMonoallelic()) throw new IllegalArgumentException();
		if (!hh2.isMonoallelic()) throw new IllegalArgumentException();
		//Passes if at least one affected has both
		// and no unaffected has both
		List<Individual> aff = getAllAffected();
		List<Individual> unaff = getAllUnaffected();
		if (aff.isEmpty()) return false;
		Variant v1 = hh1.getVariant();
		Variant v2 = hh2.getVariant();
		int a1 = hh1.getAllele();
		int a2 = hh2.getAllele();
		for (Individual u : unaff)
		{	
			//Get genotypes
			Genotype g1 = v1.getSampleGenotype(u.getName());
			Genotype g2 = v2.getSampleGenotype(u.getName());
			if (g1 == null || g2 == null) continue;
			
			//See which alleles from each variant u has...
			//Eliminate any combinations found in u...
			if (g1.hasAllele(a1) && g2.hasAllele(a2)) return false; //That won't work, right?
		}
		
		//Check affecteds to make sure at least one has combo
		for (Individual a : aff)
		{	
			Genotype g1 = v1.getSampleGenotype(a.getName());
			Genotype g2 = v2.getSampleGenotype(a.getName());
			if (g1 == null || g2 == null) continue;
			if (g1.hasAllele(a1) && g2.hasAllele(a2)) return true; //We good?
		}
		
		return false;
	}

}
