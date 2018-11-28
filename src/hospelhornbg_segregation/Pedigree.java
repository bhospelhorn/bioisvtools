package hospelhornbg_segregation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

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
	
	protected Pedigree()
	{
		
	}
	
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
	
	public Individual getProband()
	{
		return proband;
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
	
	private void annotateCandidate(Map<Individual, Genotype> genomap, Individual pb, Candidate c)
	{
		Genotype g = genomap.get(pb);
		if (!g.hasAllele(c.getAllele())) return;
		if (g.isHomozygous()) Inheritor.checkHomozygousCandidate(genomap, pb, c);
		else Inheritor.checkHeterozygousCandidate(genomap, pb, c);
		
		//Check denovo
		Inheritor.adjustInheritance(genomap, pb, c);
	}
	
	public List<Candidate> toCandidates(Variant v)
	{
		List<Candidate> clist = new LinkedList<Candidate>();
		
		//Get all alleles
		List<Individual> affected = getAllAffected();
		Set<Integer> alleles = new HashSet<Integer>();
		for (Individual aff : affected)
		{
			Genotype g = v.getSampleGenotype(aff.getName());
			int[] all = g.getAlleles();
			for (int a : all) alleles.add(a);
		}
		
		//Generate a candidate for each allele
		for (int a : alleles)
		{
			Candidate c = new Candidate(v, a);
			clist.add(c);
		}
		
		//Build genomap
		Map<Individual, Genotype> genomap = new HashMap<Individual, Genotype>();
		Collection<Individual> allindivs = indivMap.values();
		for (Individual i : allindivs)
		{
			genomap.put(i, v.getSampleGenotype(i.getName()));
		}
		
		//Get inheritance pattern for each candidate for each affected
		for (Candidate c : clist)
		{
			for (Individual i : affected)
			{
				annotateCandidate(genomap, i, c);
			}
		}
		
		return clist;
	}

	public boolean checkPair(Candidate c1, Candidate c2)
	{
		//No unaffected family member can have both
		List<Individual> unaff = getAllUnaffected();
		if (unaff == null) return true;
		for (Individual u : unaff)
		{
			boolean has1 = false;
			boolean has2 = false;
			Variant v1 = c1.getVariant();
			int a1 = c1.getAllele();
			Genotype g = v1.getSampleGenotype(u.getName());
			if (g == null) continue;
			if (g.hasAllele(a1)) has1 = true;
			
			Variant v2 = c2.getVariant();
			int a2 = c2.getAllele();
			g = v2.getSampleGenotype(u.getName());
			if (g == null) continue;
			if (g.hasAllele(a2)) has2 = true;
			if (has1 && has2) return false;
		}
		return true;
	}
	
	/* --- Relationship --- */
	
	public String getRelationshipString_ENG(Individual relative)
	{
		Relationship r = proband.getRelationship(relative);
		if (r == null) return "[Unrelated]";
		return r.toString_English();
	}
	
}
