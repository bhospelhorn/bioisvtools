package hospelhornbg_segregation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
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
import hospelhornbg_bioinformatics.Sex;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.GenomeBuild;
import hospelhornbg_genomeBuild.TwoSexChromSegModel;
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
	
	public boolean setProband(String pbID)
	{
		Individual newpb = indivMap.get(pbID);
		if (newpb == null) return false;
		if (!newpb.isAffected()) return false;
		proband = newpb;
		return true;
	}
	
	public void setFamilyName(String name)
	{
		this.familyName = name;
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
	
	private void adjustSexChromGenotypes(Collection<Variant> variants, TwoSexChromSegModel sxm)
	{
		//TODO: Rewrite this method
		//Need to move Y PAR -> X remapping per VARIANT, not indiv!
		//Also. Make copies of weird calls and store in list to return!
		List<Individual> ilist = getAllMembers();
		for(Variant v : variants)
		{
			Contig c = v.getChromosome();
			if (c.getType() == Contig.SORTCLASS_SEXCHROM)
			{
				for(Individual i : ilist)
				{
					Genotype g = v.getSampleGenotype(i.getName());
					if (g != null)
					{
						int exCN = 2;
						if (c.getUDPName().contains("X"))
						{
							//Check if PAR
							if (!sxm.inHomChromPAR(v.getPosition()))
							{
								exCN = i.getExpectedXCount();	
							}
						}
						else if (c.getUDPName().contains("Y"))
						{
							//Check PAR
							if (!sxm.inHetChromPAR(v.getPosition()))
							{
								exCN = i.getExpectedYCount();
							}
							else
							{
								//Remap to X
								exCN = 2;	
								int xcoord = sxm.mapHetPosToHom(v.getPosition());
								v.setChromosome(sxm.getHomogameticChrom());
								v.setPosition(xcoord);
							}
						}
						switch (exCN)
						{
						case 0:
							//Take out all alleles, mark CN0
							int[] all0 = {-1, -1};
							g.setAlleles(all0);
							g.setCopyNumber(0);
							break;
						case 1:
							//Determine which allele, mark CN1
							int[] myall = g.getAlleles();
							int a = -1;
							for (int ia : myall) if (ia > a) a = ia;
							int[] all1 = {a};
							g.setAlleles(all1);
							g.setCopyNumber(1);
							break;
						case 2:
							//Leave alone
							break;
						case 3:
							//Add a -1 allele, mark CN3
							int[] myall3 = g.getAlleles();
							if (myall3.length < 3)
							{
								int[] all3 = {-1, -1, -1};
								if (myall3.length >= 1) all3[0] = myall3[0];
								if (myall3.length >= 2) all3[1] = myall3[1];
							}
							g.setCopyNumber(3);
							break;
						}
					}
				}	
			}
		}
	}
	
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
		
		//TODO: Handle XY variants differently
		//Genotypes need to be edited to reflect expected copy number!
		
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

	/* --- Serialization --- */
	
	public static void writeToPED(Pedigree fam, String pedpath) throws IOException
	{
		if (fam == null) return;
		if (pedpath == null) return;
		
		FileWriter fw = new FileWriter(pedpath);
		BufferedWriter bw = new BufferedWriter(fw);
		
		bw.write("#FamID\tSampleID\tFather\tMother\tSex\tAffected\n");
		List<Individual> ilist = fam.getAllMembers();
		Collections.sort(ilist);
		for(Individual i : ilist)
		{
			bw.write(fam.getFamilyName() + "\t");
			bw.write(i.getName() + "\t");
			//Get dad
			Individual dad = i.getFather();
			if (dad != null) bw.write(dad.getName() + "\t");
			else bw.write("0\t");
			//Get mom
			Individual mom = i.getMother();
			if (mom != null) bw.write(mom.getName() + "\t");
			else bw.write("0\t");
			//The others...
			Sex isex = i.getSex();
			if (isex == Sex.MALE) bw.write("1\t");
			else if (isex == Sex.FEMALE) bw.write("2\t");
			else bw.write("0\t");
			AffectedStatus as = i.getAffectedStatus();
			if (as == AffectedStatus.UNAFFECTED) bw.write("1\n");
			else if (as == AffectedStatus.AFFECTED) bw.write("2\n");
			else bw.write("0\n");
		}
		
		bw.close();
	}
	
}
