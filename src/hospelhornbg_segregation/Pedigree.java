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
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.TwoSexChromSegModel;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class Pedigree{
	
	public static final int SC_GENO_TYPE_AUTO = 0;
	public static final int SC_GENO_TYPE_X = 1;
	public static final int SC_GENO_TYPE_Y = 2;

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
	
	public int countMembers()
	{
		if (indivMap == null) return 0;
		return indivMap.size();
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
	
	public void changeSampleName(String oldname, String newname)
	{
		Individual i = indivMap.remove(oldname);
		if (i == null) return;
		i.setSampleName(newname);
		indivMap.put(newname, i);
	}
	
	protected void setProband(Individual pb)
	{
		proband = pb;
	}
	
	public Individual removeIndividual(String name)
	{
		iUnaffected = null;
		iAffected = null;
		
		Individual i = indivMap.remove(name);
		if (i == null) return null;
		
		if(proband != null && proband.getName().equals(name)) proband = null;
		
		return i;
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
	
	public void adjustSexChromGenotypes(Collection<Variant> variants, TwoSexChromSegModel sxm)
	{
		//Get indiv list
		List<Individual> ilist = getAllMembers();
		for(Variant v : variants)
		{
			int treatAs = adjustSexChromCall(v, sxm);
			if (treatAs == Pedigree.SC_GENO_TYPE_AUTO) continue;
			for(Individual indiv : ilist)
			{
				int cn = 2;
				Genotype g = v.getSampleGenotype(indiv.getName());
				if (g == null) 
				{
					g = new Genotype();
					v.addGenotype(indiv.getName(), g);
				}
				//Determine new CN...
				if(treatAs == Pedigree.SC_GENO_TYPE_X) cn = indiv.getExpectedXCount();
				else if (treatAs == Pedigree.SC_GENO_TYPE_Y) cn = indiv.getExpectedYCount();
				//Alter genotype to match expected CN
				switch(cn)
				{
				case 0: 
					//Delete genotype from variant?
					//For now, I'll set CN0 and make one -1 allele...
					g.setCopyNumber(0);
					int[] a0 = {-1};
					g.setAlleles(a0);
					break;
				case 1: 
					//Delete all alleles but one, prioritizing alt alleles
					int[] a1 = g.getAlleles();
					g.setCopyNumber(1);
					int all1 = 0;
					for (int a : a1)
					{
						if (a > all1) all1 = a;
					}
					int[] a1_1 = {all1};
					g.setAlleles(a1_1);
					break;
				case 2: 
					//Do nothing, should be CN2 by default!
					break;
				case 3: 
					//Add a new ref allele?
					int[] a3 = g.getAlleles();
					int a3_1 = 0;
					int a3_2 = 0;
					int a3_3 = 0;
					if (a3.length >= 1) a3_1 = a3[0];
					if (a3.length >= 2) a3_2 = a3[1];
					g.setCopyNumber(3);
					int[] a3_new = {a3_1, a3_2, a3_3};
					g.setAlleles(a3_new);
					break;
				}
			}
		}
	}
	
	private int adjustSexChromCall(Variant v, TwoSexChromSegModel sxm)
	{
		//Get chrom
		int pval = Pedigree.SC_GENO_TYPE_AUTO;
		Contig c = v.getChromosome();
		if (c == sxm.getHomogameticChrom())
		{
			//Must be same object
			//See if pseudoautosomal
			int pos = v.getPosition();
			if(!sxm.inHomChromPAR(pos)) pval = Pedigree.SC_GENO_TYPE_X;
			//If not, mark as X
		}
		else if (c == sxm.getHeterogameticChrom())
		{
			//See if pseudoautosomal
			//If so, remap variant to X
			//If not, mark as Y
			int pos = v.getPosition();
			if(sxm.inHetChromPAR(pos))
			{
				int xpos = sxm.mapHetPosToHom(pos);
				v.setChromosome(sxm.getHomogameticChrom());
				v.setPosition(xpos);
			}
			else pval = Pedigree.SC_GENO_TYPE_Y;
		}
		
		//Repeat with end if variant is a structural variant...
		int eval = Pedigree.SC_GENO_TYPE_AUTO;
		if (v instanceof StructuralVariant)
		{
			StructuralVariant sv = (StructuralVariant) v;
			Contig endChrom = sv.getEndChromosome();
			if (endChrom != null && endChrom != c)
			{
				int endpos = sv.getEndPosition();
				if (endChrom == sxm.getHomogameticChrom())
				{
					if(!sxm.inHomChromPAR(endpos)) eval = Pedigree.SC_GENO_TYPE_X;
				}
				else if (endChrom == sxm.getHeterogameticChrom())
				{
					if(sxm.inHetChromPAR(endpos))
					{
						int xpos = sxm.mapHetPosToHom(endpos);
						sv.setEndChromosome(sxm.getHomogameticChrom());
						sv.setEndPosition(xpos);
					}
					else eval = Pedigree.SC_GENO_TYPE_Y;
				}	
			}
		}
		else return pval;
		
		//If equal...
		if (pval == eval) return pval;
		//X or Y over Auto
		if (pval == Pedigree.SC_GENO_TYPE_AUTO) return eval;
		if (eval == Pedigree.SC_GENO_TYPE_AUTO) return pval;
		//X vs. Y? Treat like Y?
		if (pval == Pedigree.SC_GENO_TYPE_X) return eval;
		if (eval == Pedigree.SC_GENO_TYPE_X) return pval;
		
		return Pedigree.SC_GENO_TYPE_AUTO;
	}
	
	private int getSegregationType(Variant v, TwoSexChromSegModel sxm)
	{
		Contig c = v.getChromosome();
		if (c == null) return Pedigree.SC_GENO_TYPE_AUTO;
		
		int type = Pedigree.SC_GENO_TYPE_AUTO;
		if(c == sxm.getHomogameticChrom()) {
			//Check pseudoautosomal
			if(!sxm.inHomChromPAR(v.getPosition())) type = Pedigree.SC_GENO_TYPE_X;
		}
		if(c == sxm.getHeterogameticChrom()) {
			if(!sxm.inHetChromPAR(v.getPosition())) return Pedigree.SC_GENO_TYPE_Y;
		}
		
		if(v instanceof StructuralVariant)
		{
			StructuralVariant sv = (StructuralVariant)v;
			Contig e = sv.getEndChromosome();
			if(e != null && e != c)
			{
				if(e == sxm.getHomogameticChrom()) {
					//Check pseudoautosomal
					if(!sxm.inHomChromPAR(sv.getEndPosition())) type = Pedigree.SC_GENO_TYPE_X;
				}
				if(e == sxm.getHeterogameticChrom()) {
					if(!sxm.inHetChromPAR(sv.getEndPosition())) return Pedigree.SC_GENO_TYPE_Y;
				}
			}
		}
		
		return type;
	}
	
	private void annotateCandidate(Map<Individual, Genotype> genomap, Individual pb, Candidate c, TwoSexChromSegModel sxm)
	{
		Genotype g = genomap.get(pb);
		if (!g.hasAllele(c.getAllele())) return;
		
		//Check type
		int treatAs = getSegregationType(c.getVariant(), sxm);
		switch(treatAs)
		{
		case SC_GENO_TYPE_AUTO:
			if (g.isHomozygous()) Inheritor.checkHomozygousCandidate(genomap, pb, c);
			else Inheritor.checkHeterozygousCandidate(genomap, pb, c);
			
			//Check denovo
			Inheritor.adjustInheritance(genomap, pb, c);
			break;
		case SC_GENO_TYPE_X:
			SexChromInheritor.checkMammalXCandidate(genomap, pb, c);
		case SC_GENO_TYPE_Y:
			SexChromInheritor.checkMammalYCandidate(genomap, pb, c);
		}
	}
	
	public List<Candidate> toCandidates(Variant v, TwoSexChromSegModel sxm)
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
				annotateCandidate(genomap, i, c, sxm);
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
