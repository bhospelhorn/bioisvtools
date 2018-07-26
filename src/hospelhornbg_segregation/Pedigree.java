package hospelhornbg_segregation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import hospelhornbg_bioinformatics.AffectedStatus;
import hospelhornbg_bioinformatics.Sex;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class Pedigree {

	//Can read from PED file
	
	private String familyName;
	private Individual proband;
	
	private Map<String, Individual> indivMap;
	
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

	
	
}
