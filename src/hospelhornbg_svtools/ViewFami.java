package hospelhornbg_svtools;

import java.io.IOException;
import java.util.Collection;
import java.util.List;

import hospelhornbg_segregation.Family;
import hospelhornbg_segregation.FamilyMember;
import hospelhornbg_segregation.Individual;
import hospelhornbg_segregation.Population;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class ViewFami {
	
	public static void runViewFami(String[] args)
	{
		//The last argument should be the fami file
		
		String famfile = args[args.length - 1];
		
		if(!FileBuffer.fileExists(famfile))
		{
			System.err.println("ERROR: FAMI file \"" + famfile + "\" does not appear to exist!");
			System.exit(1);
		}
		
		Family fam = null;
		try 
		{
			fam = Family.readFromFAMI(famfile);
		} 
		catch (IOException e)
		{
			System.err.println("ERROR: FAMI file \"" + famfile + "\" could not be opened!");
			e.printStackTrace();
			System.exit(1);
		} 
		catch (UnsupportedFileTypeException e) 
		{
			System.err.println("ERROR: FAMI file \"" + famfile + "\" could not be read!");
			e.printStackTrace();
			System.exit(1);
		}
		
		
		//Print family information
		System.out.println("Family Name: " + fam.getFamilyName());
		Individual pb = fam.getProband();
		System.out.println("Proband: " + pb.getName());
		
		//By member
		List<FamilyMember> members = fam.getAllFamilyMembers();
		String cond = fam.getSegregatingCondition();
		for(FamilyMember mem : members)
		{
			System.out.println();
			System.out.println(mem.getName());
			System.out.println("Full Name/Alias: " + mem.getFirstName() + " " + mem.getLastName());
			System.out.println("Internal UID: " + Integer.toHexString(mem.getUID()));
			System.out.println("Relationship to Proband: " + fam.getRelationshipString_ENG(mem));
			System.out.println("Chromosomal Sex: " + mem.getSex());
			System.out.println("Phenotypic Sex: " + mem.getPhenotypicSex());
			System.out.println("Affected Status for " + cond + ": " + mem.getAffectedStatus());
			System.out.print("Population Tags: ");
			Collection<Population> ptags = mem.getPopulationTags();
			for (Population p : ptags) System.out.print(p.getShortString() + " ");
			System.out.println();
			Individual father = mem.getFather();
			Individual mother = mem.getMother();
			if (mother != null) System.out.println("Mother: " + mother.getName());
			else System.out.println("Mother: (UNKNOWN)");
			if (father != null) System.out.println("Father: " + father.getName());
			else System.out.println("Father: (UNKNOWN)");
		}
		
	}

}
