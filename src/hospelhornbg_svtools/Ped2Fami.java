package hospelhornbg_svtools;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import hospelhornbg_segregation.Family;
import hospelhornbg_segregation.FamilyMember;
import hospelhornbg_segregation.Population;

public class Ped2Fami {
	
	public static final String OP_INPUT_PED = "-i"; 
	public static final String OP_INPUT_TABLE = "-t";
	public static final String OP_PROBAND = "-P";
	public static final String OP_OUTPUT = "-o"; 
	
	public static final String TBLOP_SAMPLENAME = "SAMPLENAME"; 
	public static final String TBLOP_FIRSTNAME = "GIVENNAME"; 
	public static final String TBLOP_LASTNAME = "FAMNAME"; 
	public static final String TBLOP_BYEAR = "BIRTHYEAR"; 
	public static final String TBLOP_DYEAR = "DEATHYEAR"; 
	public static final String TBLOP_ETHLIST = "ETHNICITY_TAGS"; 
	
	/*
	 * Input Table Information...
	 * 	It must be a tsv. The fields can be in any order.
	 * 	The only required field is SAMPLENAME.
	 */
	
	public static void printUsage()
	{
		System.out.println("--------------------------------------------------------------------------------------------------------------");
		System.out.println("BioisvTools || PED2FAMI");
		System.out.println();
		System.out.println("Purpose: Converts a standard text-based PED file to a binary FAMI (.fam) file.");
		System.out.println("FAMI files can contain more information and metadata.");
		System.out.println();
		System.out.println("Input Formats:");
		System.out.println("\tPED file for basic family data.");
		System.out.println("\tAdditionally, a tsv (tab-separated) metadata table can be provided. (See below for valid fields)");
		System.out.println();
		System.out.println("Flags:");
		System.out.println("\t-i\tFILE\t[Required]\t\tPath to input PED file.");
		System.out.println("\t-o\tFILE\t[Required]\t\tDesired FAMI output path.");
		System.out.println("\t-t\tFILE\t[Optional]\t\tPath to tsv metadata table.");
		System.out.println("\t-P\tSTRING\t[Optional]\t\tProband ID");
		System.out.println();
		System.out.println("Metadata Table Formatting:");
		System.out.println("The column header line must start with #");
		System.out.println("The following column header fields are recognized:");
		System.out.println("SAMPLENAME [Required] - Name of sample/individual as it appears in PED.");
		System.out.println("GIVENNAME [Optional] - Given (first) name of individual.");
		System.out.println("FAMNAME [Optional] - Family (last) name of individual.");
		System.out.println("BIRTHYEAR [Optional] - Birth year of individual. Must be integer.");
		System.out.println("DEATHYEAR [Optional] - Death year of individual. Must be integer or N/A.");
		System.out.println("ETHNICITY_TAGS [Optional] - Comma delimited list of ethnic groups to associate with individual.");
		System.out.println("\tThese can be used for population allele totals. The following tags are recognized at this time:");
		System.out.println("\t\tNFE - Non-Finnish European");
		System.out.println("\t\tAFR - African/African-American");
		System.out.println("\t\tAMR - Native American/Hispanic/Latino");
		System.out.println("\t\tSAS - South Asian");
		System.out.println("\t\tEAS - East Asian");
		System.out.println("\t\tASJ - Ashkenazi Jewish");
		System.out.println("\t\tFIN - Finnish");
		System.out.println("\t\tOTH - Other");
		System.out.println();
		System.out.println("Sample Usage:");
		System.out.println("java -jar bioisvtools.jar ped2fami -p myfam.ped -o myfam.fam");
		System.out.println("java -jar bioisvtools.jar ped2fami -p myfam.ped -t myfam_metatable.tsv -o myfam.fam");
		System.out.println("java -jar bioisvtools.jar ped2fami -v -p myfam.ped -o myfam.fam -P SAMPLE-PB");
		System.out.println();
		System.out.println("--------------------------------------------------------------------------------------------------------------");
	}
	
	public static void runPed2Fami(String[] args)
	{
		String pedPath = null;
		String tblPath = null;
		String outPath = null;
		String pbName = null;
		
		//Read args
		for (int i = 0; i < args.length; i++)
		{
			String s = args[i];
			if (s.equals(OP_INPUT_PED))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_INPUT_PED + " flag MUST be followed by input PED path!");
					printUsage();
					System.exit(1);
				}
				pedPath = args[i+1];
			}
			else if (s.equals(OP_INPUT_TABLE))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_INPUT_TABLE + " flag MUST be followed by input tsv path!");
					printUsage();
					System.exit(1);
				}
				tblPath = args[i+1];
			}
			else if (s.equals(OP_OUTPUT))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_OUTPUT + " flag MUST be followed by output fam path!");
					printUsage();
					System.exit(1);
				}
				outPath = args[i+1];
			}
			else if (s.equals(OP_PROBAND))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_PROBAND + " flag MUST be followed by proband sample ID!");
					printUsage();
					System.exit(1);
				}
				pbName = args[i+1];
			}
		}
		
		
		//Check Args
		if(pedPath == null || pedPath.isEmpty())
		{
			System.err.println("ERROR! Input PED file is required!");
			printUsage();
			System.exit(1);
		}
		if(outPath == null || outPath.isEmpty())
		{
			System.err.println("ERROR! Output path is required!");
			printUsage();
			System.exit(1);
		}
		
		//Read PED file
		Family fam = null;
		try
		{
			Map<String, Family> fmap = Family.readFromPED(pedPath);
			if (fmap.size() != 1)
			{
				System.err.println("ERROR! Input PED file contains more than one family! Please provide a single family PED file!");
				System.exit(1);
			}
			for(Family f : fmap.values())
			{
				fam = f;
				break;
			}
		}
		catch(IOException e)
		{
			System.err.println("ERROR! IO Error - Input PED file could not be opened!");
			e.printStackTrace();
			System.exit(1);
		}
		
		//Set Proband
		if (pbName != null && !pbName.isEmpty())
		{
			fam.setProband(pbName);
		}
		
		//Read Table
		Map<String, Map<String, String>> table = new HashMap<String, Map<String, String>>();
		if(tblPath != null &&  !tblPath.isEmpty())
		{
			System.err.println("Reading metadata table...");
			try
			{
				BufferedReader br = new BufferedReader(new FileReader(tblPath));
				String line = null;
				//Get header line
				line = br.readLine();
				//String # from header line
				line = line.substring(1);
				String[] columns = line.split("\t");
				//Find the sample name column...
				int snind = -1;
				for(int i = 0; i < columns.length; i++)
				{
					if(columns[i].equals(TBLOP_SAMPLENAME))
					{
						snind = i;
						break;
					}
				}
				if (snind < 0)
				{
					System.err.println("ERROR! Sample name field was not found in table! Table will be ignored...");
					br.close();
				}
				else
				{
					while((line = br.readLine()) != null)
					{
						String[] fields = line.split("\t");
						//Get sample name...
						String sn = fields[snind];
						//Dump other fields into map
						Map<String, String> smap = new HashMap<String, String>();
						for(int i = 0; i < fields.length; i++)
						{
							if (i == snind) continue;
							if (i >= columns.length) break;
							smap.put(columns[i], fields[i]);
						}
						table.put(sn, smap);
					}
					br.close();	
				}
				
			}
			catch(IOException e)
			{
				System.err.println("ERROR! IO Error - Input metadata table could not be opened!");
				e.printStackTrace();
				System.exit(1);
			}
		}
		
		//Apply metadata
		List<FamilyMember> members = fam.getAllFamilyMembers();
		for(FamilyMember m : members)
		{
			Map<String, String> smap = table.get(m.getName());
			if (smap == null) continue;
			//Scan through all fields...
			String val = null;
			val = smap.get(TBLOP_FIRSTNAME);
			if (val != null) m.setFirstName(val);
			val = smap.get(TBLOP_LASTNAME);
			if (val != null) m.setLastName(val);
			val = smap.get(TBLOP_BYEAR);
			if (val != null)
			{
				try {m.setBirthYear(Integer.parseInt(val));}
				catch(NumberFormatException e) {e.printStackTrace();}
			}
			val = smap.get(TBLOP_DYEAR);
			if (val != null)
			{
				try {m.setDeathYear(Integer.parseInt(val));}
				catch(NumberFormatException e) {e.printStackTrace();}
			}
			val = smap.get(TBLOP_ETHLIST);
			if(val != null)
			{
				String[] groups = val.split(",");
				if (groups.length > 0)
				{
					for(String g : groups)
					{
						Population p = Population.getPopulation(g);
						if (p != null) m.addPopulationTag(p);
					}
				}
			}
		}
		
		//Write output file
		try 
		{
			Family.writeToFAMI(fam, outPath, true);
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR! IO Error - Output file could not be written!");
			e.printStackTrace();
			System.exit(1);
		}
		
	}


}
