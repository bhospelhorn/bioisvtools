package hospelhornbg_segregation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import hospelhornbg_bioinformatics.AffectedStatus;
import hospelhornbg_bioinformatics.Sex;
import waffleoRai_Compression.huffman.Huffman;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;
import waffleoRai_Utils.SerializedString;
import waffleoRai_Utils.BinFieldSize;

public class Family extends Pedigree{
	
	//Pedigree with more information
	//Multiple conditions/phenotypes allowed. "Affected" refers to currently set condition
	
	public static final int MAX_CONDITIONS = 48;
	public static final String DEFO_CONDITION = "DISEASE";
	
	public static final String FAMI_MAGIC = "FAMI";
	public static final String FAMH_MAGIC = "FAMH";
	public static final int CURRENT_FAMI_VERSION = 1;
	
	private FamilyMember iProband;
	
	private Map<Integer, FamilyMember> iMembers;
	
	private String[] conditions;
	private int segCond;
	
	protected Family()
	{
		super();
		iProband = null;
		iMembers = new HashMap<Integer, FamilyMember>();
		segCond = 0;
		conditions = new String[MAX_CONDITIONS];
		conditions[0] = DEFO_CONDITION;
	}
	
	public Family(FamilyMember proband)
	{
		iProband = proband;
		segCond = 0;
		conditions = new String[MAX_CONDITIONS];
		conditions[0] = "DISEASE";
		iMembers = new HashMap<Integer, FamilyMember>();
		if (iProband != null)
		{
			Collection<FamilyMember> relatives = iProband.getAllRelativesAsFamilyMembers();	
			for (FamilyMember i : relatives)
			{
				iMembers.put(i.getUID(), i);
			}
		}
		
	}
	
	/* --- Getters --- */
	
	public Individual getProband()
	{
		return iProband;
	}
	
	public List<Individual> getAllMembers()
	{
		List<Individual> ilist = new LinkedList<Individual>();
		ilist.addAll(iMembers.values());
		return ilist;
	}
	
	public List<FamilyMember> getAllFamilyMembers()
	{
		int sz = iMembers.size() + 1;
		List<FamilyMember> list = new ArrayList<FamilyMember>(sz);
		list.addAll(iMembers.values());
		return list;
	}
	
	public String getSegregatingCondition()
	{
		return conditions[segCond];
	}
	
	public String getCondition(int index)
	{
		if (index < 0) return null;
		if (index >= MAX_CONDITIONS) return null;
		return conditions[index];
	}
	
	public String[] getAllConditions()
	{
		return conditions.clone();
	}
	
	public FamilyMember getMember(int UID)
	{
		return this.iMembers.get(UID);
	}
	
	public FamilyMember getMemberByName(String samplename)
	{
		//First hash the sample name and see if it's a match.
		//If not, scan all members
		if (samplename == null) return null;
		int hash = samplename.hashCode();
		FamilyMember fm = iMembers.get(hash);
		if (fm != null)
		{
			if (samplename.equals(fm.getName())) return fm;
		}
		List<FamilyMember> list = this.getAllFamilyMembers();
		for (FamilyMember m : list)
		{
			if (samplename.equals(m.getName())) return m;
		}
		return null;
	}
	
	/* --- Setters --- */
	
	public boolean setProband(String pbID)
	{
		//Scan individuals...
		List<Individual> members = this.getAllMembers();
		for (Individual i : members)
		{
			if (i.getName().equals(pbID))
			{
				if (i instanceof FamilyMember)
				{
					iProband = (FamilyMember)i;
					return true;
				}
			}
		}
		return false;
	}
	
	public boolean setProband(int pbID)
	{
		FamilyMember fm = iMembers.get(pbID);
		if (fm == null) return false;
		iProband = fm;
		return true;
	}
	
	public boolean setSegregatingCondition(int index)
	{
		if (index < 0) return false;
		if (index >= MAX_CONDITIONS) return false;
		
		String c = this.conditions[index];
		if (c == null) return false;
		segCond = index;
		updateAffectedStatus();
 		return true;
	}
	
	public boolean replaceSegregatingCondition(int index, String newCond)
	{
		if (index < 0) return false;
		if (index >= MAX_CONDITIONS) return false;
		
		String c = this.conditions[index];
		if (c != null) this.removeAffectedStatus(c);
		conditions[index] = newCond;
		refreshAffectedStatus(newCond);
		
		return true;
	}
	
	public boolean addSegregatingCondition(String newCond)
	{
		//Just finds the first null slot and adds condition
		//Returns false if no null slot
		int i = 0;
		for (i = 0; i < MAX_CONDITIONS; i++)
		{
			if (conditions[i] == null) break;
		}
		if (i >= MAX_CONDITIONS) return false;
		conditions[i] = newCond;
		this.refreshAffectedStatus(newCond);
		
		return true;
	}
	
	public void addMember(FamilyMember indiv)
	{
		iMembers.put(indiv.getUID(), indiv);
	}
	
	/* --- Conditions --- */
	
	private void updateAffectedStatus()
	{
		//Updates the primary "affected status" var in indivs for current condition
		Collection<FamilyMember> members = this.iMembers.values();
		String cond = this.getSegregatingCondition();
		for(FamilyMember m : members) m.setAffectedCondition(cond);
	}
	
	private void refreshAffectedStatus(String condition)
	{
		//Adds an entry in the affected map of all fam members for
		// this condition. Sets to unaffected for everyone.
		Collection<FamilyMember> members = this.iMembers.values();
		for(FamilyMember m : members) m.setAffectedStatus(condition, AffectedStatus.UNAFFECTED);
	}
	
	private void removeAffectedStatus(String condition)
	{
		//Removes entry in affected map for all fam members
		Collection<FamilyMember> members = this.iMembers.values();
		for(FamilyMember m : members) m.removeAffectedStatus(condition);
	}
	
	/* --- Parsing --- */
	
	private static class UnlinkedFam
	{
		public Family family;
		public Map<Integer, String> momMap;
		public Map<Integer, String> dadMap;
		
		public UnlinkedFam(Family f)
		{
			family = f;
			momMap = new HashMap<Integer, String>();
			dadMap = new HashMap<Integer, String>();
		}
	}
	
	private static void printPEDParsingError(int lineno)
	{
		System.err.println("Family.readFromPED || Parsing Error - Line " + lineno + " invalid PED record. Line ignored.");
	}
	
	public static Map<String, Family> readFromPED(String pedpath) throws IOException
	{
		//Fam Indiv Father Mother Sex Affected
		Map<String, UnlinkedFam> rawmap = new HashMap<String, UnlinkedFam>();
		Map<String, Family> fammap = new HashMap<String, Family>();
		
		//Open the file!
		FileReader fr = new FileReader(pedpath);
		BufferedReader br = new BufferedReader(fr);
		
		String line = null;
		int lcount = 0;
		while((line = br.readLine()) != null)
		{
			lcount++;
			if (line.startsWith("#")) continue; //Ignore these
			String[] fields = line.split("\t");
			if (fields.length != 6)
			{
				printPEDParsingError(lcount);
				continue;
			}
			String fam = fields[0];
			UnlinkedFam f = rawmap.get(fam);
			if (f == null)
			{
				f = new UnlinkedFam(new Family());
				f.family.setFamilyName(fam);
				rawmap.put(fam, f);
			}
			
			//Get Indiv name
			FamilyMember indiv = new FamilyMember(fields[1]);
			String sexstr = fields[4];
			String affstr = fields[5];
			if (sexstr.equals("1")) 
			{
				indiv.setSex(Sex.MALE);
				indiv.setPhenotypicSex(Sex.MALE);
			}
			else if (sexstr.equals("2")) 
			{
				indiv.setSex(Sex.FEMALE);
				indiv.setPhenotypicSex(Sex.FEMALE);
			}
			else 
			{
				indiv.setSex(Sex.UNKNOWN);
				indiv.setPhenotypicSex(Sex.UNKNOWN);
			}
			
			if (affstr.equals("1"))
			{
				indiv.setAffectedStatus(AffectedStatus.UNAFFECTED);
				indiv.setAffectedStatus(DEFO_CONDITION, AffectedStatus.UNAFFECTED);
			}
			else if (affstr.equals("2")) 
			{
				indiv.setAffectedStatus(AffectedStatus.AFFECTED);
				indiv.setAffectedStatus(DEFO_CONDITION, AffectedStatus.AFFECTED);
			}
			else 
			{
				indiv.setAffectedStatus(AffectedStatus.UNKNOWN);
				indiv.setAffectedStatus(DEFO_CONDITION, AffectedStatus.UNKNOWN);
			}
			
			//Parents
			String mname = fields[3];
			if(!mname.equals("0")) f.momMap.put(indiv.getUID(), mname);
			String fname = fields[2];
			if(!fname.equals("0")) f.dadMap.put(indiv.getUID(), fname);
		}
		
		br.close();
		
		//Link parents & dump to return map
		Collection<UnlinkedFam> allfam = rawmap.values();
		for (UnlinkedFam uf : allfam)
		{
			//Go through the mom map...
			Set<Integer> keys = uf.momMap.keySet();
			for (Integer i : keys)
			{
				//Link the individual to their mother
				FamilyMember indiv = uf.family.getMember(i);
				FamilyMember mom = uf.family.getMemberByName(uf.momMap.get(i));
				if (mom != null) indiv.setParent1(mom);
			}

			//Go through the dad map...
			keys = uf.dadMap.keySet();
			for (Integer i : keys)
			{
				//Link the individual to their father
				FamilyMember indiv = uf.family.getMember(i);
				FamilyMember dad = uf.family.getMemberByName(uf.dadMap.get(i));
				if (dad != null) {
					if(indiv.getParent1() != null) indiv.setParent2(dad);
					else indiv.setParent1(dad);
				}
			}
			
			//Dump the family in the final map
			fammap.put(uf.family.getFamilyName(), uf.family);
		}
		
		
		return fammap;
	}
	
	private static Sex getFAMISexEnum(int famiVal)
	{
		switch(famiVal)
		{
		case 0: return Sex.FEMALE;
		case 1: return Sex.MALE;
		case 2: return Sex.OTHER;
		}
		return Sex.UNKNOWN;
	}
	
	private static AffectedStatus getFAMIAffectedEnum(int asVal)
	{
		switch(asVal)
		{
		case 0: return AffectedStatus.UNAFFECTED;
		case 1: return AffectedStatus.AFFECTED;
		case 2: return AffectedStatus.PARTIALLY_AFFECTED;
		case 3: return AffectedStatus.POSSIBLY_AFFECTED;
		}
		return AffectedStatus.UNKNOWN;
	}
	
	public static Family readFromFAMI(String path) throws IOException, UnsupportedFileTypeException
	{
		FileBuffer myfile = FileBuffer.createBuffer(path, true);
		long cpos = myfile.findString(0, 0x10, FAMI_MAGIC);
		
		if (cpos != 0)
		{
			//Try decompressing
			cpos = myfile.findString(0, 0x10, FAMH_MAGIC);
			if (cpos != 0) throw new FileBuffer.UnsupportedFileTypeException();
			FileBuffer comp = myfile;
			myfile = Huffman.HuffDecodeFile(comp, FAMH_MAGIC.length());
			cpos = myfile.findString(0, 0x10, FAMI_MAGIC);
			if (cpos != 0) throw new FileBuffer.UnsupportedFileTypeException();
		}
		
		cpos += 4; //Skip magic now
		int pbid = myfile.intFromFile(cpos); cpos += 4;
		int indivCount = myfile.intFromFile(cpos); cpos += 4;
		int version = myfile.intFromFile(cpos); cpos += 4;
		if (version > CURRENT_FAMI_VERSION) throw new FileBuffer.UnsupportedFileTypeException();
		
		String famName = myfile.getASCII_string(cpos, 48); cpos += 48;
		
		//Instantiate family...
		Family f = new Family();
		f.setFamilyName(famName);
		
		//Condition name table...
		int ncond = myfile.intFromFile(cpos); cpos += 4;
		if (ncond > Family.MAX_CONDITIONS) throw new FileBuffer.UnsupportedFileTypeException();
		for (int i = 0; i < ncond; i++)
		{
			SerializedString ss = myfile.readVariableLengthString(cpos, BinFieldSize.WORD, 2);
			cpos += ss.getSizeOnDisk();
			f.replaceSegregatingCondition(i, ss.getString());
		}
		
		//Individual table...
		Map<Integer, Integer> p1map = new HashMap<Integer, Integer>();
		Map<Integer, Integer> p2map = new HashMap<Integer, Integer>();
		
		for(int i = 0; i < indivCount; i++)
		{
			int uid = myfile.intFromFile(cpos); cpos += 4;
			int p1id = myfile.intFromFile(cpos); cpos += 4;
			int p2id = myfile.intFromFile(cpos); cpos += 4;
			int csx = myfile.getValueAsInt(cpos, BinFieldSize.WORD); cpos += 2;
			int psx = myfile.getValueAsInt(cpos, BinFieldSize.WORD); cpos += 2;
			int[] aff = new int[MAX_CONDITIONS];
			for(int j = 0; j < MAX_CONDITIONS; j += 2)
			{
				int n = Byte.toUnsignedInt(myfile.getByte(cpos)); cpos++;
				aff[j] = (n >>> 4) & 0xF;
				aff[j+1] = n & 0xF;
			}
			//Skip next two bytes
			cpos += 2;
			int byear = myfile.getValueAsInt(cpos, BinFieldSize.WORD); cpos += 2;
			//Skip next two bytes
			cpos += 2;
			int dyear = myfile.getValueAsInt(cpos, BinFieldSize.WORD); cpos += 2;
			
			//Sample name
			SerializedString ss = myfile.readVariableLengthString(cpos, BinFieldSize.WORD, 2);
			cpos += ss.getSizeOnDisk();
			String samplename = ss.getString();
			
			//Last name
			ss = myfile.readVariableLengthString(cpos, BinFieldSize.WORD, 2);
			cpos += ss.getSizeOnDisk();
			String lastname = ss.getString();
			
			//First name
			ss = myfile.readVariableLengthString(cpos, BinFieldSize.WORD, 2);
			cpos += ss.getSizeOnDisk();
			String firstname = ss.getString();
			
			//Middle name(s)
			int nmn = myfile.intFromFile(cpos); cpos += 4;
			String[] midnames = null;
			if (nmn > 0)
			{
				midnames = new String[nmn];
				for (int j = 0; j < nmn; j++)
				{
					ss = myfile.readVariableLengthString(cpos, BinFieldSize.WORD, 2);
					cpos += ss.getSizeOnDisk();
					midnames[j] = ss.getString();
				}
			}
			
			//Create indiv, map parent ids, dump in fam
			FamilyMember fm = new FamilyMember(samplename);
			fm.setUID(uid);
			if(p1id != -1) p1map.put(uid, p1id);
			if(p2id != -1) p2map.put(uid, p2id);
			fm.setSex(getFAMISexEnum(csx));
			fm.setPhenotypicSex(getFAMISexEnum(psx));
			for (int j = 0; j < ncond; j++)
			{
				fm.setAffectedStatus(f.getCondition(j), Family.getFAMIAffectedEnum(aff[j]));
			}
			fm.setFirstName(firstname);
			fm.setLastName(lastname);
			fm.setBirthYear(byear);
			fm.setDeathYear(dyear);
			if(midnames != null)
			{
				for (int j = 0; j < midnames.length; j++) fm.addMiddleName(midnames[j]);
			}
			
			f.addMember(fm);
		}
		
		//Link parents
		Set<Integer> hasp1 = p1map.keySet();
		for (Integer i : hasp1)
		{
			FamilyMember p1 = f.getMember(p1map.get(i));
			if (p1 != null)
			{
				FamilyMember child = f.getMember(i);
				if(child != null) child.setParent1(p1);
			}
		}
		
		Set<Integer> hasp2 = p2map.keySet();
		for (Integer i : hasp2)
		{
			FamilyMember p2 = f.getMember(p2map.get(i));
			if (p2 != null)
			{
				FamilyMember child = f.getMember(i);
				if(child != null) child.setParent1(p2);
			}
		}
		
		//Set proband
		f.setProband(pbid);
		
		return f;
	}
	
	/* --- Serialization --- */
	
	private static short getFAMISexEnum(Sex s)
	{
		if(s == null) return -1;
		switch(s)
		{
		case FEMALE: return 0;
		case MALE: return 1;
		case OTHER: return 2;
		case UNKNOWN: return 3;
		}
		return -1;
	}
	
	private static int getFAMIAffectedEnum(AffectedStatus as)
	{
		if(as == null) return -1;
		switch(as)
		{
		case AFFECTED: return 1;
		case PARTIALLY_AFFECTED: return 2;
		case POSSIBLY_AFFECTED: return 3;
		case UNAFFECTED: return 0;
		case UNKNOWN: return 4;
		}
		return -1;
	}
	
	public static void writeToFAMI(Family fam, String path, boolean huff) throws IOException
	{
		if (fam == null) return;
		//Estimate needed size for buffer...
		int estsz = 64 + (48 * 64); //Header
		List<FamilyMember> ilist = fam.getAllFamilyMembers();
		int nindiv = ilist.size();
		estsz += nindiv * (48 + 64);
		
		FileBuffer myfile = FileBuffer.createWritableBuffer("famiwrite", estsz, true);
		
		myfile.printASCIIToFile(FAMI_MAGIC);
		if (fam.iProband != null) myfile.addToFile(fam.iProband.getUID());
		else myfile.addToFile(-1);
		myfile.addToFile(nindiv);
		String famname = fam.getFamilyName();
		if (famname.length() > 48) famname = famname.substring(0, 48);
		myfile.printASCIIToFile(famname);
		int slen = famname.length();
		if (slen < 48)
		{
			for (int i = slen; i < 48; i++)
			{
				myfile.addToFile(FileBuffer.ZERO_BYTE);
			}
		}
		myfile.addToFile(CURRENT_FAMI_VERSION);
		
		//Affected condition string table
		String[] ctbl = new String[MAX_CONDITIONS];
		int j = 0;
		for(int i = 0; i < MAX_CONDITIONS; i++)
		{
			String cname = fam.getCondition(i);
			if (cname != null)
			{
				ctbl[j] = cname;
				j++;
			}
		}
		myfile.addToFile(j);
		for(int i = 0; i < j; i++)
		{
			myfile.addVariableLengthString(ctbl[i], BinFieldSize.WORD, 2);
		}
		
		//Add individuals
		for (FamilyMember m : ilist)
		{
			myfile.addToFile(m.getUID());
			FamilyMember p1 = m.getParent1AsFamilyMember();
			if (p1 != null) myfile.addToFile(p1.getUID());
			else myfile.addToFile(-1);
			FamilyMember p2 = m.getParent2AsFamilyMember();
			if (p2 != null) myfile.addToFile(p2.getUID());
			else myfile.addToFile(-1);
			myfile.addToFile(getFAMISexEnum(m.getSex()));
			myfile.addToFile(getFAMISexEnum(m.getPhenotypicSex()));
			//Condition enums...
			for(int i = 0; i < MAX_CONDITIONS; i += 2)
			{
				//Grab two conditions...
				String c1 = ctbl[i];
				String c2 = ctbl[i+1];
				
				int wbyte = 0;
				if(c1 != null)
				{
					AffectedStatus as = m.getAffectedStatus(c1);
					int asi = getFAMIAffectedEnum(as);
					asi = asi << 4;
					wbyte |= asi;
				}
				if (c2 != null)
				{
					AffectedStatus as = m.getAffectedStatus(c2);
					int asi = getFAMIAffectedEnum(as);
					wbyte |= asi;
				}
				
				myfile.addToFile((byte)wbyte);
			}
			
			//Dates (Month and day unused, set to 1/1)
			myfile.addToFile((short)0x0101);
			myfile.addToFile(m.getBirthYear());
			myfile.addToFile((short)0x0101);
			myfile.addToFile(m.getDeathYear());
			//Sample name
			myfile.addVariableLengthString(m.getName(), BinFieldSize.WORD, 2);
			//Last name
			String name = m.getLastName();
			if (name == null) m.setLastName(fam.getFamilyName());
			myfile.addVariableLengthString(m.getLastName(), BinFieldSize.WORD, 2);
			//First name
			name = m.getFirstName();
			if (name == null) m.setFirstName(m.getName());
			myfile.addVariableLengthString(m.getFirstName(), BinFieldSize.WORD, 2);
			//Middle name(s)
			List<String> mn = m.getMiddleNames();
			if (mn == null || mn.isEmpty()) myfile.addToFile(0);
			else
			{
				myfile.addToFile(mn.size());
				for(String n : mn)
				{
					myfile.addVariableLengthString(n, BinFieldSize.WORD, 2);	
				}
			}
		}
		
		
		//Write to disk
		myfile.writeFile(path);
		
	}
	
	
}
