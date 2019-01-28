package hospelhornbg_genomeBuild;

import java.io.IOException;
import java.io.InputStream;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import waffleoRai_Utils.CompositeBuffer;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;
import waffleoRai_Compression.huffman.Huffman;

/*
 * UPDATES
 * 
 * 1.2.0 | February 26, 2018
 * 	Forgot to add 4 to offset when feeding to huff decoder!
 * 
 * 1.2.1 | March 6, 2018
 * 	Added a method that returns sorted list of chromosomes
 * 
 * 1.2.2 | July 18, 2018
 * 	Input stream constructor was not advancing input stream!!
 * 
 * 1.2.3 | January 4, 2019
 * 	Huffman class was moved. Resolved import reference.
 * 
 * 1.3.0 | January 17, 2019
 * 	!!IMPORTANT!! Fields added to file format (will have to regen existing files!)
 * 
 * 1.4.0 | January 18, 2019
 * 	!!IMPORTANT!! Fields added to file format (will have to regen existing files!)
 * 	Now stores Contig UIDs in file and GenomeBuild object
 * 
 */

/**
 * A container for information about a genome build, such as the contigs present, their
 * various aliases, and their lengths.
 * @author Blythe Hospelhorn
 * @version 1.4.0
 * @since January 18, 2019
 *
 */
public class GenomeBuild {

	public static final int CURRENT_VERSION = 2;
	
	public static final String GBLD_MAGIC = "GBLD";
	public static final String GBDH_MAGIC = "GBDH";
	
	private static final int FIELDLEN_UCSCNAME = 64;
	private static final int FIELDLEN_UDPNAME = 16;
	
	public static final String PACKAGEPATH_36 = "resources/NCBI36.gbdh";
	public static final String PACKAGEPATH_37 = "resources/GRCh37.gbdh";
	public static final String PACKAGEPATH_38 = "resources/GRCh38.gbdh";
	
	//public static final String DEBUGPATH_36 = "C:\\Users\\Blythe\\NCBI36.gbdh";
	//public static final String DEBUGPATH_37 = "C:\\Users\\Blythe\\GRCh37.gbdh";
	//public static final String DEBUGPATH_38 = "C:\\Users\\Blythe\\GRCh38.gbdh";
	
	public static final String DEBUGPATH_36 = "X:\\usr\\hospelhornbg\\Java\\db\\NCBI36.gbdh";
	public static final String DEBUGPATH_37 = "X:\\usr\\hospelhornbg\\Java\\db\\GRCh37.gbdh";
	public static final String DEBUGPATH_38 = "X:\\usr\\hospelhornbg\\Java\\db\\GRCh38.gbdh";
	
	private String species;
	private String buildName;
	private GenomeBuildUID uid_enum;
	
	private Map<String, Contig> contigMap;
	private Map<Integer, Contig> UIDMap;
	
	public GenomeBuild(String speciesID, String buildID, GenomeBuildUID uide)
	{
		contigMap = new HashMap<String, Contig>();
		UIDMap = new HashMap<Integer, Contig>();
		species = speciesID;
		buildName = buildID;
		uid_enum = uide; 
	}
	
	public GenomeBuild(String filePath) throws IOException, UnsupportedFileTypeException
	{
		contigMap = new HashMap<String, Contig>();
		UIDMap = new HashMap<Integer, Contig>();
		parseGLBD(filePath);
	}
	
	public GenomeBuild(InputStream stream) throws IOException, UnsupportedFileTypeException
	{
		FileBuffer myFile = new FileBuffer(1024 * 500); //500KB
		contigMap = new HashMap<String, Contig>();
		UIDMap = new HashMap<Integer, Contig>();
	//	int sz = 0;
		int b = stream.read();
		while (b != -1)
		{
			//sz++;
			//System.err.println("GenomeBuild.<init> || DEBUG: sz = " + sz + " b = " + b);
			byte y = (byte)b;
			myFile.addToFile(y);
			b = stream.read();
		}
		parseGLBD(myFile);
	}
	
	private void parseGLBD(String filePath) throws IOException, UnsupportedFileTypeException
	{
		FileBuffer genome = FileBuffer.createBuffer(filePath, true);
		parseGLBD(genome);
	}
	
	private void parseGLBD(FileBuffer genome) throws UnsupportedFileTypeException, IOException
	{
		//Multi byte fields are Big-Endian
		
		// ASCII "GBLD" [4]
		// Version [4] (v2+)
		// GB UID [4] (-1 if not known) (v2+)
		// Species name length [2]
		// Species name [variable, ascii]
		// Padding
		// Build name length [2]
		// Build name [variable, ascii]
		// Padding
		// Number of contigs [4]
		// (Contig blocks)
		
		// Contig block
			// Contig block size [4]
			// Contig length [8]
			// Contig type [4]
			// UCSC name [64]
			// UDP name [16]
			// Number of other names [4]
				// Name length [2]
				// Name
				// Padding
		
		long cPos = genome.findString(0, 0x10, GBLD_MAGIC);
		if (cPos < 0){
			//Check if compressed!
			cPos = genome.findString(0, 0x10, GBDH_MAGIC);
			if (cPos < 0) throw new FileBuffer.UnsupportedFileTypeException();
			FileBuffer compressed = genome;
			genome = Huffman.HuffDecodeFile(compressed, cPos + 4);
			cPos = genome.findString(0, 0x10, GBLD_MAGIC);
			if (cPos < 0) throw new FileBuffer.UnsupportedFileTypeException();
		}
		
		cPos += 4;
		
		//Version
		// This field isn't in version 1, so if the first nlen field (and following chars) in a 
		// version 1 file happens to be 2, then the file will be parsed as a 
		// v2 file (ie incorrectly)
		// It's best to just update the files!
		int version = genome.intFromFile(cPos);
		if (version < 2 || version > 2) version = 1;
		else cPos += 4; //Only advances cPos if v2+
		
		//GB UID
		if (version >= 2)
		{
			int gbuid = genome.intFromFile(cPos);
			cPos += 4;
			uid_enum = GenomeBuildUID.getByID(gbuid);
		}
		
		int nlen = (int)genome.shortFromFile(cPos); cPos += 2;
		species = genome.getASCII_string(cPos, nlen); cPos += nlen;
		if (nlen % 2 != 0) cPos++;
		
		nlen = (int)genome.shortFromFile(cPos); cPos += 2;
		buildName = genome.getASCII_string(cPos, nlen); cPos += nlen;
		if (nlen % 2 != 0) cPos++;
		
		int contigCount = genome.intFromFile(cPos); cPos += 4;
		
		for (int i = 0; i < contigCount; i++)
		{
			long sPos = cPos + 4;
			int bSz = genome.intFromFile(cPos); cPos += 4;
			long cLen = genome.longFromFile(cPos); cPos += 8;
			int cType = genome.intFromFile(cPos); cPos += 4;
			
			String UCSC = genome.getASCII_string(cPos, FIELDLEN_UCSCNAME); cPos += FIELDLEN_UCSCNAME;
			String UDP = genome.getASCII_string(cPos, FIELDLEN_UDPNAME); cPos += FIELDLEN_UDPNAME;
			
			int nameCount = genome.intFromFile(cPos); cPos += 4;
			Set<String> nameSet = new HashSet<String>();
			for (int j = 0; j < nameCount; j++)
			{
				nlen = (int)genome.shortFromFile(cPos); cPos += 2;
				String cname = genome.getASCII_string(cPos, nlen); cPos += nlen;
				if (nlen % 2 != 0) cPos++;
				nameSet.add(cname);
			}
			
			Contig c = new Contig();
			c.setLength(cLen);
			c.setType(cType);
			c.setUCSCName(UCSC);
			c.setUDPName(UDP);
			for (String n : nameSet) c.addName(n);
			
			addContig(c);
			
			cPos = sPos + bSz;
			
		}
	}
	
	public void addContig(Contig c)
	{
		if (c == null) return;
		Collection<String> allnames = c.getAllNames();
		for (String n : allnames)
		{
			contigMap.put(n, c);
		}
		UIDMap.put(c.getUDPName().hashCode(), c);
	}
	
	public Contig getContig(String contigName)
	{
		return contigMap.get(contigName);
	}
	
	public Contig getContigByUID(int uid)
	{
		return UIDMap.get(uid);
	}
	
	public void removeContig(String contigName)
	{
		Contig c = contigMap.get(contigName);
		if (c == null) return;
		Collection<String> cnames = c.getAllNames();
		for (String n : cnames) contigMap.remove(n);
		/*Set<Integer> uidkeyset = new HashSet<Integer>();
		uidkeyset.addAll(UIDMap.keySet());
		for (Integer k : uidkeyset)
		{
			Contig v = UIDMap.get(k);
			if (v == c) UIDMap.remove(k);		
		}*/
		UIDMap.remove(c.getUDPName().hashCode());
	}
	
	public void removeContig(int uid)
	{
		Contig c = UIDMap.remove(uid);
		if (c == null) return;
		Collection<String> cnames = c.getAllNames();
		for (String n : cnames) contigMap.remove(n);
	}
	
	public String getSpeciesID()
	{
		return species;
	}
	
	public String getBuildName()
	{
		return buildName;
	}
	
	public GenomeBuildUID getUIDEnum()
	{
		return this.uid_enum;
	}
	
	public List<Contig> getChromosomes()
	{
		List<Contig> clist = new LinkedList<Contig>(); //Sorts
		Set<Contig> cSet = new HashSet<Contig>(); //Removes duplicates
		Collection<Contig> vcoll = contigMap.values(); //Access
		cSet.addAll(vcoll);
		clist.addAll(cSet);
		Collections.sort(clist);
		
		return clist;
	}
	
	public void saveGLBD(String writePath, boolean compress) throws IOException
	{
		if (writePath == null) return;
		if (writePath.isEmpty()) return;
		
		Collection<Contig> allContigs = contigMap.values();
		
		FileBuffer file = new CompositeBuffer(1 + allContigs.size());
		FileBuffer header = new FileBuffer(12 + 8 + buildName.length() + 1 + species.length() + 1);
		
		header.printASCIIToFile(GBLD_MAGIC);
		header.addToFile(CURRENT_VERSION);
		if (this.uid_enum != null) header.addToFile(this.uid_enum.getUID());
		else header.addToFile(-1);
		
		short nlen = (short)species.length();
		header.addToFile(nlen);
		header.printASCIIToFile(species);
		if (nlen % 2 != 0) header.addToFile((byte)0x00);
		
		nlen = (short)buildName.length();
		header.addToFile(nlen);
		header.printASCIIToFile(buildName);
		if (nlen % 2 != 0) header.addToFile((byte)0x00);
		
		header.addToFile(allContigs.size());
		
		file.addToFile(header);
		for (Contig c : allContigs)
		{
			file.addToFile(c.serialize());
		}
		
		//System.out.println("GenomeBuild.saveGLBD || File serialized! Number of contigs: " + allContigs.size());
		
		if (compress)
		{
			FileBuffer raw = file;
			file = Huffman.HuffEncodeFile(raw, 8, GBDH_MAGIC);
		}
		
		file.writeFile(writePath);
	}
	
	private static Map<String, StandardBuild> standardBuildMap;
	private static Map<String, GenomeBuild> loadedBuildMap;
	
	private static class StandardBuild
	{
		public String packagePath;
		public String debugPath;
		public Set<String> names;
		
		public StandardBuild(String pPath, String dPath)
		{
			packagePath = pPath;
			debugPath = dPath;
			names = new HashSet<String>();
		}
	}
	
	public static void populateStandardMap()
	{
		standardBuildMap = new HashMap<String, StandardBuild>();
		
		// Human genome 36
		StandardBuild b = new StandardBuild(PACKAGEPATH_36, DEBUGPATH_36);
		b.names.add("ncbi36");
		standardBuildMap.put("ncbi36", b);
		b.names.add("hg18");
		standardBuildMap.put("hg18", b);
		
		// Human genome 37
		b = new StandardBuild(PACKAGEPATH_37, DEBUGPATH_37);
		b.names.add("grch37");
		standardBuildMap.put("grch37", b);
		b.names.add("hg19");
		standardBuildMap.put("hg19", b);
		
		// Human genome 38
		b = new StandardBuild(PACKAGEPATH_38, DEBUGPATH_38);
		b.names.add("grch38");
		standardBuildMap.put("grch38", b);
		b.names.add("hg38");
		standardBuildMap.put("hg38", b);
		
	}
	
	public static GenomeBuild loadStandardBuild(String name)
	{
		if (loadedBuildMap == null) loadedBuildMap = new HashMap<String, GenomeBuild>();
		String key = name.toLowerCase();
		
		//First, check builds that have already been loaded!
		GenomeBuild gb = loadedBuildMap.get(key);
		if (gb != null) return gb;
		
		//If not loaded, try to load it.
		if (standardBuildMap == null) populateStandardMap();
		StandardBuild b = standardBuildMap.get(key);
		if (b == null) return null;
		
		//First try jar path
		String path = b.packagePath;
		InputStream is = GenomeBuild.class.getResourceAsStream(path);
		//System.err.println("GenomeBuild.loadStandardBuild || DEBUG: path = " + path + " IS null? " + (is == null));
		if (is == null)
		{
			//If that didn't work, try debug path
			System.err.println("GenomeBuild.loadStandardBuild || Standard path for " + key + " (" + path + ") did not work. Trying debug path...");
			path = b.debugPath;	
			try {
				GenomeBuild build = new GenomeBuild(path);
				for (String n : b.names) loadedBuildMap.put(n, build);
				return build;
			} 
			catch (IOException e) 
			{
				System.err.println("GenomeBuild.loadStandardBuild || Genome \"" + name + "\" could not be found!");
				e.printStackTrace();
			} 
			catch (UnsupportedFileTypeException e) 
			{
				System.err.println("GenomeBuild.loadStandardBuild || Genome \"" + name + "\" could not be read!");
				e.printStackTrace();
			}
		}
		else
		{
			try {
				GenomeBuild build = new GenomeBuild(is);
				is.close();
				for (String n : b.names) loadedBuildMap.put(n, build);
				return build;
			} 
			catch (IOException e) 
			{
				System.err.println("GenomeBuild.loadStandardBuild || Genome \"" + name + "\" could not be found!");
				e.printStackTrace();
			} 
			catch (UnsupportedFileTypeException e) 
			{
				System.err.println("GenomeBuild.loadStandardBuild || Genome \"" + name + "\" could not be read!");
				e.printStackTrace();
			}	
		}
		
		return null;
	}
	
	public static Set<String> getAllStandardBuildKeys()
	{
		if (standardBuildMap == null) populateStandardMap();
		return standardBuildMap.keySet();
	}
	
	public static String getStandardBuildJARPath(String buildKey)
	{
		if (standardBuildMap == null) populateStandardMap();
		StandardBuild b = standardBuildMap.get(buildKey);
		if (b == null) return null;
		return b.packagePath;
	}
	
	public void printMe()
	{
		List<Contig> clist = getChromosomes();
		System.out.println("Build Name: " + buildName);
		System.out.println("Species: " + species);
		for (Contig c : clist)
		{
			System.out.println(c.printInfo());
		}
	}
	
	public void setUID(GenomeBuildUID uid)
	{
		this.uid_enum = uid;
	}
	
}
