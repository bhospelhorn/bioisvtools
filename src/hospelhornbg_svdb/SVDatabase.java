package hospelhornbg_svdb;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.time.OffsetDateTime;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import hospelhornbg_bioinformatics.Genotype;
import hospelhornbg_bioinformatics.Interval;
import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.VCF;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.Gene;
import hospelhornbg_genomeBuild.GeneFunc;
import hospelhornbg_genomeBuild.GeneSet;
import hospelhornbg_genomeBuild.GenomeBuild;
import hospelhornbg_genomeBuild.GenomeBuildUID;
import hospelhornbg_genomeBuild.TwoSexChromSegModel;
import hospelhornbg_segregation.Candidate;
import hospelhornbg_segregation.Family;
import hospelhornbg_segregation.FamilyMember;
import hospelhornbg_segregation.Individual;
import hospelhornbg_segregation.Inheritance;
import hospelhornbg_segregation.Inheritor;
import hospelhornbg_segregation.Population;
import hospelhornbg_svdb.SVDBGenotype.SVDBAllele;
import waffleoRai_Utils.BinFieldSize;
import waffleoRai_Utils.CompositeBuffer;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;
import waffleoRai_Utils.SerializedString;
import waffleoRai_Utils.StreamBuffer;

//TODO: I think the DBSample is a bit redundant with the family information?
//Maybe just get rid of the DBSample altogether and just map FamilyMember objects to UIDs/Sample IDs?

public class SVDatabase {
	
	/*
	 * Path: dir/svdb_name/(files)
	 * 
	 * 
	 */
	
	/* --- Constants --- */
	
	public static final String TABLE_EXTENSION = ".tsv";
	public static final String INDEX_EXTENSION = ".cidx";
	public static final String FULL_INDEX_EXTENSION = ".vidx";
	
	public static final String INDIV_TABLE_EXT = "sampt";
	
	public static final String DBINFO_MAGIC = "svdbINFO";
	public static final String DBINFO_EXT = "dbinfo";
	
	public static final String FAMILY_DIRECTORY = "fam";
	
	public static final double POPFREQ_CUTOFF_G2 = 0.05;
	public static final double POPFREQ_CUTOFF_G3 = 0.02;
	public static final double POPFREQ_CUTOFF_G4 = 0.01;
	
	/* --- Instance Variables --- */
	
	private String dirPath;
	private String dbName;
	private GenomeBuild genomeBuild;
	private GeneSet geneSet;
	
	private int leeway; //0-1000
	
	private int indivCount;
	//private int famCount;
	
	private Map<Population, Integer> popIndivCount;
	
	private Map<String, DBSample> sampleMap;
	private Map<Integer, DBSample> sampleIDMap;
	private Map<String, Family> familyMap;
	
	private GenotypeTable genoTable;
	
	/* --- Construction --- */
	
	private SVDatabase(String directory, String name)
	{
		dirPath = directory;
		dbName = name;
		resetPopulationCountMap();
		sampleMap = new HashMap<String, DBSample>();
		familyMap = new HashMap<String, Family>();
	}
	
	private SVDatabase(String directory, String name, GenomeBuild genome)
	{
		dirPath = directory;
		dbName = name;
		genomeBuild = genome;
		resetPopulationCountMap();
		sampleMap = new HashMap<String, DBSample>();
		familyMap = new HashMap<String, Family>();
		geneSet = GeneSet.loadRefGene(genome);
	}
	
	private void resetPopulationCountMap()
	{
		popIndivCount = new HashMap<Population, Integer>();
		Population[] allp = Population.values();
		for (Population p : allp) popIndivCount.put(p, 0);
	}
	
	public static SVDatabase readDatabase(String directory, String name) throws IOException
	{
		//Load genome build...
		SVDatabase svdb = new SVDatabase(directory, name);
		
		String dbinfo = svdb.getDBInfoFilePath();
		FileBuffer file = FileBuffer.createBuffer(dbinfo, true);
		long cpos = file.findString(0, 0x10, DBINFO_MAGIC);
		if (cpos != 0) return null;
		cpos = 8;
		SerializedString ss = file.readVariableLengthString(cpos, BinFieldSize.WORD, 2);
		cpos += ss.getSizeOnDisk();
		
		String gbName = ss.getString();
		svdb.genomeBuild = GenomeBuild.loadStandardBuild(gbName);
		svdb.geneSet = GeneSet.loadRefGene(svdb.genomeBuild);
		svdb.leeway = file.intFromFile(cpos);
		
		//Load sample table...
		svdb.readSampleTable();
		svdb.updatePopulationTotals();
		
		//Load family info...
		svdb.readFamilyInfo();
		
		//Set path for genotype table...
		svdb.genoTable = new GenotypeTable(svdb.getGenoTablePath());
		
		return svdb;
	}
	
	public static SVDatabase newDatabase(String directory, String name, GenomeBuild genome, int leeway) throws IOException
	{
		SVDatabase svdb = new SVDatabase(directory, name, genome);
		
		//Create files and folders...
		String dbdir = svdb.getDBDirPath();
		if (!FileBuffer.directoryExists(dbdir)) Files.createDirectories(Paths.get(dbdir));
		String famdir = dbdir + File.separator + FAMILY_DIRECTORY;
		if (!FileBuffer.directoryExists(famdir)) Files.createDirectory(Paths.get(famdir));
		svdb.leeway = leeway;
		
		//Generate info file
		svdb.writeInfoFile();
		
		//Generate empty genotype table
		svdb.genoTable = new GenotypeTable(svdb.getGenoTablePath());
		svdb.genoTable.generateEmptyTable();
		
		return svdb;
	}
	
	/* --- Paths --- */
	
	public String getDBDirPath()
	{
		return dirPath + File.separator + "svdb_" + dbName;
	}
	
	public String getVariantTablePath(SVType svtype)
	{
		return getDBDirPath() + File.separator + "svdb_" + dbName + "_" + svtype.getString() + TABLE_EXTENSION;
	}
	
	public String getGenoTablePath(SVType svtype)
	{
		return getDBDirPath() + File.separator + "svdb_" + dbName + "_genoTable_" + svtype.getString() + TABLE_EXTENSION;
	}
	
	public String getGenoTablePath()
	{
		return getDBDirPath() + File.separator + "svdb_" + dbName + "_genoTable_allTypes" + TABLE_EXTENSION;
	}
	
	public String getSampleTablePath()
	{
		return getDBDirPath() + File.separator + "svdb_" + dbName + "_sampleTable." + INDIV_TABLE_EXT;
	}
	
	public String getDBInfoFilePath()
	{
		return getDBDirPath() + File.separator + "svdb_" + dbName + "." + DBINFO_EXT;
	}
	
	public String getFamDirPath()
	{
		return getDBDirPath() + File.separator + FAMILY_DIRECTORY;
	}
	
	public String getFamFilePath(String famname)
	{
		return getFamDirPath() + File.separator + famname + ".fam";
	}
	
	public String getMasterVariantIndexPath()
	{
		return getDBDirPath() + File.separator + "svdb_" + dbName + "_masterIndex" + FULL_INDEX_EXTENSION;
	}
	
	/* --- Read --- */
	
	public void readSampleTable() throws IOException
	{
		String stpath = getSampleTablePath();
		if(!FileBuffer.fileExists(stpath)) return;
		FileBuffer tbl = FileBuffer.createBuffer(stpath, true);
		long fsz = tbl.getFileSize();
		long cpos = 0;
		
		while(cpos < fsz)
		{
			DBSample s = DBSample.readFromSampleTable(tbl, cpos);
			cpos += s.calculateSerializedSize();
			sampleMap.put(s.getName(), s);
		}
	}
	
	public void readFamilyInfo()
	{
		//Get families from sample info...
		Set<String> famnames = new HashSet<String>();
		for(DBSample s : sampleMap.values())
		{
			famnames.add(s.getFamilyName());
		}
		
		//Clear fam map
		familyMap.clear();
		
		//Read files...
		for (String famname : famnames)
		{
			String fpath = getFamFilePath(famname);
			if (!FileBuffer.fileExists(fpath))
			{
				System.err.println("SVDatabase.readFamilyInfo || WARNING: FAM file for family \"" + famname + "\" does not exist!");
			}
			else
			{
				try 
				{
					Family f = Family.readFromFAMI(fpath);
					familyMap.put(famname, f);
				}
				catch (Exception e)
				{
					System.err.println("SVDatabase.readFamilyInfo || WARNING: FAM file for family \"" + famname + "\" could not be read!");
				} 
			}
		}
		
	}
	
	public Map<Integer, VarIndexRecord> readMasterVariantIndex() throws IOException
	{
		//VarID [4]
		//SVType [1]
		//Line [3] (Unsigned)
		
		//Change this if it's likely files will exceed 16.7 million lines...
		Map<Integer, VarIndexRecord> map = new TreeMap<Integer, VarIndexRecord>();
		String path = getMasterVariantIndexPath();
		if (!FileBuffer.fileExists(path)) return map;
		
		long fsz = (FileBuffer.fileSize(path));
		long cpos = 0;
		StreamBuffer file = new StreamBuffer(path, true);
		while(cpos < fsz)
		{
			int varid = file.intFromFile(cpos); cpos += 4;
			byte btype = file.getByte(cpos); cpos++;
			int line = file.shortishFromFile(cpos); cpos += 3;
			
			SVType t = SVType.getTypeByID(Byte.toUnsignedInt(btype));
			VarIndexRecord r = new VarIndexRecord(t, line);
			map.put(varid, r);
		}
		
		return map;
	}
	
	public VarIndexRecord lookupVariant(int varUID) throws IOException
	{
		String path = getMasterVariantIndexPath();
		if (!FileBuffer.fileExists(path)) return null;
		
		long fsz = (FileBuffer.fileSize(path));
		long cpos = 0;
		StreamBuffer file = new StreamBuffer(path, true);
		while(cpos < fsz)
		{
			int varid = file.intFromFile(cpos); cpos += 4;
			if (varid == varUID)
			{
				byte btype = file.getByte(cpos); cpos++;
				int line = file.shortishFromFile(cpos); cpos += 3;
				
				SVType t = SVType.getTypeByID(Byte.toUnsignedInt(btype));
				VarIndexRecord r = new VarIndexRecord(t, line);
				return r;
			}
			cpos += 4;
		}
		
		return null;
	}
	
	public Set<Integer> getAllVariantUIDs() throws IOException
	{
		Set<Integer> idset = new TreeSet<Integer>();
		String path = getMasterVariantIndexPath();
		if (!FileBuffer.fileExists(path)) return idset;
		
		long fsz = (FileBuffer.fileSize(path));
		long cpos = 0;
		StreamBuffer file = new StreamBuffer(path, true);
		while(cpos < fsz)
		{
			int varid = file.intFromFile(cpos); cpos += 8;
			idset.add(varid);
		}
		
		
		return idset;
	}
	
	public List<SVDBGenotype> getGenotypesForVariant(int varUID)
	{
		if(genoTable == null) return null;
		try {
		return genoTable.getGenotypesForVariant(varUID);
		}
		catch (IOException e) {return null;}
	}
	
	public DBVariant getVariant(int varUID) throws IOException
	{
		VarIndexRecord location = lookupVariant(varUID);
		if (location == null) return null;
		
		SVType t = location.getType();
		int ln = location.getLine();
		
		String tblPath = this.getVariantTablePath(t);
		if (!FileBuffer.fileExists(tblPath)) return null;
		
		int lctr = 0;
		BufferedReader br = new BufferedReader(new FileReader(tblPath));
		String line = null;
		while((line = br.readLine()) != null)
		{
			lctr++;
			if(lctr == ln)
			{
				br.close();
				return DBVariant.getFromDBRecord(line, this.genomeBuild, this.geneSet);
			}
		}
		br.close();
		
		return null;
	}
	
	public List<DBVariant> getVariantsInRange(Contig chrom, int start, int end) throws IOException
	{
		List<DBVariant> found = new LinkedList<DBVariant>();
		for(SVType t : SVType.values())
		{
			String tblPath = this.getVariantTablePath(t);
			if(!FileBuffer.fileExists(tblPath)) continue;
			
			int firstLine = 1;
			if (t != SVType.TRA)
			{
				String idxPath = tblPath + INDEX_EXTENSION;
				
				if(FileBuffer.fileExists(idxPath))
				{
					VariantIndex idx = VariantIndex.buildIndexFromTable(tblPath, genomeBuild);
					firstLine = idx.getFirstLineOfContig(chrom);
				}	
			}
			
			//Open & Fast forward
			BufferedReader br = new BufferedReader(new FileReader(tblPath));
			String line = null;
			int ln = 0;
			while((line = br.readLine()) != null)
			{
				ln++;
				if(ln < firstLine) continue;
				if(line.isEmpty()) continue;
				if(line.startsWith("#")) continue;
				
				DBVariant var = DBVariant.getFromDBRecord(line, genomeBuild, geneSet);
				if (t == SVType.TRA)
				{
					//Check for either end
					if(chrom.equals(var.getChrom()))
					{
						Interval st = var.getStartPosition();
						if(start < st.getEnd() && end > st.getStart()) found.add(var);
					}
					else if (chrom.equals(var.getEndChrom()))
					{
						Interval ed = var.getEndPosition();
						if(start < ed.getEnd() && end > ed.getStart()) found.add(var);
					}
				}
				else
				{
					if(!chrom.equals(var.getChrom())) break; //Break if next chrom
					if(var.inRange(start, end)) found.add(var);	
				}
			}
			
			br.close();
		}
		
		return found;
	}
	
	public List<DBVariant> getVariantsInGene(Gene g) throws IOException
	{
		if (g == null) return null;
		Contig c = g.getChromosome();
		int st = g.getTranscriptStart();
		int ed = g.getTranscriptEnd();
		return getVariantsInRange(c, st, ed);
	}
	
	/* --- Write --- */
	
	public void indexVariantTable(SVType type)
	{
		String tblPath = getVariantTablePath(type);
		VariantIndex idx;
		try 
		{
			idx = VariantIndex.buildIndexFromTable(tblPath, genomeBuild);
			idx.writeIndexToDisk(tblPath + INDEX_EXTENSION);
		} 
		catch (IOException e) 
		{
			System.err.println("SVDatabase.indexSampleTable || ERROR! Variant Table could not be indexed!");
			e.printStackTrace();
		}
	}
	
	public void indexVariants() throws IOException
	{
		//It doesn't look like varIDs NEED to be sorted for now, so
		//I guess I'll just stream out the records? Hm.
		String ipath = this.getMasterVariantIndexPath();
		Path p = Paths.get(ipath);
		Files.deleteIfExists(p); //We will be appending.
		Files.createFile(p);
		for(SVType t : SVType.values())
		{
			String tblPath = this.getVariantTablePath(t);
			if (!FileBuffer.fileExists(tblPath)) continue;
			int lineNumber = 0;
			String line = null;
			BufferedReader br = new BufferedReader(new FileReader(tblPath));
			while((line = br.readLine()) != null)
			{
				lineNumber++;
				if (line.isEmpty()) continue;
				if (line.startsWith("#")) continue;
				//Nab the ID number from the record
				String[] fields = line.split("\t");
				if (fields.length < 2) continue;
				try
				{
					int id = Integer.parseUnsignedInt(fields[1], 16);
					byte[] rec = new byte[8];
					rec[0] = (byte)((id >>> 24) & 0xFF);
					rec[1] = (byte)((id >>> 16) & 0xFF);
					rec[2] = (byte)((id >>> 8) & 0xFF);
					rec[3] = (byte)(id & 0xFF);
					rec[4] = (byte)t.getID();
					rec[5] = (byte)((lineNumber >>> 16) & 0xFF);
					rec[6] = (byte)((lineNumber >>> 8) & 0xFF);
					rec[7] = (byte)(lineNumber & 0xFF);
					
					Files.write(p, rec, StandardOpenOption.APPEND);
					
				}
				catch(NumberFormatException e)
				{
					System.err.println("SVDatabase.indexVariants || Variant table line " + lineNumber + " could not be read!");
				}
			}
			br.close();
		}
		
	}
	
	public void writeInfoFile() throws IOException
	{
		String path = this.getDBInfoFilePath();
		String bname = genomeBuild.getBuildName();
		FileBuffer outfile = new FileBuffer(8 + bname.length() + 3 + 4, true);
		outfile.printASCIIToFile(DBINFO_MAGIC);
		outfile.addVariableLengthString(bname, BinFieldSize.WORD, 2);
		outfile.addToFile(leeway);
		outfile.writeFile(path);
	}
	
	public void saveDatabase() throws IOException
	{
		//Variants and genotypes are always on disk. 
		//Only samples and family data need to be updated.
		writeSampleTable();
		writeFamilyFiles();
	}
	
	public void writeSampleTable() throws IOException
	{
		//This is a simple table, no header
		FileBuffer tbl = new CompositeBuffer(sampleMap.size());
		for(DBSample s : sampleMap.values())
		{
			FileBuffer srec = s.serializeSample();
			tbl.addToFile(srec);
		}
		tbl.writeFile(this.getSampleTablePath());
	}
	
	public void writeFamilyFiles() throws IOException
	{
		for(Family f : familyMap.values())
		{
			String fpath = this.getFamFilePath(f.getFamilyName());
			Family.writeToFAMI(f, fpath, true);
		}
	}
	
	/* --- Variant Addition/Deletion --- */
	
	public boolean addFamily(String vcfpath, Family fam) throws IOException, UnsupportedFileTypeException
	{
		//Import family
		familyMap.put(fam.getFamilyName(), fam);
		if (sampleIDMap == null) mapSampleIDs();
		List<FamilyMember> members = fam.getAllFamilyMembers();
		List<DBSample> samples = new ArrayList<DBSample>(members.size() + 1);
		for(FamilyMember m : members)
		{
			DBSample s = DBSample.createNewSample(m.getName());
			m.setUID(s.getUID());
			//If needed, regenerate ID until unique...
			while(sampleIDMap.containsKey(s.getUID()))
			{
				s.regenerateUID();
				m.setUID(s.getUID());
			}
			
			//Map
			sampleMap.put(s.getName(), s);
			sampleIDMap.put(s.getUID(), s);
			
			samples.add(s);
		}
		updatePopulationTotals();
		
		//Save sample data
		saveDatabase();
		
		//Import variants
		VariantPool famvars = VCF.readVCF(vcfpath, genomeBuild, true);
		addVariants(famvars, samples);
		
		return true;
	}
	
	private void addVariants(VariantPool varpool, Collection<DBSample> samples) throws IOException
	{
		//Split up all the variants, cast as SVs and sort by type, then chrom1
		List<Variant> rawlist = varpool.getVariants();
		
		Map<Integer, Collection<SVDBGenotype>> newgenos = new TreeMap<Integer, Collection<SVDBGenotype>>();
		
		Map<SVType, List<StructuralVariant>> typeMap = new HashMap<SVType, List<StructuralVariant>>();
		for(Variant v : rawlist)
		{
			if (v instanceof StructuralVariant)
			{
				StructuralVariant sv = (StructuralVariant)v;
				SVType t = sv.getType();
				if (typeMap.get(t) == null)
				{
					typeMap.put(t, new LinkedList<StructuralVariant>());
				}
				typeMap.get(t).add(sv);
			}
		}
		
		Set<Integer> usedUIDs = getAllVariantUIDs();
		SVType[] alltypes = SVType.values();
		for (SVType t : alltypes)
		{
			List<StructuralVariant> vlist = typeMap.get(t);
			if (vlist == null || vlist.isEmpty()) continue;
			
			/*Map<Contig, List<StructuralVariant>> chromMap = new HashMap<Contig, List<StructuralVariant>>();
			for(StructuralVariant sv : vlist)
			{
				Contig c1 = sv.getChromosome();
				if (chromMap.get(c1) == null)
				{
					chromMap.put(c1, new LinkedList<StructuralVariant>());
				}
				chromMap.get(c1).add(sv);
			}	*/
			
			//Now I guess we just stream each var file...
			String tblPath = getVariantTablePath(t);
			if (!FileBuffer.fileExists(tblPath))
			{
				//Nothing to check against! Just write everything!
				BufferedWriter bw = new BufferedWriter(new FileWriter(tblPath));
				for(StructuralVariant sv : vlist)
				{
					String vname = t.name() + sv.getChromosomeName() + "-" + sv.getPosition() + "_" + OffsetDateTime.now().getNano();
					DBVariant var = DBVariant.getFromVariant(sv, vname);
					while(usedUIDs.contains(var.getIntegerID())) var.regenerateIntegerID();
					usedUIDs.add(var.getIntegerID());
					
					//Generate genotype info!
					
					List<SVDBGenotype> genolist = new LinkedList<SVDBGenotype>();
					for(DBSample s : samples)
					{
						Genotype gt = sv.getSampleGenotype(s.getName());
						if (gt != null)
						{
							//Check if has
							boolean has = gt.hasAllele(1);
							//Check if hom
							boolean hom = false;
							if(has) hom = gt.isHomozygous();
							//Note for all populations this person is in
							if(has)var.incrementTotalCount();
							if(hom)var.incrementHomozygoteCount();
							List<Population> plist = s.getPopulationFlags();
							if(plist != null && !plist.isEmpty())
							{
								for(Population p : plist)
								{
									if(has)var.incrementTotalCount(p);
									if(hom)var.incrementHomozygoteCount(p);
								}	
							}
							if(has)
							{
								//Add genotype
								SVDBGenotype geno = new SVDBGenotype(s.getUID(), 1);
								int ct = gt.countAlleleOccurrences(1);
								geno.addAllele(ct, sv.getPosition(), sv.getEndPosition());
								genolist.add(geno);
							}
						}
					}
					if(!genolist.isEmpty()) newgenos.put(var.getIntegerID(), genolist);
					var.adjustPopulationFrequency(indivCount);
					for(Population p : Population.values()) var.adjustPopulationFrequency(this.popIndivCount.get(p), p);

					//Update gene info!
					var.noteGenes(geneSet);
					
					//Update validation string!
					var.setValidationNotes("Not validated");
					
					bw.write(var.toDBRecord() + "\n");
				}
				bw.close();
			}
			else
			{
				
				String tempold = FileBuffer.generateTemporaryPath("svdb_addfam_matchedvar");
				String tempnew = FileBuffer.generateTemporaryPath("svdb_addfam_newvar");
				
				//Read from existing table, see if anything in new pool matches
				//If so, update population totals and genotypes, write to tempold
				//If not, make new var, write to tempnew
				BufferedReader br = new BufferedReader(new FileReader(tblPath));
				//BufferedWriter ntemp = new BufferedWriter(new FileWriter(tempnew));
				BufferedWriter otemp = new BufferedWriter(new FileWriter(tempold));
				String line = null;
				int stind = 0;
				double pLeeway = (double)leeway/1000.0;
				while((line = br.readLine()) != null)
				{
					if (line.isEmpty()) continue;
					if (line.startsWith("#")) {
						otemp.write(line + "\n");
						continue;
					}
					
					//Read line
					DBVariant var = DBVariant.getFromDBRecord(line, genomeBuild, geneSet);
					//Scan through list to see if matches any incoming variants
					int i = -1;
					int remi = -1;
					for(StructuralVariant sv : vlist)
					{
						i++;
						if (i < stind) continue; //Sorting allows us to skip checking over and over again
						//Compare
						if(var.svIsEquivalent(sv, pLeeway))
						{
							remi = i;
							break;
						}
					}
					//If so, remove the incoming variant from the list
					//If not, just copy the line to otemp
					if (remi >= 0)
					{
						//Found a hit
						StructuralVariant sv = vlist.remove(remi);
						stind = remi;
						//Add genotypes and update population info!
						List<SVDBGenotype> genolist = new LinkedList<SVDBGenotype>();
						for(DBSample s : samples)
						{
							Genotype gt = sv.getSampleGenotype(s.getName());
							if (gt != null)
							{
								//Check if has
								boolean has = gt.hasAllele(1);
								//Check if hom
								boolean hom = false;
								if(has) hom = gt.isHomozygous();
								//Note for all populations this person is in
								if(has)var.incrementTotalCount();
								if(hom)var.incrementHomozygoteCount();
								List<Population> plist = s.getPopulationFlags();
								if(plist != null && !plist.isEmpty())
								{
									for(Population p : plist)
									{
										if(has)var.incrementTotalCount(p);
										if(hom)var.incrementHomozygoteCount(p);
									}	
								}
								if(has)
								{
									//Add genotype
									SVDBGenotype geno = new SVDBGenotype(s.getUID(), 1);
									int ct = gt.countAlleleOccurrences(1);
									geno.addAllele(ct, sv.getPosition(), sv.getEndPosition());
									genolist.add(geno);
								}
							}
						}
						if(!genolist.isEmpty()) newgenos.put(var.getIntegerID(), genolist);
						var.adjustPopulationFrequency(indivCount);
						for(Population p : Population.values()) var.adjustPopulationFrequency(this.popIndivCount.get(p), p);
						otemp.write(var.toDBRecord() + "\n");
					}
					else
					{
						otemp.write(line + "\n");
					}
				}
				br.close();
				//ntemp.close();
				otemp.close();
				
				//Now write new variants!
				BufferedWriter ntemp = new BufferedWriter(new FileWriter(tempnew));
				//Only variants left in vlist should be new
				for(StructuralVariant sv: vlist)
				{
					String vname = t.name() + sv.getChromosomeName() + "-" + sv.getPosition() + "_" + OffsetDateTime.now().getNano();
					DBVariant var = DBVariant.getFromVariant(sv, vname);
					while(usedUIDs.contains(var.getIntegerID())) var.regenerateIntegerID();
					usedUIDs.add(var.getIntegerID());
					
					//Generate genotype info!
					
					List<SVDBGenotype> genolist = new LinkedList<SVDBGenotype>();
					for(DBSample s : samples)
					{
						Genotype gt = sv.getSampleGenotype(s.getName());
						if (gt != null)
						{
							//Check if has
							boolean has = gt.hasAllele(1);
							//Check if hom
							boolean hom = false;
							if(has) hom = gt.isHomozygous();
							//Note for all populations this person is in
							if(has)var.incrementTotalCount();
							if(hom)var.incrementHomozygoteCount();
							List<Population> plist = s.getPopulationFlags();
							if(plist != null && !plist.isEmpty())
							{
								for(Population p : plist)
								{
									if(has)var.incrementTotalCount(p);
									if(hom)var.incrementHomozygoteCount(p);
								}	
							}
							if(has)
							{
								//Add genotype
								SVDBGenotype geno = new SVDBGenotype(s.getUID(), 1);
								int ct = gt.countAlleleOccurrences(1);
								geno.addAllele(ct, sv.getPosition(), sv.getEndPosition());
								genolist.add(geno);
							}
						}
					}
					if(!genolist.isEmpty()) newgenos.put(var.getIntegerID(), genolist);
					var.adjustPopulationFrequency(indivCount);
					for(Population p : Population.values()) var.adjustPopulationFrequency(this.popIndivCount.get(p), p);

					//Update gene info!
					var.noteGenes(geneSet);
					
					//Update validation string!
					var.setValidationNotes("Not validated");
					
					ntemp.write(var.toDBRecord() + "\n");
				}
				ntemp.close();
				
				//Now merge the tables!
				combineTempTables(tempold, tempnew, tblPath);
				
			}
			//Index local table
			indexVariantTable(t);
		}
		
		//Master index var tables
		indexVariants();
		
		//Genotype dump!
		genoTable.addGenotypes(newgenos);
		
	}
	
	private void combineTempTables(String oldvars, String newvars, String targetPath) throws IOException
	{
		if (!FileBuffer.fileExists(oldvars))
		{
			Files.move(Paths.get(newvars), Paths.get(targetPath));
			return;
		}
		if (!FileBuffer.fileExists(newvars))
		{
			Files.move(Paths.get(oldvars), Paths.get(targetPath));
			return;
		}
		
		
		BufferedReader or = new BufferedReader(new FileReader(oldvars));
		BufferedReader nr = new BufferedReader(new FileReader(newvars));
		BufferedWriter bw = new BufferedWriter(new FileWriter(targetPath));
		
		//First, dump all header lines
		String oline = null;
		String nline = null;
		DBVariant ov = null;
		DBVariant nv = null;
		while((oline = or.readLine()) != null)
		{
			if (oline.startsWith("#")) bw.write(oline + "\n");
			else {
				ov = DBVariant.getFromDBRecord(oline, genomeBuild, geneSet);
				break;
			}
		}
		while((nline = nr.readLine()) != null)
		{
			if (nline.startsWith("#")) bw.write(nline + "\n");
			else {
				nv = DBVariant.getFromDBRecord(nline, genomeBuild, geneSet);
				break;
			}
		}
		
		boolean odone = (ov == null);
		boolean ndone = (nv == null);
		while(!(odone && ndone))
		{
			if(odone)
			{
				//Just write next new
				if (nv != null)
				{
					//A line leftover from last pass
					bw.write(nv.toDBRecord() + "\n");
					nv = null;
				}
				else
				{
					//No leftover lines, read next line
					nline = nr.readLine();
					if (nline == null) ndone = true;
					else bw.write(nline + "\n");	
				}
			}
			else if (ndone)
			{
				//Just write next old
				if (ov != null)
				{
					bw.write(ov.toDBRecord() + "\n");
					ov = null;
				}
				else
				{
					oline = or.readLine();
					if (oline == null) odone = true;
					else bw.write(oline + "\n");	
				}
			}
			else
			{
				//Compare
				//Neither should be null...
				int compare = ov.compareTo(nv);
				if(compare > 0)
				{
					//ov is greater, so nv goes next
					bw.write(nv.toDBRecord() + "\n");
					nline = nr.readLine();
					if(nline == null)
					{
						nv = null;
						ndone = true;
					}
					else nv = DBVariant.getFromDBRecord(nline, genomeBuild, geneSet);
				}
				else if (compare < 0)
				{
					//nv is greater, so ov goes next
					bw.write(ov.toDBRecord() + "\n");
					oline = or.readLine();
					if(oline == null)
					{
						ov = null;
						odone = true;
					}
					else ov = DBVariant.getFromDBRecord(oline, genomeBuild, geneSet);
				}
				else
				{
					//Write both
					bw.write(ov.toDBRecord() + "\n");
					bw.write(nv.toDBRecord() + "\n");
					
					nline = nr.readLine();
					if(nline == null)
					{
						nv = null;
						ndone = true;
					}
					else nv = DBVariant.getFromDBRecord(nline, genomeBuild, geneSet);
				
					oline = or.readLine();
					if(oline == null)
					{
						ov = null;
						odone = true;
					}
					else ov = DBVariant.getFromDBRecord(oline, genomeBuild, geneSet);
				}
			}
		}
		
		bw.close();
		or.close();
		nr.close();
		
	}
	
	public boolean updateFamily(Family fam, String newVCF) throws IOException, UnsupportedFileTypeException
	{
		if (fam == null) return false;
		if (newVCF == null || newVCF.isEmpty()) return false;
		
		//Remove all variants and genotypes for the family, then re-add
		
		if(!removeFamily(fam)) return false;
		if(!addFamily(newVCF, fam)) return false;
		
		return true;
	}
	
	public boolean removeFamily(Family fam) throws IOException
	{
		if (fam == null) return false;
		
		List<FamilyMember> members = fam.getAllFamilyMembers();
		Set<Integer> remVars = new TreeSet<Integer>();
		for(FamilyMember mem : members)
		{
			int uid = mem.getUID();
			List<Integer> varids = genoTable.removeSampleFromGenotypeTable(uid);
			remVars.addAll(varids);
			//removeVariants(varids);
		}
		
		removeVariants(remVars);
		return true;
	}
	
	private void removeVariants(Collection<Integer> varIDList) throws IOException
	{
		if(varIDList == null || varIDList.isEmpty()) return;
		
		SVType[] types = SVType.values();
		for(SVType t : types)
		{
			String tpath = this.getVariantTablePath(t);
			if (!FileBuffer.fileExists(tpath)) continue;
			
			String temppath = FileBuffer.generateTemporaryPath("svdbtemptbl");
			BufferedReader br = new BufferedReader(new FileReader(tpath));
			BufferedWriter temp = new BufferedWriter(new FileWriter(temppath));
			
			String line = null;
			while((line = br.readLine()) != null)
			{
				if(line.isEmpty()) continue;
				if (line.startsWith("#")) {
					temp.write(line + "\n");
					continue;
				}
				//Read the hex id, which is field 2
				String[] fields = line.split("\t");
				if (fields.length < 2) continue;
				try
				{
					int id = Integer.parseUnsignedInt(fields[1], 16);
					if(varIDList.contains(id)) continue; //This line will not be written
				}
				catch(NumberFormatException e) {continue;}
				temp.write(line + "\n");
			}
			
			temp.close();
			br.close();
			
			//Move temp back...
			Files.delete(Paths.get(tpath));
			Files.move(Paths.get(temppath), Paths.get(tpath));
			
			indexVariantTable(t);
		}
		
		indexVariants();
	}
	
	/* --- Sample Query --- */
	
	public Family getFamily(String familyName)
	{
		return familyMap.get(familyName);
	}
	
	/* --- Variant Query --- */
	
	private static double calculateUpperLimit(int hits, int total)
	{
		//http://sphweb.bumc.bu.edu/otlt/MPH-Modules/QuantCore/PH717_ConfidenceIntervals-OneSample/PH717_ConfidenceIntervals-OneSample5.html
		final double z = 1.96; //95% CI
		double pfreq = (double)hits/(double)total;
		
		double factor = pfreq * (1.0 - pfreq);
		factor /= (double)total;
		factor = Math.sqrt(factor);
		factor *= z;
		
		return pfreq + factor;
	}
	
	private boolean variantBelowPopThreshold(DBVariant var, double threshold)
	{
		//Fails if the frequency is too high for ANY population
		if (var == null) return false;
		
		//Total
		double maxFreq = calculateUpperLimit(var.getIndividualCount(), indivCount);
		if (maxFreq > threshold) return false;
		
		//By population
		for (Population p : Population.values())
		{
			int hits = var.getIndividualCount(p);
			int total = this.popIndivCount.get(p);
			maxFreq = calculateUpperLimit(hits, total);
			if (maxFreq > threshold) return false;
		}
		
		return true;
	}
	
	public List<DBVariant> queryForVariants(Collection<QueryCondition> conditions) throws IOException
	{
		return this.queryForVariants(conditions, SVType.allTypes());
	}
	
	public List<DBVariant> queryForVariants(Collection<QueryCondition> conditions, Collection<SVType> includeTypes) throws IOException
	{
		List<DBVariant> list = new LinkedList<DBVariant>();
		if(includeTypes == null) return list;
		
		for(SVType t : includeTypes)
		{
			String vtblpath = this.getVariantTablePath(t);
			if(!FileBuffer.fileExists(vtblpath)) continue;
			BufferedReader br = new BufferedReader(new FileReader(vtblpath));
			String line = null;
			while((line = br.readLine()) != null)
			{
				if (line.isEmpty()) continue;
				if (line.startsWith("#")) continue;
				DBVariant var = DBVariant.getFromDBRecord(line, genomeBuild, geneSet);
				if(var == null) continue;
				boolean passes = true;
				if(conditions != null)
				{
					for(QueryCondition c : conditions)
					{
						if(!c.passes(var))
						{
							passes = false;
							break;
						}
					}
				}
				if(passes) list.add(var);
			}
			br.close();
		}
		
		return list;
	}
	
	private Map<Integer, List<SVDBGenotype>> getFamilyGenotypes(Family fam)
	{
		if (fam == null) return null;
		
		Map<Integer, List<SVDBGenotype>> map = new TreeMap<Integer, List<SVDBGenotype>>();
		List<FamilyMember> members = fam.getAllFamilyMembers();
		
		for (FamilyMember m : members)
		{
			Map<Integer, SVDBGenotype> indivMap = genoTable.getGenotypesForSample(m.getUID());
			if (indivMap == null) continue;
			List<Integer> klist = new ArrayList<Integer>(indivMap.size() + 1);
			klist.addAll(indivMap.keySet());
			for(Integer k : klist)
			{
				List<SVDBGenotype> l = map.get(k);
				if (l == null)
				{
					l = new LinkedList<SVDBGenotype>();
					map.put(k, l);
				}
				l.add(indivMap.get(k));
			}
		}
		
		
		return map;
	}
	
	public Map<DBVariant, List<SVDBGenotype>> getFamilyVariants(Family fam) throws IOException
	{
		if (fam == null) return null;
		Map<DBVariant, List<SVDBGenotype>> map = new HashMap<DBVariant, List<SVDBGenotype>>();
		Map<Integer, List<SVDBGenotype>> idmap = getFamilyGenotypes(fam);
		
		List<Integer> klist = new ArrayList<Integer>(idmap.size() + 1);
		klist.addAll(idmap.keySet());
		for(Integer varID : klist)
		{
			DBVariant var = this.getVariant(varID);
			if(var == null) continue;
			map.put(var, idmap.get(varID));
		}
		
		return map;
	}
	
	public List<Candidate> getCandidates(Map<DBVariant, List<SVDBGenotype>> varset, Family fam)
	{
		if (fam == null) return null;
		if (varset == null) return null;
		
		List<DBVariant> keys = new ArrayList<DBVariant>(varset.size() + 1);
		keys.addAll(varset.keySet());
		
		List<FamilyMember> members = fam.getAllFamilyMembers();
		
		List<Variant> vlist = new ArrayList<Variant>(varset.size() + 1);
		for(DBVariant var : keys)
		{
			StructuralVariant sv = var.toStructuralVariant();
			//Add genotypes
			List<SVDBGenotype> glist = varset.get(var);
			for(FamilyMember m : members)
			{
				//See if there is a genotype for this person
				SVDBGenotype geno = null;
				if (glist != null && !glist.isEmpty())
				{
					for (SVDBGenotype g : glist)
					{
						if (g.getIndividualUID() == m.getUID())
						{
							geno = g;
							break;
						}
					}
				}
				
				Genotype gt = new Genotype();
				if(geno != null)
				{
					//Copy allele data
					Collection<SVDBAllele> alleles = geno.getAlleles();
					List<Integer> alist = new LinkedList<Integer>();
					int aind = 1;
					for(SVDBAllele a : alleles)
					{
						for(int i = 0; i < a.getAlleleCount(); i++) alist.add(aind);
						aind++;
					}
					while(alist.size() < 2)
					{
						//Add some refs
						alist.add(0);
					}
					//Convert to an array and dump in genotype
					int acount = alist.size();
					int[] allarr = new int[acount];
					int i = 0;
					for(Integer a : alist)
					{
						allarr[i] = a;
						i++;
					}
					gt.setAlleles(allarr);
				}
				else
				{
					//Make ref/ref
					gt.setAlleles("0/0");
				}
				sv.addGenotype(m.getName(), gt);
			}
			vlist.add(sv);
		}
		
		//Re-genotype XY variants
		fam.adjustSexChromGenotypes(vlist, new TwoSexChromSegModel(genomeBuild.getContig("X"), genomeBuild.getContig("Y"), genomeBuild));
		
		//Change to variant pool
		VariantPool pool = new VariantPool(members.size());
		pool.addVariants(vlist);
		
		//Candidates
		List<Candidate> clist = Inheritor.getCandidates(pool, fam, geneSet);
		
		return clist;
	}
	
	private boolean prioritizeCandidate(Candidate c, DBVariant dbv, List<Individual> affected)
	{
		if (c == null) return false;
		if (affected == null) return false;
		
		Variant v = c.getVariant();
		if (v == null) return false;
		if(!(v instanceof StructuralVariant)) return false;
		StructuralVariant sv = (StructuralVariant)v;
		SVType t = sv.getType();
		if (t == null) return false;
		if (t == SVType.TRA) return false;
		if (t == SVType.BND) return false;
		if (t == SVType.INV) return false;
		
		//Must segregate or be halfhet for at least one affected
		boolean pass = false;
		for(Individual a : affected)
		{
			Inheritance i = c.getInheritancePattern(a);
			if (i != null && i != Inheritance.UNRESOLVED)
			{
				pass = true;
				break;
			}
		}
		if(!pass) return false;
		
		//Position Effect
		GeneFunc eff = c.getPositionEffect();
		if (eff == null) return false;
		//Only prioritize rare vars for some effects
		/*
		 * Tier 1:
		 * 	Exonic, Splicing
		 * Tier 2:
		 * 	UTR5, ncRNA, Upstream
		 * Tier 3:
		 * 	UTR3, Intronic
		 * Tier 4:
		 * 	Downstream, Intergenic
		 */
		switch(eff)
		{
		case DOWNSTREAM:
			return variantBelowPopThreshold(dbv, POPFREQ_CUTOFF_G4);
		case EXONIC:
			return true; //Always prioritize
		case INTERGENIC:
			return variantBelowPopThreshold(dbv, POPFREQ_CUTOFF_G4);
		case INTRONIC:
			return variantBelowPopThreshold(dbv, POPFREQ_CUTOFF_G3);
		case NCRNA:
			return variantBelowPopThreshold(dbv, POPFREQ_CUTOFF_G2);
		case SPLICING:
			return true; //Always prioritize
		case UPSTREAM:
			return variantBelowPopThreshold(dbv, POPFREQ_CUTOFF_G2);
		case UTR3:
			return variantBelowPopThreshold(dbv, POPFREQ_CUTOFF_G3);
		case UTR5:
			return variantBelowPopThreshold(dbv, POPFREQ_CUTOFF_G2);
		}
		
		return true;
	}
	
	public void writeFamilyTable(String famname, String outputDir) throws IOException
	{
		//Prioritize variants that are...
		//	Exonic and splicing del/dup/ins that are half-hets or segregating (ANY)
		//	Rare del/dup/ins (HH or seg) that are ncRNA, UTR, intronic, or upstream/downstream
		// The only inv or tra that should be prioritized are those called as comp het partners
		
		/*
		 * Table Fields ----- (CSV)
		 * 		Var Name
		 * 		Var ID (Hex)
		 * 		Chrom 1
		 * 		Chrom 2 (Empty if not TRA)
		 * 		Pos (Range)
		 * 		End (Range)
		 * 		Insertion Sequence (Empty if not INS/INS:ME)
		 * 		SVType
		 * 		Size
		 * 		PosEff
		 * 		Gene Name
		 * 		Transcript ID
		 * 		OMIM Info (Blank if none found)
		 * 		DECIPHER Score (Blank if none found)
		 * 		UDP Cohort Total Indivs with variant
		 * 		UDP Cohort Hom Count
		 * 		UDP Cohort Pop Freq
		 * 		UDP Cohort Gene Hit Count (Variant)
		 * 		UDP Cohort Gene Hit Count (Indivs)
		 * 		UDP Cohort Exon Hit Count (Variant)
		 * 		UDP Cohort Exon Hit Count (Indivs)
		 * 		[Population] Allele Count
		 * 		[Population] Hom Count
		 * 		[Population] Pop Freq
		 * 			...
		 * 		[Seg Patterns]
		 * 			Per affected...
		 * 		CompHet Partners (String IDs(HexID))
		 * 		[Genotypes] #Occurrences of allele
		 * 		[Genotypes]	Start Position
		 * 		[Genotypes]	End Position
		 * 			...
		 */
		
		//Check output directory
		if(!FileBuffer.directoryExists(outputDir)) Files.createDirectories(Paths.get(outputDir));
		
		//Write db Summary
		String dbInfoPath = outputDir + File.separator + "databaseInfo.txt";
		BufferedWriter bw = new BufferedWriter(new FileWriter(dbInfoPath));
		bw.write("Database Name: " + dbName + "\n");
		GenomeBuildUID gbid = genomeBuild.getUIDEnum();
		if (gbid != null) bw.write("Genome Build: " + gbid.toString() + "\n");
		else bw.write("Genome Build: " + genomeBuild.getBuildName() + " (" + genomeBuild.getSpeciesID() + ")" + "\n");
		bw.write("Leeway Value: " + leeway + "\n");
		bw.write("Total Families: " + familyMap.size() + "\n");
		bw.write("Total Individuals: " + indivCount + "\n");
		bw.write("Individuals by Population: \n");
		for(Population p : Population.values())
		{
			int ct = this.popIndivCount.get(p);
			bw.write("\t" + p.getShortString() + " (" + p.toString() + "): " + ct + "\n");
		}
		bw.close();
		
		//Write family summary
		Family f = familyMap.get(famname);
		if (f == null)
		{
			System.err.println("SVDatabase.writeFamilyTable || ERROR: Family \"" + famname + "\" not found in database!");
			return;
		}
		String faminfopath = outputDir + File.separator + "family_" + famname + ".txt";
		bw = new BufferedWriter(new FileWriter(faminfopath));
		bw.write("Family Name: " + f.getFamilyName() + "\n");
		Individual proband = f.getProband();
		if(proband != null) bw.write("Family Proband: " + proband.getName() + "\n");
		bw.write("Number of Members: " + f.countMembers() + "\n");
		List<FamilyMember> members = f.getAllFamilyMembers();
		bw.write("SampleName\tSex\tAffectedStatus\tRelationToProband\tPopulationTags\tMother\tFather\n");
		for(FamilyMember m : members)
		{
			bw.write(m.getName() + "\t");
			bw.write(m.getSex().name() + "\t");
			bw.write(m.getAffectedStatus().name() + "\t");
			if(proband != null) bw.write(f.getRelationshipString_ENG(m) + "\t");
			else bw.write("unknown\t");
			Collection<Population> ptags = m.getPopulationTags();
			boolean first = true;
			for (Population p : ptags)
			{
				if (!first) bw.write(";");
				bw.write(p.getShortString());
				first = false;
			}
			bw.write("\t");
			Individual mom = m.getMother();
			if (mom != null) bw.write(mom.getName() + "\t");
			else bw.write("unknown\t");
			Individual dad = m.getFather();
			if (dad != null) bw.write(dad.getName() + "\t");
			else bw.write("unknown\t");
		}
		bw.close();
		
		//Now, retrieve variants....
		Map<DBVariant, List<SVDBGenotype>> vmap = this.getFamilyVariants(f);
		List<Candidate> rawlist = this.getCandidates(vmap, f);
		
		//Dump allele 0 candidates
		List<Candidate> clist = new LinkedList<Candidate>();
		for(Candidate c : rawlist)
		{
			if (c.getAllele() != 0) clist.add(c);
		}
		
		//Remap the variants for faster lookup
		Map<Integer, DBVariant> idmap = new TreeMap<Integer, DBVariant>();
		Map<Integer, List<SVDBGenotype>> gmap = new TreeMap<Integer, List<SVDBGenotype>>();
		
		List<DBVariant> alist = new ArrayList<DBVariant>(vmap.size() + 1);
		alist.addAll(vmap.keySet());
		for(DBVariant v : alist)
		{
			int id = v.getIntegerID();
			idmap.put(id, v);
			gmap.put(id, vmap.get(v));
		}
		
		List<Individual> aff = f.getAllAffected();
		
		//Prepare table headers
		StringBuilder sb = new StringBuilder(2048);
		sb.append("VariantName\t");
		sb.append("VariantUID\t");
		sb.append("Chrom\t");
		sb.append("Chrom2(TRA)\t");
		sb.append("StartPosition\t");
		sb.append("EndPosition\t");
		sb.append("InsertionSequence(INS)\t");
		sb.append("SVType\t");
		sb.append("SVLen\t");
		sb.append("PositionEffect\t");
		sb.append("Gene\t");
		sb.append("TranscriptID\t");
		sb.append("OMIM_Morbidity\t");
		sb.append("DECIPHER\t");
		sb.append("Cohort_Count\t");
		sb.append("Cohort_HomCount\t");
		sb.append("Cohort_Freq\t");
		sb.append("Cohort_GeneHits_Var\t");
		sb.append("Cohort_GeneHits_Indiv\t");
		sb.append("Cohort_GeneHitsExon_Var\t");
		sb.append("Cohort_GeneHitsExon_Indiv\t");
		Population[] pall = Population.values();
		for(Population p : pall)
		{
			sb.append(p.getShortString() + "_Count\t");
			sb.append(p.getShortString() + "_HomCount\t");
			sb.append(p.getShortString() + "_Freq\t");
		}
		for(Individual a : aff) sb.append("SEG_" + a.getName() + "\t");
		sb.append("COMPHET_PARTNERS\t");
		for(FamilyMember m : members)
		{
			sb.append(m.getName() + "_GENO_ACNT\t");
			sb.append(m.getName() + "_GENO_START\t");
			sb.append(m.getName() + "_GENO_END\t");
		}
		String hraw = sb.toString();
		String header = hraw.substring(0, hraw.length() - 1); //Chop off last tab

		//Open output streams
		String fullPath = outputDir + File.separator + f.getFamilyName() + "_fullTable.tsv";
		String priPath = outputDir + File.separator + f.getFamilyName() + "_prioritized.tsv";
		BufferedWriter ftbl = new BufferedWriter(new FileWriter(fullPath));
		BufferedWriter ptbl = new BufferedWriter(new FileWriter(priPath));
		ftbl.write(header + "\n");
		ptbl.write(header + "\n");
		
		//Go down list of candidates
		Collections.sort(clist);
		for(Candidate c : clist)
		{
			//Prepare record
			Variant v = c.getVariant();
			if(v == null)
			{
				System.err.println("SVDatabase.writeFamilyTable || ERROR: Candidate lacks linked variant! Skipping...");
				continue;
			}
			
			int id = v.getSingleIntInfoEntry(DBVariant.ID_INFO_KEY);
			DBVariant dbv = idmap.get(id);
			if(dbv == null)
			{
				System.err.println("SVDatabase.writeFamilyTable || ERROR: Variant database UID is invalid: " + Integer.toHexString(id));
				continue;
			}
			List<SVDBGenotype> dbgenos = gmap.get(id);
			
			
			sb = new StringBuilder(4096);
			//IDs
			sb.append(dbv.getName() + "\t");
			sb.append(Integer.toHexString(dbv.getIntegerID()) + "\t");
			//Chroms
			Contig c1 = dbv.getChrom();
			if (c1 != null) sb.append(c1.getUDPName() + "\t");
			else sb.append("*\t");
			if(dbv.getType() == SVType.TRA)
			{
				Contig c2 = dbv.getEndChrom();
				if (c2 == null)
				{
					if (c1 != null) sb.append(c1.getUDPName() + "\t");
					else sb.append("*\t");
				}
				else sb.append(c2.getUDPName() + "\t");
			}
			else sb.append("\t"); //Empty
			//Positions
			Interval stpos = dbv.getStartPosition();
			Interval edpos = dbv.getEndPosition();
			int st1 = stpos.getStart();
			int st2 = stpos.getEnd();
			int ed1 = edpos.getStart();
			int ed2 = edpos.getEnd();
			if (st1 == st2) sb.append(st1 + "\t");
			else sb.append(st1 + "-" + st2 + "\t");
			if (ed1 == ed2) sb.append(ed1 + "\t");
			else sb.append(ed1 + "-" + ed2 + "\t");
			//Insertion Sequence
			SVType t = dbv.getType();
			if (t == SVType.INS || t == SVType.INSME)
			{
				String iseq = dbv.getAltAlleleString();
				if (iseq == null || iseq.isEmpty()) sb.append("<UNKNOWN>\t");
				else sb.append(iseq + "\t");
			}
			else sb.append("<N/A>\t");
			//Other SV Data
			sb.append(t.getString() + "\t");
			int len = ed2 - st1;
			sb.append(len + "bp\t");
			//Gene stuff
			GeneFunc poseff = c.getPositionEffect();
			sb.append(poseff.toString() + "\t");
			Gene gene = c.getGene();
			if (gene == null) sb.append("<N/A>\t<N/A>\t");
			else
			{
				sb.append(gene.getName() + "\t");
				sb.append(gene.getID() + "\t");
			}
			//TODO: Insert OMIM annotation here when ready!
			sb.append("(TBI)\t");
			//TODO: Insert DECIPHER annotation here when ready!
			sb.append("(TBI)\t");
			//Cohort Counts
			sb.append(dbv.getIndividualCount() + "\t");
			sb.append(dbv.getHomozygoteCount() + "\t");
			sb.append(dbv.getCohortFreq() + "\t");
			//Cohort Gene Hits
			if (gene != null)
			{
				List<DBVariant> vhits = this.getVariantsInGene(gene);	
				if (vhits == null) sb.append("0\t0\t0\t0\t");
				else
				{
					sb.append(vhits.size() + "\t");
					int tot = 0;
					for (DBVariant h : vhits) tot += h.getIndividualCount();
					sb.append(tot + "\t");
					//Isolate exonic hits
					int vcount = 0;
					tot = 0;
					for (DBVariant h : vhits)
					{
						GeneFunc eff = null;
						if (h.getType() == SVType.INV)
						{
							//Check ends
							eff = gene.getRelativeRegionLocationEffect(h.getStartPosition().getStart(), h.getStartPosition().getEnd());
							GeneFunc e2 = gene.getRelativeRegionLocationEffect(h.getEndPosition().getStart(), h.getEndPosition().getEnd());
							if(e2.getPriority() < eff.getPriority()) eff = e2;
						}
						else if (h.getType() == SVType.INS || h.getType() == SVType.INSME)
						{
							//Only check the position
							eff = gene.getRelativeLocationEffect(h.getStartPosition().getCenter());
						}
						else if(h.getType() == SVType.TRA)
						{
							//Have to check which chrom it is
							eff = GeneFunc.INTERGENIC;
							Contig hc1 = h.getChrom();
							if (gene.getChromosome().equals(hc1))
							{
								eff = gene.getRelativeRegionLocationEffect(h.getStartPosition().getStart(), h.getStartPosition().getEnd());
							}
							Contig hc2 = h.getEndChrom();
							if (hc2 == null) hc2 = hc1;
							if (gene.getChromosome().equals(hc2))
							{
								GeneFunc e2 = gene.getRelativeRegionLocationEffect(h.getEndPosition().getStart(), h.getEndPosition().getEnd());
								if(e2.getPriority() < eff.getPriority()) eff = e2;
							}
						}
						else
						{
							eff = gene.getRelativeRegionLocationEffect(h.getStartPosition().getStart(), h.getEndPosition().getEnd());
						}
					
						//Determine if exonic
						if (eff == GeneFunc.EXONIC)
						{
							vcount++;
							tot += h.getIndividualCount();
						}
					}//End for
					sb.append(vcount + "\t");
					sb.append(tot + "\t");
				}
			}
			else sb.append("<N/A>\t<N/A>\t<N/A>\t<N/A>\t");
			//Population Counts
			for(Population p : pall)
			{
				sb.append(dbv.getIndividualCount(p) + "\t");
				sb.append(dbv.getHomozygoteCount(p) + "\t");
				sb.append(dbv.getCohortFreq(p) + "\t");
			}
			//Segregation
			for(Individual ai : aff)
			{
				Inheritance inh = c.getInheritancePattern(ai);
				if (inh != null) sb.append(inh.toString() + "\t");
				else sb.append("Undetermined\t");
			}
			//Comphet Partners
			Collection<Variant> plist = c.getAllPartnerVariants();
			if (plist != null && !plist.isEmpty())
			{
				boolean first = true;
				for(Variant p : plist)
				{
					if (!first) sb.append(";");
					sb.append(p.getVarID());
					int pid = p.getSingleIntInfoEntry(DBVariant.ID_INFO_KEY);
					sb.append("(" + Integer.toHexString(pid) + ")");
					first = false;
				}
				sb.append("\t");
			}
			else sb.append("<N/A>\t");
			//Genotypes
			for(FamilyMember mem : members)
			{
				SVDBGenotype mg = null;
				//Search
				for(SVDBGenotype sg : dbgenos)
				{
					if (sg.getIndividualUID() == mem.getUID())
					{
						mg = sg;
						break;
					}
				}
				if (mg == null) sb.append("0\t<N/A>\t<N/A>\t");
				else
				{
					//Just first allele...
					Collection<SVDBAllele> gall = mg.getAlleles();
					for(SVDBAllele a : gall)
					{
						sb.append(a.getAlleleCount() + "\t");
						sb.append(a.getAllele().getStart() + "\t");
						sb.append(a.getAllele().getEnd() + "\t");
					}
					if (gall.size() > 1)
					{
						System.err.println("SVDatabase.writeFamilyTable || WARNING: Family Member " + mem.getName() + " has multiple unique alleles for variant " + dbv.getName() + ". Only first allele info will be tabled.");
					}
				}
			}
			int rlen = sb.length();
			//Trim last tab
			sb.deleteCharAt(rlen-1);

			//Write!
			String recLine = sb.toString();
			ftbl.write(recLine + "\n");
			//Prioritize variant?
			if(prioritizeCandidate(c, dbv, aff)) ptbl.write(recLine + "\n");
		}
		
		//Close streams
		ftbl.close();
		ptbl.close();
	}
	
	public void writeBED(List<DBVariant> vars, String path)
	{
		//TODO: Write
	}
	
	/* --- Management of Totals --- */
	
	public boolean updatePopulationCounts()
	{
		updatePopulationTotals();
		
		//Update for each variant
		try 
		{
			String temppath = FileBuffer.generateTemporaryPath("svdb_tempvartbl");	
			
			SVType[] types = SVType.values();
			for(SVType t : types)
			{
				String tpath = this.getVariantTablePath(t);
				//Check if exists...
				if (!FileBuffer.fileExists(tpath)) continue;
				
				//Else, we need to update counts!
				BufferedReader in = new BufferedReader(new FileReader(tpath));
				BufferedWriter temp = new BufferedWriter(new FileWriter(temppath));
				
				String line = null;
				while((line = in.readLine()) != null)
				{
					if (line.isEmpty()) continue;
					if (line.charAt(0) == '#') {
						temp.write(line + "\n");
						continue;
					}
					DBVariant var = DBVariant.getFromDBRecord(line, genomeBuild, geneSet);
					List<SVDBGenotype> genos = genoTable.getGenotypesForVariant(var.getIntegerID());
					if (genos != null && !genos.isEmpty())
					{
						//TODO: X and Y??? Hemizygote counts???
						int ngenos = genos.size();
						//Set total
						var.setTotalCount(ngenos);
						var.adjustPopulationFrequency(indivCount);
						//Set homozygotes
						int homcount = 0;
						for(SVDBGenotype g : genos)
						{
							for(SVDBAllele a : g.getAlleles())
							{
								if (a.getAlleleCount() > 1) {
									homcount++;
									break;
								}
							}
						}
						var.setHomozygoteCount(homcount);
						//Set by population
						Population[] allpop = Population.values();
						for(Population p : allpop)
						{
							int tcount = 0;
							homcount = 0;
							for(SVDBGenotype g : genos)
							{
								DBSample s = getSample(g.getIndividualUID());
								if (s == null) continue;
								if(s.isInPopulation(p))
								{
									tcount++;
									for(SVDBAllele a : g.getAlleles())
									{
										if (a.getAlleleCount() > 1) {
											homcount++;
											break;
										}
									}
								}
							}
							var.setTotalCount(tcount, p);
							var.adjustPopulationFrequency(this.popIndivCount.get(p), p);
							var.setHomozygoteCount(homcount, p);
						}
						//Write var back to table
						temp.write(var.toDBRecord() + "\n");
					}
					//Otherwise the var doesn't get copied to the new file!
				}
				
				temp.close();
				in.close();
				
				//Copy back & Reindex
				Files.move(Paths.get(temppath), Paths.get(tpath));
				indexVariantTable(t);
			}
			
		}
		catch(IOException e)
		{
			e.printStackTrace();
			return false;
		}
		return true;
	}
	
	public void updatePopulationTotals()
	{
		indivCount = sampleMap.size();
		Collection<DBSample> allsamps = sampleMap.values();
		resetPopulationCountMap();
		
		for(DBSample s : allsamps)
		{
			List<Population> ptags = s.getPopulationFlags();
			if(ptags == null) continue;
			for(Population p : ptags)
			{
				Integer i = popIndivCount.get(p);
				if (i == null)
				{
					popIndivCount.put(p, 0);
					i = 0;
				}
				popIndivCount.put(p, i+1);
			}
		}
	}
	
	/* --- Sample Information --- */
	
	public DBSample getSample(String sampleName)
	{
		return sampleMap.get(sampleName);
	}
	
	private void mapSampleIDs()
	{
		sampleIDMap = new TreeMap<Integer, DBSample>();
		for(DBSample s : sampleMap.values())
		{
			sampleIDMap.put(s.getUID(), s);
		}
	}
	
	public DBSample getSample(int UID)
	{
		if (sampleIDMap == null) mapSampleIDs();
		return sampleIDMap.get(UID);
	}
	
	public void updateFamilyName(String oldname, String newname)
	{
		Family f = familyMap.remove(oldname);
		if (f == null) return;
		f.setFamilyName(newname);
		familyMap.put(newname, f);
		
		//Update individuals
		List<FamilyMember> fmlist = f.getAllFamilyMembers();
		for(FamilyMember m : fmlist)
		{
			DBSample s = sampleMap.get(m.getName());
			s.setFamilyName(newname);
		}
	}
	
	public void updateSampleName(String oldname, String newname)
	{
		DBSample s = sampleMap.remove(oldname);
		if (s == null) return;
		s.setSampleName(newname);
		
		//Update in family, if there is one
		Family f = familyMap.get(s.getFamilyName());
		if (f == null) return;
		FamilyMember mem = f.getMemberByName(oldname);
		if (mem == null) return;
		f.changeSampleName(oldname, newname);
	}
	
	public void clearSampleEthnicities(String sampleName)
	{
		DBSample s = sampleMap.get(sampleName);
		s.clearPopulationFlags();
		
		Family f = familyMap.get(s.getFamilyName());
		if (f == null) return;
		FamilyMember mem = f.getMemberByName(sampleName);
		if (mem == null) return;
		mem.clearPopulationTags();
		
		updatePopulationCounts();
	}
	
	public void addSampleEthnicities(String sampleName, Collection<Population> pflags)
	{
		DBSample s = sampleMap.get(sampleName);
		
		Family f = familyMap.get(s.getFamilyName());
		if (f == null) return;
		FamilyMember mem = f.getMemberByName(sampleName);
		if (mem == null) return;
		for(Population p : pflags)
		{
			s.addPopulationFlag(p);
			mem.addPopulationTag(p);
		}
		
		updatePopulationCounts();
	}

}
