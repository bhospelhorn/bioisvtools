package hospelhornbg_svdb;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.time.OffsetDateTime;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import hospelhornbg_bioinformatics.Genotype;
import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.VCF;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.GeneSet;
import hospelhornbg_genomeBuild.GenomeBuild;
import hospelhornbg_segregation.Family;
import hospelhornbg_segregation.FamilyMember;
import hospelhornbg_segregation.Population;
import hospelhornbg_svdb.SVDBGenotype.SVDBAllele;
import waffleoRai_Utils.BinFieldSize;
import waffleoRai_Utils.CompositeBuffer;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;
import waffleoRai_Utils.SerializedString;
import waffleoRai_Utils.StreamBuffer;

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
	
	/* --- Instance Variables --- */
	
	private String dirPath;
	private String dbName;
	private GenomeBuild genomeBuild;
	private GeneSet geneSet;
	
	private int bpLeeway;
	
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
		svdb.bpLeeway = file.intFromFile(cpos);
		
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
		svdb.bpLeeway = leeway;
		
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
	
	public void indexVariants()
	{
		//TODO: Write
	}
	
	public void writeInfoFile() throws IOException
	{
		String path = this.getDBInfoFilePath();
		String bname = genomeBuild.getBuildName();
		FileBuffer outfile = new FileBuffer(8 + bname.length() + 3 + 4, true);
		outfile.printASCIIToFile(DBINFO_MAGIC);
		outfile.addVariableLengthString(bname, BinFieldSize.WORD, 2);
		outfile.addToFile(bpLeeway);
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
	
	/* --- Variant Addition --- */
	
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
						if(var.svIsEquivalent(sv, bpLeeway))
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
	
	private void combineTempTables(String oldvars, String newvars, String targetPath)
	{
		//TODO: Write
	}
	
	/* --- Variant Query --- */
	
	public List<String> queryForRecords(Collection<QueryCondition> conditions)
	{
		//TODO: Write
		return null;
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
