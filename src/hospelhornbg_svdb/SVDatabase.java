package hospelhornbg_svdb;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import hospelhornbg_bioinformatics.Genotype;
import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_genomeBuild.Gene;
import hospelhornbg_genomeBuild.GeneFunc;
import hospelhornbg_genomeBuild.GeneSet;
import hospelhornbg_genomeBuild.GenomeBuild;
import hospelhornbg_genomeBuild.GenomeBuildUID;
import hospelhornbg_genomeBuild.OMIMGeneMapImporter;
import hospelhornbg_segregation.Candidate;
import hospelhornbg_segregation.Family;
import hospelhornbg_segregation.FamilyMember;
import hospelhornbg_segregation.Individual;
import hospelhornbg_segregation.Inheritance;
import hospelhornbg_segregation.Inheritor;
import hospelhornbg_segregation.Population;
import hospelhornbg_svdb.DBVariantTable.GeneHitCounter;
import hospelhornbg_svdb.SVDBGenotype.SVDBAllele;
import waffleoRai_Utils.BinFieldSize;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.SerializedString;

public class SVDatabase {
	
	/* ----- Constants ----- */
	
	public static final double POPFREQ_CUTOFF_G2 = 0.05;
	public static final double POPFREQ_CUTOFF_G3 = 0.02;
	public static final double POPFREQ_CUTOFF_G4 = 0.01;
	
	public static final String SETTINGS_FILE = "SVDB_settings.bin";
	public static final int SETTINGS_VERSION = 2;
	public static final String SETTINGS_MAGIC = "svDB";
	
	public static final int VARTBL_TYPE_DEFO = 0;
	public static final int VARTBL_TYPE_SQL = 1;
	
	/* ----- Instance Variables ----- */
	
	private DBSampleTable sampleTable;
	private VariantTable variantTable;
	
	private String directory;
	private String dbName;
	
	private int mergeFactor;
	private GenomeBuild genome;
	private GeneSet genes;
	
	private String omim_table_path;
	
	private int varTableType;
	private String sqlURL;
	private String sqlUser;
	private String sqlPassword;
	
	/* ----- Private Constructors ----- */
	
	private SVDatabase(String dir, boolean indexLimiter) throws IOException, SQLException
	{
		//Load existing
		directory = dir;
		loadSettings();
		
		//Load tables
		sampleTable = new DBSampleTable(directory);
		//variantTable = new DBVariantTable(genome, genes, directory, indexLimiter);
	
		switch(varTableType)
		{
		case VARTBL_TYPE_DEFO:
			variantTable = new DBVariantTable(genome, genes, directory, indexLimiter);
			break;
		case VARTBL_TYPE_SQL:
			variantTable = new SQLVariantTable(sqlURL, sqlUser, sqlPassword, genome, genes);
			break;
		}
		
	}
	
	private SVDatabase(String dir, String name, int leewayValue, GenomeBuild gb, boolean indexLimiter) throws IOException
	{
		//New DB
		directory = dir;
		dbName = name;
		mergeFactor = leewayValue;
		genome = gb;
		genes = GeneSet.loadRefGene(genome);
		
		sampleTable = new DBSampleTable(directory);
		variantTable = new DBVariantTable(genome, genes, directory, indexLimiter);
		varTableType = VARTBL_TYPE_DEFO;
		
		saveDatabase();
	}
	
	private SVDatabase(String dir, String name, int leewayValue, GenomeBuild gb, String sqlPath) throws IOException, SQLException
	{
		//New DB
		directory = dir;
		dbName = name;
		mergeFactor = leewayValue;
		genome = gb;
		genes = GeneSet.loadRefGene(genome);
		
		sqlURL = sqlPath;
		sqlUser = "svdb_bioi";
		sqlPassword = "bIoI%pAsS";
		
		sampleTable = new DBSampleTable(directory);
		variantTable = new SQLVariantTable(sqlURL, sqlUser, sqlPassword, genome, genes);
		varTableType = VARTBL_TYPE_SQL;
		
		saveDatabase();
	}
	
	private void loadSettings() throws IOException
	{
		//Settings bin format:
		//	Magic [4] "svDB"
		// 	Version [4]
		//	GB UID[4]
		//	Merge Factor[4]
		//	Variant Table Type [2]
		//	DB Name [VLS 2x2]
		// 	OMIM Table Path [VLS 2x2]
		//	SQL DB Path [VLS 2x2]
		//	SQL DB Username [VLS 2x2]
		//	SQL DB Password [VLS 2x2]
		
		String settings_path = directory + File.separator + SETTINGS_FILE;
		FileBuffer settings = new FileBuffer(settings_path, true);
		
		long cpos = 8; //Skip magic and version
		int gbuid = settings.intFromFile(cpos); cpos += 4;
		mergeFactor = settings.intFromFile(cpos); cpos += 4;
		
		varTableType = Short.toUnsignedInt(settings.shortFromFile(cpos)); cpos += 2;
		
		SerializedString ss = settings.readVariableLengthString(cpos, BinFieldSize.WORD, 2);
		dbName = ss.getString();
		cpos += ss.getSizeOnDisk();
		
		ss = settings.readVariableLengthString(cpos, BinFieldSize.WORD, 2);	
		omim_table_path = ss.getString();
		cpos += ss.getSizeOnDisk();
		if(omim_table_path.equals("null")) omim_table_path = null;
		
		//Load genome build...
		GenomeBuildUID en = GenomeBuildUID.getByID(gbuid);
		if (en != null)
		{
			genome = GenomeBuild.loadStandardBuild(en.getName());
			genes = GeneSet.loadRefGene(genome);
		}
		
		//Get the SQL stuff if needed.
		if(varTableType == VARTBL_TYPE_SQL)
		{
			ss = settings.readVariableLengthString(cpos, BinFieldSize.WORD, 2);
			sqlURL = ss.getString();
			cpos += ss.getSizeOnDisk();
			
			ss = settings.readVariableLengthString(cpos, BinFieldSize.WORD, 2);
			sqlUser = ss.getString();
			cpos += ss.getSizeOnDisk();
			
			ss = settings.readVariableLengthString(cpos, BinFieldSize.WORD, 2);
			sqlPassword = ss.getString();
		}
		
		
	}
	
	public void saveDatabase() throws IOException
	{
		//Settings bin format:
		//	Magic [4] "svDB"
		// 	Version [4]
		//	GB UID[4]
		//	Merge Factor[4]
		//	Variant Table Type [2]
		//	DB Name [VLS 2x2]
		// 	OMIM Table Path [VLS 2x2]
		//	SQL DB Path [VLS 2x2]
		//	SQL DB Username [VLS 2x2]
		//	SQL DB Password [VLS 2x2]
		
		//Write settings
		String settings_path = directory + File.separator + SETTINGS_FILE;
		int sz = 18;
		sz += this.dbName.length() + 4;
		if(this.omim_table_path != null) sz += this.omim_table_path.length() + 4;
		else sz += 8;
		if(this.sqlURL != null) sz += this.sqlURL.length() + 4;
		if(this.sqlUser != null) sz += this.sqlUser.length() + 4;
		if(this.sqlPassword != null) sz += this.sqlPassword.length() + 4;
		FileBuffer file = new FileBuffer(sz, true);
		file.printASCIIToFile(SETTINGS_MAGIC);
		file.addToFile(SETTINGS_VERSION);
		if(genome != null)
		{
			GenomeBuildUID gbuid = genome.getUIDEnum();
			if(gbuid != null) file.addToFile(gbuid.getUID());
			else file.addToFile(-1);
		}
		else file.addToFile(-1);
		file.addToFile(mergeFactor);
		file.addToFile((short)varTableType);
		file.addVariableLengthString(dbName, BinFieldSize.WORD, 2);
		if(this.omim_table_path != null) file.addVariableLengthString(this.omim_table_path, BinFieldSize.WORD, 2);
		else file.addVariableLengthString("null", BinFieldSize.WORD, 2);
		
		if(varTableType == VARTBL_TYPE_SQL)
		{
			file.addVariableLengthString(sqlURL, BinFieldSize.WORD, 2);
			file.addVariableLengthString(sqlUser, BinFieldSize.WORD, 2);
			file.addVariableLengthString(sqlPassword, BinFieldSize.WORD, 2);
		}
		
		file.writeFile(settings_path);
		
		//Save samples
		sampleTable.saveTable();
		
		//Save variants
		variantTable.save();
		
	}
	
	/* ----- Sample Management ----- */
	
	public void addFamily(Family fam)
	{
		sampleTable.addOrReplaceFamily(fam);
	}

	public boolean removeFamily(int famUID)
	{
		//Don't forget to remove from variant table...
		Family f = sampleTable.removeFamily(famUID);
		if(f == null) return false;
		
		return variantTable.removeFamily(f);
	}
	
	public boolean removeFamily(String famName)
	{
		Family f = sampleTable.removeFamily(famName);
		if(f == null) return false;
		
		return variantTable.removeFamily(f);
	}
	
	public boolean removeSample(int sampleUID)
	{
		FamilyMember s = sampleTable.removeSample(sampleUID);
		if(s == null) return false;
		return variantTable.removeSample(s);
	}
	
	public Family getFamily(int famUID)
	{
		return sampleTable.getFamily(famUID);
	}
	
	public Family getFamily(String famName)
	{
		return sampleTable.getFamily(famName);
	}
	
	public FamilyMember getSample(int sampleUID)
	{
		return sampleTable.getSample(sampleUID);
	}
	
	public FamilyMember getSample(String sampleName)
	{
		return sampleTable.getSample(sampleName);
	}
	
	/* ----- Variant Management ----- */
	
	private List<String> getVCFSampleNames(String vcfpath) throws IOException
	{
		List<String> slist = new LinkedList<String>();
		BufferedReader br = new BufferedReader(new FileReader(vcfpath));
		String line = null;
		while((line = br.readLine()) != null)
		{
			if(line.isEmpty()) continue;
			if(line.startsWith("#CHROM"))
			{
				String[] fields = line.split("\t");
				if(fields.length > 9)
				{
					for(int i = 9; i < fields.length; i++)
					{
						slist.add(fields[i]);
					}
				}
				break;
			}
		}
		
		br.close();
		
		return slist;
	}
	
	public boolean addVCF(String vcfpath, boolean verbose) throws IOException
	{
		//See if VCF even exists
		if(!FileBuffer.fileExists(vcfpath)) return false;
		
		//Get sample names and pull sample info
		List<String> snames = getVCFSampleNames(vcfpath);
		Map<String, FamilyMember> smap = new HashMap<String, FamilyMember>();
		for(String s : snames)
		{
			FamilyMember fm = sampleTable.getSample(s);
			if (fm == null) 
			{
				if(verbose) System.err.println("SVDatabase.addVCF || WARNING: Sample \"" + s + "\" was not found in database! Skipping...");
				continue;
			}
			smap.put(s, fm);
			sampleTable.markSample(fm.getUID());
		}
		
		//Add to variant table...
		return variantTable.addVCF(vcfpath, smap, mergeFactor);
	}
	
	/* ----- Family Data Dumping ----- */
	
	private Map<String, GeneHitCounter> geneHitMap;
	
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
		double maxFreq = calculateUpperLimit(var.getIndividualCount(), sampleTable.countSamples());
		if (maxFreq > threshold) return false;
		
		//By population
		for (Population p : Population.values())
		{
			int hits = var.getIndividualCount(p);
			int total = sampleTable.countSamplesInPopulation(p);
			maxFreq = calculateUpperLimit(hits, total);
			if (maxFreq > threshold) return false;
		}
		
		return true;
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
	
	public void loadOMIMAnno()
	{
		if(omim_table_path != null && !omim_table_path.isEmpty());
		OMIMGeneMapImporter o = new OMIMGeneMapImporter(omim_table_path);
		o.importTable(genes);
	}
	
	public void dumpFamily(Family f, String dumpDir, boolean verbose) throws IOException
	{
		/*
		 * Table Fields ----- (TSV)
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
		 * 		[Genotypes] VCF Style
		 * 		[Genotypes] #Occurrences of allele
		 * 		[Genotypes]	Start Position
		 * 		[Genotypes]	End Position
		 * 			...
		 */
		
		//This might be VERY SLOW
		if(verbose) System.err.println("Checking family variants...");
		List<Long> varIDList = variantTable.getVariantIDsForFamily(f);
		if(verbose) System.err.println("Generating gene hit map...");
		if(geneHitMap == null) geneHitMap = variantTable.generateGeneHitMap();
		if(verbose) System.err.println("Loading OMIM annotations...");
		loadOMIMAnno();
		if(verbose) System.err.println("Counting samples...");
		int indivCount = sampleTable.countSamples();
		Map<Population, Integer> popIndivCount = new HashMap<Population, Integer>();
		for (Population p : Population.values()) popIndivCount.put(p, sampleTable.countSamplesInPopulation(p));
		
		//Load all vars (eeeekkk) and do seg analysis
		if(verbose) System.err.println("Loading variants...");
		List<FamilyMember> members = f.getAllFamilyMembers();
		VariantPool pool = new VariantPool(f.countMembers());
		Map<String, Long> vidMap = new HashMap<String, Long>();
		for(Long vid : varIDList)
		{
			DBVariant v = variantTable.getVariant(vid);
			VariantGenotype vg = variantTable.getGenotype(vid);
			vidMap.put(v.getName(), vid);
			
			StructuralVariant sv = v.toStructuralVariant();
			for(FamilyMember m : members)
			{
				SVDBGenotype dbg = vg.getGenotype(m.getUID());
				Genotype g = new Genotype();
				if(dbg == null)
				{
					//Ref/Ref or ./.
					if(sampleTable.sampleHasData(m.getUID())) g.setAlleles("0/0");
					else g.setAlleles("./.");	
					g.setCopyNumber(2);
				}
				else
				{
					Collection<SVDBAllele> alist = dbg.getAlleles();
					List<Integer> alleles = new LinkedList<Integer>();
					int j = 0;
					for(SVDBAllele a : alist)
					{
						j++;
						for(int i = 0; i < a.getAlleleCount(); i++) alleles.add(j);
					}
					j = 0;
					int acount = alleles.size();
					int[] aarr = new int[acount];
					if (acount < 2) {
						aarr[j] = 0;
						acount++; //Put a ref in there
						j++;
					}
					for(Integer i : alleles)
					{
						aarr[j] = i;
						j++;
					}
					g.setCopyNumber(acount);
					g.setAlleles(aarr);
				}
				sv.addGenotype(m.getName(), g);
			}
			pool.addVariant(sv);
		}
		
		//Seg analysis
		if(verbose) System.err.println("Running segregation analysis...");
		List<Candidate> candidates = Inheritor.getCandidates(pool, f, genes);
		
		//Prepare for output
		if(verbose) System.err.println("Preparing table outputs...");
		String full_table_path = dumpDir + File.separator + f.getFamilyName() + "_fullTable.tsv";
		String pri_table_path = dumpDir + File.separator + f.getFamilyName() + "_priorizedTable.tsv";
		String header = "Var Name\tVar ID (Hex)\tChrom\tEnd Chrom (TRA)\tPos\tEnd\t";
		header += "INS Seq\tSV Type\tSize\tPos Effect\tGene Name\tTranscriptID\tOMIM\t";
		header += "DECIPHER\tCohort Indiv Count\tCohort Hom Count\tCohort Pop Freq\tCohort Variant Hit Count\t";
		header += "Cohort Indiv Hit Count\tCohort Variant Exon Hit Count\tCohort Indiv Exon Hit Count\t";
		Population[] plist = Population.values();
		for(Population p : plist)
		{
			header += p + "_ALLELE_COUNT\t";
			header += p + "_HOM_COUNT\t";
			header += p + "_POP_FREQ\t";
		}
		List<Individual> aff = f.getAllAffected();
		for(Individual a : aff) header += a.getName() + "_SEG\t";
		for(Individual a : aff) header += a.getName() + "_HET_PARTNERS\t";
		for(FamilyMember m : members)
		{
			header += m.getName() + " Genotype (VCF)\t";
		}
		for(FamilyMember m : members)
		{
			header += m.getName() + " Genotype (DB)\t";
		}
		
		//Open output streams and write tables
		if(verbose) System.err.println("Writing tables...");
		BufferedWriter fullTable = new BufferedWriter(new FileWriter(full_table_path));
		BufferedWriter priTable = new BufferedWriter(new FileWriter(pri_table_path));
		fullTable.write(header + "\n");
		priTable.write(header + "\n");
		int cCount = 0;
		for(Candidate c : candidates)
		{
			if(c.getAllele() == 0) continue; //Skip ref allele candidates
			
			//Get variant & geno data
			Long vid = vidMap.get(c.getVariant().getVarID());
			DBVariant dbv = variantTable.getVariant(vid);
			VariantGenotype vg = variantTable.getGenotype(vid);
			
			//Generate record
			StringBuilder sb = new StringBuilder();
			sb.append(dbv.getName() + "\t");
			sb.append(Long.toHexString(vid) + "\t");
			sb.append(dbv.getChrom().getUDPName() + "\t");
			if(dbv.getType() == SVType.TRA) sb.append(dbv.getEndChrom().getUDPName() + "\t");
			else sb.append("[N/A]\t");
			if(dbv.getStartPosition().getStart() == dbv.getStartPosition().getStart()) sb.append(dbv.getStartPosition().getStart() + "\t");
			else sb.append(dbv.getStartPosition().getStart() + "-" + dbv.getStartPosition().getEnd() + "\t");
			if(dbv.getEndPosition().getStart() == dbv.getEndPosition().getStart()) sb.append(dbv.getEndPosition().getStart() + "\t");
			else sb.append(dbv.getEndPosition().getStart() + "-" + dbv.getEndPosition().getEnd() + "\t");
			if(dbv.getType() == SVType.INS || dbv.getType() == SVType.INSME) sb.append(dbv.getAltAlleleString() + "\t");
			else sb.append("[N/A]\t");
			sb.append(dbv.getType().getString() + "\t");
			if(dbv.getType() != SVType.TRA) sb.append(dbv.getEndPosition().getEnd() - dbv.getStartPosition().getStart() + "bp\t");
			else sb.append("[N/A]\t");
			sb.append(c.getPositionEffect()+ "\t");
			Gene g = c.getGene();
			if(g != null)
			{
				sb.append(g.getName() + "\t");
				sb.append(g.getID() + "\t");
				String omim = g.getAnnotation(OMIMGeneMapImporter.ANNO_KEY);
				if(omim != null) sb.append(omim + "\t");
				else sb.append("[NONE]\t");
				sb.append("TBI\t");
			}
			else sb.append("[N/A]\t[N/A]\t[N/A]\t[N/A]\t");
			sb.append(dbv.getIndividualCount() + "\t");
			sb.append(dbv.getHomozygoteCount() + "\t");
			double pfreq = (double)dbv.getIndividualCount() / (double)indivCount;
			sb.append(String.format("%.3f", pfreq) + "\t");
			
			if(g != null)
			{
				GeneHitCounter ghc = geneHitMap.get(g.getID());
				if(ghc != null)
				{
					sb.append(ghc.total_hits_var + "\t");
					sb.append(ghc.total_hits_indiv.size() + "\t");
					sb.append(ghc.exon_hits_var + "\t");
					sb.append(ghc.exon_hits_indiv.size() + "\t");
				}
				else sb.append("0\t0\t0\t0\t");
			}
			else sb.append("[N/A]\t[N/A]\t[N/A]\t[N/A]\t");
			
			for(Population p : plist)
			{
				sb.append(dbv.getIndividualCount(p) + "\t");
				sb.append(dbv.getHomozygoteCount(p) + "\t");
				pfreq = (double)dbv.getIndividualCount(p) / (double)popIndivCount.get(p);
				sb.append(String.format("%.3f", pfreq) + "\t");
			}
			 
			for(Individual a : aff) sb.append(c.getInheritancePattern(a)+ "\t");
			for(Individual a : aff) 
			{
				List<Candidate> partners = c.getAllPartners(a);
				if(partners == null || partners.isEmpty()) sb.append(c.getInheritancePattern(a)+ "N/A\t");
				else
				{
					boolean first = true;
					for(Candidate p : partners)
					{
						if(!first)sb.append(";");
						String vname = p.getVariant().getVarID();
						long pid = vidMap.get(vname);
						sb.append(Long.toHexString(pid));
						first = false;
					}
					sb.append("\t");
				}
			}
			
			for(FamilyMember m : members)
			{
				//VCF style genotypes
				Genotype geno = c.getVariant().getSampleGenotype(m.getName());
				if(geno == null) sb.append("0/0\t");
				else
				{
					int[] alleles = geno.getAlleles();
					boolean first = true;
					for(int a : alleles)
					{
						if(!first) sb.append("/");
						first = false;
						sb.append(a);
					}
					sb.append("\t");
				}
			}
			for(FamilyMember m : members)
			{
				//DB style genotypes
				SVDBGenotype geno = vg.getGenotype(m.getUID());
				if(geno == null) sb.append("REF\t");
				else
				{
					Collection<SVDBAllele> alist = geno.getAlleles();
					boolean first = true;
					for(SVDBAllele a : alist)
					{
						if(!first) sb.append(";");
						sb.append(a.getAllele().getStart() + "-" + a.getAllele().getEnd() + "(" + a.getAlleleCount() + ")");
						first = false;
					}
					sb.append("\t");
				}
			}
			
			//Write to full table
			String record = sb.toString();
			fullTable.write(record + "\n");
			
			//Determine whether to prioritize and write to p table if so
			if(prioritizeCandidate(c, dbv, aff))
			{
				priTable.write(record + "\n");
			}
			cCount++;
			if(verbose && cCount % 200 == 0)  System.err.println("Records processed: " + cCount);
		}
		fullTable.close();
		priTable.close();
		if(verbose) System.err.println("Family table dump complete!");
	}
	
	public void setOMIMTablePath(String path)
	{
		omim_table_path = path;
	}
	
	public void dumpToBED(int sampleUID, SVType type, String outpath) throws IOException
	{
		List<Long> idList = variantTable.getVariantIDsForSampleOfType(sampleUID, type);
		if(idList != null && !idList.isEmpty())
		{
			BufferedWriter bed = new BufferedWriter(new FileWriter(outpath));
			bed.write("#CHROM\tPOS\tEND\n");
			for(Long id : idList)
			{
				DBVariant dbv = variantTable.getVariant(id);
				if(dbv != null)
				{
					bed.write(dbv.getChrom().getUDPName() + "\t");
					bed.write(dbv.getStartPosition().getStart() + "\t");
					bed.write(dbv.getEndPosition().getEnd() + "\n");
				}
			}
			bed.close();
		}
	}
	
	/* ----- Database Loaders ----- */
	
	public static SVDatabase loadDatabase(String dir) throws IOException, SQLException
	{
		SVDatabase db = new SVDatabase(dir, true);
		return db;
	}
	
	public static SVDatabase newDatabase(String dir, String name, int mergeFactor, GenomeBuildUID gbid) throws IOException
	{
		//Load genome...
		GenomeBuild gb = GenomeBuild.loadStandardBuild(gbid.getName());
		
		//Generate db
		SVDatabase db = new SVDatabase(dir, name, mergeFactor, gb, true);
		return db;
	}

	public static SVDatabase newDatabase(String dir, String name, int mergeFactor, GenomeBuild gb) throws IOException
	{
		SVDatabase db = new SVDatabase(dir, name, mergeFactor, gb, true);
		return db;
	}
	
	public static SVDatabase newDatabase(String dir, String name, int mergeFactor, GenomeBuild gb, String sqlPath) throws IOException, SQLException
	{
		SVDatabase db = new SVDatabase(dir, name, mergeFactor, gb, sqlPath);
		return db;
	}

	
}
