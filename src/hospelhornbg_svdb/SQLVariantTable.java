package hospelhornbg_svdb;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.sql.Blob;
import java.sql.Connection;
import java.sql.DatabaseMetaData;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ConcurrentMap;

import hospelhornbg_bioinformatics.Genotype;
import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.Gene;
import hospelhornbg_genomeBuild.Gene.Exon;
import hospelhornbg_genomeBuild.GeneFunc;
import hospelhornbg_genomeBuild.GeneSet;
import hospelhornbg_genomeBuild.GenomeBuild;
import hospelhornbg_segregation.Family;
import hospelhornbg_segregation.FamilyMember;
import hospelhornbg_segregation.Population;
import hospelhornbg_svdb.DBVariantTable.GeneHitCounter;
import waffleoRai_Utils.FileBuffer;

public class SQLVariantTable implements VariantTable{
	
	/* ----- Constants ----- */
	
	public static final int MAX_VARS_INQUERY = 256;
	
	public static final String TABLENAME_VARIANTS = "VARIANTS";
	public static final String TABLENAME_SAMPLEGENO = "SAMPLEGENO";
	public static final String TABLENAME_GENEHITS = "GENES";
	
	public static final String FIELDNAME_VARUID = "VARUID";
	public static final String FIELDNAME_CTG1 = "CONTIG1";
	public static final String FIELDNAME_START1 = "START1";
	public static final String FIELDNAME_START2 = "START2";
	public static final String FIELDNAME_END1 = "END1";
	public static final String FIELDNAME_END2 = "END2";
	public static final String FIELDNAME_SVTYPE = "SVTYPE";
	public static final String FIELDNAME_POSEFF = "POSEFF";
	public static final String FIELDNAME_VARNAME = "VARNAME";
	
	public static final String FIELDNAME_ACOUNT_TOT = "ALLELECOUNT_TOTAL";
	public static final String FIELDNAME_HCOUNT_TOT = "HOMOZYGOTECOUNT_TOTAL";
	public static final String FIELDNAME_ACOUNT_NFE = "ALLELECOUNT_NFE";
	public static final String FIELDNAME_HCOUNT_NFE = "HOMOZYGOTECOUNT_NFE";
	public static final String FIELDNAME_ACOUNT_AFR = "ALLELECOUNT_AFR";
	public static final String FIELDNAME_HCOUNT_AFR = "HOMOZYGOTECOUNT_AFR";
	public static final String FIELDNAME_ACOUNT_AMR = "ALLELECOUNT_AMR";
	public static final String FIELDNAME_HCOUNT_AMR = "HOMOZYGOTECOUNT_AMR";
	public static final String FIELDNAME_ACOUNT_FIN = "ALLELECOUNT_FIN";
	public static final String FIELDNAME_HCOUNT_FIN = "HOMOZYGOTECOUNT_FIN";
	public static final String FIELDNAME_ACOUNT_EAS = "ALLELECOUNT_EAS";
	public static final String FIELDNAME_HCOUNT_EAS = "HOMOZYGOTECOUNT_EAS";
	public static final String FIELDNAME_ACOUNT_SAS = "ALLELECOUNT_SAS";
	public static final String FIELDNAME_HCOUNT_SAS = "HOMOZYGOTECOUNT_SAS";
	public static final String FIELDNAME_ACOUNT_ASJ = "ALLELECOUNT_ASJ";
	public static final String FIELDNAME_HCOUNT_ASJ = "HOMOZYGOTECOUNT_ASJ";
	public static final String FIELDNAME_ACOUNT_OTH = "ALLELECOUNT_OTH";
	public static final String FIELDNAME_HCOUNT_OTH = "HOMOZYGOTECOUNT_OTH";
	
	public static final String FIELDNAME_GENELIST = "GENES";
	public static final String FIELDNAME_VALNOTES = "VALIDATIONNOTES";
	public static final String FIELDNAME_CTG2 = "CONTIG2";
	public static final String FIELDNAME_INSSEQ = "INSERTION_SEQ";
	public static final String FIELDNAME_GENOTYPES = "GENOTYPES";
	
	public static final String FIELDNAME_SAMPLEUID = "SAMPLEUID";
	public static final String FIELDNAME_SVARLIST_HOM = "HOMVAR";
	public static final String FIELDNAME_SVARLIST_HET = "HETVAR";
	public static final String FIELDNAME_SVARLIST_OTH = "OTHVAR";
	
	public static final String FIELDNAME_GH_GENEUID = "TRANSCRIPT_UID"; //int
	public static final String FIELDNAME_GH_HITS_T = "TOTAL_VARS"; //int
	public static final String FIELDNAME_GH_HITS_E = "EXON_VARS"; //int
	public static final String FIELDNAME_GH_HITS_TI = "TOTAL_INDIVS"; //blob
	public static final String FIELDNAME_GH_HITS_EI = "EXON_INDIVS"; //blob
	
	public static final String REGIDX_FILESTEM = "regidx_";
	
	public static final String[][] VAR_COLUMNS = {
			{FIELDNAME_VARUID, "BIGINT"},{FIELDNAME_CTG1, "INTEGER"},{FIELDNAME_START1, "INTEGER"},
			{FIELDNAME_START2, "INTEGER"},{FIELDNAME_END1, "INTEGER"},{FIELDNAME_END2, "INTEGER"},
			{FIELDNAME_SVTYPE, "SMALLINT"},{FIELDNAME_POSEFF, "SMALLINT"},{FIELDNAME_VARNAME, "VARCHAR(255)"},
			{FIELDNAME_ACOUNT_TOT, "INTEGER"},{FIELDNAME_HCOUNT_TOT, "INTEGER"},
			{FIELDNAME_ACOUNT_NFE, "INTEGER"},{FIELDNAME_HCOUNT_NFE, "INTEGER"},
			{FIELDNAME_ACOUNT_AFR, "INTEGER"},{FIELDNAME_HCOUNT_AFR, "INTEGER"},
			{FIELDNAME_ACOUNT_AMR, "INTEGER"},{FIELDNAME_HCOUNT_AMR, "INTEGER"},
			{FIELDNAME_ACOUNT_FIN, "INTEGER"},{FIELDNAME_HCOUNT_FIN, "INTEGER"},
			{FIELDNAME_ACOUNT_EAS, "INTEGER"},{FIELDNAME_HCOUNT_EAS, "INTEGER"},
			{FIELDNAME_ACOUNT_SAS, "INTEGER"},{FIELDNAME_HCOUNT_SAS, "INTEGER"},
			{FIELDNAME_ACOUNT_ASJ, "INTEGER"},{FIELDNAME_HCOUNT_ASJ, "INTEGER"},
			{FIELDNAME_ACOUNT_OTH, "INTEGER"},{FIELDNAME_HCOUNT_OTH, "INTEGER"},
			{FIELDNAME_GENELIST, "BLOB"},{FIELDNAME_VALNOTES, "VARCHAR(20000)"},{FIELDNAME_CTG2, "INTEGER"},
			{FIELDNAME_INSSEQ, "BLOB"},{FIELDNAME_GENOTYPES, "BLOB"}};
	
	public static final String[][] SAMPLEGENO_COLUMNS = {{FIELDNAME_SAMPLEUID, "INTEGER"},
			{FIELDNAME_SVARLIST_HOM, "BLOB"},
			{FIELDNAME_SVARLIST_HET, "BLOB"},
			{FIELDNAME_SVARLIST_OTH, "BLOB"},};
	
	public static final String[][] GENEHITS_COLUMNS = {{FIELDNAME_GH_GENEUID, "INTEGER"},
													   {FIELDNAME_GH_HITS_T, "INTEGER"},
													   {FIELDNAME_GH_HITS_E, "INTEGER"},
													   {FIELDNAME_GH_HITS_TI, "BLOB"},
													   {FIELDNAME_GH_HITS_EI, "BLOB"},};
	
	/* ----- Instance Variables ----- */
	
	private String dbDir;
	private String dbURL;
	private String username;
	private String password;
	
	private Connection connection;
	private SQLManager sqlManager;
	//private StatementPrepper sprepper;
	
	private int mergeFactor;
	private double percLeeway;

	private GenomeBuild genome;
	private GeneSet genes;
	private GenomeIndex uidIndex;
	
	private boolean threadlock;
	private ReadCache read_cache;
	private RegionSearchCache rs_cache;
	private RegionIndex regidx;
	
	//private List<String> tempBlobFiles;
	
	/* ----- Construction ----- */
	
	public SQLVariantTable(String dbdir, String url, String user, String pw, GenomeBuild gb, GeneSet gs, int mf) throws SQLException
	{
		dbDir = dbdir;
		genome = gb;
		genes = gs;
		uidIndex = new GenomeIndex(genome);
		//Attempt to connect
		dbURL = url;
		username = user.toUpperCase();
		password = pw;
		read_cache = new ReadCache(false);
		rs_cache = new RegionSearchCache();
		System.err.println("Now connecting...");
		connect();
		System.err.println("Connection successful!");
		//Check for tables, create if not there
		if(!varTableExists()) createVarTable();
		if(!sampleGenoTableExists()) createSampleGenoTable();
		if(!geneHitTableExists()) createGeneHitTable();
		
		sqlManager = new SQLManager(connection);
		mergeFactor = mf;
		
		percLeeway = (double)mergeFactor / 1000.0;
		loadGeneHitTable();
	}
	
	private void connect() throws SQLException
	{
		//System.err.println("URL: " + dbURL + " | Username: " + username + " | PW: " + password);
		connection = DriverManager.getConnection(dbURL, username, password);
		//cstat = connection.createStatement();
		connection.setAutoCommit(false);
	}
	
	/* ----- SQL Table Management ----- */
	
	private boolean varTableExists() throws SQLException
	{
		DatabaseMetaData meta = connection.getMetaData();
		String[] ttypes = {"TABLE"};
		ResultSet rs = meta.getTables(null, null, TABLENAME_VARIANTS.toUpperCase(), ttypes);
		boolean b = rs.next();
		rs.close();
		return b;
	}
	
	private void createVarTable() throws SQLException
	{
		System.err.println("Variant table not found. Creating variant table...");
		String sqlcmd = "CREATE TABLE " + SQLVariantTable.TABLENAME_VARIANTS + "(";
		boolean first = true;
		for(String[] field : VAR_COLUMNS)
		{
			if(!first) sqlcmd += ",";
			sqlcmd += field[0] + " " + field[1];
			first = false;
		}
		sqlcmd += ")";
		//System.err.println("-DEBUG- SQLCOMMAND: " + sqlcmd);
		//System.err.println("-DEBUG- SQLCOMMAND (Short): " + sqlcmd.substring(0, 157));
		//System.exit(1);
		Statement cstat = connection.createStatement();
		cstat.executeUpdate(sqlcmd);
		cstat.closeOnCompletion();
		
		/*Index*/
		sqlcmd = "CREATE UNIQUE INDEX VARUID_INDEX ON " + SQLVariantTable.TABLENAME_VARIANTS + "(" + SQLVariantTable.FIELDNAME_VARUID + ")";
		cstat = connection.createStatement();
		cstat.executeUpdate(sqlcmd);
		cstat.closeOnCompletion();
	}
	
	private boolean sampleGenoTableExists() throws SQLException
	{
		DatabaseMetaData meta = connection.getMetaData();
		String[] ttypes = {"TABLE"};
		ResultSet rs = meta.getTables(null, null, TABLENAME_SAMPLEGENO.toUpperCase(), ttypes);
		boolean b = rs.next();
		rs.close();
		return b;
	}
	
	private void createSampleGenoTable() throws SQLException
	{
		System.err.println("Sample geno table not found. Creating table...");
		String sqlcmd = "CREATE TABLE " + SQLVariantTable.TABLENAME_SAMPLEGENO + "(";
		boolean first = true;
		for(String[] field : SAMPLEGENO_COLUMNS)
		{
			if(!first) sqlcmd += ",";
			sqlcmd += field[0] + " " + field[1];
			first = false;
		}
		sqlcmd += ")";
		Statement cstat = connection.createStatement();
		cstat.executeUpdate(sqlcmd);
		cstat.closeOnCompletion();
		
		/*Index*/
		sqlcmd = "CREATE UNIQUE INDEX SGENOUID_INDEX ON " + SQLVariantTable.TABLENAME_SAMPLEGENO + "(" + SQLVariantTable.FIELDNAME_SAMPLEUID + ")";
		cstat = connection.createStatement();
		cstat.executeUpdate(sqlcmd);
		cstat.closeOnCompletion();
	}
	
	private boolean geneHitTableExists() throws SQLException
	{
		DatabaseMetaData meta = connection.getMetaData();
		String[] ttypes = {"TABLE"};
		ResultSet rs = meta.getTables(null, null, TABLENAME_GENEHITS.toUpperCase(), ttypes);
		boolean b = rs.next();
		rs.close();
		return b;
	}
	
	private void createGeneHitTable() throws SQLException
	{
		System.err.println("Gene hit table not found. Creating table...");
		String sqlcmd = "CREATE TABLE " + SQLVariantTable.TABLENAME_GENEHITS + "(";
		boolean first = true;
		for(String[] field : GENEHITS_COLUMNS)
		{
			if(!first) sqlcmd += ",";
			sqlcmd += field[0] + " " + field[1];
			first = false;
		}
		sqlcmd += ")";
		Statement cstat = connection.createStatement();
		cstat.executeUpdate(sqlcmd);
		cstat.closeOnCompletion();
		
		//Add all genes with 0's
		zeroGeneHitTable();
		
		/*Index*/
		sqlcmd = "CREATE UNIQUE INDEX GENEUID_INDEX ON " + SQLVariantTable.TABLENAME_GENEHITS + "(" + SQLVariantTable.FIELDNAME_GH_GENEUID + ")";
		cstat = connection.createStatement();
		cstat.executeUpdate(sqlcmd);
		cstat.closeOnCompletion();
	}
	
	private void zeroGeneHitTable() throws SQLException
	{
		//
		List<Gene> glist = genes.getAllGenes();
		//StatementPrepper sprepper = sqlManager.getStatementGenerator();
		StatementPrepper sprepper = new StatementPrepper(connection);
		PreparedStatement gh_insert = sprepper.getGeneHitInsertStatement();
		int dbc = 0;
		Set<Integer> uidset = new TreeSet<Integer>();
		for(Gene g : glist)
		{
			gh_insert.setInt(StatementPrepper.GENEHIT_UID, g.getGUID());
			gh_insert.setInt(StatementPrepper.GENEHIT_TOT, 0);
			gh_insert.setInt(StatementPrepper.GENEHIT_EXON, 0);
			
			byte[] neg1 = {-1, -1, -1, -1};
			Blob b1 = sprepper.wrapInBlob(neg1);
			gh_insert.setBlob(StatementPrepper.GENEHIT_TOT_INDIV, b1);
			Blob b2 = sprepper.wrapInBlob(neg1);
			gh_insert.setBlob(StatementPrepper.GENEHIT_EXON_INDIV, b2);
			
			gh_insert.executeUpdate();
			dbc++;
			uidset.add(g.getGUID());
		}
		System.err.println("Gene hit table zeroing complete! " + dbc + " records written!");
		System.err.println("Unique IDs: " + uidset.size());
	}
	
	/* ----- Parsing Retrieved ----- */
	
	private DBVariant readFromResultSet(ResultSet rs) throws SQLException, IOException
	{
		if(rs == null) return null;
		DBVariant dbv = DBVariant.getEmptyVariant();
		dbv.setLongUID(rs.getLong(FIELDNAME_VARUID));
		int c1id = rs.getInt(FIELDNAME_CTG1);
		Contig c1 = genome.getContigByUID(c1id);
		dbv.setContig1(c1);
		
		int st1 = rs.getInt(FIELDNAME_START1);
		int st2 = rs.getInt(FIELDNAME_START2);
		int ed1 = rs.getInt(FIELDNAME_END1);
		int ed2 = rs.getInt(FIELDNAME_END2);
		dbv.setStartInterval(st1, st2);
		dbv.setEndInterval(ed1, ed2);
		dbv.setSVType(rs.getInt(FIELDNAME_SVTYPE));
		dbv.setPosEff(rs.getInt(FIELDNAME_POSEFF));
		dbv.setVariantName(rs.getString(FIELDNAME_VARNAME));
		
		//Counts
		dbv.setTotalCount(rs.getInt(FIELDNAME_ACOUNT_TOT));
		dbv.setTotalCount(rs.getInt(FIELDNAME_ACOUNT_NFE), Population.NFE);
		dbv.setTotalCount(rs.getInt(FIELDNAME_ACOUNT_AFR), Population.AFR);
		dbv.setTotalCount(rs.getInt(FIELDNAME_ACOUNT_AMR), Population.AMR);
		dbv.setTotalCount(rs.getInt(FIELDNAME_ACOUNT_FIN), Population.FIN);
		dbv.setTotalCount(rs.getInt(FIELDNAME_ACOUNT_SAS), Population.SAS);
		dbv.setTotalCount(rs.getInt(FIELDNAME_ACOUNT_EAS), Population.EAS);
		dbv.setTotalCount(rs.getInt(FIELDNAME_ACOUNT_ASJ), Population.ASJ);
		dbv.setTotalCount(rs.getInt(FIELDNAME_ACOUNT_OTH), Population.OTH);
		
		dbv.setHomozygoteCount(rs.getInt(FIELDNAME_HCOUNT_TOT));
		dbv.setHomozygoteCount(rs.getInt(FIELDNAME_HCOUNT_NFE), Population.NFE);
		dbv.setHomozygoteCount(rs.getInt(FIELDNAME_HCOUNT_AFR), Population.AFR);
		dbv.setHomozygoteCount(rs.getInt(FIELDNAME_HCOUNT_AMR), Population.AMR);
		dbv.setHomozygoteCount(rs.getInt(FIELDNAME_HCOUNT_FIN), Population.FIN);
		dbv.setHomozygoteCount(rs.getInt(FIELDNAME_HCOUNT_SAS), Population.SAS);
		dbv.setHomozygoteCount(rs.getInt(FIELDNAME_HCOUNT_EAS), Population.EAS);
		dbv.setHomozygoteCount(rs.getInt(FIELDNAME_HCOUNT_ASJ), Population.ASJ);
		dbv.setHomozygoteCount(rs.getInt(FIELDNAME_HCOUNT_OTH), Population.OTH);
		
		Blob geneblob = rs.getBlob(FIELDNAME_GENELIST);
		dbv.loadGeneListFromBLOB(geneblob, genes);
		dbv.setValidationNotes(rs.getString(FIELDNAME_VALNOTES));
		
		if(dbv.getType() == SVType.TRA || dbv.getType() == SVType.BND)
		{
			int c2id = rs.getInt(FIELDNAME_CTG2);
			Contig c2 = genome.getContigByUID(c2id);
			dbv.setContig2(c2);
		}
		
		if(dbv.getType() == SVType.INS || dbv.getType() == SVType.INSME)
		{
			//Blob, apparently.
			Blob insblob = rs.getBlob(FIELDNAME_INSSEQ);
			InputStream is = insblob.getBinaryStream();
			StringBuilder sb = new StringBuilder(4096);
			int b = -1;
			try 
			{
				while((b = is.read()) != -1) {
					if(b == 0) break;
					sb.append((char)b);
				}
				dbv.setInsSeq(sb.toString());
				is.close();
			} 
			catch (IOException e) 
			{
				e.printStackTrace();
				return null;
			}
		}
		
		return dbv;
	}
	
	private List<Long> readVarUIDListBlob(Blob blob) throws SQLException, IOException
	{
		List<Long> list = new LinkedList<Long>();
		
		if(blob == null) return list;
		FileBuffer fblob = new FileBuffer((int)blob.length(), true);
		InputStream is = blob.getBinaryStream();
		int i = -1;
		while((i = is.read()) != -1)
		{
			fblob.addToFile((byte)i);
		}
		is.close();
		
		long fsz = fblob.getFileSize();
		long cpos = 0;
		while(cpos < fsz)
		{
			long vid = fblob.longFromFile(cpos); cpos += 8;
			if(vid == -1L) break;
			list.add(vid);
		}
		
		return list;
	}
	
	/* ----- Storage Prep ----- */
	
	private PreparedStatement generateFullVarInsertStatement(DBVariant var, VariantGenotype vgeno, boolean threadlock) throws IOException, SQLException
	{
		PreparedStatement pstat = null;
		StatementPrepper sprepper = sqlManager.getStatementGenerator();
		if(threadlock) pstat = sqlManager.requestStatement(SQLManager.SKEY_VAR_INSERT, true);
		else pstat = sqlManager.getStatementGenerator().getFullInsertStatement();
		
		Contig ctg = var.getChrom();
		
		pstat.setLong(StatementPrepper.FULLINS_VARUID, var.getLongID());
		pstat.setInt(StatementPrepper.FULLINS_CHR1, ctg.getUID());
		pstat.setInt(StatementPrepper.FULLINS_ST1, var.getStartPosition().getStart());
		pstat.setInt(StatementPrepper.FULLINS_ST2, var.getStartPosition().getEnd());
		pstat.setInt(StatementPrepper.FULLINS_ED1, var.getEndPosition().getStart());
		pstat.setInt(StatementPrepper.FULLINS_ED2, var.getEndPosition().getEnd());
		pstat.setInt(StatementPrepper.FULLINS_SVTYPE, var.getType().getID());
		pstat.setInt(StatementPrepper.FULLINS_POSEFF, var.getPositionEffect().getPriority());
		pstat.setString(StatementPrepper.FULLINS_VARNAME, var.getName());
		
		pstat.setInt(StatementPrepper.FULLINS_ACOUNT_TOT, var.getIndividualCount()); //System.err.println("Total Indiv Count: " + var.getIndividualCount());
		pstat.setInt(StatementPrepper.FULLINS_HCOUNT_TOT, var.getHomozygoteCount()); //System.err.println("Homozygote Count: " + var.getHomozygoteCount());
		pstat.setInt(StatementPrepper.FULLINS_ACOUNT_NFE, var.getIndividualCount(Population.NFE));
		pstat.setInt(StatementPrepper.FULLINS_HCOUNT_NFE, var.getHomozygoteCount(Population.NFE));
		pstat.setInt(StatementPrepper.FULLINS_ACOUNT_AFR, var.getIndividualCount(Population.AFR));
		pstat.setInt(StatementPrepper.FULLINS_HCOUNT_AFR, var.getHomozygoteCount(Population.AFR));
		pstat.setInt(StatementPrepper.FULLINS_ACOUNT_AMR, var.getIndividualCount(Population.AMR));
		pstat.setInt(StatementPrepper.FULLINS_HCOUNT_AMR, var.getHomozygoteCount(Population.AMR));
		pstat.setInt(StatementPrepper.FULLINS_ACOUNT_FIN, var.getIndividualCount(Population.FIN));
		pstat.setInt(StatementPrepper.FULLINS_HCOUNT_FIN, var.getHomozygoteCount(Population.FIN));
		pstat.setInt(StatementPrepper.FULLINS_ACOUNT_EAS, var.getIndividualCount(Population.EAS));
		pstat.setInt(StatementPrepper.FULLINS_HCOUNT_EAS, var.getHomozygoteCount(Population.EAS));
		pstat.setInt(StatementPrepper.FULLINS_ACOUNT_SAS, var.getIndividualCount(Population.SAS));
		pstat.setInt(StatementPrepper.FULLINS_HCOUNT_SAS, var.getHomozygoteCount(Population.SAS));
		pstat.setInt(StatementPrepper.FULLINS_ACOUNT_ASJ, var.getIndividualCount(Population.ASJ));
		pstat.setInt(StatementPrepper.FULLINS_HCOUNT_ASJ, var.getHomozygoteCount(Population.ASJ));
		pstat.setInt(StatementPrepper.FULLINS_ACOUNT_OTH, var.getIndividualCount(Population.OTH));
		pstat.setInt(StatementPrepper.FULLINS_HCOUNT_OTH, var.getHomozygoteCount(Population.OTH));
		
		//Genelist blob
		if(threadlock) sqlManager.requestStatementExecution(true);
		Blob blob = sprepper.wrapInBlob(var.getGeneListAsBLOBBytes());
		pstat.setBlob(StatementPrepper.FULLINS_GENELIST, blob);
		//System.err.println("Gene list blob is null: " + (blob == null));
		
		String valnotes = var.getValidationNotes();
		if(valnotes == null) valnotes = "N/A";
		pstat.setString(StatementPrepper.FULLINS_VALNOTES, valnotes);
		
		Contig ctg2 = var.getEndChrom();
		if(ctg2 == null) pstat.setInt(StatementPrepper.FULLINS_CTG2, ctg.getUDPName().hashCode());
		else pstat.setInt(StatementPrepper.FULLINS_CTG2, ctg2.getUID());
		
		//String insseq = var.getAltAlleleString();
		//if(insseq == null) insseq = "N/A";
		//pstat.setString(StatementPrepper.FULLINS_INSSEQ, insseq);
		blob = sprepper.wrapInBlob(var.getInsseqAsBLOBBytes());
		pstat.setBlob(StatementPrepper.FULLINS_INSSEQ, blob);
		//System.err.println("Insseq blob is null: " + (blob == null));
		
		//Genotype blob
		blob = sprepper.wrapInBlob(vgeno.getGenotypesAsBLOBBytes());
		pstat.setBlob(StatementPrepper.FULLINS_GENOTYPES, blob);
		//System.err.println("Geno blob is null: " + (blob == null));
		
		return pstat;
	}
	
	private PreparedStatement generateAbridgedSetVarUpdateStatement(DBVariant var, VariantGenotype vgeno, boolean threadlock) throws IOException, SQLException
	{
		//Updates:
		//	Start, End
		//	Allele Counts
		//	Genes
		//  PosEff
		//	Validation notes
		//	Genotype
		
		//PreparedStatement pstat = sprepper.getShortUpdateStatement();
		
		PreparedStatement pstat = null;
		StatementPrepper sprepper = sqlManager.getStatementGenerator();
		if(threadlock) pstat = sqlManager.requestStatement(SQLManager.SKEY_VAR_UPDATE_SHORT, true);
		else pstat = sprepper.getShortUpdateStatement();
		
		pstat.setInt(StatementPrepper.SHORTUD_ST1, var.getStartPosition().getStart()); //System.err.println("Start1 = " + var.getStartPosition().getStart());
		pstat.setInt(StatementPrepper.SHORTUD_ST2, var.getStartPosition().getEnd()); //System.err.println("Start2 = " + var.getStartPosition().getEnd());
		pstat.setInt(StatementPrepper.SHORTUD_ED1, var.getEndPosition().getStart()); //System.err.println("End1 = " + var.getEndPosition().getStart());
		pstat.setInt(StatementPrepper.SHORTUD_ED2, var.getEndPosition().getEnd()); //System.err.println("End2 = " + var.getEndPosition().getEnd());
		pstat.setInt(StatementPrepper.SHORTUD_POSEFF, var.getPositionEffect().getPriority()); //System.err.println("PosEff = " + var.getPositionEffect().getPriority());
		
		pstat.setInt(StatementPrepper.SHORTUD_ACOUNT_TOT, var.getIndividualCount()); //System.err.println("Total Allele = " + var.getIndividualCount());
		pstat.setInt(StatementPrepper.SHORTUD_HCOUNT_TOT, var.getHomozygoteCount()); //System.err.println("Total Hom = " + var.getHomozygoteCount());
		pstat.setInt(StatementPrepper.SHORTUD_ACOUNT_NFE, var.getIndividualCount(Population.NFE)); //System.err.println("NFE Count = " + var.getIndividualCount(Population.NFE));
		pstat.setInt(StatementPrepper.SHORTUD_HCOUNT_NFE, var.getHomozygoteCount(Population.NFE)); //System.err.println("NFE Hom = " + var.getHomozygoteCount(Population.NFE));
		pstat.setInt(StatementPrepper.SHORTUD_ACOUNT_AFR, var.getIndividualCount(Population.AFR)); //System.err.println("AFR Count = " + var.getIndividualCount(Population.AFR));
		pstat.setInt(StatementPrepper.SHORTUD_HCOUNT_AFR, var.getHomozygoteCount(Population.AFR)); //System.err.println("AFR Hom = " + var.getHomozygoteCount(Population.AFR));
		pstat.setInt(StatementPrepper.SHORTUD_ACOUNT_AMR, var.getIndividualCount(Population.AMR)); //System.err.println("AMR Count = " + var.getIndividualCount(Population.AMR));
		pstat.setInt(StatementPrepper.SHORTUD_HCOUNT_AMR, var.getHomozygoteCount(Population.AMR)); //System.err.println("AMR Hom = " + var.getHomozygoteCount(Population.AMR));
		pstat.setInt(StatementPrepper.SHORTUD_ACOUNT_FIN, var.getIndividualCount(Population.FIN)); //System.err.println("FIN Count = " + var.getIndividualCount(Population.FIN));
		pstat.setInt(StatementPrepper.SHORTUD_HCOUNT_FIN, var.getHomozygoteCount(Population.FIN)); //System.err.println("FIN Hom = " + var.getHomozygoteCount(Population.FIN));
		pstat.setInt(StatementPrepper.SHORTUD_ACOUNT_EAS, var.getIndividualCount(Population.EAS)); //System.err.println("EAS Count = " + var.getIndividualCount(Population.EAS));
		pstat.setInt(StatementPrepper.SHORTUD_HCOUNT_EAS, var.getHomozygoteCount(Population.EAS)); //System.err.println("EAS Hom = " + var.getHomozygoteCount(Population.EAS));
		pstat.setInt(StatementPrepper.SHORTUD_ACOUNT_SAS, var.getIndividualCount(Population.SAS)); //System.err.println("SAS Count = " + var.getIndividualCount(Population.SAS));
		pstat.setInt(StatementPrepper.SHORTUD_HCOUNT_SAS, var.getHomozygoteCount(Population.SAS)); //System.err.println("SAS Hom = " + var.getHomozygoteCount(Population.SAS));
		pstat.setInt(StatementPrepper.SHORTUD_ACOUNT_ASJ, var.getIndividualCount(Population.ASJ)); //System.err.println("ASJ Count = " + var.getIndividualCount(Population.ASJ));
		pstat.setInt(StatementPrepper.SHORTUD_HCOUNT_ASJ, var.getHomozygoteCount(Population.ASJ)); //System.err.println("ASJ Hom = " + var.getHomozygoteCount(Population.ASJ));
		pstat.setInt(StatementPrepper.SHORTUD_ACOUNT_OTH, var.getIndividualCount(Population.OTH)); //System.err.println("OTH Count = " + var.getIndividualCount(Population.OTH));
		pstat.setInt(StatementPrepper.SHORTUD_HCOUNT_OTH, var.getHomozygoteCount(Population.OTH)); //System.err.println("OTH Hom = " + var.getHomozygoteCount(Population.OTH));
		
		if(threadlock) sqlManager.requestStatementExecution(true);
		Blob blob = sprepper.wrapInBlob(var.getGeneListAsBLOBBytes());
		pstat.setBlob(StatementPrepper.SHORTUD_GENELIST, blob);
		
		String valnotes = var.getValidationNotes();
		if(valnotes == null) valnotes = "N/A";
		//System.err.println("ValNotes = " + valnotes);
		pstat.setString(StatementPrepper.SHORTUD_VALNOTES, valnotes);
		
		blob = sprepper.wrapInBlob(vgeno.getGenotypesAsBLOBBytes());
		pstat.setBlob(StatementPrepper.SHORTUD_GENOTYPES, blob);
		
		pstat.setLong(StatementPrepper.SHORTUD_QUERYID, var.getLongID()); //System.err.println("VarUID = " + Long.toHexString(var.getLongID()));
		
		return pstat;
	}
	
	private PreparedStatement generatePopulationSetVarUpdateStatement(DBVariant dbv, boolean threadlock) throws SQLException
	{
		//PreparedStatement pstat = sprepper.getPopUpdateStatement();
		
		PreparedStatement pstat = null;
		StatementPrepper sprepper = sqlManager.getStatementGenerator();
		if(threadlock) pstat = sqlManager.requestStatement(SQLManager.SKEY_VAR_UPDATE_POP, true);
		else pstat = sprepper.getPopUpdateStatement();
		
		pstat.setInt(StatementPrepper.POPUD_ACOUNT_TOT, dbv.getIndividualCount());
		pstat.setInt(StatementPrepper.POPUD_HCOUNT_TOT, dbv.getHomozygoteCount());
		pstat.setInt(StatementPrepper.POPUD_ACOUNT_NFE, dbv.getIndividualCount(Population.NFE));
		pstat.setInt(StatementPrepper.POPUD_HCOUNT_NFE, dbv.getHomozygoteCount(Population.NFE));
		pstat.setInt(StatementPrepper.POPUD_ACOUNT_AFR, dbv.getIndividualCount(Population.AFR));
		pstat.setInt(StatementPrepper.POPUD_HCOUNT_AFR, dbv.getHomozygoteCount(Population.AFR));
		pstat.setInt(StatementPrepper.POPUD_ACOUNT_AMR, dbv.getIndividualCount(Population.AMR));
		pstat.setInt(StatementPrepper.POPUD_HCOUNT_AMR, dbv.getHomozygoteCount(Population.AMR));
		pstat.setInt(StatementPrepper.POPUD_ACOUNT_FIN, dbv.getIndividualCount(Population.FIN));
		pstat.setInt(StatementPrepper.POPUD_HCOUNT_FIN, dbv.getHomozygoteCount(Population.FIN));
		pstat.setInt(StatementPrepper.POPUD_ACOUNT_EAS, dbv.getIndividualCount(Population.EAS));
		pstat.setInt(StatementPrepper.POPUD_HCOUNT_EAS, dbv.getHomozygoteCount(Population.EAS));
		pstat.setInt(StatementPrepper.POPUD_ACOUNT_SAS, dbv.getIndividualCount(Population.SAS));
		pstat.setInt(StatementPrepper.POPUD_HCOUNT_SAS, dbv.getHomozygoteCount(Population.SAS));
		pstat.setInt(StatementPrepper.POPUD_ACOUNT_ASJ, dbv.getIndividualCount(Population.ASJ));
		pstat.setInt(StatementPrepper.POPUD_HCOUNT_ASJ, dbv.getHomozygoteCount(Population.ASJ));
		pstat.setInt(StatementPrepper.POPUD_ACOUNT_OTH, dbv.getIndividualCount(Population.OTH));
		pstat.setInt(StatementPrepper.POPUD_HCOUNT_OTH, dbv.getHomozygoteCount(Population.OTH));
		
		pstat.setLong(StatementPrepper.POPUD_QUERYID, dbv.getLongID());
		
		return pstat;
	}
	
	/* ----- Caching ----- */
	
	private void loadRegionIndex() throws IOException
	{
		String pathstem = dbDir + File.separator + REGIDX_FILESTEM;
		regidx = RegionIndex.loadIndex(pathstem, genome, true);
	}
	
	private void saveRegionIndex() throws IOException
	{
		String pathstem = dbDir + File.separator + REGIDX_FILESTEM;
		regidx.saveIndex(pathstem);
	}
	
	protected class ReadCache
	{
		public static final int CACHE_SIZE = 2048;
		
		private ConcurrentMap<Long, DBVariant> v_map;
		private ConcurrentMap<Long, VariantGenotype> g_map;
		
		private ConcurrentLinkedQueue<Long> v_queue;
		private ConcurrentLinkedQueue<Long> g_queue;
		
		private boolean thread_locked;
		
		public ReadCache(boolean threadlock)
		{
			v_map = new ConcurrentHashMap<Long, DBVariant>();
			g_map = new ConcurrentHashMap<Long, VariantGenotype>();
			v_queue = new ConcurrentLinkedQueue<Long>();
			g_queue = new ConcurrentLinkedQueue<Long>();
			thread_locked = threadlock;
		}
		
		public DBVariant getVariant(long uid)
		{
			DBVariant v = v_map.get(uid);
			if(v != null) 
			{
				v_queue.remove(uid);
				v_queue.add(uid);
				return v;
			}
			
			//Cache miss
			try 
			{
				//PreparedStatement pstat = sprepper.getVarGetterStatement();
				PreparedStatement pstat = null;
				StatementPrepper sprepper = sqlManager.getStatementGenerator();
				if(thread_locked) pstat = sqlManager.requestStatement(SQLManager.SKEY_GETVAR, true);
				else pstat = sprepper.getVarGetterStatement();
				
				pstat.setLong(StatementPrepper.VARGET_VARUID, uid);
				if(thread_locked) sqlManager.requestStatementExecution(true);
				ResultSet rs = pstat.executeQuery();
				if(!rs.next()) {rs.close(); return null;}
				DBVariant var = readFromResultSet(rs);
				rs.close();
				if(thread_locked) {
					sqlManager.acknowledgeStatementExecution();
					sqlManager.releaseStatement(SQLManager.SKEY_GETVAR, pstat);
				}
				
				//Cache
				if(v_map.size() >= CACHE_SIZE)
				{
					long id = v_queue.poll();
					v_map.remove(id);
				}
				v_map.put(var.getLongID(), var);
				v_queue.add(uid);
				
				//Return
				return var;
			} 
			catch (Exception e) 
			{
				e.printStackTrace();
				return null;
			}
		}
		
		public List<DBVariant> getVariants(Collection<Long> varUIDs)
		{
			Set<Long> misses = new TreeSet<Long>();
			List<DBVariant> out = new LinkedList<DBVariant>();
			
			for(Long vid : varUIDs)
			{
				DBVariant v = v_map.get(vid);
				if(v != null) 
				{
					v_queue.remove(vid);
					v_queue.add(vid);
					out.add(v);
				}
				else misses.add(vid);
			}
			
			//Get misses
			if(!misses.isEmpty())
			{
				int len = misses.size();
				//Reorder into queue
				LinkedList<Long> q = new LinkedList<Long>();
				for(Long vid : misses) q.add(vid);
				try 
				{
					//I'm also gonna try doing chunks at a time in case there are too many...
					int remaining = len;
					while(remaining > 0)
					{
						int grab = MAX_VARS_INQUERY;
						if(MAX_VARS_INQUERY > remaining) grab = remaining;
						
						if(thread_locked) sqlManager.requestStatementExecution(true);
						PreparedStatement pstat = sqlManager.getStatementGenerator().generateMultiVarGetterStatement(grab);
						//int i = 1;
						/*for(Long vid : misses)
						{
							pstat.setLong(i, vid);
							i++;
						}*/
						for(int i = 1; i <= grab; i++)
						{
							long pop = q.pop();
							pstat.setLong(i, pop);
						}
						ResultSet rs = pstat.executeQuery();
						while(rs.next())
						{
							DBVariant v = readFromResultSet(rs);
							
							if(v_map.size() >= CACHE_SIZE)
							{
								long id = v_queue.poll();
								v_map.remove(id);
							}
							
							v_map.put(v.getLongID(), v);
							v_queue.add(v.getLongID());
							out.add(v);
						}
						rs.close();
						if(thread_locked) sqlManager.acknowledgeStatementExecution();
						remaining -= grab;
					}
				} 
				catch (Exception e) 
				{
					e.printStackTrace();
					return out;
				}	
			}
			
			
			return out;
		}
		
		public VariantGenotype getGenotype(long var_uid)
		{
			VariantGenotype vg = g_map.get(var_uid);
			if(vg != null) 
			{
				g_queue.remove(var_uid);
				g_queue.add(var_uid);
				return vg;
			}
			
			//Cache miss
			try 
			{
				//PreparedStatement pstat = sprepper.getGenoGetterStatement();
				PreparedStatement pstat = null;
				StatementPrepper sprepper = sqlManager.getStatementGenerator();
				if(thread_locked) pstat = sqlManager.requestStatement(SQLManager.SKEY_GETGENO, true);
				else pstat = sprepper.getGenoGetterStatement();
				
				pstat.setLong(StatementPrepper.GENOGET_VARUID, var_uid);
				
				if(thread_locked) sqlManager.requestStatementExecution(true);
				ResultSet rs = pstat.executeQuery();
				if(!rs.next()) {rs.close(); return null;}
				Blob genoblob = rs.getBlob(FIELDNAME_GENOTYPES);
				vg = new VariantGenotype(var_uid);
				vg.readDataFromBLOB(genoblob);
				rs.close();
				if(thread_locked) {
					sqlManager.acknowledgeStatementExecution();
					sqlManager.releaseStatement(SQLManager.SKEY_GETGENO, pstat);
				}
				
				//Cache
				if(g_map.size() >= CACHE_SIZE)
				{
					long id = g_queue.poll();
					g_map.remove(id);
				}
				g_map.put(vg.getVariantUID(), vg);
				g_queue.add(vg.getVariantUID());
				
				return vg;
			} 
			catch (Exception e) 
			{
				e.printStackTrace();
				return null;
			}
		}
		
		public List<VariantGenotype> getGenotypes(Collection<Long> varUIDs)
		{
			Set<Long> misses = new TreeSet<Long>();
			List<VariantGenotype> out = new LinkedList<VariantGenotype>();
			
			for(Long vid : varUIDs)
			{
				VariantGenotype vg = g_map.get(vid);
				if(vg != null) 
				{
					g_queue.remove(vid);
					g_queue.add(vid);
					out.add(vg);
				}
				else misses.add(vid);
			}
			
			//Get misses
			if(!misses.isEmpty())
			{
				int len = misses.size();
				try 
				{
					PreparedStatement pstat = sqlManager.getStatementGenerator().generateMultiGenoGetterStatement(len);
					int i = 1;
					for(Long vid : varUIDs)
					{
						pstat.setLong(i, vid);
						i++;
					}
					
					if(thread_locked) sqlManager.requestStatementExecution(true);
					ResultSet rs = pstat.executeQuery();
					while(rs.next())
					{
						long varUID = rs.getLong(FIELDNAME_VARUID);
						Blob genoblob = rs.getBlob(FIELDNAME_GENOTYPES);
						VariantGenotype vg = new VariantGenotype(varUID);
						vg.readDataFromBLOB(genoblob);
						
						if(g_map.size() >= CACHE_SIZE)
						{
							long id = g_queue.poll();
							g_map.remove(id);
						}
						
						g_map.put(vg.getVariantUID(), vg);
						g_queue.add(vg.getVariantUID());
						
						out.add(vg);
					}
					rs.close();
					if(thread_locked) sqlManager.acknowledgeStatementExecution();
				} 
				catch (Exception e) 
				{
					e.printStackTrace();
					return out;
				}
			}
			return out;
		}
		
		public void clear()
		{
			v_map.clear();
			g_map.clear();
			v_queue.clear();
			g_queue.clear();
		}
		
		public void setThreadLocked(boolean threadlock)
		{
			thread_locked = threadlock;
		}
		
	}
	
	protected class RegionSearchCache
	{
		//public static final int MAX_REG_SIZE_HELD = 10000000;
		
		//private boolean threadlock;
		
		private Contig lastq_contig;
		private int lastq_start;
		private int lastq_end;
		
		private boolean lastq_notra;
		private SVType lastq_type;
		
		//private volatile Collection<Long> last_ids;
		private Collection<DBVariant> last_vars;
		
		private class DiskGopher implements Runnable
		{

			private Contig chr;
			private int start;
			private int end;
			
			private SVType type;
			private boolean ignoreTRA;
			
			private Collection<DBVariant> output;
			
			public DiskGopher(Contig c, int s, int e, boolean noTRA, SVType svtype, Collection<DBVariant> targ)
			{
				chr = c;
				start = s;
				end = e;
				output = targ;
				ignoreTRA = noTRA;
				type = svtype;
			}
			
			@Override
			public void run() 
			{
				try{
					Collection<DBVariant> col = null;
					if(type == null) col = getFromDisk(chr, start, end, ignoreTRA); 
					else col = getFromDisk(chr, start, end, type); 
				output.addAll(col);
				}
				catch(Exception e)
				{
					e.printStackTrace();
					return;
				}
			}
			
		}
		
		private Collection<DBVariant> scanLastQuery(int start, int end)
		{
			Collection<DBVariant> out = new LinkedList<DBVariant>();
			if(lastq_contig == null) return out;
			if(last_vars == null || last_vars.isEmpty()) return out;
			
			for(DBVariant v : last_vars)
			{
				if(v.getType() != SVType.TRA && v.getType() != SVType.BND)
				{
					int vst = v.getStartPosition().getStart();
					int ved = v.getEndPosition().getEnd();
					
					if(end <= vst) continue;
					if(start >= ved) continue;
					
					out.add(v);
				}
				else
				{
					Contig c1 = v.getChrom();
					if(lastq_contig.equals(c1))
					{
						int vst = v.getStartPosition().getStart();
						int ved = v.getStartPosition().getEnd();
						
						if(end > vst && start < ved) {out.add(v); continue;}
					}
					
					Contig c2 = v.getEndChrom();
					if(c2 == null) c2 = c1;
					if(lastq_contig.equals(c2))
					{
						int vst = v.getEndPosition().getStart();
						int ved = v.getEndPosition().getEnd();
						
						if(end > vst && start < ved) {out.add(v); continue;}
					}
				}
				
			}
			
			return out;
		}
		
		private Collection<DBVariant> getFromDisk(Contig c, int start, int end, boolean ignoreTRA)
		{
			List<DBVariant> varlist = new LinkedList<DBVariant>();
			if(c == null) return varlist;
			//System.err.println("-DEBUG- Disk Query-- " + c.getUDPName() + ":" + start + "-" + end);

			try 
			{
				//Index!
				if(regidx == null) loadRegionIndex();
				
				//Get variant UIDs from index...
				Collection<Long> vuids = regidx.getVariantIDsInApproximateRegion(c, start, end, ignoreTRA);
				
				//Load variants from SQL DB
				Collection<DBVariant> vars = getVariants(vuids);
				for(DBVariant v : vars)
				{
					//Double check to see if it's really in the requested region...
					if(!ignoreTRA && (v.getType() == SVType.TRA || v.getType() == SVType.BND))
					{
						if(v.getChrom().equals(c))
						{
							if(!((v.getStartPosition().getStart() >= end)&&(v.getStartPosition().getEnd() < start)))
							{
								varlist.add(v);
								continue;
							}
						}
						if (v.getEndChrom().equals(c))
						{
							if(!((v.getEndPosition().getStart() >= end)&&(v.getEndPosition().getEnd() < start)))
							{
								varlist.add(v);
								continue;
							}
						}
					}
					else
					{
						if(v.getStartPosition().getStart() >= end) continue;
						if(v.getEndPosition().getEnd() < start) continue;
						varlist.add(v);
					}
				}
			} 
			catch (Exception e) 
			{
				e.printStackTrace();
			}
			
			return varlist;
		}
		
		private Collection<DBVariant> getFromDisk(Contig c, int start, int end, SVType type)
		{
			List<DBVariant> varlist = new LinkedList<DBVariant>();
			if(c == null) return varlist;

			try 
			{
				if(regidx == null) loadRegionIndex();
				
				//Get variant UIDs from index...
				Collection<Long> vuids = regidx.getVariantIDsInApproximateRegion(type, c, start, end);
				
				//Load variants from SQL DB
				Collection<DBVariant> vars = getVariants(vuids);
				boolean typetra = type == SVType.TRA || type == SVType.BND;
				for(DBVariant v : vars)
				{
					//Double check to see if it's really in the requested region...
					if(typetra)
					{
						if(v.getChrom().equals(c))
						{
							if(!((v.getStartPosition().getStart() >= end)&&(v.getStartPosition().getEnd() < start)))
							{
								varlist.add(v);
								continue;
							}
						}
						if (v.getEndChrom().equals(c))
						{
							if(!((v.getEndPosition().getStart() >= end)&&(v.getEndPosition().getEnd() < start)))
							{
								varlist.add(v);
								continue;
							}
						}
					}
					else
					{
						if(v.getStartPosition().getStart() >= end) continue;
						if(v.getEndPosition().getEnd() < start) continue;
						varlist.add(v);
					}
				}
			} 
			catch (Exception e) 
			{
				e.printStackTrace();
			}
			
			return varlist;
		}
		
		public Collection<DBVariant> getVariantsInRegion(Contig c, int start, int end, boolean ignoreTRA)
		{
			//System.err.println("-DEBUG- Searching region-- " + c.getUDPName() + ":" + start + "-" + end);
			if(lastq_contig != null && lastq_contig.equals(c) && (lastq_notra == ignoreTRA) && (lastq_type == null))
			{
				//System.err.println("-DEBUG- Contig matches last search!");
				//Check for overlap
				if(start == lastq_start && end == lastq_end)
				{
					//Exact same query
					//if(lastq_end - lastq_start <= MAX_REG_SIZE_HELD) return last_vars;
					//else return getVariants(last_ids);
					return last_vars;
				}
				
				int in_s = -1;
				int in_e = -1;
				int ex_s = -1;
				int ex_e = -1;
				
				Collection<DBVariant> out = null;
				
				if(start <= lastq_start)
				{
					if(end >= lastq_end)
					{
						//Includes full last query
						out = new LinkedList<DBVariant>();
						out.addAll(last_vars);
						//Query front
						if(start < lastq_start) out.addAll(getFromDisk(c, start, lastq_start, ignoreTRA));
						//Query back
						if(end > lastq_end) out.addAll(getFromDisk(c, lastq_end, end, ignoreTRA));
					}
					else
					{
						//Includes first part
						in_s = lastq_start;
						in_e = end;
						if(start < lastq_start)
						{
							ex_s = start;
							ex_e = lastq_start;
						}
					}
				}
				else
				{
					if(start < lastq_end)
					{
						//Start is in between last query start and end
						in_s = start;
						in_e = lastq_end;
						if(end > lastq_end)
						{
							ex_s = lastq_end;
							ex_e = end;
						}
					}
					else
					{
						//Completely outside
						ex_s = start;
						ex_e = end;
					}
				}
				
				boolean doin = (in_s > -1 && in_e > -1);
				boolean doex = (ex_s > -1 && ex_e > -1);
				
				if(doin && doex)
				{
					//Do the DB fetching and cache scan in a different threads
					out = new ConcurrentLinkedQueue<DBVariant>();
					DiskGopher runner = new DiskGopher(c, ex_s, ex_e, ignoreTRA, null, out);
					Thread t = new Thread(runner);
					t.start();
					out.addAll(scanLastQuery(in_s, in_e));
					
					try {t.join();} 
					catch (InterruptedException e) 
					{
						//Shouldn't happen, maybe?
						e.printStackTrace();
					}
				}
				else if (doin && !doex)
				{
					out = scanLastQuery(in_s, in_e);
				}
				else if (doex && !doin)
				{
					out = getFromDisk(c, ex_s, ex_e, ignoreTRA);
				}
				
				lastq_start = start;
				lastq_end = end;
				last_vars = out;
				lastq_type = null;
				lastq_notra = ignoreTRA;
				return out;
			}
			
			//System.err.println("-DEBUG- New search!");
			//Completely new query
			Collection<DBVariant> out = getFromDisk(c, start, end, ignoreTRA);
			lastq_contig = c;
			lastq_start = start;
			lastq_end = end;
			last_vars = out;
			lastq_type = null;
			lastq_notra = ignoreTRA;
			return out;
		}
		
		public Collection<DBVariant> getVariantsInRegion(Contig c, int start, int end, SVType svtype)
		{
			//System.err.println("-DEBUG- Searching region-- " + c.getUDPName() + ":" + start + "-" + end);
			if(lastq_contig != null && lastq_contig.equals(c) && (lastq_type == svtype))
			{
				//System.err.println("-DEBUG- Contig matches last search!");
				//Check for overlap
				if(start == lastq_start && end == lastq_end)
				{
					//Exact same query
					//if(lastq_end - lastq_start <= MAX_REG_SIZE_HELD) return last_vars;
					//else return getVariants(last_ids);
					return last_vars;
				}
				
				int in_s = -1;
				int in_e = -1;
				int ex_s = -1;
				int ex_e = -1;
				
				Collection<DBVariant> out = null;
				
				if(start <= lastq_start)
				{
					if(end >= lastq_end)
					{
						//Includes full last query
						out = new LinkedList<DBVariant>();
						out.addAll(last_vars);
						//Query front
						if(start < lastq_start) out.addAll(getFromDisk(c, start, lastq_start, svtype));
						//Query back
						if(end > lastq_end) out.addAll(getFromDisk(c, lastq_end, end, svtype));
					}
					else
					{
						//Includes first part
						in_s = lastq_start;
						in_e = end;
						if(start < lastq_start)
						{
							ex_s = start;
							ex_e = lastq_start;
						}
					}
				}
				else
				{
					if(start < lastq_end)
					{
						//Start is in between last query start and end
						in_s = start;
						in_e = lastq_end;
						if(end > lastq_end)
						{
							ex_s = lastq_end;
							ex_e = end;
						}
					}
					else
					{
						//Completely outside
						ex_s = start;
						ex_e = end;
					}
				}
				
				boolean doin = (in_s > -1 && in_e > -1);
				boolean doex = (ex_s > -1 && ex_e > -1);
				
				if(doin && doex)
				{
					//Do the DB fetching and cache scan in a different threads
					out = new ConcurrentLinkedQueue<DBVariant>();
					DiskGopher runner = new DiskGopher(c, ex_s, ex_e, false, svtype, out);
					Thread t = new Thread(runner);
					t.start();
					out.addAll(scanLastQuery(in_s, in_e));
					
					try {t.join();} 
					catch (InterruptedException e) 
					{
						//Shouldn't happen, maybe?
						e.printStackTrace();
					}
				}
				else if (doin && !doex)
				{
					out = scanLastQuery(in_s, in_e);
				}
				else if (doex && !doin)
				{
					out = getFromDisk(c, ex_s, ex_e, svtype);
				}
				
				lastq_start = start;
				lastq_end = end;
				last_vars = out;
				lastq_notra = false;
				lastq_type = svtype;
				return out;
			}
			
			//System.err.println("-DEBUG- New search!");
			//Completely new query
			Collection<DBVariant> out = getFromDisk(c, start, end, svtype);
			lastq_contig = c;
			lastq_start = start;
			lastq_end = end;
			last_vars = out;
			lastq_notra = false;
			lastq_type = svtype;
			return out;
		}
		
		public void flush()
		{
			lastq_contig = null;
			lastq_start = -1;
			lastq_end = -1;
			last_vars = null;
		}
		
		public void setThreadLock(boolean b)
		{
			threadlock = b;
		}
		
	}
	
	public void indexByRegion() throws IOException
	{
		regidx = new RegionIndex(true);
		
		//Go through all variants...
		try
		{
			PreparedStatement statement = sqlManager.getStatementGenerator().getVariantGetAllStatement();
			ResultSet rs = statement.executeQuery();
			while(rs.next())
			{
				DBVariant var = readFromResultSet(rs);
				regidx.indexVariant(var);
			}
		}
		catch(SQLException x)
		{
			x.printStackTrace();
			throw new IOException();
		}
		
		saveRegionIndex();
	}
	
	/* ----- Multithreading ----- */
	
	/* ----- Gene Hit Cache ----- */
	
	private Map<Integer, GeneHitCounter> ghc_cache;
	private boolean ghc_cache_dirty;
	
	private boolean loadGeneHitTable()
	{
		ghc_cache = new ConcurrentHashMap<Integer, GeneHitCounter>();
		try
		{
			//int dbc = 0;
			PreparedStatement s = sqlManager.getStatementGenerator().getGeneHitGetAllStatement();
			ResultSet rs = s.executeQuery();
			while(rs.next())
			{
				int guid = rs.getInt(FIELDNAME_GH_GENEUID);
				GeneHitCounter ghc = readFromTableRecord(rs);
				ghc_cache.put(guid, ghc);
				//dbc++;
			}
			rs.close();
		//	System.err.println("-DEBUG: Gene hit records read: " + dbc);
		}
		catch(Exception e)
		{
			e.printStackTrace();
			return false;
		}
		
		return true;
	}
	
	private GeneHitCounter readFromTableRecord(ResultSet rs) throws SQLException
	{
		GeneHitCounter ghc = new GeneHitCounter();
		//ghc.total_hits_var = rs.getInt(FIELDNAME_GH_HITS_T);
		ghc.setTotalHits(rs.getInt(FIELDNAME_GH_HITS_T));
		//ghc.exon_hits_var = rs.getInt(FIELDNAME_GH_HITS_E);
		ghc.setExonHits(rs.getInt(FIELDNAME_GH_HITS_E));
		
		//Blobs
		Blob b = rs.getBlob(FIELDNAME_GH_HITS_TI);
		ByteBuffer bb = ByteBuffer.wrap(b.getBytes(1, (int)b.length()));
		bb.order(ByteOrder.BIG_ENDIAN);
		while(bb.hasRemaining())
		{
			int i = bb.getInt();
			if(i == -1) break;
			//ghc.total_hits_indiv.add(i);
			ghc.addIndivTotal(i);
		}
		
		b = rs.getBlob(FIELDNAME_GH_HITS_EI);
		bb = ByteBuffer.wrap(b.getBytes(1, (int)b.length()));
		bb.order(ByteOrder.BIG_ENDIAN);
		while(bb.hasRemaining())
		{
			int i = bb.getInt();
			if(i == -1) break;
			//ghc.total_hits_indiv.add(i);
			ghc.addIndivExon(i);
		}
		
		return ghc;
	}
	
	private boolean saveGeneHitTable()
	{
		if(ghc_cache == null) return false;
		
		int dbc = 0;
		boolean b = true;
		Integer badid = null;
		try 
		{	
			//Wipe table
			//It seems to be faster to just rewrite the table from scratch
			//(That way it doesn't have to look up the records to edit...)
			PreparedStatement delstatement = sqlManager.getStatementGenerator().getGeneHitTableWipeStatement();
			delstatement.executeUpdate();
			
			
			List<Integer> ids = new ArrayList<Integer>(ghc_cache.size() + 1);
			ids.addAll(ghc_cache.keySet());
			for(Integer id : ids)
			{
				GeneHitCounter ghc = ghc_cache.get(id);
				b = b && updateGeneHitRecord(id, ghc);
				dbc++;
				badid = id;
				//System.err.println("-DEBUG: Gene hit records written: " + dbc);
			}
			ghc_cache_dirty = false;
		}
		catch(Exception e)
		{
			e.printStackTrace();
			System.err.println("-DEBUG: Gene hit records written: " + dbc);
			System.err.println("-- DEBUG INFO --");
			System.err.println("Bad Gene UID Is Null?: " + (badid == null));
			if(badid != null)
			{
				System.err.println("Bad Gene UID: 0x" + Integer.toHexString(badid));
				Gene g = genes.getGeneByTranscriptGUID(badid);
				System.err.println("Bad Gene Is Null?: " + (g == null));
				System.err.println("Bad Gene: " + g.getName() + " (" + g.getID() + ")");
			}
			return false;
		}
		System.err.println("-DEBUG: Gene hit records written: " + dbc);
		
		return b;
	}
	
	private boolean updateGeneHitRecord(int guid, GeneHitCounter ghc) throws SQLException
	{
		//Prepare blobs
		StatementPrepper sprepper = sqlManager.getStatementGenerator();
		//int ti_count = ghc.total_hits_indiv.size();
		int ti_count = ghc.countIndivTotal();
		Blob ti_blob = null;
		if(ti_count < 1)
		{
			byte[] neg1 = {-1, -1, -1, -1};
			ti_blob = sprepper.wrapInBlob(neg1);
		}
		else
		{
			FileBuffer buff = new FileBuffer(ti_count * 4, true);
			Set<Integer> total_hits_indiv = ghc.getTotalHitsIndiv_setref();
			//for(Integer i : ghc.total_hits_indiv) buff.addToFile(i);
			for(Integer i : total_hits_indiv) buff.addToFile(i);
			ti_blob = sprepper.wrapInBlob(buff.getBytes());
		}
		
		//int ei_count = ghc.exon_hits_indiv.size();
		int ei_count = ghc.countIndivExon();
		Blob ei_blob = null;
		if(ei_count < 1)
		{
			byte[] neg1 = {-1, -1, -1, -1};
			ei_blob = sprepper.wrapInBlob(neg1);
		}
		else
		{
			FileBuffer buff = new FileBuffer(ei_count * 4, true);
			Set<Integer> exon_hits_indiv = ghc.getExonHitsIndiv_setref();
			//for(Integer i : ghc.exon_hits_indiv) buff.addToFile(i);
			for(Integer i : exon_hits_indiv) buff.addToFile(i);
			ei_blob = sprepper.wrapInBlob(buff.getBytes());
		}
		
		/*PreparedStatement ps = sprepper.getGeneHitUpdateStatement();
		
		ps.setInt(StatementPrepper.GENEHIT_UD_TOT, ghc.total_hits_var);
		ps.setInt(StatementPrepper.GENEHIT_UD_EXON, ghc.exon_hits_var);
		ps.setBlob(StatementPrepper.GENEHIT_UD_TOT_INDIV, ti_blob);
		ps.setBlob(StatementPrepper.GENEHIT_UD_EXON_INDIV, ei_blob);
		ps.setInt(StatementPrepper.GENEHIT_UD_UID, guid);*/
		
		PreparedStatement ps = sprepper.getGeneHitInsertStatement();
		
		ps.setInt(StatementPrepper.GENEHIT_UID, guid);
		ps.setInt(StatementPrepper.GENEHIT_TOT, ghc.getTotalHitsVar());
		ps.setInt(StatementPrepper.GENEHIT_EXON, ghc.getExonHitsVar());
		ps.setBlob(StatementPrepper.GENEHIT_TOT_INDIV, ti_blob);
		ps.setBlob(StatementPrepper.GENEHIT_EXON_INDIV, ei_blob);
		
		int count = ps.executeUpdate();
		
		return (count == 1);
	}
	
	public static boolean variantExonic(DBVariant var, Gene g)
	{
		if(var == null || g == null) return false;
		if(var.getType() == SVType.TRA || var.getType() == SVType.BND)
		{
			if(g.getRelativeRegionLocationEffect(var.getStartPosition().getStart(), var.getStartPosition().getEnd()) == GeneFunc.EXONIC) return true;
			if(g.getRelativeRegionLocationEffect(var.getEndPosition().getStart(), var.getEndPosition().getEnd()) == GeneFunc.EXONIC) return true;
		}
		else
		{
			if(g.getRelativeRegionLocationEffect(var.getStartPosition().getStart(), var.getEndPosition().getEnd()) == GeneFunc.EXONIC) return true;
		}
		
		return false;
	}
	
	protected GeneHitCounter getGeneHitRecord(int guid)
	{
		return ghc_cache.get(guid);
	}
	
	protected synchronized void tickGHCCDirty()
	{
		ghc_cache_dirty = true;
	}
	
	/* ----- Getters ----- */
	
	public boolean variantExists(long varUID)
	{
		try 
		{
			//PreparedStatement pstat = sprepper.getVarUIDCheckStatement();
			PreparedStatement pstat = null;
			if(threadlock) pstat = sqlManager.requestStatement(SQLManager.SKEY_CHECKVARUID, true);
			else pstat = sqlManager.getStatementGenerator().getVarUIDCheckStatement();
			pstat.setLong(StatementPrepper.VARUIDCHECK_VARUID, varUID);
			if(threadlock) sqlManager.requestStatementExecution(true);
			ResultSet rs = pstat.executeQuery();
			boolean b = rs.next();
			rs.close();
			
			if(threadlock)
			{
				sqlManager.acknowledgeStatementExecution();
				sqlManager.releaseStatement(SQLManager.SKEY_CHECKVARUID, pstat);
			}
			
			return b;
		} 
		catch (SQLException e) 
		{
			e.printStackTrace();
			return false;
		}
	}
	
	public DBVariant getVariant(long varUID) 
	{	
		return read_cache.getVariant(varUID);
	}

	public List<DBVariant> getVariants(Collection<Long> varUIDs)
	{
		if(varUIDs == null || varUIDs.isEmpty()) return new LinkedList<DBVariant>();
		return read_cache.getVariants(varUIDs);
	}
	
	@Override
	public VariantGenotype getGenotype(long varUID) 
	{
		return read_cache.getGenotype(varUID);
	}

	public List<VariantGenotype> getGenotypes(Collection<Long> varUIDs)
	{
		if(varUIDs == null || varUIDs.isEmpty()) return null;
		return read_cache.getGenotypes(varUIDs);
	}
	
	public List<Long> getVariantIDsForSample(int sampleUID) 
	{
		try 
		{
			//PreparedStatement pstat = sprepper.getSampleVarGetterStatement();
			PreparedStatement pstat = null;
			if(threadlock) pstat = sqlManager.requestStatement(SQLManager.SKEY_GETSAMPLEVAR, true);
			else pstat = sqlManager.getStatementGenerator().getSampleVarGetterStatement();
			
			pstat.setInt(StatementPrepper.SVARGET_SAMPUID, sampleUID);
			if(threadlock) sqlManager.requestStatementExecution(true);
			ResultSet rs = pstat.executeQuery();
			if(!rs.next()) {rs.close(); return null;}
			
			Set<Long> vars = new HashSet<Long>();
			vars.addAll(readVarUIDListBlob(rs.getBlob(FIELDNAME_SVARLIST_HET)));
			vars.addAll(readVarUIDListBlob(rs.getBlob(FIELDNAME_SVARLIST_HOM)));
			vars.addAll(readVarUIDListBlob(rs.getBlob(FIELDNAME_SVARLIST_OTH)));
			
			List<Long> list = new ArrayList<Long>(vars.size() + 1);
			list.addAll(vars);
			Collections.sort(list);
			rs.close();
			
			if(threadlock)
			{
				sqlManager.acknowledgeStatementExecution();
				sqlManager.releaseStatement(SQLManager.SKEY_GETSAMPLEVAR, pstat);
			}
			
			return list;
		}
		catch(SQLException e)
		{
			e.printStackTrace();
			return null;
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
			return null;
		}

	}

	@Override
	public List<Long> getVariantIDsForFamily(Family fam) 
	{
		if(fam == null) return null;
		
		Set<Long> set = new HashSet<Long>();
		List<FamilyMember> members = fam.getAllFamilyMembers();
		
		for(FamilyMember mem : members)
		{
			List<Long> indivVars = this.getVariantIDsForSample(mem.getUID());
			if(indivVars != null) set.addAll(indivVars);
		}
		
		List<Long> list = new ArrayList<Long>(set.size() + 1);
		list.addAll(set);
		Collections.sort(list);
		
		return list;
	}

	@Override
	public List<Long> getVariantIDsForSampleOfType(int sampleUID, SVType type) 
	{
		List<Long> list1 = this.getVariantIDsForSample(sampleUID);
		if(list1 == null) return null;
		
		List<Long> list2 = new LinkedList<Long>();
		for(Long varid : list1)
		{
			//Very clunky, but whatever.
			DBVariant dbv = this.getVariant(varid);
			if(dbv.getType() == type) list2.add(varid);
		}
		
		return list2;
	}

	public Collection<DBVariant> getVariantsInRegion(Contig c, int start, int end) 
	{
		return rs_cache.getVariantsInRegion(c, start, end, false);
	}
	
	public Collection<DBVariant> getVariantsInRegion(Contig c, int start, int end, boolean ignoreTRA) 
	{
		return rs_cache.getVariantsInRegion(c, start, end, ignoreTRA);
	}
	
	public Collection<DBVariant> getVariantsInRegion(Contig c, int start, int end, SVType type) 
	{
		return rs_cache.getVariantsInRegion(c, start, end, type);
	}

	public GeneHitCounter getGeneHitCounts(int transcriptUID)
	{
		if(ghc_cache == null) loadGeneHitTable();
		return ghc_cache.get(transcriptUID);
	}

	public GenomeBuild getGenomeBuild()
	{
		return this.genome;
	}
	
	public GeneSet getTranscriptSet()
	{
		return this.genes;
	}
	
	public double getPercentLeeway()
	{
		return percLeeway;
	}
	
	/* ----- Setters ----- */
	
	public void setThreadlock(boolean b)
	{
		threadlock = b;
		read_cache.setThreadLocked(b);
		rs_cache.setThreadLock(b);
	}
	
	/* ----- Analysis ----- */
	
	public boolean variantInRegion(DBVariant var, Contig c, int start, int end)
	{
		if(var == null) return false;
		if(c == null) return false;
		
		boolean tra = var.getType() == SVType.TRA || var.getType() == SVType.BND;
		if(tra)
		{
			//Check start
			if(c.equals(var.getChrom()))
			{
				boolean b = 
						(var.getStartPosition().getStart() < end) &&
						(var.getStartPosition().getEnd() > start);
				if(b) return true;
			}
			//Check end
			if(c.equals(var.getEndChrom()))
			{
				boolean b = 
						(var.getEndPosition().getStart() < end) &&
						(var.getEndPosition().getEnd() > start);
				if(b) return true;
			}
			return false;
		}
		else
		{
			if(!c.equals(var.getChrom())) return false;
			if(var.getEndPosition().getEnd() <= start) return false;
			if(var.getStartPosition().getStart() >= end) return false;
		}
		
		return true;
	}
	
	public boolean sampleHasVariantInGene(int sampleUID, Gene g)
	{
		if(g == null) return false;
		
		List<Long> svars = getVariantIDsForSample(sampleUID);
		if(svars == null) return false;
		
		Collection<DBVariant> regvars = this.getVariantsInRegion(g.getChromosome(), g.getTranscriptStart(), g.getTranscriptEnd());
		if(regvars == null || regvars.isEmpty()) return false;
		
		for(DBVariant v : regvars)
		{
			if(svars.contains(v.getLongID())) return true;
		}
		
		return false;
	}
	
	public boolean sampleHasExonicVariantInGene(int sampleUID, Gene g)
	{
		if(g == null) return false;
		
		List<Long> svars = getVariantIDsForSample(sampleUID);
		if(svars == null) return false;
		
		//Now do by exon...
		int ecount = g.getExonCount();
		for(int i = 0; i < ecount; i++)
		{
			Exon e = g.getExon(i);
			Collection<DBVariant> regvars = this.getVariantsInRegion(g.getChromosome(), e.getStart(), e.getEnd());
			if(regvars == null || regvars.isEmpty()) return false;
			
			for(DBVariant v : regvars)
			{
				if(svars.contains(v.getLongID())) return true;
			}
		}
		
		return false;
	}
	
	public boolean regenerateGeneHitMap()
	{
		//Analyzes variants to regenerate hit map anew
		
		//Clear cache
		if(ghc_cache == null)
		{
			ghc_cache = new ConcurrentHashMap<Integer, GeneHitCounter>();
		}
		else ghc_cache.clear();
		
		List<Gene> glist = genes.getAllGenes();
		int debugctr = 0;
		System.err.println("Analyzing genes... ");
		//Map<String, GeneHitCounter> ghmap = new HashMap<String, GeneHitCounter>();
		for(Gene g : glist)
		{
			//System.err.println("SQLVariantTable.generateGeneHitMap || -DEBUG- Loop: " + debugctr);
			GeneHitCounter ghc = new GeneHitCounter();
			//ghmap.put(g.getID(), ghc);
			ghc_cache.put(g.getGUID(), ghc);
			
			Collection<DBVariant> hits = getVariantsInRegion(g.getChromosome(), g.getTranscriptStart(), g.getTranscriptEnd());
			if(hits == null) continue;
			//ghc.total_hits_var = hits.size();
			ghc.setTotalHits(hits.size());
			
			List<DBVariant> exonhits = new LinkedList<DBVariant>();
			for(DBVariant var : hits)
			{
				if(g.getRelativeRegionLocationEffect(var.getStartPosition().getStart(), var.getEndPosition().getEnd()) == GeneFunc.EXONIC)
				{
					exonhits.add(var);
				}
			}
			//ghc.exon_hits_var = exonhits.size();
			ghc.setExonHits(exonhits.size());
			
			Map<Long, VariantGenotype> vgmap = new HashMap<Long, VariantGenotype>();
			for(DBVariant var : hits)
			{
				VariantGenotype vg = getGenotype(var.getLongID());
				//if(vg != null) ghc.total_hits_indiv.addAll(vg.getAllIndividuals());
				if(vg != null) ghc.addIndividuals_Total(vg.getAllIndividuals());
				vgmap.put(vg.getVariantUID(), vg);
			}
			
			for(DBVariant var : exonhits)
			{
				//VariantGenotype vg = getGenotype(var.getLongID());
				VariantGenotype vg = vgmap.get(var.getLongID());
				//if(vg != null) ghc.exon_hits_indiv.addAll(vg.getAllIndividuals());
				if(vg != null) ghc.addIndividuals_Exon(vg.getAllIndividuals());
			}
			debugctr++;
			if(debugctr%1000 == 0) System.err.println(debugctr + " transcripts analyzed...");
		}
		
		//save
		return saveGeneHitTable();
	}
	
	public Map<String, GeneHitCounter> generateGeneHitMap() 
	{
		if(ghc_cache == null) loadGeneHitTable();
		
		//Remap to string transcript IDs
		Map<String, GeneHitCounter> smap = new HashMap<String, GeneHitCounter>();
		List<Gene> glist = genes.getAllGenes();
		for(Gene g : glist)
		{
			smap.put(g.getID(), ghc_cache.get(g.getGUID()));
		}
		
		return smap;
	}
	
	public void dumpTable(String directory)
	{
		//TODO
		//To see what's inside
		
		//Doesn't need to thread lock, doesn't check
		String vtpath = directory + File.separator + "vartbl.csv";
		String stpath = directory + File.separator + "sgenotbl.csv";
		String gtpath = directory + File.separator + "genehittbl.csv";
		
		try
		{
			System.err.println("Dumping variants...");
			//PreparedStatement ps = sprepper.getVariantGetAllStatement();
			PreparedStatement ps = sqlManager.getStatementGenerator().getVariantGetAllStatement();
			ResultSet rs = ps.executeQuery();
			
			BufferedWriter bw = new BufferedWriter(new FileWriter(vtpath));
			//Header
			boolean first = true;
			for(String[] spair : VAR_COLUMNS) {
				if(!first)bw.write(",");
				bw.write(spair[0]); 
				first = false;
			}
			bw.write("\n");
			//Records
			while(rs.next())
			{
				DBVariant v = readFromResultSet(rs);
				bw.write("0x" + Long.toHexString(v.getLongID()) + ",");
				bw.write(v.getChrom().getUDPName() + ",");
				bw.write(v.getStartPosition().getStart() + ",");
				bw.write(v.getStartPosition().getEnd() + ",");
				bw.write(v.getEndPosition().getStart() + ",");
				bw.write(v.getEndPosition().getEnd() + ",");
				bw.write(v.getType().getString() + ",");
				bw.write(v.getPositionEffect().name() + ",");
				bw.write(v.getName() + ",");
				
				bw.write(v.getIndividualCount() + ",");
				bw.write(v.getHomozygoteCount() + ",");
				
				bw.write(v.getIndividualCount(Population.NFE) + ",");
				bw.write(v.getHomozygoteCount(Population.NFE) + ",");
				bw.write(v.getIndividualCount(Population.AFR) + ",");
				bw.write(v.getHomozygoteCount(Population.AFR) + ",");
				bw.write(v.getIndividualCount(Population.AMR) + ",");
				bw.write(v.getHomozygoteCount(Population.AMR) + ",");
				bw.write(v.getIndividualCount(Population.FIN) + ",");
				bw.write(v.getHomozygoteCount(Population.FIN) + ",");
				bw.write(v.getIndividualCount(Population.EAS) + ",");
				bw.write(v.getHomozygoteCount(Population.EAS) + ",");
				bw.write(v.getIndividualCount(Population.SAS) + ",");
				bw.write(v.getHomozygoteCount(Population.SAS) + ",");
				bw.write(v.getIndividualCount(Population.ASJ) + ",");
				bw.write(v.getHomozygoteCount(Population.ASJ) + ",");
				bw.write(v.getIndividualCount(Population.OTH) + ",");
				bw.write(v.getHomozygoteCount(Population.OTH) + ",");
				
				//Gene list...
				List<Gene> glist = v.getGeneListReference();
				if(glist == null || glist.isEmpty()) bw.write("[NONE],");
				else
				{
					first = true;
					Set<String> nset = new HashSet<String>();
					for(Gene g : glist)
					{
						//if(!first)bw.write(";");
						//first = false;
						//bw.write(g.getName() + "(" + g.getID() + ")");
						//bw.write(g.getName());
						nset.add(g.getName());
					}
					List<String> gnames = new ArrayList<String>(nset.size()+1);
					gnames.addAll(nset);
					Collections.sort(gnames);
					for(String n : gnames)
					{
						if(!first)bw.write(";");
						first = false;
						bw.write(n);
					}
					bw.write(",");
				}
				
				bw.write(v.getValidationNotes() + ",");
				if(v.getType() == SVType.TRA || v.getType() == SVType.BND)
				{
					Contig c2 = v.getEndChrom();
					if(c2 == null) c2 = v.getChrom();
					bw.write(c2.getUDPName() + ",");
				}
				else
				{
					bw.write("N/A,");
				}
				
				if(v.getType() == SVType.INS || v.getType() == SVType.INSME)
				{
					bw.write(v.getAltAlleleString() + ",");
				}
				else
				{
					bw.write("N/A,");
				}
				
				//Genotypes
				VariantGenotype vg = getGenotype(v.getLongID());
				if(vg.isCorrupted()) bw.write("CORRUPTED|");
				Collection<Integer> sids = vg.getAllIndividuals();
				first = true;
				for(Integer sid : sids)
				{
					Genotype g = vg.getAsGenotypeObject(sid, v.getType(), true);
					if(!first)bw.write(";");
					first = false;
					bw.write(Integer.toHexString(sid) + "=" + g.getField("GT"));
				}
				
				bw.write("\n");
			}
			rs.close();
			bw.close();
			
			System.err.println("Dumping sample genotype table...");
			//ps = sprepper.getSGenoGetAllStatement();
			ps = sqlManager.getStatementGenerator().getSGenoGetAllStatement();
			rs = ps.executeQuery();
			
			bw = new BufferedWriter(new FileWriter(stpath));
			//Header
			first = true;
			for(String[] spair : SAMPLEGENO_COLUMNS) {
				if(!first)bw.write(",");
				bw.write(spair[0]); 
				first = false;
			}
			bw.write("\n");
			
			while(rs.next())
			{
				bw.write("0x" + Integer.toHexString(rs.getInt(FIELDNAME_SAMPLEUID)) + ",");
		
				//vars.addAll(readVarUIDListBlob(rs.getBlob(FIELDNAME_SVARLIST_HET)));
				List<Long> vlist = this.readVarUIDListBlob(rs.getBlob(FIELDNAME_SVARLIST_HOM));
				if(vlist == null || vlist.isEmpty()) bw.write("[NONE]");
				else
				{
					first = true;
					for(Long l : vlist)
					{
						if(!first)bw.write(";");
						first = false;
						bw.write(Long.toHexString(l));
					}
				}
				bw.write(",");
				
				vlist = this.readVarUIDListBlob(rs.getBlob(FIELDNAME_SVARLIST_HET));
				if(vlist == null || vlist.isEmpty()) bw.write("[NONE]");
				else
				{
					first = true;
					for(Long l : vlist)
					{
						if(!first)bw.write(";");
						first = false;
						bw.write(Long.toHexString(l));
					}
				}
				bw.write(",");
				
				vlist = this.readVarUIDListBlob(rs.getBlob(FIELDNAME_SVARLIST_OTH));
				if(vlist == null || vlist.isEmpty()) bw.write("[NONE]");
				else
				{
					first = true;
					for(Long l : vlist)
					{
						if(!first)bw.write(";");
						first = false;
						bw.write(Long.toHexString(l));
					}
				}
				bw.write("\n");
				
			}
			rs.close();
			bw.close();
			
			System.err.println("Dumping gene hit table...");
			//ps = sprepper.getGeneHitGetAllStatement();
			ps = sqlManager.getStatementGenerator().getGeneHitGetAllStatement();
			rs = ps.executeQuery();
			
			bw = new BufferedWriter(new FileWriter(gtpath));
			//Header
			first = true;
			for(String[] spair : GENEHITS_COLUMNS) {
				if(!first)bw.write(",");
				bw.write(spair[0]); 
				first = false;
			}
			bw.write("GENE_NAME,");
			bw.write("TRANSCRIPT_ID\n");
			
			while(rs.next())
			{
				int guid = rs.getInt(FIELDNAME_GH_GENEUID);
				Gene g = genes.getGeneByTranscriptGUID(guid);
				
				bw.write("0x" + Integer.toHexString(guid) + ",");
				
				//bw.write(rs.getInt(FIELDNAME_GH_HITS_T) + ",");
				//bw.write(rs.getInt(FIELDNAME_GH_HITS_E) + ",");
				GeneHitCounter ghc = this.readFromTableRecord(rs);
				//bw.write(ghc.total_hits_var + ",");
				//bw.write(ghc.exon_hits_var + ",");
				
				int total_hits_var = ghc.getTotalHitsVar();
				int exon_hits_var = ghc.getExonHitsVar();
				bw.write(total_hits_var + ",");
				bw.write(exon_hits_var + ",");
				
				Set<Integer> total_hits_indiv = ghc.getTotalHitsIndiv_setref();
				Set<Integer> exon_hits_indiv = ghc.getExonHitsIndiv_setref();
				
				if(total_hits_indiv == null || total_hits_indiv.isEmpty()) bw.write("[NONE],");
				else
				{
					first = true;
					for(Integer i : total_hits_indiv)
					{
						if(!first) bw.write(";");
						first = false;
						bw.write(Integer.toHexString(i));
					}
					bw.write(",");
				}
				
				if(exon_hits_indiv == null || exon_hits_indiv.isEmpty()) bw.write("[NONE],");
				else
				{
					first = true;
					for(Integer i : exon_hits_indiv)
					{
						if(!first) bw.write(";");
						first = false;
						bw.write(Integer.toHexString(i));
					}
					bw.write(",");
				}
				
				bw.write(g.getName() + ",");
				bw.write(g.getID() + "\n");
				
			}
			
			rs.close();
			bw.close();
			
		}
		catch(Exception e)
		{
			e.printStackTrace();
			return;
		}
	}
	
	/* ----- Variant Addition ----- */
	
	protected long generateUID(Contig c, int pos)
	{
		long prefix = uidIndex.getVarIDPrefix(c, pos);
		//Generate random suffixes until the ID is unique...
		Random r = new Random();
		long suffix = r.nextLong() & 0xFFFFFFFFL;
		long id = prefix | suffix;
		while(variantExists(id) || id == -1L)
		{
			suffix = r.nextLong() & 0xFFFFFFFFL;
			id = prefix | suffix;
		}
		
		return id;
	}
	
	protected boolean addOrUpdateVariant(DBVariant var, VariantGenotype vgeno, boolean isnew)
	{
		try
		{
			PreparedStatement baseStatement = null;
			if(isnew) baseStatement = generateFullVarInsertStatement(var, vgeno, threadlock);
			else {
				var.noteGenes(genes);
				baseStatement = generateAbridgedSetVarUpdateStatement(var, vgeno, threadlock);
			}
			
			//Statement cstat = connection.createStatement();
			//int count = cstat.executeUpdate(baseStatement);
			//if(threadlock) sqlManager.requestStatementExecution(true);
			int count = baseStatement.executeUpdate();
			if(threadlock)
			{
				sqlManager.acknowledgeStatementExecution();
				String skey = null;
				if(isnew) skey = SQLManager.SKEY_VAR_INSERT;
				else skey = SQLManager.SKEY_VAR_UPDATE_SHORT;
				sqlManager.releaseStatement(skey, baseStatement);
			}
			if(count != 1) return false;
			
			if(regidx == null) this.loadRegionIndex();
			regidx.indexVariant(var);
		}
		catch(Exception e)
		{
			e.printStackTrace();
			return false;
		}
		
		//Update gene hit counts...
		/*List<Gene> glist = var.getGeneListReference();
		if(glist != null && !glist.isEmpty())
		{
			for(Gene g : glist)
			{
				GeneHitCounter ghc = this.getGeneHitCounts(g.getGUID());
				if(ghc == null) continue;
				ghc.total_hits_indiv.addAll(vgeno.getAllIndividuals());
				if(isnew) ghc.total_hits_var++;
				//If exonic...
				boolean exonic = false;
				if(var.getType() == SVType.TRA || var.getType() == SVType.BND)
				{
					exonic = g.getRelativeRegionLocationEffect(var.getStartPosition().getStart(), var.getStartPosition().getEnd()) == GeneFunc.EXONIC;
					if(!exonic) exonic = g.getRelativeRegionLocationEffect(var.getEndPosition().getStart(), var.getEndPosition().getEnd()) == GeneFunc.EXONIC;
				}
				else
				{
					exonic = g.getRelativeRegionLocationEffect(var.getStartPosition().getStart(), var.getEndPosition().getEnd()) == GeneFunc.EXONIC;
				}
				
				if(exonic) 
				{
					ghc.exon_hits_indiv.addAll(vgeno.getAllIndividuals());
					if(isnew) ghc.exon_hits_var++;
				}
				ghc_cache_dirty = true;
			}
		}*/
		
		return true;
	}
	
	@Override
	public boolean addVCF(String vcfpath, Map<String, FamilyMember> sampleMap, int mergeFactor, boolean ignoreTRA, int threads) throws IOException 
	{
		if(vcfpath == null) return false;
		if(sampleMap == null) return false;
		if (mergeFactor < 0) return false;
		
		SQLVTVarAdder add_engine = new SQLVTVarAdder(this, sampleMap, true);
		this.mergeFactor = mergeFactor;
		this.percLeeway = (double)mergeFactor / 1000.0;
		
		long time = 0;
		try
		{
			add_engine.start(vcfpath, ignoreTRA);
			while(!add_engine.isDone())
			{
				Thread.sleep(1000);
				time++;
				//add_engine.printQueueLoads();
				//add_engine.printCounts();
			}
		}
		catch(Exception e)
		{
			e.printStackTrace();
			return false;
		}
		
		System.err.println("VCF Addition completed! App. time elapsed: " + time + " seconds");
		
		return true;
	}
	
	/* ----- Variant Deletion ----- */
	
	private boolean removeSampleGenoRecord(int sampleID)
	{
		try 
		{
			//Statement cstat = connection.createStatement();
			//PreparedStatement pstat = sprepper.getSampleGenoDeleteStatment();
			PreparedStatement pstat = null;
			if(threadlock) pstat = sqlManager.requestStatement(SQLManager.SKEY_SGENO_DELETE, true);
			else pstat = sqlManager.getStatementGenerator().getSampleGenoDeleteStatment();
			
			pstat.setInt(1, sampleID);
			if(threadlock) sqlManager.requestStatementExecution(true);
			int count = pstat.executeUpdate();
			if(threadlock)
			{
				sqlManager.acknowledgeStatementExecution();
				sqlManager.releaseStatement(SQLManager.SKEY_SGENO_DELETE, pstat);
			}
			return (count == 1);
		} 
		catch (SQLException e) 
		{
			e.printStackTrace();
			return false;
		}
	}
	
	public boolean removeVariant(long varUID)
	{
		//Get variant so can remove gene hit info
		DBVariant var = getVariant(varUID);
		VariantGenotype vg = getGenotype(varUID);
		
		try 
		{
			//Statement cstat = connection.createStatement();
			//int count = cstat.executeUpdate(sqlQuery);
			//PreparedStatement pstat = sprepper.getVarDeleteStatment();
			PreparedStatement pstat = null;
			if(threadlock) pstat = sqlManager.requestStatement(SQLManager.SKEY_VAR_DELETE, true);
			else pstat = sqlManager.getStatementGenerator().getVarDeleteStatment();
			
			pstat.setLong(1, varUID);
			
			if(threadlock) sqlManager.requestStatementExecution(true);
			int count = pstat.executeUpdate();
			if(threadlock)
			{
				sqlManager.acknowledgeStatementExecution();
				sqlManager.releaseStatement(SQLManager.SKEY_VAR_DELETE, pstat);
			}
			
			//Update sample geno table...
			if(vg != null)
			{
				Set<Integer> samps = vg.getAllIndividuals();
				Map<Integer, List<Long>> removes = new HashMap<Integer, List<Long>>();
				for(Integer sid : samps)
				{
					List<Long> list = new LinkedList<Long>();
					list.add(varUID);
					removes.put(sid, list);
				}
				updateSampleGenoTable(null, removes);
			}
			
			List<Gene> vgenes = var.getGeneListReference();
			if(vgenes != null)
			{
				for(Gene g : vgenes)
				{
					ghc_cache_dirty = true;
					GeneHitCounter ghc = getGeneHitCounts(g.getGUID());
					boolean exonic = variantExonic(var, g);
					//ghc.total_hits_var--;
					//if(exonic) ghc.exon_hits_var--;
					
					ghc.decrementTotalHits();
					if(exonic) ghc.decrementExonHits();
					
					for(Integer sid : vg.getAllIndividuals())
					{
						//if(!this.sampleHasVariantInGene(sid, g)) ghc.total_hits_indiv.remove(sid);
						//if(exonic && !this.sampleHasExonicVariantInGene(sid, g)) ghc.exon_hits_indiv.remove(sid);
						if(!this.sampleHasVariantInGene(sid, g)) ghc.removeIndividual_Total(sid);
						if(exonic && !this.sampleHasExonicVariantInGene(sid, g)) ghc.removeIndividual_Exon(sid);
					}
				}
			}
			
			return (count == 1);
		} 
		catch (SQLException e) 
		{
			e.printStackTrace();
			return false;
		}
	}
	
	public boolean removeVariants(Collection<Long> varUIDs)
	{
		if(varUIDs == null || varUIDs.isEmpty()) return false;
		int sz = varUIDs.size();
		
		Collection<VariantGenotype> vglist = getGenotypes(varUIDs);
		Collection<DBVariant> varlist = getVariants(varUIDs);
		
		int len = varUIDs.size();
		
		try 
		{
			//PreparedStatement pstat = sprepper.generateMultiVarDeleteStatement(len);
			PreparedStatement pstat = sqlManager.getStatementGenerator().generateMultiVarDeleteStatement(len);
			int i = 1;
			for(Long vid : varUIDs) 
			{
				pstat.setLong(i, vid);
				i++;
			}
			if(threadlock) sqlManager.requestStatementExecution(true);
			int count = pstat.executeUpdate();
			if(threadlock) sqlManager.acknowledgeStatementExecution();

			//Update sample geno table...
			if(vglist != null)
			{
				Map<Integer, List<Long>> removes = new HashMap<Integer, List<Long>>();
				for(VariantGenotype vg : vglist)
				{
					Set<Integer> samps = vg.getAllIndividuals();
					for(Integer sid : samps)
					{
						List<Long> list = new LinkedList<Long>();
						list.add(vg.getVariantUID());
						removes.put(sid, list);
					}
				}
				updateSampleGenoTable(null, removes);	
			}
			
			for(DBVariant v : varlist)
			{
				List<Gene> vgenes = v.getGeneListReference();
				if(vgenes != null && !vgenes.isEmpty())
				{
					//Gene genotype
					VariantGenotype vg = null;
					for(VariantGenotype itr : vglist)
					{
						if(itr.getVariantUID() == v.getLongID())
						{
							vg = itr;
							break;
						}
					}
					for(Gene g : vgenes)
					{
						ghc_cache_dirty = true;
						GeneHitCounter ghc = getGeneHitCounts(g.getGUID());
						boolean exonic = variantExonic(v, g);
						//ghc.total_hits_var--;
						//if(exonic) ghc.exon_hits_var--;
						
						ghc.decrementTotalHits();
						if(exonic) ghc.decrementExonHits();
						
						for(Integer sid : vg.getAllIndividuals())
						{
							//if(!this.sampleHasVariantInGene(sid, g)) ghc.total_hits_indiv.remove(sid);
							//if(exonic && !this.sampleHasExonicVariantInGene(sid, g)) ghc.exon_hits_indiv.remove(sid);
							if(!this.sampleHasVariantInGene(sid, g)) ghc.removeIndividual_Total(sid);
							if(exonic && !this.sampleHasExonicVariantInGene(sid, g)) ghc.removeIndividual_Exon(sid);
						}
					}
				}
			}
			
			return (count == sz);
		} 
		catch (SQLException e) 
		{
			e.printStackTrace();
			return false;
		}
		
	}
	
	public boolean removeSample(FamilyMember sample) 
	{
		if(sample == null) return true;
		int uid = sample.getUID();
		
		boolean b = true;
		List<Long> sampleVars = getVariantIDsForSample(uid);
		//See which ones only have this sample under genotypes
		List<VariantGenotype> vglist = getGenotypes(sampleVars);
		
		//Remove sample geno record
		b = b && removeSampleGenoRecord(uid);
		
		//Remove individual from gene hit counts
		ghc_cache_dirty = true;
		List<Gene> glist = genes.getAllGenes();
		for(Gene g : glist)
		{
			GeneHitCounter ghc = getGeneHitCounts(g.getGUID());
			//ghc.exon_hits_indiv.remove(uid);
			//ghc.total_hits_indiv.remove(uid);
			ghc.removeIndividual_Total(uid);
			ghc.removeIndividual_Exon(uid);
		}
		
		//Remove variants or geno records
		if(vglist != null)
		{
			for(VariantGenotype vg : vglist)
			{
				SVDBGenotype gt = vg.removeGenotype(uid);
				if(vg.isEmpty()) 
				{
					b = b && removeVariant(vg.getVariantUID());
				}
				else
				{
					//Pull variant and update counts
					DBVariant dbv = getVariant(vg.getVariantUID());
					boolean hom = gt.isHomozygous();
					dbv.decrementTotalCount();
					if(hom) dbv.decrementHomozygoteCount();
					Collection<Population> poplist = sample.getPopulationTags();
					for(Population p : poplist)
					{
						dbv.decrementTotalCount(p);
						if(hom)dbv.decrementHomozygoteCount(p);
					}
					
					b = b && addOrUpdateVariant(dbv, vg, false);
				}
			}
		}
		
		return b;
	}

	@Override
	public boolean removeFamily(Family fam) 
	{
		if(fam == null) return false;
		boolean allgood = true;
		List<FamilyMember> members = fam.getAllFamilyMembers();
		for(FamilyMember mem : members)
		{
			allgood = allgood && removeSample(mem);
		}
		return allgood;
	}

	/* ----- Variant Modification ----- */
	
	private boolean updateSampleGenoTable(Map<Integer, List<Long>> additions, Map<Integer, List<Long>> removes)
	{
		if(additions == null && removes == null) return true;
		
		Set<Integer> sSet = new HashSet<Integer>();
		if(additions != null) sSet.addAll(additions.keySet());
		if(removes != null) sSet.addAll(removes.keySet());
		
		try 
		{
			for(Integer sid : sSet)
			{
				List<Long> a = null;
				List<Long> d = null;
				if(additions != null) a = additions.get(sid);
				if(removes != null) d = removes.get(sid);
			
				//Fetch record, if there
				List<Long> nowvars = this.getVariantIDsForSample(sid);
				if(nowvars == null) 
				{
					Set<Long> vidset = new HashSet<Long>();
					if(a != null && !a.isEmpty())
					{
						vidset.addAll(a);
						//Create record and add all to blob!
						FileBuffer hetblob = new FileBuffer(8, true);
						FileBuffer homblob = new FileBuffer(8, true);
						FileBuffer othblob = new FileBuffer(vidset.size()*8, true);
						
						hetblob.addToFile(-1L);
						homblob.addToFile(-1L);
						for(Long l : vidset) othblob.addToFile(l);
						
						//Statement cstat = connection.createStatement();
						//cstat.executeUpdate(sqlCmd);
						StatementPrepper sprepper = sqlManager.getStatementGenerator();
						PreparedStatement pstat = sprepper.getSGenoInsertStatement();
						pstat.setInt(StatementPrepper.SGENOINS_SID, sid);
						
						Blob b = sprepper.wrapInBlob(homblob.getBytes());
						pstat.setBlob(StatementPrepper.SGENOINS_HOM, b);
						
						b = sprepper.wrapInBlob(hetblob.getBytes());
						pstat.setBlob(StatementPrepper.SGENOINS_HET, b);
						
						b = sprepper.wrapInBlob(othblob.getBytes());
						pstat.setBlob(StatementPrepper.SGENOINS_OTH, b);
						
						pstat.executeUpdate();
					}
				}
				else
				{
					//Make a set by combining all three...
					Set<Long> vidset = new HashSet<Long>();
					vidset.addAll(nowvars);
					if(a != null && !a.isEmpty()) vidset.addAll(a);
					if(d != null && !d.isEmpty()) vidset.removeAll(d);
					
					FileBuffer hetblob = new FileBuffer(8, true);
					FileBuffer homblob = new FileBuffer(8, true);
					FileBuffer othblob = new FileBuffer(vidset.size()*8, true);

					hetblob.addToFile(-1L);
					homblob.addToFile(-1L);
					for(Long l : vidset) othblob.addToFile(l);
					
					//Statement cstat = connection.createStatement();
					//cstat.executeUpdate(sqlCmd);
					StatementPrepper sprepper = sqlManager.getStatementGenerator();
					PreparedStatement pstat = sprepper.getSGenoUpdateStatement();
					pstat.setInt(StatementPrepper.SGENOUD_SID, sid);
					
					Blob b = sprepper.wrapInBlob(homblob.getBytes());
					pstat.setBlob(StatementPrepper.SGENOUD_HOM, b);
					
					b = sprepper.wrapInBlob(hetblob.getBytes());
					pstat.setBlob(StatementPrepper.SGENOUD_HET, b);
					
					b = sprepper.wrapInBlob(othblob.getBytes());
					pstat.setBlob(StatementPrepper.SGENOUD_OTH, b);
					
					pstat.executeUpdate();
				}
				
			
			}
		}
		catch(Exception e)
		{
			e.printStackTrace();
			return false;
		}
		
		return true;
	}
	
	private boolean updateSampleCountsForVariant(DBVariant dbv, DBSampleTable sampleTable)
	{
		boolean allgood = true;
		
		dbv.setTotalCount(0);
		dbv.setHomozygoteCount(0);
		
		Population[] allpops = Population.values();
		for(Population p : allpops)
		{
			dbv.setTotalCount(0, p);
			dbv.setHomozygoteCount(0, p);
		}
		
		VariantGenotype vg = getGenotype(dbv.getLongID());
		Set<Integer> samples = vg.getAllIndividuals();
		for(Integer sid : samples)
		{
			SVDBGenotype gt = vg.getGenotype(sid);
			FamilyMember mem = sampleTable.getSample(sid);
			if(mem == null || gt == null) continue;
			
			Collection<Population> pTags = mem.getPopulationTags();
			boolean hom = gt.isHomozygous();
			
			dbv.incrementTotalCount();
			if(hom)dbv.incrementHomozygoteCount();
			for(Population p : pTags)
			{
				dbv.incrementTotalCount(p);
				if(hom)dbv.incrementHomozygoteCount(p);
			}
		}
		
		try
		{
			//Statement cstat = connection.createStatement();
			//int count = cstat.executeUpdate(sqlcmd);
			PreparedStatement pstat = generatePopulationSetVarUpdateStatement(dbv, threadlock);
			if(threadlock) sqlManager.requestStatementExecution(true);
			int count = pstat.executeUpdate();
			if(threadlock)
			{
				sqlManager.acknowledgeStatementExecution();
				sqlManager.releaseStatement(SQLManager.SKEY_VAR_UPDATE_POP, pstat);
			}
			allgood = allgood && (count == 1);
		}
		catch(Exception e)
		{
			e.printStackTrace();
			return false;
		}
	
		return allgood;
	}
	
	@Override
	public boolean updateSampleCounts(DBSampleTable sampleTable) 
	{
		String sqlQuery = "SELECT * FROM " + TABLENAME_VARIANTS;
		boolean allgood = true;
		
		try 
		{
			Statement cstat = connection.createStatement();
			ResultSet rs = cstat.executeQuery(sqlQuery);
			while(rs.next())
			{
				DBVariant dbv = readFromResultSet(rs);
				allgood = allgood && updateSampleCountsForVariant(dbv, sampleTable);
			}

		} 
		catch (Exception e) 
		{
			e.printStackTrace();
			return false;
		}
		
		return true;
	}

	@Override
	public boolean updateSamplePopulationFlags(FamilyMember sample, Collection<Population> oldFlags) 
	{
		if(sample == null) return true;
		
		boolean allgood = true;
		
		List<Long> varlist = this.getVariantIDsForSample(sample.getUID());
		for(Long vid : varlist)
		{
			DBVariant dbv = getVariant(vid);
			VariantGenotype vg = getGenotype(vid);
			SVDBGenotype gt = vg.getGenotype(sample.getUID());
			boolean hom = gt.isHomozygous();
			
			if(oldFlags != null)
			{
				for(Population p : oldFlags)
				{
					dbv.decrementTotalCount(p);
					if(hom) dbv.decrementHomozygoteCount(p);
				}
			}
			
			//New flags
			Collection<Population> tags = sample.getPopulationTags();
			for(Population p : tags)
			{
				dbv.incrementTotalCount(p);
				if(hom) dbv.incrementHomozygoteCount(p);
			}
			
			try
			{
				PreparedStatement pstat = generatePopulationSetVarUpdateStatement(dbv, true);
				
				if(threadlock) sqlManager.requestStatementExecution(true);
				int count = pstat.executeUpdate();
				if(threadlock)
				{
					sqlManager.acknowledgeStatementExecution();
					sqlManager.releaseStatement(SQLManager.SKEY_VAR_UPDATE_POP, pstat);
				}
				allgood = allgood && (count == 1);
			}
			catch(Exception e)
			{
				e.printStackTrace();
				return false;
			}
		}
		
		return allgood;
	}
	
	public boolean updateSampleGenotypeTable() throws SQLException
	{
		Map<Integer, List<Long>> map = new TreeMap<Integer, List<Long>>();
		PreparedStatement statement = sqlManager.getStatementGenerator().getVariantGetAllStatement();
		
		ResultSet rs = statement.executeQuery();
		while(rs.next())
		{
			long vuid = rs.getLong(FIELDNAME_VARUID);
			Blob genoblob = rs.getBlob(FIELDNAME_GENOTYPES);
			VariantGenotype vg = new VariantGenotype(vuid);
			try {vg.readDataFromBLOB(genoblob);} 
			catch (IOException e) {e.printStackTrace(); return false;}
			if(vg.isCorrupted())continue;
			Set<Integer> indivs = vg.getAllIndividuals();
			for(int sid : indivs)
			{
				List<Long> list = map.get(sid);
				if(list == null)
				{
					list = new LinkedList<Long>();
					map.put(sid, list);
				}
				list.add(vuid);
			}
		}
		
		//Rewrite sample geno table...
		statement = sqlManager.getStatementGenerator().getSampleGenoTableWipeStatement();
		statement.executeUpdate();
		Set<Integer> indivs = map.keySet();
		for(int sid : indivs)
		{
			Set<Long> vidset = new TreeSet<Long>();
			List<Long> list = map.get(sid);
			if(list != null && !list.isEmpty())
			{
				vidset.addAll(list);
				//Create record and add all to blob!
				FileBuffer hetblob = new FileBuffer(8, true);
				FileBuffer homblob = new FileBuffer(8, true);
				FileBuffer othblob = new FileBuffer(vidset.size()*8, true);
				
				hetblob.addToFile(-1L);
				homblob.addToFile(-1L);
				for(Long l : vidset) othblob.addToFile(l);
				
				//Statement cstat = connection.createStatement();
				//cstat.executeUpdate(sqlCmd);
				StatementPrepper sprepper = sqlManager.getStatementGenerator();
				PreparedStatement pstat = sprepper.getSGenoInsertStatement();
				pstat.setInt(StatementPrepper.SGENOINS_SID, sid);
				
				Blob b = sprepper.wrapInBlob(homblob.getBytes());
				pstat.setBlob(StatementPrepper.SGENOINS_HOM, b);
				
				b = sprepper.wrapInBlob(hetblob.getBytes());
				pstat.setBlob(StatementPrepper.SGENOINS_HET, b);
				
				b = sprepper.wrapInBlob(othblob.getBytes());
				pstat.setBlob(StatementPrepper.SGENOINS_OTH, b);
				
				pstat.executeUpdate();
			}
		}
		
		return true;
	}
	
	/* ----- Cleanup ----- */
	
	public void commitUpdates() throws SQLException
	{
		if(threadlock) sqlManager.requestStatementExecution(true);
		connection.commit();
		if(threadlock) sqlManager.acknowledgeStatementExecution();
	}
	
	public void flushCache()
	{
		read_cache.clear();
		rs_cache.flush();
	}
	
	@Override
	public void save() throws IOException 
	{
		System.err.println("Saving database updates...");
		if(ghc_cache_dirty && ghc_cache != null)
		{
			saveGeneHitTable();
		}
		if(regidx != null) saveRegionIndex();
		try {
			commitUpdates();
		} catch (SQLException e) {
			e.printStackTrace();
			throw new IOException();
		}
	}

	public void clearVariantTable()
	{
		try
		{
			StatementPrepper sprepper = sqlManager.getStatementGenerator();
			System.err.println("Deleting variants...");
			PreparedStatement ps = sprepper.getVarTableWipeStatement();
			ps.executeUpdate();
			ps.close();
			
			System.err.println("Deleting sample genotype mapping data...");
			ps = sprepper.getSampleGenoTableWipeStatement();
			ps.executeUpdate();
			ps.close();
			
			System.err.println("Deleting gene hit data...");
			ps = sprepper.getGeneHitTableWipeStatement();
			ps.executeUpdate();
			ps.close();
			
			this.zeroGeneHitTable();
			System.err.println("Variant table reset complete!");
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
		
		
	}
	
	public void close() throws SQLException
	{
		//connection.close();
	}
	
}
