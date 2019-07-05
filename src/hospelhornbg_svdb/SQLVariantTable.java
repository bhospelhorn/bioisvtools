package hospelhornbg_svdb;

import java.io.IOException;
import java.io.InputStream;
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
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;

import hospelhornbg_bioinformatics.Genotype;
import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.VCFReadStreamer;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.Gene;
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
	
	public static final String TABLENAME_VARIANTS = "Variants";
	public static final String TABLENAME_SAMPLEGENO = "SampleGeno";
	
	public static final String FIELDNAME_VARUID = "VarUID";
	public static final String FIELDNAME_CTG1 = "Contig1";
	public static final String FIELDNAME_START1 = "Start1";
	public static final String FIELDNAME_START2 = "Start2";
	public static final String FIELDNAME_END1 = "End1";
	public static final String FIELDNAME_END2 = "End2";
	public static final String FIELDNAME_SVTYPE = "SVType";
	public static final String FIELDNAME_POSEFF = "PosEff";
	public static final String FIELDNAME_VARNAME = "VarName";
	
	public static final String FIELDNAME_ACOUNT_TOT = "AlleleCount_Total";
	public static final String FIELDNAME_HCOUNT_TOT = "HomozygoteCount_Total";
	public static final String FIELDNAME_ACOUNT_NFE = "AlleleCount_NFE";
	public static final String FIELDNAME_HCOUNT_NFE = "HomozygoteCount_NFE";
	public static final String FIELDNAME_ACOUNT_AFR = "AlleleCount_AFR";
	public static final String FIELDNAME_HCOUNT_AFR = "HomozygoteCount_AFR";
	public static final String FIELDNAME_ACOUNT_AMR = "AlleleCount_AMR";
	public static final String FIELDNAME_HCOUNT_AMR = "HomozygoteCount_AMR";
	public static final String FIELDNAME_ACOUNT_FIN = "AlleleCount_FIN";
	public static final String FIELDNAME_HCOUNT_FIN = "HomozygoteCount_FIN";
	public static final String FIELDNAME_ACOUNT_EAS = "AlleleCount_EAS";
	public static final String FIELDNAME_HCOUNT_EAS = "HomozygoteCount_EAS";
	public static final String FIELDNAME_ACOUNT_SAS = "AlleleCount_SAS";
	public static final String FIELDNAME_HCOUNT_SAS = "HomozygoteCount_SAS";
	public static final String FIELDNAME_ACOUNT_ASJ = "AlleleCount_ASJ";
	public static final String FIELDNAME_HCOUNT_ASJ = "HomozygoteCount_ASJ";
	public static final String FIELDNAME_ACOUNT_OTH = "AlleleCount_OTH";
	public static final String FIELDNAME_HCOUNT_OTH = "HomozygoteCount_OTH";
	
	public static final String FIELDNAME_GENELIST = "Genes";
	public static final String FIELDNAME_VALNOTES = "ValidationNotes";
	public static final String FIELDNAME_CTG2 = "Contig2";
	public static final String FIELDNAME_INSSEQ = "InsertionSequence";
	public static final String FIELDNAME_GENOTYPES = "Genotypes";
	
	public static final String FIELDNAME_SAMPLEUID = "SampleUID";
	public static final String FIELDNAME_SVARLIST_HOM = "HomVar";
	public static final String FIELDNAME_SVARLIST_HET = "HetVar";
	public static final String FIELDNAME_SVARLIST_OTH = "OthVar";
	
	
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
	
	/* ----- Instance Variables ----- */
	
	private String dbURL;
	private String username;
	private String password;
	
	private Connection connection;
	private StatementPrepper sprepper;

	private GenomeBuild genome;
	private GeneSet genes;
	private GenomeIndex uidIndex;
	
	//private List<String> tempBlobFiles;
	
	/* ----- Construction ----- */
	
	public SQLVariantTable(String url, String user, String pw, GenomeBuild gb, GeneSet gs) throws SQLException
	{
		genome = gb;
		genes = gs;
		uidIndex = new GenomeIndex(genome);
		//Attempt to connect
		dbURL = url;
		username = user;
		password = pw;
		System.err.println("Now connecting...");
		connect();
		System.err.println("Connection successful!");
		sprepper = new StatementPrepper(connection);
		//Check for tables, create if not there
		if(!varTableExists()) createVarTable();
		if(!sampleGenoTableExists()) createSampleGenoTable();
		//tempBlobFiles = new LinkedList<String>();
	}
	
	private void connect() throws SQLException
	{
		connection = DriverManager.getConnection(dbURL, username, password);
		//cstat = connection.createStatement();
	}
	
	/* ----- SQL Table Management ----- */
	
	private boolean varTableExists() throws SQLException
	{
		DatabaseMetaData meta = connection.getMetaData();
		String[] ttypes = {"TABLE"};
		ResultSet rs = meta.getTables(null, null, TABLENAME_VARIANTS.toUpperCase(), ttypes);
		return rs.next();
	}
	
	private void createVarTable() throws SQLException
	{
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
	}
	
	private boolean sampleGenoTableExists() throws SQLException
	{
		DatabaseMetaData meta = connection.getMetaData();
		String[] ttypes = {"TABLE"};
		ResultSet rs = meta.getTables(null, null, TABLENAME_SAMPLEGENO.toUpperCase(), ttypes);
		return rs.next();
	}
	
	private void createSampleGenoTable() throws SQLException
	{
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
	
	private PreparedStatement generateFullVarInsertStatement(DBVariant var, VariantGenotype vgeno) throws IOException, SQLException
	{
		PreparedStatement pstat = sprepper.getFullInsertStatement();
		
		Contig ctg = var.getChrom();
		
		pstat.setLong(StatementPrepper.FULLINS_VARUID, var.getLongID());
		pstat.setInt(StatementPrepper.FULLINS_CHR1, ctg.getUDPName().hashCode());
		pstat.setInt(StatementPrepper.FULLINS_ST1, var.getStartPosition().getStart());
		pstat.setInt(StatementPrepper.FULLINS_ST2, var.getStartPosition().getEnd());
		pstat.setInt(StatementPrepper.FULLINS_ED1, var.getEndPosition().getStart());
		pstat.setInt(StatementPrepper.FULLINS_ED2, var.getEndPosition().getEnd());
		pstat.setInt(StatementPrepper.FULLINS_SVTYPE, var.getType().getID());
		pstat.setInt(StatementPrepper.FULLINS_POSEFF, var.getPositionEffect().getPriority());
		pstat.setString(StatementPrepper.FULLINS_VARNAME, var.getName());
		
		pstat.setInt(StatementPrepper.FULLINS_ACOUNT_TOT, var.getIndividualCount());
		pstat.setInt(StatementPrepper.FULLINS_HCOUNT_TOT, var.getHomozygoteCount());
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
		Blob blob = sprepper.toBlob(var.getGeneListAsBLOBBytes());
		pstat.setBlob(StatementPrepper.FULLINS_GENELIST, blob);
		
		String valnotes = var.getValidationNotes();
		if(valnotes == null) valnotes = "N/A";
		pstat.setString(StatementPrepper.FULLINS_VALNOTES, valnotes);
		
		Contig ctg2 = var.getEndChrom();
		if(ctg2 == null) pstat.setInt(StatementPrepper.FULLINS_CTG2, ctg.getUDPName().hashCode());
		else pstat.setInt(StatementPrepper.FULLINS_CTG2, ctg2.getUDPName().hashCode());
		
		//String insseq = var.getAltAlleleString();
		//if(insseq == null) insseq = "N/A";
		//pstat.setString(StatementPrepper.FULLINS_INSSEQ, insseq);
		blob = sprepper.toBlob(var.getInsseqAsBLOBBytes());
		pstat.setBlob(StatementPrepper.FULLINS_INSSEQ, blob);
		
		//Genotype blob
		blob = sprepper.toBlob(vgeno.getGenotypesAsBLOBBytes());
		pstat.setBlob(StatementPrepper.FULLINS_GENOTYPES, blob);
		
		return pstat;
	}
	
	private PreparedStatement generateAbridgedSetVarUpdateStatement(DBVariant var, VariantGenotype vgeno) throws IOException, SQLException
	{
		//Updates:
		//	Start, End
		//	Allele Counts
		//	Genes
		//  PosEff
		//	Validation notes
		//	Genotype
		
		PreparedStatement pstat = sprepper.getShortUpdateStatement();
		
		pstat.setInt(StatementPrepper.SHORTUD_ST1, var.getStartPosition().getStart()); System.err.println("Start1 = " + var.getStartPosition().getStart());
		pstat.setInt(StatementPrepper.SHORTUD_ST2, var.getStartPosition().getEnd()); System.err.println("Start2 = " + var.getStartPosition().getEnd());
		pstat.setInt(StatementPrepper.SHORTUD_ED1, var.getEndPosition().getStart()); System.err.println("End1 = " + var.getEndPosition().getStart());
		pstat.setInt(StatementPrepper.SHORTUD_ED2, var.getEndPosition().getEnd()); System.err.println("End2 = " + var.getEndPosition().getEnd());
		pstat.setInt(StatementPrepper.SHORTUD_POSEFF, var.getPositionEffect().getPriority()); System.err.println("PosEff = " + var.getPositionEffect().getPriority());
		
		pstat.setInt(StatementPrepper.SHORTUD_ACOUNT_TOT, var.getIndividualCount());
		pstat.setInt(StatementPrepper.SHORTUD_HCOUNT_TOT, var.getHomozygoteCount());
		pstat.setInt(StatementPrepper.SHORTUD_ACOUNT_NFE, var.getIndividualCount(Population.NFE));
		pstat.setInt(StatementPrepper.SHORTUD_HCOUNT_NFE, var.getHomozygoteCount(Population.NFE));
		pstat.setInt(StatementPrepper.SHORTUD_ACOUNT_AFR, var.getIndividualCount(Population.AFR));
		pstat.setInt(StatementPrepper.SHORTUD_HCOUNT_AFR, var.getHomozygoteCount(Population.AFR));
		pstat.setInt(StatementPrepper.SHORTUD_ACOUNT_AMR, var.getIndividualCount(Population.AMR));
		pstat.setInt(StatementPrepper.SHORTUD_HCOUNT_AMR, var.getHomozygoteCount(Population.AMR));
		pstat.setInt(StatementPrepper.SHORTUD_ACOUNT_FIN, var.getIndividualCount(Population.FIN));
		pstat.setInt(StatementPrepper.SHORTUD_HCOUNT_FIN, var.getHomozygoteCount(Population.FIN));
		pstat.setInt(StatementPrepper.SHORTUD_ACOUNT_EAS, var.getIndividualCount(Population.EAS));
		pstat.setInt(StatementPrepper.SHORTUD_HCOUNT_EAS, var.getHomozygoteCount(Population.EAS));
		pstat.setInt(StatementPrepper.SHORTUD_ACOUNT_SAS, var.getIndividualCount(Population.SAS));
		pstat.setInt(StatementPrepper.SHORTUD_HCOUNT_SAS, var.getHomozygoteCount(Population.SAS));
		pstat.setInt(StatementPrepper.SHORTUD_ACOUNT_ASJ, var.getIndividualCount(Population.ASJ));
		pstat.setInt(StatementPrepper.SHORTUD_HCOUNT_ASJ, var.getHomozygoteCount(Population.ASJ));
		pstat.setInt(StatementPrepper.SHORTUD_ACOUNT_OTH, var.getIndividualCount(Population.OTH));
		pstat.setInt(StatementPrepper.SHORTUD_HCOUNT_OTH, var.getHomozygoteCount(Population.OTH));
		
		Blob blob = sprepper.toBlob(var.getGeneListAsBLOBBytes());
		pstat.setBlob(StatementPrepper.SHORTUD_GENELIST, blob);
		
		String valnotes = var.getValidationNotes();
		if(valnotes == null) valnotes = "N/A";
		System.err.println("ValNotes = " + valnotes);
		pstat.setString(StatementPrepper.SHORTUD_VALNOTES, valnotes);
		
		blob = sprepper.toBlob(vgeno.getGenotypesAsBLOBBytes());
		pstat.setBlob(StatementPrepper.SHORTUD_GENOTYPES, blob);
		
		pstat.setLong(StatementPrepper.SHORTUD_QUERYID, var.getLongID()); System.err.println("VarUID = " + Long.toHexString(var.getLongID()));
		
		return pstat;
	}
	
	private PreparedStatement generatePopulationSetVarUpdateStatement(DBVariant dbv) throws SQLException
	{
		PreparedStatement pstat = sprepper.getPopUpdateStatement();
		
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
	
	/* ----- Getters ----- */
	
	public boolean variantExists(long varUID)
	{
		try 
		{
			PreparedStatement pstat = sprepper.getVarUIDCheckStatement();
			pstat.setLong(StatementPrepper.VARUIDCHECK_VARUID, varUID);
			ResultSet rs = pstat.executeQuery();
			return rs.next();
		} 
		catch (SQLException e) 
		{
			e.printStackTrace();
			return false;
		}
	}
	
	public DBVariant getVariant(long varUID) 
	{	
		try 
		{
			PreparedStatement pstat = sprepper.getVarGetterStatement();
			pstat.setLong(StatementPrepper.VARGET_VARUID, varUID);
			ResultSet rs = pstat.executeQuery();
			if(!rs.next()) return null;
			DBVariant var = readFromResultSet(rs);
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
		if(varUIDs == null || varUIDs.isEmpty()) return null;
		List<DBVariant> vlist = new LinkedList<DBVariant>();
		int len = varUIDs.size();
		
		try 
		{
			PreparedStatement pstat = sprepper.generateMultiVarGetterStatement(len);
			int i = 1;
			for(Long vid : varUIDs)
			{
				pstat.setLong(i, vid);
				i++;
			}
			ResultSet rs = pstat.executeQuery();
			while(rs.next())
			{
				vlist.add(readFromResultSet(rs));
			}
		} 
		catch (Exception e) 
		{
			e.printStackTrace();
			return null;
		}
		
		return vlist;
	}
	
	@Override
	public VariantGenotype getGenotype(long varUID) 
	{
		try 
		{
			PreparedStatement pstat = sprepper.getGenoGetterStatement();
			pstat.setLong(StatementPrepper.GENOGET_VARUID, varUID);
			ResultSet rs = pstat.executeQuery();
			if(!rs.next()) return null;
			Blob genoblob = rs.getBlob(FIELDNAME_GENOTYPES);
			VariantGenotype vg = new VariantGenotype(varUID);
			vg.readDataFromBLOB(genoblob);
			
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
		if(varUIDs == null || varUIDs.isEmpty()) return null;
		List<VariantGenotype> glist = new LinkedList<VariantGenotype>();
		
		int len = varUIDs.size();
		
		try 
		{
			PreparedStatement pstat = sprepper.generateMultiGenoGetterStatement(len);
			int i = 1;
			for(Long vid : varUIDs)
			{
				pstat.setLong(i, vid);
				i++;
			}
			ResultSet rs = pstat.executeQuery();
			while(rs.next())
			{
				long varUID = rs.getLong(FIELDNAME_VARUID);
				Blob genoblob = rs.getBlob(FIELDNAME_GENOTYPES);
				VariantGenotype vg = new VariantGenotype(varUID);
				vg.readDataFromBLOB(genoblob);
				glist.add(vg);
			}
		} 
		catch (Exception e) 
		{
			e.printStackTrace();
			return null;
		}
		
		return glist;
	}
	
	public List<Long> getVariantIDsForSample(int sampleUID) 
	{
		try 
		{
			PreparedStatement pstat = sprepper.getSampleVarGetterStatement();
			pstat.setInt(StatementPrepper.SVARGET_SAMPUID, sampleUID);
			ResultSet rs = pstat.executeQuery();
			if(!rs.next()) return null;
			
			Set<Long> vars = new HashSet<Long>();
			vars.addAll(readVarUIDListBlob(rs.getBlob(FIELDNAME_SVARLIST_HET)));
			vars.addAll(readVarUIDListBlob(rs.getBlob(FIELDNAME_SVARLIST_HOM)));
			vars.addAll(readVarUIDListBlob(rs.getBlob(FIELDNAME_SVARLIST_OTH)));
			
			List<Long> list = new ArrayList<Long>(vars.size() + 1);
			list.addAll(vars);
			Collections.sort(list);
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
		List<DBVariant> varlist = new LinkedList<DBVariant>();
		if(c == null) return varlist;
		
		//Don't forget to check the other end of TRAs!
		int cuid = c.getUDPName().hashCode();
		/*String sqlQuery = "SELECT * FROM " + TABLENAME_VARIANTS;
		sqlQuery += " WHERE ";
		
		String statement_type = FIELDNAME_SVTYPE + "=" + SVType.TRA.getID() + " OR " + FIELDNAME_SVTYPE + "=" + SVType.BND.getID();
		String statement_ctg1 = FIELDNAME_CTG1 + "=" + cuid;
		String statement_ctg2 = FIELDNAME_CTG2 + "=" + cuid;
		String statement_reg1 = FIELDNAME_START1 + "<" + end;
		String statement_reg2 = FIELDNAME_END2 + ">=" + start;
		
		String statement_trareg1 = FIELDNAME_START1 + "<" + end;
		String statement_trareg2 = FIELDNAME_START2 + ">=" + start;
		String statement_trareg3 = FIELDNAME_END1 + "<" + end;
		String statement_trareg4 = FIELDNAME_END2 + ">=" + start;
		
		//Combine for full statement
		*
		 *  (!TRA) && ctg1 && reg1 && reg2
		 *  or
		 *  TRA && ((ctg1 && treg1 && treg2) || (ctg2 && treg3 && treg4))
		 
		
		String nottra = "(NOT " + statement_type + ")";
		String normQuery = nottra + " AND " + statement_ctg1 + " AND " + statement_reg1 + " AND " + statement_reg2;
		String tQuery1 = statement_ctg1 + " AND " + statement_trareg1 + " AND " + statement_trareg2;
		String tQuery2 = statement_ctg2 + " AND " + statement_trareg3 + " AND " + statement_trareg4;
		String traQuery = statement_type + " AND ((" + tQuery1 + ") OR (" + tQuery2 + "))";
		sqlQuery += "(" + normQuery + ") OR (" + traQuery + ")";*/
		//System.err.println("-DEBUG- sqlQuery: " + sqlQuery);
				
		try 
		{
			PreparedStatement pstat = sprepper.getRegionVarGetterStatement();
			for(int i : StatementPrepper.VARS_REG_CUID) pstat.setInt(i, cuid);
			for(int i : StatementPrepper.VARS_REG_START) pstat.setInt(i, start);
			for(int i : StatementPrepper.VARS_REG_END) pstat.setInt(i, end);
			ResultSet rs = pstat.executeQuery();
			while(rs.next())
			{
				DBVariant var = readFromResultSet(rs);
				if(var != null) varlist.add(var);
			}
		} 
		catch (Exception e) 
		{
			e.printStackTrace();
		}
		
		return varlist;
	}

	/* ----- Analysis ----- */
	
	public Map<String, GeneHitCounter> generateGeneHitMap() 
	{
		List<Gene> glist = genes.getAllGenes();
		Map<String, GeneHitCounter> ghmap = new HashMap<String, GeneHitCounter>();
		for(Gene g : glist)
		{
			GeneHitCounter ghc = new GeneHitCounter();
			ghmap.put(g.getID(), ghc);
			
			Collection<DBVariant> hits = getVariantsInRegion(g.getChromosome(), g.getTranscriptStart(), g.getTranscriptEnd());
			if(hits == null) continue;
			ghc.total_hits_var = hits.size();
			
			List<DBVariant> exonhits = new LinkedList<DBVariant>();
			for(DBVariant var : hits)
			{
				if(g.getRelativeRegionLocationEffect(var.getStartPosition().getStart(), var.getEndPosition().getEnd()) == GeneFunc.EXONIC)
				{
					exonhits.add(var);
				}
			}
			ghc.exon_hits_var = exonhits.size();
			
			for(DBVariant var : hits)
			{
				VariantGenotype vg = getGenotype(var.getLongID());
				if(vg != null) ghc.total_hits_indiv.addAll(vg.getAllIndividuals());
			}
			
			for(DBVariant var : exonhits)
			{
				VariantGenotype vg = getGenotype(var.getLongID());
				if(vg != null) ghc.exon_hits_indiv.addAll(vg.getAllIndividuals());
			}
			
		}
		
		return ghmap;
	}
	
	/* ----- Variant Addition ----- */
	
	private long generateUID(Contig c, int pos)
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
	
	private boolean addOrUpdateVariant(DBVariant var, VariantGenotype vgeno, boolean isnew)
	{
		try
		{
			PreparedStatement baseStatement = null;
			if(isnew) baseStatement = generateFullVarInsertStatement(var, vgeno);
			else {
				var.noteGenes(genes);
				baseStatement = generateAbridgedSetVarUpdateStatement(var, vgeno);
			}
			
			//Statement cstat = connection.createStatement();
			//int count = cstat.executeUpdate(baseStatement);
			int count = baseStatement.executeUpdate();
			if(count != 1) return false;
		}
		catch(Exception e)
		{
			e.printStackTrace();
			return false;
		}
		
		return true;
	}
	
	@Override
	public boolean addVCF(String vcfpath, Map<String, FamilyMember> sampleMap, int mergeFactor) throws IOException 
	{
		if(vcfpath == null) return false;
		if(sampleMap == null) return false;
		if (mergeFactor < 0) return false;
		
		Map<Integer, List<Long>> sgAdd = new TreeMap<Integer, List<Long>>();
		Map<Integer, List<Long>> sgDel = new TreeMap<Integer, List<Long>>();
		
		List<String> samples = new ArrayList<String>(sampleMap.size() + 1);
		samples.addAll(sampleMap.keySet());
		
		VCFReadStreamer vcfReader = new VCFReadStreamer(vcfpath, genome);
		vcfReader.open();
		
		int debugctr = 0;
		
		Iterator<StructuralVariant> sviter = vcfReader.getSVIterator();
		while(sviter.hasNext())
		{
			StructuralVariant sv = sviter.next();
			//Look for merge candidates...
			long id = -1L;
			//Calculate merge distance in bp...
			double mfac = (double)mergeFactor / 1000.0;
			int svsz = sv.getCIPosition(true, false, true) - sv.getCIPosition(false, false, false);
			int bp = (int)Math.round(mfac * (double)svsz);
			
			//Grab variants in region to look up
			int search_start = sv.getCIPosition(false, false, false) - bp;
			int search_end = sv.getCIPosition(true, false, true) + bp;
			Collection<DBVariant> nearvars = this.getVariantsInRegion(sv.getChromosome(), search_start, search_end);
			if(nearvars != null && !nearvars.isEmpty())
			{
				//Look for matches
				for(DBVariant dbv : nearvars)
				{
					if(dbv.svIsEquivalent(sv, mfac))
					{
						id = dbv.getLongID();
						break;
					}
				}
			}
			
			//If matched, need to update record.
			//If not matched, need new record.
			boolean newvar = (id == -1L);
			if(!newvar)
			{
				DBVariant dbv = this.getVariant(id);
				//Update start/end
				int npos = sv.getCIPosition(false, false, false);
				int opos = dbv.getStartPosition().getStart();
				if(npos < opos) dbv.getStartPosition().setStart(npos);
				npos = sv.getCIPosition(false, false, true);
				opos = dbv.getStartPosition().getEnd();
				if(npos > opos) dbv.getStartPosition().setEnd(npos);
				
				npos = sv.getCIPosition(true, false, false);
				opos = dbv.getEndPosition().getStart();
				if(npos < opos) dbv.getEndPosition().setStart(npos);
				npos = sv.getCIPosition(true, false, true);
				opos = dbv.getEndPosition().getEnd();
				if(npos > opos) dbv.getEndPosition().setEnd(npos);
				
				//Update genotype & allele counts
				VariantGenotype vg = this.getGenotype(id);
				for(String sample : samples)
				{
					FamilyMember mem = sampleMap.get(sample);
					if(mem == null) continue;
					//See if already has genotype...
					SVDBGenotype g = vg.getGenotype(mem.getUID());
					Genotype geno = sv.getSampleGenotype(sample);
					if(geno != null)
					{
						//See if homref.
						//If so remove if there or continue
						if(geno.isGenotypeUnknown() || (geno.isHomozygous() && geno.hasAllele(0)))
						{
							if(g != null)
							{
								vg.removeGenotype(mem.getUID());
								boolean homalt = g.isHomozygous();
								Collection<Population> pflags = mem.getPopulationTags();
								for(Population p : pflags)
								{
									dbv.decrementTotalCount(p);
									if(homalt)dbv.decrementHomozygoteCount(p);
								}
								dbv.decrementTotalCount();
								if(homalt)dbv.decrementHomozygoteCount();
								List<Long> dlist = sgDel.get(mem.getUID());
								if(dlist == null) {dlist = new LinkedList<Long>(); sgDel.put(mem.getUID(), dlist);}
								dlist.add(id);
							}
						}
						else
						{
							//Add/update
							SVDBGenotype newgeno = SVDBGenotype.generateGenotype(mem.getUID(), geno, sv);
							vg.addGenotype(newgeno);
							if(g != null)
							{
								boolean homalt = newgeno.isHomozygous();
								boolean washom = g.isHomozygous();
								if(homalt && !washom)
								{
									//Increment hom count
									Collection<Population> pflags = mem.getPopulationTags();
									dbv.incrementHomozygoteCount();
									for(Population p : pflags) dbv.incrementHomozygoteCount(p);
									List<Long> alist = sgAdd.get(mem.getUID());
									if(alist == null) {alist = new LinkedList<Long>(); sgAdd.put(mem.getUID(), alist);}
									alist.add(id);
								}
								else if (!homalt && washom)
								{
									//Decrement hom count
									Collection<Population> pflags = mem.getPopulationTags();
									dbv.decrementHomozygoteCount();
									for(Population p : pflags) dbv.decrementHomozygoteCount(p);
									List<Long> dlist = sgDel.get(mem.getUID());
									if(dlist == null) {dlist = new LinkedList<Long>(); sgDel.put(mem.getUID(), dlist);}
									dlist.add(id);
								}
							}
							else
							{
								boolean homalt = newgeno.isHomozygous();
								Collection<Population> pflags = mem.getPopulationTags();
								for(Population p : pflags)
								{
									dbv.incrementTotalCount(p);
									if(homalt)dbv.incrementHomozygoteCount(p);
								}
								dbv.incrementTotalCount();
								if(homalt)dbv.incrementHomozygoteCount();
								List<Long> alist = sgAdd.get(mem.getUID());
								if(alist == null) {alist = new LinkedList<Long>(); sgAdd.put(mem.getUID(), alist);}
								alist.add(id);
							}
						}
					}
					else
					{
						//Remove if there.
						//Otherwise do nothing
						if(g != null)
						{
							vg.removeGenotype(mem.getUID());
							boolean homalt = g.isHomozygous();
							Collection<Population> pflags = mem.getPopulationTags();
							for(Population p : pflags)
							{
								dbv.decrementTotalCount(p);
								if(homalt)dbv.decrementHomozygoteCount(p);
							}
							dbv.decrementTotalCount();
							if(homalt)dbv.decrementHomozygoteCount();
							List<Long> dlist = sgDel.get(mem.getUID());
							if(dlist == null) {dlist = new LinkedList<Long>(); sgDel.put(mem.getUID(), dlist);}
							dlist.add(id);
						}
					}
				}
				this.addOrUpdateVariant(dbv, vg, false);
			}
			else
			{
				id = generateUID(sv.getChromosome(), sv.getPosition());
				DBVariant dbv = DBVariant.getFromVariant(sv, sv.getVarID());
				dbv.setLongUID(id);
				dbv.noteGenes(genes);
				VariantGenotype vg = new VariantGenotype(id);
				
				for(String sample : samples)
				{
					FamilyMember mem = sampleMap.get(sample);
					if(mem == null) continue;
					Genotype gt = sv.getSampleGenotype(sample);
					if(gt == null) continue;
					if(gt.isGenotypeUnknown()) continue;
					if(gt.isHomozygous() && gt.hasAllele(0)) continue;
					SVDBGenotype g = SVDBGenotype.generateGenotype(mem.getUID(), gt, sv);
					vg.addGenotype(g);
					
					boolean homalt = g.isHomozygous();
					Collection<Population> pflags = mem.getPopulationTags();
					for(Population p : pflags)
					{
						dbv.incrementTotalCount(p);
						if(homalt)dbv.incrementHomozygoteCount(p);
					}
					dbv.incrementTotalCount();
					if(homalt)dbv.incrementHomozygoteCount();
					
					List<Long> alist = sgAdd.get(mem.getUID());
					if(alist == null) {alist = new LinkedList<Long>(); sgAdd.put(mem.getUID(), alist);}
					alist.add(id);
				}
				
				addOrUpdateVariant(dbv, vg, true);
			}
			//System.exit(2);
			debugctr++;
			if(debugctr >= 1) System.exit(2);
		}
		
		vcfReader.close();
		return this.updateSampleGenoTable(sgAdd, sgDel);
	}
	
	/* ----- Variant Deletion ----- */
	
	private boolean removeSampleGenoRecord(int sampleID)
	{
		try 
		{
			//Statement cstat = connection.createStatement();
			PreparedStatement pstat = sprepper.getSampleGenoDeleteStatment();
			pstat.setInt(1, sampleID);
			int count = pstat.executeUpdate();
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
		//Get genotype
		VariantGenotype vg = getGenotype(varUID);
		
		try 
		{
			//Statement cstat = connection.createStatement();
			//int count = cstat.executeUpdate(sqlQuery);
			PreparedStatement pstat = sprepper.getVarDeleteStatment();
			pstat.setLong(1, varUID);
			int count = pstat.executeUpdate();
			
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
		
		int len = varUIDs.size();
		
		try 
		{
			PreparedStatement pstat = sprepper.generateMultiVarDeleteStatement(len);
			int i = 1;
			for(Long vid : varUIDs) 
			{
				pstat.setLong(i, vid);
				i++;
			}
			int count = pstat.executeUpdate();

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
						PreparedStatement pstat = sprepper.getSGenoInsertStatement();
						pstat.setInt(StatementPrepper.SGENOINS_SID, sid);
						
						Blob b = sprepper.toBlob(homblob.getBytes());
						pstat.setBlob(StatementPrepper.SGENOINS_HOM, b);
						
						b = sprepper.toBlob(hetblob.getBytes());
						pstat.setBlob(StatementPrepper.SGENOINS_HET, b);
						
						b = sprepper.toBlob(othblob.getBytes());
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
					PreparedStatement pstat = sprepper.getSGenoUpdateStatement();
					pstat.setInt(StatementPrepper.SGENOUD_SID, sid);
					
					Blob b = sprepper.toBlob(homblob.getBytes());
					pstat.setBlob(StatementPrepper.SGENOUD_HOM, b);
					
					b = sprepper.toBlob(hetblob.getBytes());
					pstat.setBlob(StatementPrepper.SGENOUD_HET, b);
					
					b = sprepper.toBlob(othblob.getBytes());
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
			PreparedStatement pstat = generatePopulationSetVarUpdateStatement(dbv);
			int count = pstat.executeUpdate();
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
				PreparedStatement pstat = generatePopulationSetVarUpdateStatement(dbv);
				int count = pstat.executeUpdate();
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
	
	/* ----- Cleanup ----- */
	
	@Override
	public void save() throws IOException 
	{
		//Shouldn't need to do anything???
	}

}
