package hospelhornbg_svdb;

import java.sql.Blob;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.SQLException;

import hospelhornbg_bioinformatics.SVType;

public class StatementPrepper {
	
	public static final int VARUIDCHECK_VARUID = 1;
	
	public static final int VARGET_VARUID = 1;
	public static final int GENOGET_VARUID = 1;
	public static final int SVARGET_SAMPUID = 1;
	
	public static final int[] VARS_REG_CUID = {1, 4, 7};
	public static final int[] VARS_REG_START = {3, 6, 9};
	public static final int[] VARS_REG_END = {2, 5, 8};
	
	public static final int VARS_REG_NOTRA_CUID = 1;
	public static final int VARS_REG_NOTRA_START = 2;
	public static final int VARS_REG_NOTRA_END = 3;
	
	public static final int VARS_REG_TYPE_TYPE = 1;
	public static final int VARS_REG_TYPE_CUID = 2;
	public static final int VARS_REG_TYPE_START = 3;
	public static final int VARS_REG_TYPE_END = 4;
	
	public static final int FULLINS_VARUID = 1;
	public static final int FULLINS_CHR1 = 2;
	public static final int FULLINS_ST1 = 3;
	public static final int FULLINS_ST2 = 4;
	public static final int FULLINS_ED1 = 5;
	public static final int FULLINS_ED2 = 6;
	public static final int FULLINS_SVTYPE = 7;
	public static final int FULLINS_POSEFF = 8;
	public static final int FULLINS_VARNAME = 9;
	public static final int FULLINS_ACOUNT_TOT = 10;
	public static final int FULLINS_HCOUNT_TOT = 11;
	public static final int FULLINS_ACOUNT_NFE = 12;
	public static final int FULLINS_HCOUNT_NFE = 13;
	public static final int FULLINS_ACOUNT_AFR = 14;
	public static final int FULLINS_HCOUNT_AFR = 15;
	public static final int FULLINS_ACOUNT_AMR = 16;
	public static final int FULLINS_HCOUNT_AMR = 17;
	public static final int FULLINS_ACOUNT_FIN = 18;
	public static final int FULLINS_HCOUNT_FIN = 19;
	public static final int FULLINS_ACOUNT_EAS = 20;
	public static final int FULLINS_HCOUNT_EAS = 21;
	public static final int FULLINS_ACOUNT_SAS = 22;
	public static final int FULLINS_HCOUNT_SAS = 23;
	public static final int FULLINS_ACOUNT_ASJ = 24;
	public static final int FULLINS_HCOUNT_ASJ = 25;
	public static final int FULLINS_ACOUNT_OTH = 26;
	public static final int FULLINS_HCOUNT_OTH = 27;
	public static final int FULLINS_GENELIST = 28;
	public static final int FULLINS_VALNOTES = 29;
	public static final int FULLINS_CTG2 = 30;
	public static final int FULLINS_INSSEQ = 31;
	public static final int FULLINS_GENOTYPES = 32;
	
	public static final int SHORTUD_ST1 = 1;
	public static final int SHORTUD_ST2 = 2;
	public static final int SHORTUD_ED1 = 3;
	public static final int SHORTUD_ED2 = 4;
	public static final int SHORTUD_POSEFF = 5;
	public static final int SHORTUD_ACOUNT_TOT = 6;
	public static final int SHORTUD_HCOUNT_TOT = 7;
	public static final int SHORTUD_ACOUNT_NFE = 8;
	public static final int SHORTUD_HCOUNT_NFE = 9;
	public static final int SHORTUD_ACOUNT_AFR = 10;
	public static final int SHORTUD_HCOUNT_AFR = 11;
	public static final int SHORTUD_ACOUNT_AMR = 12;
	public static final int SHORTUD_HCOUNT_AMR = 13;
	public static final int SHORTUD_ACOUNT_FIN = 14;
	public static final int SHORTUD_HCOUNT_FIN = 15;
	public static final int SHORTUD_ACOUNT_EAS = 16;
	public static final int SHORTUD_HCOUNT_EAS = 17;
	public static final int SHORTUD_ACOUNT_SAS = 18;
	public static final int SHORTUD_HCOUNT_SAS = 19;
	public static final int SHORTUD_ACOUNT_ASJ = 20;
	public static final int SHORTUD_HCOUNT_ASJ = 21;
	public static final int SHORTUD_ACOUNT_OTH = 22;
	public static final int SHORTUD_HCOUNT_OTH = 23;
	public static final int SHORTUD_GENELIST = 24;
	public static final int SHORTUD_VALNOTES = 25;
	public static final int SHORTUD_GENOTYPES = 26;
	public static final int SHORTUD_QUERYID = 27;
	
	public static final int POPUD_ACOUNT_TOT = 1;
	public static final int POPUD_HCOUNT_TOT = 2;
	public static final int POPUD_ACOUNT_NFE = 3;
	public static final int POPUD_HCOUNT_NFE = 4;
	public static final int POPUD_ACOUNT_AFR = 5;
	public static final int POPUD_HCOUNT_AFR = 6;
	public static final int POPUD_ACOUNT_AMR = 7;
	public static final int POPUD_HCOUNT_AMR = 8;
	public static final int POPUD_ACOUNT_FIN = 9;
	public static final int POPUD_HCOUNT_FIN = 10;
	public static final int POPUD_ACOUNT_EAS = 11;
	public static final int POPUD_HCOUNT_EAS = 12;
	public static final int POPUD_ACOUNT_SAS = 13;
	public static final int POPUD_HCOUNT_SAS = 14;
	public static final int POPUD_ACOUNT_ASJ = 15;
	public static final int POPUD_HCOUNT_ASJ = 16;
	public static final int POPUD_ACOUNT_OTH = 17;
	public static final int POPUD_HCOUNT_OTH = 18;
	public static final int POPUD_QUERYID = 19;
	
	public static final int SGENOINS_SID = 1;
	public static final int SGENOINS_HOM = 2;
	public static final int SGENOINS_HET = 3;
	public static final int SGENOINS_OTH = 4;
	
	public static final int SGENOUD_SID = 4;
	public static final int SGENOUD_HOM = 1;
	public static final int SGENOUD_HET = 2;
	public static final int SGENOUD_OTH = 3;
	
	public static final int GENEHIT_UID = 1;
	public static final int GENEHIT_TOT = 2;
	public static final int GENEHIT_EXON = 3;
	public static final int GENEHIT_TOT_INDIV = 4;
	public static final int GENEHIT_EXON_INDIV = 5;
	
	public static final int GENEHIT_UD_UID = 5;
	public static final int GENEHIT_UD_TOT = 1;
	public static final int GENEHIT_UD_EXON = 2;
	public static final int GENEHIT_UD_TOT_INDIV = 3;
	public static final int GENEHIT_UD_EXON_INDIV = 4;
	
	private Connection connection;
	
	private PreparedStatement varuid_check;
	private PreparedStatement var_getter;
	private PreparedStatement geno_getter;
	private PreparedStatement sampvar_getter;
	
	private PreparedStatement vars_in_reg;
	private PreparedStatement vars_in_reg_notra;
	private PreparedStatement vars_in_reg_type;
	
	private PreparedStatement insert_full;
	private PreparedStatement short_update;
	private PreparedStatement pop_update;
	
	private PreparedStatement sgeno_delete;
	private PreparedStatement var_delete;
	
	private PreparedStatement sgeno_insert;
	private PreparedStatement sgeno_update;
	
	private PreparedStatement gh_get_one;
	private PreparedStatement gh_insert;
	private PreparedStatement gh_update;

	public StatementPrepper(Connection c)
	{
		connection = c;
	}
	
	public PreparedStatement getVarUIDCheckStatement() throws SQLException
	{
		if(varuid_check == null)
		{
			String sqlQuery = "SELECT " + SQLVariantTable.FIELDNAME_VARUID + " FROM " + SQLVariantTable.TABLENAME_VARIANTS;
			sqlQuery += " WHERE " + SQLVariantTable.FIELDNAME_VARUID + " = ?";
			//System.err.println(sqlQuery);
			varuid_check = connection.prepareStatement(sqlQuery);
		}
		return varuid_check;
	}
	
	public PreparedStatement getVarGetterStatement() throws SQLException
	{
		if(var_getter == null)
		{
			String sqlQuery = "SELECT * FROM " + SQLVariantTable.TABLENAME_VARIANTS;
			sqlQuery += " WHERE " + SQLVariantTable.FIELDNAME_VARUID + " = ?";
			var_getter = connection.prepareStatement(sqlQuery);
		}
		return var_getter;
	}
	
	public PreparedStatement generateMultiVarGetterStatement(int count) throws SQLException
	{
		String sqlQuery = "SELECT * FROM " + SQLVariantTable.TABLENAME_VARIANTS + " WHERE";
		for(int i = 0; i < count; i++)
		{
			sqlQuery += " " + SQLVariantTable.FIELDNAME_VARUID + " = ?";
			if(i < (count - 1)) sqlQuery += " OR";
		}
		return connection.prepareStatement(sqlQuery);
	}
	
	public PreparedStatement getGenoGetterStatement() throws SQLException
	{
		if(geno_getter == null)
		{
			String sqlQuery = "SELECT " + SQLVariantTable.FIELDNAME_VARUID + ", " + SQLVariantTable.FIELDNAME_GENOTYPES + " FROM " + SQLVariantTable.TABLENAME_VARIANTS;
			sqlQuery += " WHERE " + SQLVariantTable.FIELDNAME_VARUID + " = ?";
			geno_getter = connection.prepareStatement(sqlQuery);
		}
		return geno_getter;
	}
	
	public PreparedStatement generateMultiGenoGetterStatement(int count) throws SQLException
	{
		String sqlQuery = "SELECT " + SQLVariantTable.FIELDNAME_VARUID + ", " + SQLVariantTable.FIELDNAME_GENOTYPES + " FROM " + SQLVariantTable.TABLENAME_VARIANTS + " WHERE";
		for(int i = 0; i < count; i++)
		{
			sqlQuery += " " + SQLVariantTable.FIELDNAME_VARUID + " = ?";
			if(i < (count - 1)) sqlQuery += " OR";
		}
		return connection.prepareStatement(sqlQuery);
	}
	
	public PreparedStatement getSampleVarGetterStatement() throws SQLException
	{
		if(sampvar_getter == null)
		{
			String sqlQuery = "SELECT * FROM " + SQLVariantTable.TABLENAME_SAMPLEGENO;
			sqlQuery += " WHERE " + SQLVariantTable.FIELDNAME_SAMPLEUID + " = ?";
			sampvar_getter = connection.prepareStatement(sqlQuery);
		}
		return sampvar_getter;
	}
	
	public PreparedStatement getRegionVarGetterStatement() throws SQLException
	{
		if(vars_in_reg == null)
		{
			String sqlQuery = "SELECT * FROM " + SQLVariantTable.TABLENAME_VARIANTS;
			sqlQuery += " WHERE ";
			
			String statement_type = SQLVariantTable.FIELDNAME_SVTYPE + "=" + SVType.TRA.getID() + " OR " + SQLVariantTable.FIELDNAME_SVTYPE + "=" + SVType.BND.getID();
			String statement_ctg1 = SQLVariantTable.FIELDNAME_CTG1 + " = ?";// + cuid; 0 3
			String statement_ctg2 = SQLVariantTable.FIELDNAME_CTG2 + " = ?";// + cuid; 6
			String statement_reg1 = SQLVariantTable.FIELDNAME_START1 + " < ?";// + end; 1
			String statement_reg2 = SQLVariantTable.FIELDNAME_END2 + " >= ?";// + start; 2
			
			String statement_trareg1 = SQLVariantTable.FIELDNAME_START1 + " < ?";// + end; 4
			String statement_trareg2 = SQLVariantTable.FIELDNAME_START2 + " >= ?";// + start; 5
			String statement_trareg3 = SQLVariantTable.FIELDNAME_END1 + " < ?";// + end; 7
			String statement_trareg4 = SQLVariantTable.FIELDNAME_END2 + " >= ?";// + start; 8
			
			//Combine for full statement
			/*
			 *  (!TRA) && ctg1 && reg1 && reg2
			 *  or
			 *  TRA && ((ctg1 && treg1 && treg2) || (ctg2 && treg3 && treg4))
			 */
			
			String nottra = "(NOT " + statement_type + ")";
			String normQuery = nottra + " AND " + statement_ctg1 + " AND " + statement_reg1 + " AND " + statement_reg2;
			String tQuery1 = statement_ctg1 + " AND " + statement_trareg1 + " AND " + statement_trareg2;
			String tQuery2 = statement_ctg2 + " AND " + statement_trareg3 + " AND " + statement_trareg4;
			String traQuery = statement_type + " AND ((" + tQuery1 + ") OR (" + tQuery2 + "))";
			sqlQuery += "(" + normQuery + ") OR (" + traQuery + ")";
			
			vars_in_reg = connection.prepareStatement(sqlQuery);
			//System.err.println("StatementPrepper.getRegionVarGetterStatement || Statement is null? " + (vars_in_reg == null));
		}
		return vars_in_reg;
	}
	
	public PreparedStatement getRegionNoTRAVarGetterStatement() throws SQLException
	{
		if(vars_in_reg_notra == null)
		{
			String sqlQuery = "SELECT * FROM " + SQLVariantTable.TABLENAME_VARIANTS;
			sqlQuery += " WHERE ";
			
			String statement_type = SQLVariantTable.FIELDNAME_SVTYPE + "=" + SVType.TRA.getID() + " OR " + SQLVariantTable.FIELDNAME_SVTYPE + "=" + SVType.BND.getID();
			String statement_ctg1 = SQLVariantTable.FIELDNAME_CTG1 + " = ?";// + cuid; 0 3
			String statement_reg1 = SQLVariantTable.FIELDNAME_START1 + " < ?";// + end; 1
			String statement_reg2 = SQLVariantTable.FIELDNAME_END2 + " >= ?";// + start; 2
			
			//Combine for full statement
			/*
			 *  (!TRA) && ctg1 && reg1 && reg2
			 */
			
			String nottra = "(NOT " + statement_type + ")";
			String normQuery = nottra + " AND " + statement_ctg1 + " AND " + statement_reg1 + " AND " + statement_reg2;
			sqlQuery += normQuery;
			
			vars_in_reg_notra = connection.prepareStatement(sqlQuery);
			//System.err.println("StatementPrepper.getRegionVarGetterStatement || Statement is null? " + (vars_in_reg == null));
		}
		return vars_in_reg_notra;
	}
	
	public PreparedStatement getRegionNoTRAVarGetterStatement_ofType() throws SQLException
	{
		if(vars_in_reg_type == null)
		{
			String sqlQuery = "SELECT * FROM " + SQLVariantTable.TABLENAME_VARIANTS;
			sqlQuery += " WHERE ";
			
			String statement_type = SQLVariantTable.FIELDNAME_SVTYPE + "= ?";
			String statement_ctg1 = SQLVariantTable.FIELDNAME_CTG1 + " = ?";// + cuid; 0 3
			String statement_reg1 = SQLVariantTable.FIELDNAME_START1 + " < ?";// + end; 1
			String statement_reg2 = SQLVariantTable.FIELDNAME_END2 + " >= ?";// + start; 2
			
			//Combine for full statement
			/*
			 *  (!TRA) && ctg1 && reg1 && reg2
			 */
			
			String normQuery = statement_type + " AND " + statement_ctg1 + " AND " + statement_reg1 + " AND " + statement_reg2;
			sqlQuery += normQuery;
			
			vars_in_reg_type = connection.prepareStatement(sqlQuery);
			//System.err.println("StatementPrepper.getRegionVarGetterStatement || Statement is null? " + (vars_in_reg == null));
		}
		return vars_in_reg_type;
	}
	
	public PreparedStatement getFullInsertStatement() throws SQLException
	{
		if(insert_full == null)
		{
			String valStatement = "";
			for(int i = 0; i < SQLVariantTable.VAR_COLUMNS.length; i++)
			{
				valStatement += "?";
				if(i < SQLVariantTable.VAR_COLUMNS.length-1) valStatement += ", ";
			}
			
			String sqlQuery = "INSERT INTO " + SQLVariantTable.TABLENAME_VARIANTS + " VALUES (" + valStatement + ")";
			//System.err.println("-DEBUG- Raw Query: " + sqlQuery);
			insert_full = connection.prepareStatement(sqlQuery);
		}
		return insert_full;
	}
	
	public PreparedStatement getShortUpdateStatement() throws SQLException
	{
		if(short_update == null)
		{
			String valStatement = "";
			
			valStatement += SQLVariantTable.FIELDNAME_START1 + " = ?, ";// + var.getStartPosition().getStart() + ", ";
			valStatement += SQLVariantTable.FIELDNAME_START2 + " = ?, ";// + var.getStartPosition().getEnd() + ", ";
			valStatement += SQLVariantTable.FIELDNAME_END1 + " = ?, ";// + var.getEndPosition().getStart() + ", ";
			valStatement += SQLVariantTable.FIELDNAME_END2 + " = ?, ";// + var.getEndPosition().getEnd() + ", ";
			valStatement += SQLVariantTable.FIELDNAME_POSEFF + " = ?, ";// + var.getPositionEffect().getPriority() + ", ";
		
			valStatement += SQLVariantTable.FIELDNAME_ACOUNT_TOT + " = ?, ";// + var.getIndividualCount() + ", ";
			valStatement += SQLVariantTable.FIELDNAME_HCOUNT_TOT + " = ?, ";// + var.getHomozygoteCount() + ", ";
			valStatement += SQLVariantTable.FIELDNAME_ACOUNT_NFE + " = ?, ";// + var.getIndividualCount(Population.NFE) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_HCOUNT_NFE + " = ?, ";// + var.getHomozygoteCount(Population.NFE) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_ACOUNT_AFR + " = ?, ";// + var.getIndividualCount(Population.AFR) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_HCOUNT_AFR + " = ?, ";// + var.getHomozygoteCount(Population.AFR) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_ACOUNT_AMR + " = ?, ";// + var.getIndividualCount(Population.AMR) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_HCOUNT_AMR + " = ?, ";// + var.getHomozygoteCount(Population.AMR) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_ACOUNT_FIN + " = ?, ";// + var.getIndividualCount(Population.FIN) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_HCOUNT_FIN + " = ?, ";// + var.getHomozygoteCount(Population.FIN) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_ACOUNT_EAS + " = ?, ";// + var.getIndividualCount(Population.EAS) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_HCOUNT_EAS + " = ?, ";// + var.getHomozygoteCount(Population.EAS) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_ACOUNT_SAS + " = ?, ";// + var.getIndividualCount(Population.SAS) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_HCOUNT_SAS + " = ?, ";// + var.getHomozygoteCount(Population.SAS) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_ACOUNT_ASJ + " = ?, ";// + var.getIndividualCount(Population.ASJ) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_HCOUNT_ASJ + " = ?, ";// + var.getHomozygoteCount(Population.ASJ) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_ACOUNT_OTH + " = ?, ";// + var.getIndividualCount(Population.OTH) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_HCOUNT_OTH + " = ?, ";// + var.getHomozygoteCount(Population.OTH) + ", ";
			
			valStatement += SQLVariantTable.FIELDNAME_GENELIST + " = ?, ";
			valStatement += SQLVariantTable.FIELDNAME_VALNOTES + " = ?, ";
			valStatement += SQLVariantTable.FIELDNAME_GENOTYPES + " = ?";
			
			String sqlQuery = "UPDATE " + SQLVariantTable.TABLENAME_VARIANTS + " SET " + valStatement + " WHERE " + SQLVariantTable.FIELDNAME_VARUID + " = ?";
			//String sqlQuery = "UPDATE " + SQLVariantTable.TABLENAME_VARIANTS + " SET " + valStatement;
			//System.err.println("-DEBUG- Raw Query: " + sqlQuery);
			short_update = connection.prepareStatement(sqlQuery);
		}
		return short_update;
	}
	
	public PreparedStatement getPopUpdateStatement() throws SQLException
	{
		if(pop_update == null)
		{
			String valStatement = "";
			
			valStatement += SQLVariantTable.FIELDNAME_ACOUNT_TOT + " = ?, ";// + var.getIndividualCount() + ", ";
			valStatement += SQLVariantTable.FIELDNAME_HCOUNT_TOT + " = ?, ";// + var.getHomozygoteCount() + ", ";
			valStatement += SQLVariantTable.FIELDNAME_ACOUNT_NFE + " = ?, ";// + var.getIndividualCount(Population.NFE) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_HCOUNT_NFE + " = ?, ";// + var.getHomozygoteCount(Population.NFE) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_ACOUNT_AFR + " = ?, ";// + var.getIndividualCount(Population.AFR) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_HCOUNT_AFR + " = ?, ";// + var.getHomozygoteCount(Population.AFR) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_ACOUNT_AMR + " = ?, ";// + var.getIndividualCount(Population.AMR) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_HCOUNT_AMR + " = ?, ";// + var.getHomozygoteCount(Population.AMR) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_ACOUNT_FIN + " = ?, ";// + var.getIndividualCount(Population.FIN) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_HCOUNT_FIN + " = ?, ";// + var.getHomozygoteCount(Population.FIN) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_ACOUNT_EAS + " = ?, ";// + var.getIndividualCount(Population.EAS) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_HCOUNT_EAS + " = ?, ";// + var.getHomozygoteCount(Population.EAS) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_ACOUNT_SAS + " = ?, ";// + var.getIndividualCount(Population.SAS) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_HCOUNT_SAS + " = ?, ";// + var.getHomozygoteCount(Population.SAS) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_ACOUNT_ASJ + " = ?, ";// + var.getIndividualCount(Population.ASJ) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_HCOUNT_ASJ + " = ?, ";// + var.getHomozygoteCount(Population.ASJ) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_ACOUNT_OTH + " = ?, ";// + var.getIndividualCount(Population.OTH) + ", ";
			valStatement += SQLVariantTable.FIELDNAME_HCOUNT_OTH + " = ?";// + var.getHomozygoteCount(Population.OTH) + ", ";
			
			
			String sqlQuery = "UPDATE " + SQLVariantTable.TABLENAME_VARIANTS + " SET " + valStatement + " WHERE " + SQLVariantTable.FIELDNAME_VARUID + " = ?";
			pop_update = connection.prepareStatement(sqlQuery);
		}
		return pop_update;
	}

	public PreparedStatement getSampleGenoDeleteStatment() throws SQLException
	{
		if(sgeno_delete == null)
		{
			String sqlQuery = "DELETE FROM " + SQLVariantTable.TABLENAME_SAMPLEGENO;
			sqlQuery += " WHERE " + SQLVariantTable.FIELDNAME_SAMPLEUID + " = ?";
			sgeno_delete = connection.prepareStatement(sqlQuery);
		}
		return sgeno_delete;
	}
	
	public PreparedStatement getVarDeleteStatment() throws SQLException
	{
		if(var_delete == null)
		{
			String sqlQuery = "DELETE FROM " + SQLVariantTable.TABLENAME_VARIANTS;
			sqlQuery += " WHERE " + SQLVariantTable.FIELDNAME_VARUID + " = ?";
			var_delete = connection.prepareStatement(sqlQuery);
		}
		return var_delete;
	}
	
	public PreparedStatement generateMultiVarDeleteStatement(int count) throws SQLException
	{
		String sqlQuery = "DELETE FROM " + SQLVariantTable.TABLENAME_VARIANTS + " WHERE";
		for(int i = 0; i < count; i++)
		{
			sqlQuery += " " + SQLVariantTable.FIELDNAME_VARUID + " = ?";
			if(i < (count - 1)) sqlQuery += " OR";
		}
		return connection.prepareStatement(sqlQuery);
	}
	
	public PreparedStatement getSGenoInsertStatement() throws SQLException
	{
		if(sgeno_insert == null)
		{
			String sqlQuery = "INSERT INTO " + SQLVariantTable.TABLENAME_SAMPLEGENO;
			sqlQuery += " VALUES (?, ?, ?, ?)";
			
			sgeno_insert = connection.prepareStatement(sqlQuery);
		}
		return sgeno_insert;
	}
	
	public PreparedStatement getSGenoUpdateStatement() throws SQLException
	{
		if(sgeno_update == null)
		{
			String sqlCmd = "UPDATE " + SQLVariantTable.TABLENAME_SAMPLEGENO;
			sqlCmd += " SET ";
			sqlCmd += SQLVariantTable.FIELDNAME_SVARLIST_HOM + " = ?, ";
			sqlCmd += SQLVariantTable.FIELDNAME_SVARLIST_HET + " = ?, ";
			sqlCmd += SQLVariantTable.FIELDNAME_SVARLIST_OTH + " = ?";
			sqlCmd += " WHERE " + SQLVariantTable.FIELDNAME_SAMPLEUID + " = ?";
			
			
			sgeno_update = connection.prepareStatement(sqlCmd);
		}
		return sgeno_update;
	}
	
	public PreparedStatement getGeneHitGetterStatement() throws SQLException
	{
		if(gh_get_one == null)
		{
			String sqlQuery = "SELECT * FROM " + SQLVariantTable.TABLENAME_GENEHITS;
			sqlQuery += " WHERE " + SQLVariantTable.FIELDNAME_GH_GENEUID + " = ?";
			gh_get_one = connection.prepareStatement(sqlQuery);
		}
		return gh_get_one;
	}
	
	public PreparedStatement getGeneHitGetAllStatement() throws SQLException
	{
		String sqlQuery = "SELECT * FROM " + SQLVariantTable.TABLENAME_GENEHITS;
		return connection.prepareStatement(sqlQuery);
	}
	
	public PreparedStatement getGeneHitInsertStatement() throws SQLException
	{
		if(gh_insert == null)
		{
			String valStatement = "";
			for(int i = 0; i < SQLVariantTable.GENEHITS_COLUMNS.length; i++)
			{
				valStatement += "?";
				if(i < SQLVariantTable.GENEHITS_COLUMNS.length-1) valStatement += ", ";
			}
			
			String sqlQuery = "INSERT INTO " + SQLVariantTable.TABLENAME_GENEHITS + " VALUES (" + valStatement + ")";
			gh_insert = connection.prepareStatement(sqlQuery);
		}
		return gh_insert;
	}
	
	public PreparedStatement getGeneHitUpdateStatement() throws SQLException
	{
		if(gh_update == null)
		{
			String valStatement = "";
			
			valStatement += SQLVariantTable.FIELDNAME_GH_HITS_T + " = ?, ";
			valStatement += SQLVariantTable.FIELDNAME_GH_HITS_E + " = ?, ";
			valStatement += SQLVariantTable.FIELDNAME_GH_HITS_TI + " = ?, ";
			valStatement += SQLVariantTable.FIELDNAME_GH_HITS_EI + " = ?";
			
			String sqlQuery = "UPDATE " + SQLVariantTable.TABLENAME_GENEHITS + " SET " + valStatement + " WHERE " + SQLVariantTable.FIELDNAME_GH_GENEUID + " = ?";
			//System.err.println("SQL Query: " + sqlQuery);
			gh_update = connection.prepareStatement(sqlQuery);
		}
		return gh_update;
	}
	
	public PreparedStatement getGeneHitTableWipeStatement() throws SQLException
	{
		String sqlQuery = "DELETE FROM " + SQLVariantTable.TABLENAME_GENEHITS;
		return connection.prepareStatement(sqlQuery);
	}
	
	public PreparedStatement getVarTableWipeStatement() throws SQLException
	{
		String sqlQuery = "DELETE FROM " + SQLVariantTable.TABLENAME_VARIANTS;
		return connection.prepareStatement(sqlQuery);
	}
	
	public PreparedStatement getSampleGenoTableWipeStatement() throws SQLException
	{
		String sqlQuery = "DELETE FROM " + SQLVariantTable.TABLENAME_SAMPLEGENO;
		return connection.prepareStatement(sqlQuery);
	}
	
	public PreparedStatement getVariantGetAllStatement() throws SQLException
	{
		String sqlQuery = "SELECT * FROM " + SQLVariantTable.TABLENAME_VARIANTS;
		return connection.prepareStatement(sqlQuery);
	}
	
	public PreparedStatement getVariantGetAllIDsStatement() throws SQLException
	{
		String sqlQuery = "SELECT " + SQLVariantTable.FIELDNAME_VARUID + " FROM " + SQLVariantTable.TABLENAME_VARIANTS;
		return connection.prepareStatement(sqlQuery);
	}
	
	public PreparedStatement getSGenoGetAllStatement() throws SQLException
	{
		String sqlQuery = "SELECT * FROM " + SQLVariantTable.TABLENAME_SAMPLEGENO;
		return connection.prepareStatement(sqlQuery);
	}
	
	public Blob wrapInBlob(byte[] bytes) throws SQLException
	{
		//System.err.println("Blob length: " + bytes.length);
		Blob b = connection.createBlob();
		//int written = b.setBytes(1, bytes);
		//System.err.println("Blob written: " + written);
		b.setBytes(1, bytes);
		return b;
	}
	
}
