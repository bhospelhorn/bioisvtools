package hospelhornbg_svdb;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;

public class SQLManager {
	
	/*--- Key Constants ---*/
	
	public static final String SKEY_CHECKVARUID = "check_varUID";
	
	public static final String SKEY_GETVAR = "get_var";
	//public static final String SKEY_GETVAR_MULTI = "get_var_multi";
	public static final String SKEY_GETGENO = "get_geno";
	//public static final String SKEY_GETGENO_MULTI = "get_geno_multi";
	public static final String SKEY_GETSAMPLEVAR = "get_smplvar";
	
	public static final String SKEY_GETVAR_REG = "vars_in_reg";
	public static final String SKEY_GETVAR_REG_NOTRA = "vars_in_reg_no_TRA";
	public static final String SKEY_GETVAR_REG_OFTYPE = "vars_in_reg_ofType";
	
	public static final String SKEY_VAR_INSERT = "insertNewVar";
	public static final String SKEY_VAR_UPDATE_SHORT = "updateVarShort";
	public static final String SKEY_VAR_UPDATE_POP = "updateVarPopulations";
	public static final String SKEY_VAR_DELETE = "deleteVar";
	
	public static final String SKEY_SGENO_INSERT = "insert_new_sampleGeno";
	public static final String SKEY_SGENO_UPDATE = "update_sampleGeno";
	public static final String SKEY_SGENO_DELETE = "delete_sampleGeno";
	
	public static final String SKEY_GENEHITS_GET = "get_genehit_record";
	public static final String SKEY_GENEHITS_INSERT = "insert_genehit_record";
	public static final String SKEY_GENEHITS_UPDATE = "update_genehit_record";
	
	public static final String[] ALL_KEYS = {SKEY_CHECKVARUID, SKEY_GETVAR, SKEY_GETGENO, SKEY_GETSAMPLEVAR,
											 SKEY_GETVAR_REG, SKEY_GETVAR_REG_NOTRA, SKEY_GETVAR_REG_OFTYPE,
											 SKEY_VAR_INSERT, SKEY_VAR_UPDATE_SHORT, SKEY_VAR_UPDATE_POP,
											 SKEY_VAR_DELETE, SKEY_SGENO_INSERT, SKEY_SGENO_UPDATE,
											 SKEY_SGENO_DELETE, SKEY_GENEHITS_GET, SKEY_GENEHITS_INSERT,
											 SKEY_GENEHITS_UPDATE};
	
	/*--- Instance Variables ---*/
	
	private StatementPrepper statement_prepper;
	
	private ConcurrentHashMap<String, PreparedStatement> psLib;
	private ConcurrentHashMap<String, Thread> psCleared;
	private ConcurrentHashMap<String, ConcurrentLinkedQueue<Thread>> psWaitList;
	
	private volatile Thread exeReady;
	private ConcurrentLinkedQueue<Thread> exeQueue;
	
	/*--- Initialization ---*/
	
	public SQLManager(Connection c) throws SQLException
	{
		statement_prepper = new StatementPrepper(c);
		psLib = new ConcurrentHashMap<String, PreparedStatement>();
		psCleared = new ConcurrentHashMap<String, Thread>();
		psWaitList = new ConcurrentHashMap<String, ConcurrentLinkedQueue<Thread>>();
		exeQueue = new ConcurrentLinkedQueue<Thread>();
		
		populateLibrary();
	}
	
	private void populateLibrary() throws SQLException
	{
		psLib.put(SKEY_CHECKVARUID, statement_prepper.getVarUIDCheckStatement());
		psLib.put(SKEY_GETVAR, statement_prepper.getVarGetterStatement());
		psLib.put(SKEY_GETGENO, statement_prepper.getGenoGetterStatement());
		psLib.put(SKEY_GETSAMPLEVAR, statement_prepper.getSampleVarGetterStatement());
		
		psLib.put(SKEY_GETVAR_REG, statement_prepper.getRegionVarGetterStatement());
		psLib.put(SKEY_GETVAR_REG_NOTRA, statement_prepper.getRegionNoTRAVarGetterStatement());
		psLib.put(SKEY_GETVAR_REG_OFTYPE, statement_prepper.getRegionNoTRAVarGetterStatement_ofType());
	
		psLib.put(SKEY_VAR_INSERT, statement_prepper.getFullInsertStatement());
		psLib.put(SKEY_VAR_UPDATE_SHORT, statement_prepper.getShortUpdateStatement());
		psLib.put(SKEY_VAR_UPDATE_POP, statement_prepper.getPopUpdateStatement());
		psLib.put(SKEY_VAR_DELETE, statement_prepper.getVarDeleteStatment());
		
		psLib.put(SKEY_SGENO_INSERT, statement_prepper.getSGenoInsertStatement());
		psLib.put(SKEY_SGENO_UPDATE, statement_prepper.getSGenoUpdateStatement());
		psLib.put(SKEY_SGENO_DELETE, statement_prepper.getSampleGenoDeleteStatment());
		
		psLib.put(SKEY_GENEHITS_GET, statement_prepper.getGeneHitGetterStatement());
		psLib.put(SKEY_GENEHITS_INSERT, statement_prepper.getGeneHitInsertStatement());
		psLib.put(SKEY_GENEHITS_UPDATE, statement_prepper.getGeneHitUpdateStatement());
		
		for(String k : ALL_KEYS)
		{
			psWaitList.put(k, new ConcurrentLinkedQueue<Thread>());
		}
	}
	
	/*--- Statement Access Concurrency ---*/
	
	public PreparedStatement requestStatement(String key, boolean ignoreExternal)
	{
		//System.err.println(Thread.currentThread().getName() + " requested " + key);
		if(psCleared.get(key) == null)
		{
			//Free pickin'
			//System.err.println("No threads ahead in queue. " + Thread.currentThread().getName() + " checking for " + key);
			PreparedStatement ps = psLib.remove(key);
			if(ps != null)
			{
				//System.err.println(Thread.currentThread().getName() + " granted " + key);
				return ps;	
			}
		}
		
		//Put on wait list and block until cleared
		ConcurrentLinkedQueue<Thread> waitlist = psWaitList.get(key);
		if(waitlist == null) return null; //Key is invalid
		Thread me = Thread.currentThread();
		waitlist.add(me);
		while(psCleared.get(key) != me && (psLib.containsKey(key)))
		{
			try 
			{
				//We don't want it checking synchronously very often
				Thread.sleep(10000);
			} 
			catch (InterruptedException e) 
			{
				//Check again!
				if(psCleared.get(key) != me)
				{
					//We assume it was an external signal.
					if(!ignoreExternal)
					{
						//System.err.println(Thread.currentThread().getName() + " received interrupt. Request for " + key + " withdrawn.");
						waitlist.remove(me);
						return null;	
					}
				}
			}
		}
		
		//System.err.println(Thread.currentThread().getName() + " granted " + key);
		psCleared.remove(key);
		return psLib.remove(key);
	}
	
	public void releaseStatement(String key, PreparedStatement ps)
	{
		//System.err.println(Thread.currentThread().getName() + " requesting release of " + key);
		ConcurrentLinkedQueue<Thread> waitlist = psWaitList.get(key);
		if(waitlist == null) return; //Key is invalid
		Thread cleared = null;
		if(!waitlist.isEmpty())
		{
			cleared = waitlist.poll();
			psCleared.put(key, cleared);
			//System.err.println(cleared.getName() + " queued next for " + key);
		}
		//If the waiting thread happens to check after its marked as cleared
		// but before the statement is put back, it'll grab a null statement!
		
		psLib.put(key, ps);
		if(cleared != null)
		{
			synchronized(cleared) {cleared.interrupt();}
		}
		//System.err.println(Thread.currentThread().getName() + " released " + key);
	}
	
	/*--- Database Access Concurrency ---*/
	
	public boolean requestStatementExecution(boolean ignoreExternalInterrupts)
	{
		//Blocks until cleared for statement execution
		System.err.println(Thread.currentThread().getName() + " requesting statement execution!");
		synchronized(this)
		{
			if(exeReady == null)
			{
				//Put me!
				exeReady = Thread.currentThread();
				System.err.println(Thread.currentThread().getName() + " granted statement execution!");
				return true;
			}
		}
		
		Thread me = Thread.currentThread();
		exeQueue.add(me);	
		
		while(exeReady != me)
		{
			try 
			{
				Thread.sleep(10000);
			} 
			catch (InterruptedException e) 
			{
				if(exeReady != me)
				{
					if(!ignoreExternalInterrupts)
					{
						System.err.println(Thread.currentThread().getName() + " received external interrupt. Withdrawing execution request...");
						exeQueue.remove(me);
						return false;	
					}
				}
				break;
			}
		}
		
		System.err.println(Thread.currentThread().getName() + " granted statement execution!");
		return true;
	}
	
	public boolean acknowledgeStatementExecution()
	{
		//If this thread was the only one cleared to execute statement, this
		//	allows queue to move forward
		Thread me = Thread.currentThread();
		if(exeReady != me) return false;
		System.err.println(Thread.currentThread().getName() + " requesting acknowledgement of statement execution");
		
		if(!exeQueue.isEmpty())
		{
			synchronized(this) {exeReady = exeQueue.poll();}
			synchronized(exeReady) {exeReady.interrupt();}
		}
		else exeReady = null;
		
		System.err.println(Thread.currentThread().getName() + " acknowledgement of execution received");
		return true;
	}

	/*--- Generated Statements ---*/
	
	public PreparedStatement generateMultiVarGetterStatement(int count) throws SQLException
	{
		return statement_prepper.generateMultiVarGetterStatement(count);
	}
	
	public PreparedStatement generateMultiVarDeleteStatement(int count) throws SQLException
	{
		return statement_prepper.generateMultiVarDeleteStatement(count);
	}
	
	public PreparedStatement generateMultiGenoGetterStatement(int count) throws SQLException
	{
		return statement_prepper.generateMultiGenoGetterStatement(count);
	}
	
	/*--- Direct Access ---*/
	
	protected StatementPrepper getStatementGenerator()
	{
		return statement_prepper;
	}
	
}
