package hospelhornbg_svtools;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ConcurrentLinkedQueue;

import hospelhornbg_bioinformatics.SAMHeaderLine;
import hospelhornbg_bioinformatics.SAMRecord;
import hospelhornbg_bioinformatics.SAMRecord.FailFlags;
import hospelhornbg_bioinformatics.SAMRecord.InvalidSAMRecordException;
import hospelhornbg_bioinformatics.SAMRecord.ParsedRecord;
import hospelhornbg_bioinformatics.SAMRecord.WarningFlags;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class SamScanner {
	
	public static final String OP_INPUT = "-i"; 
	public static final String OP_TEMPDIR = "-T"; 
	public static final String OP_THREADS = "-t"; 
	
	public static final int MAX_QUEUESIZE = 1000000000;
	
	private static class Counts
	{
		private volatile long total;
		
		private volatile int err_syntax_general;
		private volatile int err_nullseq_nnqual;
		//private volatile int err_bad_customFieldType;
		private HashMap<String, Integer> bad_customFieldTypes;
		
		//private volatile int err_invalid_RNAME;
		private volatile int err_invalid_POS;
		//private volatile int err_invalid_RNEXT;
		private volatile int err_invalid_PNEXT;
		private HashMap<String, Long> bad_RNAME;
		private HashMap<String, Long> bad_RNEXT;
		
		private HashMap<String, Integer> all_RNAME;
		private HashMap<String, Integer> all_RNEXT;
		
		private volatile int err_qualstr_len_bad;
		
		public Counts()
		{
			total = 0;
			err_syntax_general = 0;
			err_nullseq_nnqual = 0;
			//err_bad_customFieldType = 0;
			bad_customFieldTypes = new HashMap<String, Integer>();
			
			//err_invalid_RNAME = 0;
			err_invalid_POS = 0;
			//err_invalid_RNEXT = 0;
			err_invalid_PNEXT = 0;
			bad_RNAME = new HashMap<String, Long>();
			bad_RNEXT = new HashMap<String, Long>();
			all_RNAME = new HashMap<String, Integer>();
			all_RNEXT = new HashMap<String, Integer>();
			
			err_qualstr_len_bad = 0;
		}
		
		public synchronized void incrementTotal()
		{
			total++;
		}
		
		public synchronized void incrementSyntaxErrCount()
		{
			err_syntax_general++;
		}
		
		public synchronized void incrementNullSeqErrCount()
		{
			err_nullseq_nnqual++;
		}
		
		public synchronized void addBadCustomFieldErr(String badtype)
		{
			//err_bad_customFieldType++;
			Integer count = bad_customFieldTypes.get(badtype);
			if (count == null) bad_customFieldTypes.put(badtype, 1);
			else bad_customFieldTypes.put(badtype, count+1);
		}
		
		public synchronized void addInvalidRefContig(String c)
		{
			if(!bad_RNAME.containsKey(c)) bad_RNAME.put(c, 1L);
			else
			{
				Long count = bad_RNAME.get(c);
				bad_RNAME.put(c, count+1L);
			}
		}
		
		public synchronized void addInvalidRNextContig(String c)
		{
			Long count = bad_RNEXT.get(c);
			if (count == null) bad_RNEXT.put(c, 1L);
			else bad_RNEXT.put(c, count+1L);
		}
		
		public synchronized void incrementInvalidPosErrCount()
		{
			err_invalid_POS++;
		}
		
		public synchronized void incrementInvalidPNextCount()
		{
			err_invalid_PNEXT++;
		}
		
		public synchronized void incrementQualStrLengthErrCount()
		{
			err_qualstr_len_bad++;
		}
		
		public synchronized void add_RNAME(String rname)
		{
			Integer count = all_RNAME.get(rname);
			if (count == null) all_RNAME.put(rname, 1);
			else all_RNAME.put(rname, count+1);
		}
		
		public synchronized void add_RNEXT(String rname)
		{
			Integer count = all_RNEXT.get(rname);
			if (count == null) all_RNEXT.put(rname, 1);
			else all_RNEXT.put(rname, count+1);
		}
		
		public long getTotal()
		{
			return total;
		}
		
		public int getSyntaxErrCount()
		{
			return err_syntax_general;
		}
		
		public int getNullSeqErrCount()
		{
			return err_nullseq_nnqual;
		}
		
		public int getInvalidPosErrCount()
		{
			return err_invalid_POS;
		}
		
		public int getInvalidPNextCount()
		{
			return err_invalid_PNEXT;
		}
		
		public int getQualStrLengthErrCount()
		{
			return err_qualstr_len_bad;
		}
		
		public synchronized int getBadCustomFieldErrCount()
		{
			return bad_customFieldTypes.size();
		}
		
		public synchronized long getInvalidRefContigErrCount()
		{
			long tot = 0;
			Set<String> kset = bad_RNAME.keySet();
			for (String k : kset)
			{
				Long l = bad_RNAME.get(k);
				if (l != null) tot += l;
			}
			return tot;
			//return bad_RNAME.size();
		}
		
		public synchronized long getInvalidRNextContigErrCount()
		{
			long tot = 0;
			Set<String> kset = bad_RNEXT.keySet();
			for (String k : kset)
			{
				Long l = bad_RNEXT.get(k);
				if (l != null) tot += l;
			}
			return tot;
		}
		
		public synchronized List<String> getBadCustomFieldTypes()
		{
			int lsz = bad_customFieldTypes.size() + 1;
			List<String> list = new ArrayList<String>(lsz);
			list.addAll(bad_customFieldTypes.keySet());
			Collections.sort(list);
			return list;
		}
		
		public synchronized List<String> getBadRefContigs()
		{
			int lsz = bad_RNAME.size() + 1;
			List<String> list = new ArrayList<String>(lsz);
			list.addAll(bad_RNAME.keySet());
			Collections.sort(list);
			return list;
		}
		
		public synchronized List<String> getBadRNextContigs()
		{
			int lsz = bad_RNEXT.size() + 1;
			List<String> list = new ArrayList<String>(lsz);
			list.addAll(bad_RNEXT.keySet());
			Collections.sort(list);
			return list;
		}
		
		public synchronized List<String> getAll_RNAME_Strings()
		{
			int lsz = all_RNAME.size() + 1;
			List<String> list = new ArrayList<String>(lsz);
			list.addAll(all_RNAME.keySet());
			Collections.sort(list);
			return list;
		}
		
		public synchronized List<String> getAll_RNEXT_Strings()
		{
			int lsz = all_RNEXT.size() + 1;
			List<String> list = new ArrayList<String>(lsz);
			list.addAll(all_RNEXT.keySet());
			Collections.sort(list);
			return list;
		}
		
	}
	
	public static void printUsage()
	{
		System.out.println("--------------------------------------------------------------------------------");
		System.out.println("BioisvTools || SAM Scanner");
		System.out.println();
		System.out.println("Purpose: Scans SAM files and outputs information about syntax errors that are present");
		System.out.println("or oversights that may make the SAM unreadable to some tools.");
		System.out.println();
		System.out.println("Input Formats:");
		System.out.println("\tInput must be a SAM file");
		System.out.println();
		System.out.println("Flags:");
		System.out.println("\t-i\tFILE\t[Optional]\t\tInput file path. Defaults to stdin stream if not provided.");
		System.out.println("\t-t\tINT\t[Optional]\t\tNumber of threads to use.");
		System.out.println();
		System.out.println("Note:");
		System.out.println("You may specify a genome build if you wish to check the @SQ, REF, and RNEXT contigs");
		System.out.println("against a known build to see if there are any contig strings other tools may not");
		System.out.println("recognize.");
		System.out.println("If you do not specify a genome build, the scanner uses the @SQ lines for reference.");
		System.out.println();
		System.out.println("Sample Usage:");
		System.out.println("java -jar bioisvtools.jar scansam -g hg38");
		System.out.println("java -jar bioisvtools.jar scansam -v -t 16");
		System.out.println("java -jar bioisvtools.jar scansam -g GRCh37 -i sample.sam -t 24");
		System.out.println();
		System.out.println("--------------------------------------------------------------------------------");
	}
	
	private static Contig generateContig(SAMHeaderLine sq)
	{
		if (sq.getFieldCount("SN") > 1)
		{
			//Duplication! Bad!
			System.out.println("ERROR   || @SQ line has multiple sequence name records!");
			System.out.println("\t" + sq);
			return null;
		}
		String cname = sq.getValue("SN");
		if (sq.getFieldCount("LN") > 1)
		{
			//Duplication! Bad!
			System.out.println("ERROR   || @SQ line has multiple sequence length records!");
			System.out.println("\t" + sq);
			return null;
		}
		String clenstr = sq.getValue("LN");
		long clen = -1;
		if(clenstr != null)
		{
			try
			{
				clen = Long.parseUnsignedLong(clenstr);
			}
			catch(NumberFormatException e)
			{
				System.out.println("ERROR   || @SQ LN value invalid (" + clenstr + ")");
				System.out.println("\t" + sq);
				return null;
			}
		}
		if (cname == null)
		{
			System.out.println("ERROR   || @SQ line has no SN field!");
			System.out.println("\t" + sq);
			return null;
		}
		else
		{
			if (clen < 0)
			{
				System.out.println("ERROR   || @SQ line has no LN field!");
				System.out.println("\t" + sq);
				return null;
			}
			else
			{
					Contig c = new Contig();	
					c.setUDPName(cname);
					c.setUCSCName(cname);
					c.setLength(clen);
					return c;
			}
		}
	}
	
	private static GenomeBuild scanSamHeader(BufferedReader reader, GenomeBuild gb) throws IOException
	{
		//Read header
		//Check for HD line (which is mandatory)
		//If there is a provided GenomeBuild, report any SQ lines that have a contig not in genome build
		//	and any genome build contigs not in SQ (remove these from local gb copy so that warnings will be thrown if these contigs show up
		//	in records
		//If there is no provided GenomeBuild, build one from the SQ lines and regurgitate for user
		//Check for RG line
		reader.mark(1000000);
		System.out.println("===================== HEADER =====================");
		String line = reader.readLine();
		//This line should be the @HD. If not, report.
		if (!(line.startsWith("@HD")))
		{
			System.out.println("ERROR   || @HD line was NOT found!");
			if (!line.startsWith("@")) {
				System.out.println("ERROR   || SAM lacks header!");
				reader.reset();
			}
		}
		else
		{
			System.out.println("INFO   || @HD line found!");
			String[] fields = line.split("\t");
			int fc = fields.length;
			for (int i = 1; i < fc; i++) System.out.println("\t" + fields[i]);
		}
		
		Set<SAMHeaderLine> sqlines = new HashSet<SAMHeaderLine>();
		boolean rgfound = false;
		boolean rgsmfound = false;
		line = reader.readLine();
		while(line.startsWith("@"))
		{
			SAMHeaderLine hl = null;
			try {hl = new SAMHeaderLine(line);}
			catch(FileBuffer.UnsupportedFileTypeException e) {
				System.out.println("ERROR   || Syntax Error - Improperly formatted header line:");
				System.out.println("\t" + line);
				continue;
			};
			
			if (hl.getKey().equals("SQ")) sqlines.add(hl);
			else if (hl.getKey().equals("RG"))
			{
				rgfound = true;
				String sm = hl.getValue("SM");
				if (sm != null)
				{
					System.out.println("INFO   || @RG line found! Sample Name: " + sm);
					rgsmfound = true;
				}
			}
			line = reader.readLine();
		}
		
		//Handle RG
		if (!rgfound) System.out.println("WARNING || No @RG line was found! Some packages require this line to read a SAM/BAM file.");
		if (rgfound && !rgsmfound) System.out.println("WARNING || @RG line(s) found, but there were no SM (sample) fields. Some packages require this to obtain the sample name!");
		
		//Handle SQ set
		
		if (gb != null)
		{
			Set<Contig> matched = new HashSet<Contig>();
			for (SAMHeaderLine sq : sqlines)
			{
				Contig c = generateContig(sq);
				if (c != null)
				{
					String cname = c.getUDPName();
					Contig o = gb.getContig(cname);
					if (o != null)
					{
						//Lengths must match!
						long mylen = c.getLength();
						long olen = o.getLength();
						if (mylen != olen)
						{
							System.out.println("ERROR   || Contig found in build, but length is mismatched: " + cname + " [Build:" + olen + "][Local:" + mylen + "]");
						}
						else matched.add(o);
					}
					else
					{
						//No match
						System.out.println("ERROR   || Contig not found in build: " + cname);
					}
				}
			}
			
			//Find contigs in build that are NOT in SAM
			List<Contig> allc = gb.getChromosomes();
			for (Contig c : allc)
			{
				if(!matched.contains(c))
				{
					System.out.println("WARNING || Build contig " + c.getUDPName() + " was not found in SAM file header!");
					gb.removeContig(c.getUDPName());
				}
			}
			
			reader.reset();
			return gb;
		}
		else
		{
			GenomeBuild ngb = new GenomeBuild("UNKNOWN", "UNKNOWN");
			for (SAMHeaderLine sq : sqlines)
			{
				Contig c = generateContig(sq);
				if (c != null)
				{
					String cname = c.getUDPName();
					Contig o = ngb.getContig(cname);
					if (o != null)
					{
						System.out.println("ERROR   || Duplicate contig (" + cname + ")");
					}
					else ngb.addContig(c);
				}
			}
			
			reader.reset();
			return ngb;
		}
	}
	
	private static Counts getCounts(BufferedReader reader, GenomeBuild gb, boolean verbose) throws IOException
	{
		//Skip header lines
		String line = reader.readLine();
		while(line.startsWith("@")) line = reader.readLine();
		
		//Counts
		Counts c = new Counts();
		
		//Reads
		while (line != null)
		{
			c.incrementTotal();
			try 
			{
				ParsedRecord pr = SAMRecord.parseSAMRecord(line, gb, verbose);
				SAMRecord rec = pr.getRecord();
				WarningFlags wf = rec.getParserWarnings();
				if(wf.err_qualstr_len_bad) c.incrementQualStrLengthErrCount();
				if(wf.err_invalid_RNAME != null) c.addInvalidRefContig(wf.err_invalid_RNAME);
				if(wf.err_invalid_RNEXT != null) c.addInvalidRNextContig(wf.err_invalid_RNEXT);
				if(wf.err_invalid_POS != 0) c.incrementInvalidPosErrCount();
				if(wf.err_invalid_PNEXT != 0) c.incrementInvalidPNextCount();
				c.add_RNAME(pr.raw_rname);
				c.add_RNEXT(pr.raw_rnext);
			} 
			catch (UnsupportedFileTypeException e)
			{
				c.incrementSyntaxErrCount();
			} 
			catch (InvalidSAMRecordException e) 
			{
				FailFlags ff = e.getFlags();
				if (ff.err_syntax) c.incrementSyntaxErrCount();
				if (ff.err_nullseq_nnqual) c.incrementNullSeqErrCount();
				if (ff.err_bad_customFieldType != null) c.addBadCustomFieldErr(ff.err_bad_customFieldType);
			}
			line = reader.readLine();
		}
		
		return c;
	}
	
	private static Counts getCounts(BufferedReader reader, GenomeBuild gb, boolean verbose, int threads) throws IOException
	{
		//Counts
		Counts c = new Counts();
		
		//Threads
		int tcount = threads - 1;
		Thread[] tarr = new Thread[tcount];
		ConcurrentLinkedQueue<String> linequeue = new ConcurrentLinkedQueue<String>();
		for (int i = 0; i < tcount; i++)
		{
			tarr[i] = new Thread(new Runnable() {
				@Override
				public void run() 
				{
					int pcount = 0;
					while(!Thread.interrupted())
					{
						if(linequeue.isEmpty()) 
						{
							try 
							{
								Thread.sleep(10);
							} 
							catch (InterruptedException e) 
							{
								//If it's getting interrupted, that should only be a kill signal...
								if(verbose)System.err.println("DEBUG: Lines Parsed " + pcount + " (Worker: " + Thread.currentThread().getName() + ") - Worker thread returning... (From queue empty sleep)");
								return;
							}
						}
						else
						{
							String line = linequeue.poll();
							if (line == null) continue;
							pcount++;
							c.incrementTotal();
							try 
							{
								ParsedRecord pr = SAMRecord.parseSAMRecord(line, gb, verbose);
								SAMRecord rec = pr.getRecord();
								WarningFlags wf = rec.getParserWarnings();
								if(wf.err_qualstr_len_bad) c.incrementQualStrLengthErrCount();
								if(wf.err_invalid_RNAME != null) c.addInvalidRefContig(wf.err_invalid_RNAME);
								if(wf.err_invalid_RNEXT != null) c.addInvalidRNextContig(wf.err_invalid_RNEXT);
								if(wf.err_invalid_POS != 0) c.incrementInvalidPosErrCount();
								if(wf.err_invalid_PNEXT != 0) c.incrementInvalidPNextCount();
								c.add_RNAME(pr.raw_rname);
								c.add_RNEXT(pr.raw_rnext);
							} 
							catch (UnsupportedFileTypeException e)
							{
								c.incrementSyntaxErrCount();
							} 
							catch (InvalidSAMRecordException e) 
							{
								FailFlags ff = e.getFlags();
								if (ff.err_syntax) c.incrementSyntaxErrCount();
								if (ff.err_nullseq_nnqual) c.incrementNullSeqErrCount();
								if (ff.err_bad_customFieldType != null) c.addBadCustomFieldErr(ff.err_bad_customFieldType);
							}
							//if(pcount % 1000000 == 0) System.err.println("DEBUG: Lines Parsed " + pcount + " (Worker: " + Thread.currentThread().getName() + ")");
						}		
					}
					//Run until the queue is empty...
					while(!linequeue.isEmpty())
					{
						String line = linequeue.poll();
						if (line == null) break;
						pcount++;
						c.incrementTotal();
						try 
						{
							ParsedRecord pr = SAMRecord.parseSAMRecord(line, gb, verbose);
							SAMRecord rec = pr.getRecord();
							WarningFlags wf = rec.getParserWarnings();
							if(wf.err_qualstr_len_bad) c.incrementQualStrLengthErrCount();
							if(wf.err_invalid_RNAME != null) c.addInvalidRefContig(wf.err_invalid_RNAME);
							if(wf.err_invalid_RNEXT != null) c.addInvalidRNextContig(wf.err_invalid_RNEXT);
							if(wf.err_invalid_POS != 0) c.incrementInvalidPosErrCount();
							if(wf.err_invalid_PNEXT != 0) c.incrementInvalidPNextCount();
							c.add_RNAME(pr.raw_rname);
							c.add_RNEXT(pr.raw_rnext);
						} 
						catch (UnsupportedFileTypeException e)
						{
							c.incrementSyntaxErrCount();
						} 
						catch (InvalidSAMRecordException e) 
						{
							FailFlags ff = e.getFlags();
							if (ff.err_syntax) c.incrementSyntaxErrCount();
							if (ff.err_nullseq_nnqual) c.incrementNullSeqErrCount();
							if (ff.err_bad_customFieldType != null) c.addBadCustomFieldErr(ff.err_bad_customFieldType);
						}
					}
					if(verbose)System.err.println("DEBUG: Lines Parsed " + pcount + " (Worker: " + Thread.currentThread().getName() + ") - Worker thread returning...");
				}
			});
			tarr[i].setName("SamScannerWorkerThread" + i);
		}
		
		//Start worker threads
		for (Thread t : tarr) {
			t.start();
			if(verbose)System.err.println("DEBUG: Worker thread " + t.getName() + " started!");
		}
		
		//Start line reading
		String line = null;
		long lcount = 0;
		int hcount = 0;
		while((line = reader.readLine()) != null)
		{
			if (line.startsWith("@")) {
				hcount++;
				continue; //Skip header
			}
			lcount++;
			//pause if linequeue is too full
			while (linequeue.size() >= MAX_QUEUESIZE)
			{
				try 
				{
					Thread.sleep(10);
				} 
				catch (InterruptedException e) 
				{
					//Eat
					e.printStackTrace();
				}
			}
			linequeue.add(line);
			if(verbose && (lcount % 10000000 == 0)) System.err.println("DEBUG: Lines Read: " + lcount);
		}
		if(verbose)System.err.println("DEBUG: SAM readin complete. Total records read: " + lcount);
		if(verbose)System.err.println("DEBUG: Total header lines read: " + hcount);
		
		//Kill worker threads
		for (Thread t : tarr){
			synchronized(t) { t.interrupt();}
			if(verbose)System.err.println("DEBUG: Worker thread " + t.getName() + " interrupted!");
		}
		
		//Wait for worker threads to die
		if(verbose)System.err.println("DEBUG: Waiting for worker threads to die...");
		boolean anyalive = true;
		while(anyalive)
		{
			anyalive = false;
			for (Thread t : tarr)
			{
				if (t.isAlive())
				{
					anyalive = true;
					break;
				}
			}	
		}
		if(verbose)System.err.println("DEBUG: Worker threads dead.");
		
		return c;
	}
	
	public static void scanSam(BufferedReader reader, GenomeBuild gb, boolean verbose, int threads) throws IOException
	{
		//Read header (and reset stream)
		gb = scanSamHeader(reader, gb);
		
		//Dump contig list to stdout before looking at reads
		System.out.println("--------- Generated Genome Build --------- ");
		gb.printMe();
		
		//Scan reads to look for major syntactical errors
		Counts c = null;
		if (threads < 2) c = getCounts(reader, gb, verbose);
		else c = getCounts(reader, gb, verbose, threads);
		
		//Print stats on errors found
		System.out.println();
		System.out.println("===================== RECORDS =====================");
		System.out.println("Total Reads: " + c.getTotal());
		System.out.println("Records containing serious syntax errors: " + c.getSyntaxErrCount());
		System.out.println("Records containing a null SEQ field and non-null QUAL field: " + c.getNullSeqErrCount());
		System.out.println("Records where the QUAL string length doesn't match the SEQ length: " + c.getQualStrLengthErrCount());
		System.out.println("Records containing an illegal custom field type: " + c.getBadCustomFieldErrCount());
		if (c.getBadCustomFieldErrCount() > 0)
		{
			System.out.println("--- Error-Causing Custom Field Type Strings ---");
			List<String> slist = c.getBadCustomFieldTypes();
			for (String s : slist) System.out.println("\t" + s);
		}
		List<String> slist = c.getBadRefContigs();
		long count = c.getInvalidRefContigErrCount();
		System.out.println("Records containing an illegal REF contig: " + count);
		if (count > 0)
		{
			System.out.println("--- Illegal REF Contig Strings ---");
			for (String s : slist) System.out.println("\t" + s);
		}
		slist = c.getBadRNextContigs();
		count = c.getInvalidRNextContigErrCount();
		System.out.println("Records containing an illegal RNEXT contig: " + count);
		if (count > 0)
		{
			System.out.println("--- Illegal RNEXT Contig Strings ---");
			for (String s : slist) System.out.println("\t" + s);
		}
		System.out.println("Records containing an illegal POS value: " + c.getInvalidPosErrCount());
		System.out.println("Records containing an illegal PNEXT value: " + c.getInvalidPNextCount());
		
		slist = c.getAll_RNAME_Strings();
		System.out.println("--- All REF Contig Strings In Readable Records ---");
		for (String s : slist) System.out.println("\t" + s);
		
		slist = c.getAll_RNEXT_Strings();
		System.out.println("--- All RNEXT Contig Strings In Readable Records ---");
		for (String s : slist) System.out.println("\t" + s);
		
	}
	
	public static void runSamScanner(String[] args, GenomeBuild gb, boolean verbose)
	{
		//Genome build may be null. If it is null, genome build is generated from SAM header
		
		String inPath = null;
		int threads = 1;
		
		for (int i = 0; i < args.length; i++)
		{
			String s = args[i];
			if (s.equals(OP_INPUT))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_INPUT + " flag MUST be followed by input SAM path!");
					printUsage();
					System.exit(1);
				}
				inPath = args[i+1];
			}
			else if (s.equals(OP_THREADS))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_THREADS + " flag MUST be followed by postive integer greater than or equal to 1!");
					printUsage();
					System.exit(1);
				}
				String tstr = args[i+1];
				try {threads = Integer.parseInt(tstr);}
				catch(NumberFormatException e) {
					System.err.println("ERROR: " + OP_THREADS + " flag MUST be followed by postive integer greater than or equal to 1!");
					printUsage();
					System.exit(1);
				}
			}
		}
		
		//Check validity of thread argument
		if (threads < 1)
		{
			System.err.println("ERROR: " + OP_THREADS + " flag MUST be followed by postive integer greater than or equal to 1!");
			printUsage();
			System.exit(1);
		}
		
		//Input - default to stdin if there's no path
		BufferedReader br = null;
		if (inPath == null || inPath.isEmpty())
		{
			br = new BufferedReader(new InputStreamReader(System.in));
		}
		else
		{
			try 
			{
				FileReader fr = new FileReader(inPath);
				br = new BufferedReader(fr);
			} 
			catch (FileNotFoundException e)
			{
				System.err.println("ERROR: File " + inPath + " could not be opened!");
				e.printStackTrace();
				System.exit(1);
			}
		}
		
		
		try 
		{
			scanSam(br, gb, verbose, threads);
		} 
		catch (IOException e1) 
		{
			System.err.println("ERROR: Error reading input stream. Exiting...");
			e1.printStackTrace();
			try 
			{
				br.close();
			} 
			catch (IOException e) 
			{
				System.err.println("ERROR: BufferedReader could not be closed... Exiting anyway.");
				e.printStackTrace();
				System.exit(1);
			}
			System.exit(1);
		}
		
		
		try 
		{
			br.close();
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: BufferedReader could not be closed... Exiting anyway.");
			e.printStackTrace();
			System.exit(1);
		}
		
	}

}
