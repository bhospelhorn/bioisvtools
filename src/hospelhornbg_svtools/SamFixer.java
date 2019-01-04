package hospelhornbg_svtools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentLinkedQueue;

import hospelhornbg_bioinformatics.SAMHeaderLine;
import hospelhornbg_bioinformatics.SAMRecord;
import hospelhornbg_bioinformatics.SAMRecord.InvalidSAMRecordException;
import hospelhornbg_bioinformatics.SAMRecord.WarningFlags;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class SamFixer {
	
	public static final String OP_INPUT = "-i"; 
	public static final String OP_OUTPUT = "-o"; 

	public static final String OP_THREADS = "-t"; 
	public static final String OP_UCSC_NAMES = "-u"; 
	
	public static final String OP_SAMPLE_NAME = "-s"; 
	
	public static final int MAX_QUEUESIZE = 1000000000;
	
	public static void printUsage()
	{
		System.out.println("--------------------------------------------------------------------------------------------------------------");
		System.out.println("BioisvTools || SAM Fixer");
		System.out.println();
		System.out.println("Purpose: Outputs an unsorted version of the input SAM with contig names adjusted");
		System.out.println("to match requested build.");
		System.out.println();
		System.out.println("Input Formats:");
		System.out.println("\tInput stream must be a SAM file. BAM cannot be read at this time.");
		System.out.println();
		System.out.println("Flags:");
		System.out.println("\t-i\tFILE\t[Optional]\t\tInput file path. Defaults to stdin stream if not provided.");
		System.out.println("\t-o\tFILE\t[Optional]\t\tOutput file path. Defaults to stdout stream if not provided.");
		System.out.println("\t-t\tINT\t[Optional]\t\tNumber of threads to use.");
		System.out.println("\t-u\tFLAG\t[Optional]\t\tUse UCSC standard names (eg. chr1, chrX, etc.) instead of NCBI (eg. 1, X, etc.)");
		System.out.println("\t-s\tFILE\t[Optional]\t\tName of sample the SAM holds data for in case an RG line is not present.");
		System.out.println();
		System.out.println("Note:");
		System.out.println("You maust specify a genome build.");
		System.out.println();
		System.out.println("Sample Usage:");
		System.out.println("java -jar bioisvtools.jar fixsam -g hg38");
		System.out.println("java -jar bioisvtools.jar fixsam -g ncbi36 -v");
		System.out.println("java -jar bioisvtools.jar fixsam -g grch37 -v -t 24 -u");
		System.out.println("java -jar bioisvtools.jar fixsam -g grch38 -v -i sample.sam -o samplechromfix.sam -t 8 -s MySample");
		System.out.println();
		System.out.println("--------------------------------------------------------------------------------------------------------------");
	}
	
	private static class BadCounter
	{
		private long bad_RNAME;
		private long bad_RNEXT;
		
		private long bad_format;
		private long bad_record;
		
		public BadCounter()
		{
			bad_RNAME = 0;
			bad_RNEXT = 0;
			bad_format = 0;
			bad_record = 0;
		}
		
		public synchronized void increment_RNAME()
		{
			bad_RNAME++;
		}
		
		public synchronized void increment_RNEXT()
		{
			bad_RNEXT++;
		}
		
		public synchronized long get_RNAME()
		{
			return bad_RNAME;
		}
		
		public synchronized long get_RNEXT()
		{
			return bad_RNEXT;
		}
	
		public synchronized void increment_BadFormat()
		{
			bad_format++;
		}
		
		public synchronized void increment_BadRecord()
		{
			bad_record++;
		}
		
		public synchronized long get_BadFormat()
		{
			return bad_format;
		}
		
		public synchronized long get_BadRecord()
		{
			return bad_record;
		}
	
		
	
	}
	
	private static class SyncedInt
	{
		private volatile int i;
		
		public SyncedInt()
		{
			i = 0;
		}
		
		public synchronized void increment()
		{
			i++;
		}
		
		public int get()
		{
			return i;
		}
	}
	
	private static List<SAMHeaderLine> parseHeader(List<String> rawlines) throws UnsupportedFileTypeException
	{
		List<SAMHeaderLine> headerlist = new LinkedList<SAMHeaderLine>();
		for(String line : rawlines)
		{
			headerlist.add(new SAMHeaderLine(line));
		}
		
		return headerlist;
	}

	private static List<SAMHeaderLine> fixHeader(List<SAMHeaderLine> header, GenomeBuild gb, SAMHeaderLine newRG, boolean ucsc)
	{
		//Sort by key...
		Map<String, List<SAMHeaderLine>> keymap = new HashMap<String, List<SAMHeaderLine>>();
		for(SAMHeaderLine l : header)
		{
			List<SAMHeaderLine> list = keymap.get(l.getKey());
			if (list == null)
			{
				list = new LinkedList<SAMHeaderLine>();
				keymap.put(l.getKey(), list);
			}
			list.add(l);
		}
		
		List<SAMHeaderLine> outlist = new LinkedList<SAMHeaderLine>();
		
		//Look for @HD, make sure sort order is "unsorted". All but the first @HD line should be tossed.
		List<SAMHeaderLine> allhd = keymap.get("HD");
		SAMHeaderLine hdline = null;
		if(allhd != null && !allhd.isEmpty()) {
			hdline = allhd.get(0);
			hdline.setField("SO", "unsorted");
		}
		else {
			hdline = new SAMHeaderLine();
			hdline.setKey("HD");
			hdline.addField("VN", "1.5");
			hdline.addField("SO", "unsorted");
		}
		outlist.add(hdline);
		
		//Look for @SQ
			//Strip contigs in genome build that are not in the SQ listings
			//Replace all SQ lines with those that remain in genome build
		List<SAMHeaderLine> allsq = keymap.get("SQ");
		List<Contig> allchrom = gb.getChromosomes();
		Set<Contig> foundchrom = new HashSet<Contig>();
		for (SAMHeaderLine sq : allsq)
		{
			Contig c = gb.getContig(sq.getValue("SN"));
			if (c == null)
			{
				System.err.println("SamFixer.fixHeader || SQ contig found in SAM does not match any build contigs: " + sq.getValue("SN"));
				System.err.println("\tAll reads mapped to this contig will be reset to unmapped.");
				continue;
			}
			foundchrom.add(c);
		}
		for (Contig c : allchrom)
		{
			if (!foundchrom.contains(c))
			{
				//Remove it
				System.err.println("SamFixer.fixHeader || Build contig " + c.getUDPName() + " not found in SAM header.");
				System.err.println("\tThis contig will be henceforth ignored.");
				gb.removeContig(c.getUDPName());
			}
		}
		allchrom = gb.getChromosomes(); //Re-get chrom list
		for (Contig c : allchrom)
		{
			//Generate the new SQ lines
			SAMHeaderLine sq = new SAMHeaderLine();
			sq.setKey("SQ");
			if(!ucsc) sq.addField("SN", c.getUDPName());
			else sq.addField("SN", c.getUCSCName());
			sq.addField("LN", Long.toUnsignedString(c.getLength()));
			sq.addField("AS", gb.getBuildName());
			outlist.add(sq);
		}
		
		//Add @RG if not there
		List<SAMHeaderLine> allrg = keymap.get("RG");
		if(allrg == null || allrg.isEmpty()) outlist.add(newRG);
		else
		{
			for (SAMHeaderLine rg : allrg) outlist.add(rg);
		}
		
		//Add back everything else
		for(SAMHeaderLine hl : header)
		{
			if (hl.getKey().equals("HD")) continue;
			if (hl.getKey().equals("SQ")) continue;
			if (hl.getKey().equals("RG")) continue;
			outlist.add(hl); //Crude, but easiest way to keep in order
		}
		
		return outlist;
	}
	
	private static void generateOutputLine(String input, GenomeBuild gb, BadCounter counter, ConcurrentLinkedQueue<String> writequeue, boolean verbose, boolean ucsc) throws InterruptedException
	{
		//Tosses lines with bad contigs (commented out lines for unmapping bad contig lines)
		try 
		{
			SAMRecord sr = SAMRecord.parseSAMRecord(input, gb, verbose).getRecord();
			WarningFlags wf = sr.getParserWarnings();
			if (wf.err_invalid_RNAME != null) {
				counter.increment_RNAME();
				//sr.flagSegmentUnmapped(true);
				return;
			}
			if (wf.err_invalid_RNEXT != null) {
				counter.increment_RNEXT();
				//sr.flagNextSegmentUnmapped(true);
				return;
			}
			String outline = sr.writeSAMRecord(ucsc);
			//Wait until queue has space...
			while(writequeue.size() >= SamFixer.MAX_QUEUESIZE)
			{
				try 
				{
					Thread.sleep(10);
				} 
				catch (InterruptedException e) 
				{
					//If it's getting interrupted, that should only be a kill signal...
					//Dump the record in the queue regardless of queue size, and throw the exception.
					writequeue.add(outline);
					throw e;
				}
			}
			writequeue.add(outline);
		}
		catch (UnsupportedFileTypeException e) 
		{
			if(verbose) System.err.println(Thread.currentThread().getName() + " || SamFixer.generateOutputLine || Record could not be parsed!");
			counter.increment_BadFormat();
		} 
		catch (InvalidSAMRecordException e) 
		{
			if(verbose) System.err.println(Thread.currentThread().getName() + " || SamFixer.generateOutputLine || Record could not be parsed!");
			counter.increment_BadRecord();
		}
	}
	
	private static int finishParsingQueue(ConcurrentLinkedQueue<String> linequeue, GenomeBuild gb, BadCounter counter, ConcurrentLinkedQueue<String> writequeue, boolean verbose, boolean ucsc)
	{
		int lcount = 0;
		while(!linequeue.isEmpty())
		{
			String line = linequeue.poll();
			if (line == null) continue;
			try 
			{
				generateOutputLine(line, gb, counter, writequeue, verbose, ucsc);
			} 
			catch (InterruptedException e) 
			{
				//IGNORE!!! Make sure this method is used KNOWING that it ignores interrupt signals!!
				//e.printStackTrace();
			}
			lcount++;
		}
		return lcount;
	}
	
	public static BadCounter fixSam(BufferedReader input, BufferedWriter output, GenomeBuild gb, SAMHeaderLine newRG, int threads, boolean verbose, boolean ucsc) throws IOException, UnsupportedFileTypeException
	{
		System.err.println("DEBUG: fixSam method entered!");
		//Prepare Counter
		BadCounter counter = new BadCounter();
		
		//Prepare and start record parsing/fixing threads
		ConcurrentLinkedQueue<String> linequeue = new ConcurrentLinkedQueue<String>();
		ConcurrentLinkedQueue<String> writequeue = new ConcurrentLinkedQueue<String>();
		
		int tcount = threads - 2;
		Thread[] tarr = new Thread[tcount];
		SyncedInt donecount = new SyncedInt();
		for (int i = 0; i < tcount; i++)
		{
			tarr[i] = new Thread(new Runnable() {
				@Override
				public void run() 
				{
					System.err.println("DEBUG: Worker thread " + Thread.currentThread().getName() + " started!");
					long pcount = 0;
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
								donecount.increment();
								return;
							}
						}
						else
						{
							String line = linequeue.poll();
							if (line == null) continue;
							pcount++;
							
							//Try to parse
							try 
							{
								generateOutputLine(line, gb, counter, writequeue, verbose, ucsc);
							} 
							catch (InterruptedException e) 
							{
								//Interruption should only be a kill signal.
								//Run until read queue is empty, eating any more interruptions, and return.
								pcount += finishParsingQueue(linequeue, gb, counter, writequeue, verbose, ucsc);
								donecount.increment();
								if(verbose)System.err.println("DEBUG: Lines Parsed " + pcount + " (Worker: " + Thread.currentThread().getName() + ") - Worker thread returning... (From write queue full sleep)");
								return;
							}
							
						}		
					}
					//Run until the queue is empty...
					pcount += finishParsingQueue(linequeue, gb, counter, writequeue, verbose, ucsc);
					donecount.increment();
					if(verbose)System.err.println("DEBUG: Lines Parsed " + pcount + " (Worker: " + Thread.currentThread().getName() + ") - Worker thread returning...");
				}
			});
			tarr[i].setName("SamFixerWorkerThread" + i);
		}
		
		//Prepare writer thread
		Thread writerthread = new Thread(new Runnable() {

			@Override
			public void run() {
				System.err.println("DEBUG: Writer thread " + Thread.currentThread().getName() + " started!");
				long lcount = 0;
				while(!Thread.interrupted())
				{
					
					if(writequeue.isEmpty()) 
					{
						try 
						{
							Thread.sleep(10);
						} 
						catch (InterruptedException e) 
						{
							//If it's getting interrupted, that should only be a kill signal...
							break;
						}
					}
					else
					{
						//Write the line to the output stream
						String line = writequeue.poll();
						if (line == null) continue;
						try 
						{
							//output.write("\n" + line);
							//output.write(line + "\n");
							if (lcount != 0) output.write("\n" + line);
							else output.write(line);
							lcount++;
							//if(lcount > 5942000000L && lcount < 5943000000L) System.err.println("R" + lcount + "\t" + line); //DEBUG
						} 
						catch (IOException e) 
						{
							if(verbose)System.err.println("ERROR: (" + Thread.currentThread().getName() + ") Fatal IO Exception!");
							if(verbose)System.err.println("DEBUG: Lines Written: " + lcount + " (" + Thread.currentThread().getName() + ") - Writer thread returning... (IO Fail)");
							e.printStackTrace();
							return;
						}
					}
					if (verbose && (lcount % 100000000L == 0)) System.err.println("DEBUG: Lines Written: " + lcount);
				}
				//Finish writing any lines still in queue or that are waiting to be processed...
				while(!writequeue.isEmpty() || donecount.get() < tcount)
				{
					String line = writequeue.poll();
					if (line == null) continue;
					try 
					{
						output.write(line);
						lcount++;
					} 
					catch (IOException e) 
					{
						if(verbose)System.err.println("ERROR: (" + Thread.currentThread().getName() + ") Fatal IO Exception!");
						if(verbose)System.err.println("DEBUG: Lines Written: " + lcount + " (" + Thread.currentThread().getName() + ") - Writer thread returning... (IO Fail)");
						e.printStackTrace();
						return;
					}
				}
				if(verbose)System.err.println("DEBUG: Lines Written: " + lcount + " (" + Thread.currentThread().getName() + ") - Writer thread returning...");	
			}
			
		});
		writerthread.setName("SamFixerWriterThread");
		
		//Read and re-write header
		//Remove any contigs from the genome build not in header
		List<String> rawheader = new LinkedList<String>();
		String line = null;
		while((line = input.readLine()) != null)
		{
			if (line.isEmpty() || line.charAt(0) != '@') break; //We'll save that line for records in a minute..
			//Otherwise...
			rawheader.add(line);
		}
		
		//Interpret header lines
		List<SAMHeaderLine> hlist = parseHeader(rawheader);
		
		//Modify header/ genome build
		List<SAMHeaderLine> nhlist = fixHeader(hlist, gb, newRG, ucsc);
		
		//Output new header
		//boolean first = true;
		int hlcount = 0;
		for(SAMHeaderLine hl : nhlist)
		{
			//if (!first) output.write("\n");
			//first = false;
			output.write(hl.toString() + "\n");
			hlcount++;
		}
		if(verbose) System.err.println("DEBUG: New header written (" + hlcount + " lines)... Starting record processing...");
		
		//Now, the records!
		//Start the workers and writer
		for (int i = 0; i < tcount; i++) tarr[i].start();
		writerthread.start();
		linequeue.add(line); //Dump the first non-header line we read before...
		long readcount = 1;
		while((line = input.readLine()) != null)
		{
			//Check line queue size and wait if full...
			while(linequeue.size() >= SamFixer.MAX_QUEUESIZE)
			{
				try 
				{
					Thread.sleep(10);
				} 
				catch (InterruptedException e) 
				{
					System.err.println("SamFixer.fixSam || Unexpected interruption to reader thread. Ignoring...");
				}
			}
			linequeue.add(line);
			readcount++;
			if ((verbose) && (readcount % 10000000L == 0)) System.err.println("DEBUG: " + readcount + " lines read!");
		}
		
		if(verbose) System.err.println("DEBUG: Read complete!");
		//Now, the reader is done. Interrupt the workers and writer...
		for (int i = 0; i < tcount; i++) {
			synchronized(tarr[i]) {tarr[i].interrupt();}
		}
		synchronized(writerthread) {writerthread.interrupt();}
		
		//Wait for the other threads to die...
		boolean anyalive = true;
		while(anyalive)
		{
			try {
				Thread.sleep(10);
			} catch (InterruptedException e) {
				System.err.println("SamFixer.fixSam || Unexpected interruption to reader thread. Ignoring...");
			}
			
			anyalive = false;
			if (writerthread.isAlive()) {
				anyalive = true;
				continue;
			}
			for (int i = 0; i < tcount; i++) {
				if(tarr[i].isAlive())
				{
					anyalive = true;
					break;
				}
			}
		}
		
		return counter;
	}

	public static void runSamFixer(String[] args, GenomeBuild gb, boolean verbose)
	{
		if (gb == null)
		{
			System.err.println("ERROR: Genome Build must be non-null!");
			printUsage();
			System.exit(1);
		}
		
		String inPath = null;
		String outPath = null;
		int threads = 3;
		boolean ucsc = false;
		String samplename = "MySAM";
		
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
			else if (s.equals(OP_OUTPUT))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_OUTPUT + " flag MUST be followed by output SAM path!");
					printUsage();
					System.exit(1);
				}
				outPath = args[i+1];
			}
			else if (s.equals(OP_SAMPLE_NAME))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_SAMPLE_NAME + " flag MUST be followed by sample name string!");
					printUsage();
					System.exit(1);
				}
				samplename = args[i+1];
			}
			else if (s.equals(OP_UCSC_NAMES))
			{
				ucsc = true;
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
		
		//Check threads
		if (threads < 1)
		{
			System.err.println("ERROR: " + OP_THREADS + " flag MUST be followed by postive integer greater than or equal to 1!");
			printUsage();
			System.exit(1);
		}
		if (threads < 3) threads = 3; //Because I'm lazy. The 1 and 2 are a lie. Sorry.
		
		//Generate an RG line in case needed
		SAMHeaderLine rgline = new SAMHeaderLine();
		rgline.setKey("RG");
		rgline.addField("ID", "1");
		rgline.addField("SM", samplename);
		

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
		
		//Output - default to stdout if there's no path
		BufferedWriter bw = null;
		if (outPath == null || outPath.isEmpty())
		{
			bw = new BufferedWriter(new OutputStreamWriter(System.out));
		}
		else
		{
			try 
			{
				FileWriter fw = new FileWriter(outPath);
				bw = new BufferedWriter(fw);
			} 
			catch (IOException e)
			{
				System.err.println("ERROR: File " + outPath + " could not be opened!");
				e.printStackTrace();
				System.exit(1);
			}
		}
		
		//Do the thing
		BadCounter counter = null;
		try 
		{
			counter = fixSam(br, bw, gb, rgline, threads, verbose, ucsc);
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: There was an error reading or writing!");
			e.printStackTrace();
			System.exit(1);
		} 
		catch (UnsupportedFileTypeException e) 
		{
			System.err.println("ERROR: There was an error parsing the input!");
			e.printStackTrace();
			System.exit(1);
		}
		
		//Close streams
		try 
		{
			br.close();
			bw.flush();
			bw.close();
		}
		catch (IOException e) 
		{
			System.err.println("ERROR: There was an error closing one or more streams!");
			e.printStackTrace();
			System.exit(1);
		} 
		
		//Print some info
		System.err.println("============== SAM Scan Successful ==============");
		System.err.println("Unreadable Records: " + counter.get_BadFormat());
		System.err.println("Invalid Records: " + counter.get_BadRecord());
		System.err.println("Records with Illegal RNAME: " + counter.get_RNAME());
		System.err.println("Records with Illegal RNEXT: " + counter.get_RNEXT());
		
		
	}
	
	
}
