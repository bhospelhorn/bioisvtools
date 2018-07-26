package hospelhornbg_svtools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Deque;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;

import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.GenomeBuild;

public class Stdchrom {
	
	public static final String OP_INPUT = "-i";
	public static final String OP_OUTPUT = "-o";
	public static final String OP_THREADS = "-t";
	
	public static void printUsage()
	{
		System.out.println("--------------------------------------------------------------------------------");
		System.out.println("BioisvTools || stdchrom");
		System.out.println();
		System.out.println("Purpose: For standardizing the chromosome coordinates in a SAM file.");
		System.out.println();
		System.out.println("Input Formats:");
		System.out.println("\tAligned sequence in [SAM] format. Must be [SAM], NOT BAM.");
		System.out.println();
		System.out.println("Output Formats:");
		System.out.println("\tPlain text sequence alignment map [SAM] format. Use samtools to compress to BAM.");
		System.out.println();
		System.out.println("Flags:");
		System.out.println("\t-i\tFILE\t[Required]\t\tInput file path");
		System.out.println("\t-o\tFILE\t[Required]\t\tOutput file path");
		System.out.println("\t-t\tINT\t[Optional]\t\tThreads to use (Default: 1)");
		System.out.println();
		System.out.println("Sample Usage:");
		System.out.println("java -jar bioisvtools.jar stdchrom -g GRCh37 -i aligneroutput.sam -o pipelineinput.sam");
		System.out.println();
		System.out.println("--------------------------------------------------------------------------------");
	}
	
	public static String generateSQBlock(GenomeBuild g)
	{
		String s = "\n";
		List<Contig> clist = g.getChromosomes();
		int sz = clist.size();
		for (int i = 0; i < sz; i++)
		{
			Contig c = clist.get(i);
			s += "@SQ\tSN:" + c.getUDPName() + "\tLN:" + c.getLength();
			if (i < (sz - 1)) s += "\n";
		}
		return s;
	}
	
	public static class ConvertThread extends Thread
	{
		private static final String killline = "@KILL";
		
		private Deque<String> inqueue;
		private SynchronizedWriter writeBuffer;
		private WrappedGenome genome;
		private int counter;
		
		private boolean killme;
		
		public ConvertThread(SynchronizedWriter writer, WrappedGenome g)
		{
			inqueue = new LinkedList<String>();
			writeBuffer = writer;
			killme = false;
			genome = g;
			this.setDaemon(true);
			Random r = new Random();
			this.setName("Bioisvtools_stdchrom_ConvertThread_" + Integer.toHexString(r.nextInt()));
		}
		
		public void run()
		{
			//System.err.println("Thread " + Thread.currentThread().getName() + " starting...");
			while (!killme)
			{
				while (!queueEmpty())
				{
					String line = popQueue();
					if (line.equals(killline))
					{
						killme = true;
						//System.err.println("Thread " + Thread.currentThread().getName() + " kill signal received. Now terminating...");
						return;
					}
					
					String[] fields = line.split("\t");
					if (fields.length < 3)
					{
						System.err.println("ERROR: Input file is not properly formatted!");
						System.err.println("Terminating conversion...");
						System.exit(1);
					}
					String chrom = fields[2];
					if (!chrom.equals("*"))
					{
						Contig c = genome.getContig(chrom);
						/*if (counter < 10)
						{
							System.err.println(Thread.currentThread().getName() + "|| stdchrom.ConvertThread.run || chrom = " + chrom);
							System.err.println(Thread.currentThread().getName() + "|| stdchrom.ConvertThread.run || chrom output = " + c.getUDPName());
						}*/
						if (c == null)
						{
							System.err.println("ERROR: Contig \"" + chrom + "\" was not recognized!");
							System.err.println("Terminating conversion...");
							System.exit(1);
						}
						fields[2] = c.getUDPName();	
					}
					
					if (fields.length >= 7){
						String rnext = fields[6];
						if (!(rnext.equals("=") || rnext.equals("*")))
						{
							Contig cnext = genome.getContig(rnext);
							if (cnext == null)
							{
								System.err.println("ERROR: Contig \"" + rnext + "\" was not recognized!");
								System.err.println("Terminating conversion...");
								System.exit(1);
							}
							fields[6] = cnext.getUDPName();	

							//System.err.println(Thread.currentThread().getName() + "|| stdchrom.ConvertThread.run || chrom = " + chrom);
							//System.err.println(Thread.currentThread().getName() + "|| stdchrom.ConvertThread.run || chrom output = " + c.getUDPName());
							//System.err.println(Thread.currentThread().getName() + "|| stdchrom.ConvertThread.run || chrom saved = " + fields[2]);
							//System.err.println(Thread.currentThread().getName() + "|| stdchrom.ConvertThread.run || rnext = " + rnext);
							//System.err.println(Thread.currentThread().getName() + "|| stdchrom.ConvertThread.run || rnext output = " + cnext.getUDPName());
							//System.err.println(Thread.currentThread().getName() + "|| stdchrom.ConvertThread.run || rnext saved = " + fields[6]);
						}
						else 
						{
							/*if (counter < 10)
							{
								System.err.println(Thread.currentThread().getName() + "|| stdchrom.ConvertThread.run || rnext = " + rnext);
								System.err.println(Thread.currentThread().getName() + "|| stdchrom.ConvertThread.run || No conversion needed.");
							}	*/
						}
					}
					
					
					String s = "";
					for (int i = 0; i < fields.length; i++)
					{
						s += fields[i];
						if (i < (fields.length - 1)) s += "\t";
					}
					try 
					{
						writeBuffer.writeString("\n" + s);
						incrementCounter();
					} 
					catch (IOException e) 
					{
						System.err.println("ERROR: " + Thread.currentThread().getName() + " || I/O error - cannot write to output!");
						System.err.println("Thread " + Thread.currentThread().getName() + " now terminating...");
						e.printStackTrace();
						return;
					}
				}
				
			}
			System.err.println("Thread " + Thread.currentThread().getName() + " killed. Now terminating...");
		}
		
		public void kill()
		{
			killme = true;
		}
		
		public void killWhenDone()
		{
			addToQueue(killline);
		}
		
		public synchronized void pushQueue(String line)
		{
			inqueue.push(line);
		}
		
		public synchronized void addToQueue(String line)
		{
			inqueue.add(line);
		}
		
		private synchronized String popQueue()
		{
			return inqueue.pop();
		}
		
		private synchronized boolean queueEmpty()
		{
			return inqueue.isEmpty();
		}
		
		private synchronized void incrementCounter()
		{
			counter++;
		}
		
		public synchronized int countProcessed()
		{
			return counter;
		}
		
		
		
	}
	
	public static class SynchronizedWriter
	{
		private BufferedWriter writer;
		
		public SynchronizedWriter(FileWriter filewriter)
		{
			writer = new BufferedWriter(filewriter);
		}
		
		public synchronized void writeString(String s) throws IOException
		{
			writer.write(s);
		}
		
		public void close() throws IOException
		{
			writer.close();
		}
		
	}
	
	public static class WrappedGenome
	{
		private GenomeBuild genome;
		
		public WrappedGenome(GenomeBuild g)
		{
			genome = g;
		}
		
		public synchronized Contig getContig(String cName)
		{
			return genome.getContig(cName);
		}
	
	}
	
	public static boolean liveThreads(Thread[] tArr)
	{
		for (int i = 0; i < tArr.length; i++)
		{
			if (tArr[i] != null)
			{
				if (tArr[i].isAlive()) return true;
			}
		}
		return false;
	}
	
	public static void convertSAM(String inpath, String outpath, GenomeBuild g, int threads)
	{
		int cthreads = threads - 1;
		ConvertThread[] converters = new ConvertThread[cthreads];
		
		try
		{
			FileReader fr = new FileReader(inpath);
			BufferedReader br = new BufferedReader(fr);
			
			FileWriter fw = new FileWriter(outpath);
			SynchronizedWriter sw = new SynchronizedWriter(fw);
			
			int counter = 0;
			int tPtr = 0;
			
			WrappedGenome wg = new WrappedGenome(g);
			
			for (int i = 0; i < cthreads; i++)
			{
				converters[i] = new ConvertThread(sw, wg);
				converters[i].start();
			}
			
			String line = null;
			boolean sqfound = false;
			boolean skipnl = true;
			while ((line = br.readLine()) != null)
			{
				if (line.isEmpty()) continue;
				if (line.charAt(0) == '@')
				{
					if (line.startsWith("@SQ"))
					{
						if (!sqfound) sw.writeString(generateSQBlock(g));
						sqfound = true;
					}
					else if (line.startsWith("@HD"))
					{
						if (!skipnl) sw.writeString("\n");
						String[] fields = line.split("\t");
						String outline = "";
						boolean hasso = false;
						for (int i = 0; i < fields.length; i++)
						{
							if (fields[i].startsWith("SO:"))
							{
								fields[i] = "SO:unsorted";
								hasso = true;
							}
							outline += fields[i];
							if (i < fields.length - 1) outline += "\t";
						}
						if (!hasso) outline += "\tSO:unsorted";
						sw.writeString(outline);
					}
					else {
						if (!skipnl) sw.writeString("\n");
						sw.writeString(line);
					}
				}
				else
				{
					converters[tPtr].addToQueue(line);
					tPtr++;
					if (tPtr >= cthreads) tPtr = 0;
				}
				skipnl = false;
				counter++;
				if (counter % 1000000 == 0) {
					System.err.println("Update -- " + counter + " lines read...");
					int pcount = 0;
					for (int i = 0; i < cthreads; i++) pcount += converters[i].countProcessed();
					System.err.println("\t(" + pcount + " records processed)");
					//for (int i = 0; i < cthreads; i++) System.err.println("\tRecords processed by " + converters[i].getName() + " : " + converters[i].countProcessed());
				}
			}
			
			br.close();
			fr.close();
			
			for (int i = 0; i < cthreads; i++) converters[i].killWhenDone();
			while (liveThreads(converters))
			{
				try 
				{
					Thread.sleep(10000);
				} 
				catch (InterruptedException e)
				{
					System.err.println("Main thread wait sleep interrupted. Going back to sleep...");
					e.printStackTrace();
				}
				int pcount = 0;
				for (int i = 0; i < cthreads; i++) pcount += converters[i].countProcessed();
				System.err.println("Update -- " + pcount + " records processed");
			}
			
			int pcount = 0;
			for (int i = 0; i < cthreads; i++) pcount += converters[i].countProcessed();
			System.err.println("Final -- " + pcount + " records processed");
			
			sw.close();
			fw.close();	
		}
		catch (IOException e)
		{
			System.err.println("ERROR: There was an error with file I/O...");
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	public static void convertSAM(String inpath, String outpath, GenomeBuild g)
	{
		try
		{
			FileReader fr = new FileReader(inpath);
			BufferedReader br = new BufferedReader(fr);
			
			FileWriter fw = new FileWriter(outpath);
			BufferedWriter bw = new BufferedWriter(fw);
			
			int counter = 0;
			
			String line = null;
			boolean sqfound = false;
			boolean skipnl = true;
			while ((line = br.readLine()) != null)
			{
				if (line.isEmpty()) continue;
				if (!skipnl) bw.write("\n");
				if (line.charAt(0) == '@')
				{
					if (line.startsWith("@SQ"))
					{
						if (!sqfound) bw.write(generateSQBlock(g));
						sqfound = true;
					}
					else bw.write(line);
				}
				else
				{
					String[] fields = line.split("\t");
					if (fields.length < 3)
					{
						System.err.println("ERROR: Input file is not properly formatted!");
						System.err.println("Terminating conversion...");
						System.exit(1);
					}
					String chrom = fields[2];
					if (!chrom.equals("*"))
					{
						Contig c = g.getContig(chrom);
						if (c == null)
						{
							System.err.println("ERROR: Contig \"" + chrom + "\" was not recognized!");
							System.err.println("Terminating conversion...");
							System.exit(1);
						}
						fields[2] = c.getUDPName();	
					}
					
					if (fields.length >= 7){
						String rnext = fields[6];
						if (!(rnext.equals("=") || rnext.equals("*")))
						{
							Contig cnext = g.getContig(rnext);
							if (cnext == null)
							{
								System.err.println("ERROR: Contig \"" + rnext + "\" was not recognized!");
								System.err.println("Terminating conversion...");
								System.exit(1);
							}
							fields[6] = cnext.getUDPName();	
						}
					}
					
					String s = "";
					for (int i = 0; i < fields.length; i++)
					{
						s += fields[i];
						if (i < (fields.length - 1)) s += "\t";
					}
					bw.write(s);
				}
				skipnl = false;
				counter++;
				if (counter % 1000000 == 0) System.err.println("Update -- " + counter + " lines processed...");
			}
			
			br.close();
			fr.close();
			
			bw.close();
			fw.close();	
		}
		catch (IOException e)
		{
			System.err.println("ERROR: There was an error with file I/O...");
			e.printStackTrace();
			System.exit(1);
		}
		
	}
	
	public static void standardizeChrom(String[] args, GenomeBuild g)
	{
		String inPath = null;
		String outPath = null;
		int threads = 1;
		
		for (int i = 0; i < args.length; i++)
		{
			String s = args[i];
			if (s.equals(OP_INPUT))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_INPUT + " flag MUST be followed by input VCF path!");
					printUsage();
					System.exit(1);
				}
				inPath = args[i+1];
			}
			if (s.equals(OP_OUTPUT))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_OUTPUT + " flag MUST be followed by output VCF path!");
					printUsage();
					System.exit(1);
				}
				outPath = args[i+1];
			}
			if (s.equals(OP_THREADS))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_THREADS + " flag MUST be followed by number of requested threads!");
					printUsage();
					System.exit(1);
				}
				try 
				{
					threads = Integer.parseInt(args[i+1]);
				}
				catch (NumberFormatException e)
				{
					System.err.println("ERROR: " + OP_THREADS + " flag MUST be followed by a valid integer!");
					printUsage();
					System.exit(1);
				}
			}
		}
		
		
		boolean pass = true;
		if (inPath == null || inPath.isEmpty())
		{
			System.err.println("ERROR: Input path is required!");
			pass = false;
		}
		if (outPath == null || outPath.isEmpty())
		{
			System.err.println("ERROR: Output path is required!");
			pass = false;
		}
		if (!pass)
		{
			printUsage();
			System.exit(1);
		}
		
		if (g == null)
		{
			System.err.println("ERROR: Reference genome build is required!");
			printUsage();
			System.exit(1);
		}
		System.err.println("stdchrom || Genome Selected: " + g.getBuildName());
		
		if (threads < 2) convertSAM(inPath, outPath, g);
		else convertSAM(inPath, outPath, g, threads);
		
	}

}
