package hospelhornbg_svtools;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Deque;
import java.util.LinkedList;
import java.util.List;

import waffleoRai_Utils.TallyMap;

public class VarSizes {
	
	public static final String OP_VCFIN = "-i"; 
	public static final String OP_THREADS = "-t"; 
	
	public static void printUsage()
	{
		System.out.println("--------------------------------------------------------------------------------");
		System.out.println("BioisvTools || varsztally");
		System.out.println();
		System.out.println("Purpose: For determining the size of each variant in a vcf");
		System.out.println("and printing the tallies to stdout in tsv format.");
		System.out.println();
		System.out.println("Input Formats:");
		System.out.println("\tInput callset must be in [vcf] format");
		System.out.println();
		System.out.println("Output Formats:");
		System.out.println("\tTab delimited table [tsv]");
		System.out.println();
		System.out.println("Flags:");
		System.out.println("\t-i\tFILE\t[Required]\t\tInput vcf path.");
		System.out.println("\t-t\tINT\t[Optional]\t\tNumber of threads. Defaults to 1.");
		System.out.println();
		System.out.println("Sample Usage:");
		System.out.println("java -jar bioisvtools.jar varsztally -g GRCh37 -v -n NA12878 -i NA12878_vars.vcf");
		System.out.println("java -jar bioisvtools.jar varsztally -g hg38 -n NA12878 -i NA12878_vars.vcf -t 24");
		System.out.println();
		System.out.println("--------------------------------------------------------------------------------");
	}
	
	private static class ParserThread extends Thread
	{
		private static final String killsignal = "#KILL";
		
		private Deque<String> lineQueue;
		private TallyMap mapRef;
		
		public ParserThread(TallyMap map, int number)
		{
			lineQueue = new LinkedList<String>();
			this.setName("VarSizes_ParserThread_" + number);
			this.setDaemon(true);
			mapRef = map;
		}
		
		public void run()
		{
			boolean killMe = false;
			while (!killMe)
			{
				if (!emptyQueue())
				{
					String line = popLine();
					if (line != null && !line.isEmpty())
					{
						if (line.equals(killsignal)){
							killMe = true;
							break;
						}
						VarSizes.parseLine(line, mapRef);	
					}
				}
				else
				{
					try 
					{
						Thread.sleep(100);
					}
					catch (InterruptedException e) 
					{
						Thread.interrupted();
						//e.printStackTrace();
					}	
				}
			}
		}
		
		public void killWhenDone()
		{
			addLine(killsignal);
		}
		
		public synchronized void addLine(String line)
		{
			lineQueue.addLast(line);
		}
		
		public synchronized boolean emptyQueue()
		{
			return lineQueue.isEmpty();
		}
		
		public synchronized String popLine()
		{
			return lineQueue.pop();
		}
		
		
	}
	
	public static void parseLine(String line, TallyMap map)
	{
		if (line != null && !line.isEmpty())
		{
			if (line.charAt(0) != '#')
			{
				String[] fields = line.split("\t");
				if (fields.length < 8) return;
				int sv = -1;
				String[] infofields = null;
				if (!fields[7].isEmpty() && !fields[7].equals("."))
				{
					infofields = fields[7].split(";");
					for (int i = 0; i < infofields.length; i++)
					{
						if (infofields[i] != null && infofields[i].startsWith("SVTYPE="))
						{
							String type = infofields[i].substring(infofields[i].indexOf('='));
							if (type != null && !type.isEmpty())
							{
								if (type.equals("BND") || type.equals("TRA"))
								{
									//DO NOT LOG THIS VARIANT SIZE, IT DOESN'T REALLY HAVE ONE!
									return;
								}
							}
						}
						if (infofields[i] != null && infofields[i].startsWith("SVLEN="))
						{
							sv = i;
							break;
						}
					}
				}
				
				if (sv < 0)
				{
					String refAllele = fields[3];
					String[] altArr = fields[4].split(",");
					//Length is the string length difference between ref and alt(s)
					for (int i = 0; i < altArr.length; i++)
					{
						int len = Math.abs(altArr[i].length() - refAllele.length());
						map.increment(len);
					}
				}
				else
				{
					//Length is abs of SVLEN
					String SVLEN = infofields[sv];
					SVLEN = SVLEN.substring(SVLEN.indexOf('=') + 1);
					String[] svlens = SVLEN.split(",");
					for (int i = 0; i < svlens.length; i++)
					{
						try
						{
							int len = Math.abs(Integer.parseInt(svlens[i]));
							map.increment(len);
						}
						catch (NumberFormatException e)
						{
							//Skip this one, I guess?
							continue;
						}
					}
				}
				
			}	
		}
		
	}
	
	public static boolean anyAlive(Thread[] tarr)
	{
		if (tarr == null) return false;
		if (tarr.length < 1) return false;
		
		for (int i = 0; i < tarr.length; i++)
		{
			if (tarr[i] != null)
			{
				if (tarr[i].isAlive()) return true;
			}
		}
		return false;
	}
	
	public static TallyMap tallySizes(String vcfpath, int threads)
	{
		//Stream variants one by one...
		TallyMap myMap = new TallyMap();
		
		try 
		{
			FileReader fr = new FileReader(vcfpath);
			BufferedReader br = new BufferedReader(fr);
			
			if (threads > 1)
			{
				//int pthreads = threads - 1;
				ParserThread[] tarr = new ParserThread[threads - 1];
				int ind = 0;
				for (int i = 0; i < tarr.length; i++)
				{
					tarr[i] = new ParserThread(myMap, i);
					tarr[i].start();
				}
				String line = null;
				while((line = br.readLine()) != null)
				{
					tarr[ind].addLine(line);
					ind++;
					if (ind >= tarr.length) ind = 0;
				}
				
				//Request termination of other threads, and wait for them to terminate...
				for (int i = 0; i < tarr.length; i++)
				{
					tarr[i].killWhenDone();
				}
				while(anyAlive(tarr))
				{
					try 
					{
						Thread.sleep(1000);
					} 
					catch (InterruptedException e) 
					{
						System.err.println("Main thread wait sleep interrupted. Checking parser thread states...");
						Thread.interrupted();
						//e.printStackTrace();
					}
				}
				
			}
			else if (threads == 1)
			{
				String line = null;
				while((line = br.readLine()) != null)
				{
					parseLine(line, myMap);
				}
			}
			
			br.close();
			fr.close();
		} 
		catch (IOException e) 
		{
			System.err.println("VarSizes.tallySizes || File " + vcfpath + " could not be read!");
			e.printStackTrace();
		}
		
		return myMap;
	}
	
	public static void varsizes(String[] args)
	{
		String inFile = null;
		int threads = 1;
		
		for (int i = 0; i < args.length; i++)
		{
			String s = args[i];
			if (s.equals(OP_VCFIN))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_VCFIN + " flag MUST be followed by input VCF path!");
					printUsage();
					System.exit(1);
				}
				inFile = args[i+1];
			}
			else if (s.equals(OP_THREADS))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_THREADS + " flag MUST be followed by a valid integer!");
					printUsage();
					System.exit(1);
				}
				String str = args[i+1];
				try
				{
					threads = Integer.parseInt(str);
				}
				catch(NumberFormatException e)
				{
					System.err.println("ERROR: " + OP_THREADS + " flag MUST be followed by a valid integer!");
					printUsage();
					System.exit(1);
				}
			}
		}
		
		//
		if (inFile == null || inFile.isEmpty())
		{
			System.err.println("ERROR: Input path is required!");
			printUsage();
			System.exit(1);
		}
		
		//
		TallyMap myMap = tallySizes(inFile, threads);
		
		//Print results to stdout in tsv format
		System.out.println("SIZE\tHITS");
		List<Integer> vals = myMap.getAllValues();
		
		for(Integer i : vals)
		{
			System.out.println(i + "\t" + myMap.getCount(i));
		}
		
	}

}
