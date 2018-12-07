package hospelhornbg_svtools;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.Gene;
import hospelhornbg_genomeBuild.GeneSet;
import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.FileBuffer;

public class LimsReportGenes {
	
	public static final String OP_INPUT = "-i"; 
	public static final String OP_GBUILD = "-g";
	
	public static final int CHR_FIELD = 2;
	public static final int START_FIELD = 3;
	public static final int END_FIELD = 4;
	public static final int GENE_FIELD = 6;
	
	public static final int MIN_FIELDS = 8;

	public static void main(String[] args) {
		//All this does is fetch the genes for a certain region reading it out of a LIMS csv file.
		
		String inPath = null;
		String gbname = null;
		GenomeBuild gb = null;
		GeneSet genes = null;
		
		for (int i = 0; i < args.length; i++)
		{
			String s = args[i];
			if (s.equals(OP_INPUT))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_INPUT + " flag MUST be followed by input csv path!");
					System.exit(1);
				}
				inPath = args[i+1];
			}
			else if (s.equals(OP_GBUILD))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_GBUILD + " flag MUST be followed by genome build name!");
					System.exit(1);
				}
				gbname = args[i+1];
			}
		}
		
		//Check input path validity
		if(inPath == null || inPath.isEmpty())
		{
			System.err.println("ERROR: " + OP_INPUT + " flag followed by path to LIMS report csv is required!");
			System.exit(1);
		}
		if(!FileBuffer.fileExists(inPath))
		{
			System.err.println("ERROR: File \"" + inPath + "\" does not appear to exist!");
			System.exit(1);
		}
		
		//Check genome build validity
		String homedir = ConsoleMain.getUserHome(true);
		gb = ConsoleMain.loadBuild(homedir, gbname, true);
		if(gb == null)
		{
			System.err.println("ERROR: Genome build \"" + gbname + "\" was not recognized.");
			System.exit(1);
		}
		genes = ConsoleMain.loadRefSeq(homedir, gb, true);
		if(genes == null)
		{
			System.err.println("ERROR: RefSeq for \"" + gbname + "\" could not be loaded!");
			System.exit(1);
		}
		
		//Open input file
		try
		{
			FileReader fr = new FileReader(inPath);
			BufferedReader br = new BufferedReader(fr);
			
			String line = br.readLine();
			System.out.println(line); //First line just gets printed back out (header)
			while ((line = br.readLine()) != null)
			{
				String[] fields = line.split(","); 
				//WARNING: This DOES split fields that have internal commas, but since those come
					// after the fields we're interested in, I'm going to ignore that for now.
				if (fields.length < MIN_FIELDS)
				{
					System.out.println(line);
					continue;
				}
				//Get chrom
				Contig c = gb.getContig(fields[CHR_FIELD]);
				if (c == null)
				{
					System.out.println(line);
					continue;
				}
				//Get start
				int st = -1;
				try {st = Integer.parseInt(fields[START_FIELD]);}
				catch(NumberFormatException e)
				{
					System.out.println(line);
					continue;
				}
				//Get end
				int ed = -1;
				try {ed = Integer.parseInt(fields[END_FIELD]);}
				catch(NumberFormatException e)
				{
					System.out.println(line);
					continue;
				}
				
				//Get genes
				List<Gene> glist = genes.getGenesInRegion(c, st, ed);
				
				//If list is non-null and not empty, populate genes field
				if (glist != null && !glist.isEmpty())
				{
					//Print fields before genes field...
					for (int i = 0; i < GENE_FIELD; i++) System.out.print(fields[i] + ",");
					
					//Remove redundant names (same gene, different transcripts)
					List<String> gnames = new LinkedList<String>();
					for(Gene g : glist)
					{
						String gname = g.getName();
						if(!gnames.contains(gname)) gnames.add(gname);
					}
					
					//Generate field string
					boolean first = true;
					for (String g : gnames)
					{
						if(!first) System.out.print(";");
						first = false;
						System.out.print(g);
					}
					

					//Print fields after genes field...
					for (int i = GENE_FIELD + 1; i < fields.length; i++) System.out.print("," + fields[i]);
					System.out.print("\n");
				}
				else System.out.println(line);
			}
			
			br.close();
		}
		catch(IOException e)
		{
			System.err.println("ERROR: Error reading input stream!");
			System.exit(1);
		}
		
	}

}
