package hospelhornbg_svtools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_genomeBuild.Gene;
import hospelhornbg_genomeBuild.GeneFunc;
import hospelhornbg_genomeBuild.GeneSet;
import hospelhornbg_genomeBuild.GenomeBuild;
import hospelhornbg_genomeBuild.OMIMGeneMapImporter;

public class GeneTally {
	
	public static final String OP_INFILE = "-n"; 
	
	public static final String OP_OMIM_TABLE = "-a"; 
	public static final String OP_TYPELIST = "-t"; 
	public static final String OP_INTRONIC = "-i"; 
	public static final String OP_INTERGENIC = "-I"; 
	public static final String OP_EXONIC_ONLY = "-e"; 
	
	public static void printUsage()
	{
		System.out.println("--------------------------------------------------------------------------------");
		System.out.println("BioisvTools || svght");
		System.out.println();
		System.out.println("Purpose: For analyzing gene lists output by svreport and counting the number of individuals ");
		System.out.println("there is a hit for each gene in.");
		System.out.println();
		System.out.println("Input Formats:");
		System.out.println("\tInput table is a simple tsv (tab delimited) two-column list of samples and svreport output directories");
		System.out.println();
		System.out.println("Output Formats:");
		System.out.println("\tTab-separated table [tsv] of genes found, number of hits, and list of samples hit was found in. Printed to stdout.");
		System.out.println();
		System.out.println("Flags:");
		System.out.println("\t-n\tFILE\t[Required]\t\tInput table path.");
		System.out.println("\t-t\tSTRING\t[Optional]\t\tComma delimited list of SV types to include. (Default: DEL,DUP)");
		System.out.println("\t-a\tFILE\t[Optional]\t\tOMIM genemap2.txt file for OMIM annotation.");
		System.out.println("\t-i\tFLAG\t[Optional]\t\tInclude intronic hits.");
		System.out.println("\t-I\tFLAG\t[Optional]\t\tInclude downstream and intergenic hits.");
		System.out.println("\t-e\tFLAG\t[Optional]\t\tInclude only exonic hits.");
		System.out.println();
		System.out.println("Sample Usage:");
		System.out.println("java -jar bioisvtools.jar svght -n svreport_paths.tsv");
		System.out.println("java -jar bioisvtools.jar svght -v -n svreport_paths.tsv -t DEL,INV -i");
		System.out.println();
		System.out.println("--------------------------------------------------------------------------------");
	}
	
	
	private static class Hit
	{
		//public String samplename;
		public boolean seg;
		public SVType type;
		public GeneFunc poseff;
		
	}
	
	private static class GeneHit implements Comparable<GeneHit>
	{
		public String gene;
		public int hits;

		public int compareTo(GeneHit o) 
		{
			if (o == null) return 1;
			if (this.hits != o.hits) return this.hits - o.hits;
			
			return this.gene.compareTo(o.gene);
		}
		
		
	}
	
	public static Map<String, Map<String, Hit>> generateHitMap(Map<String,String> pathmap, Collection<SVType> types, boolean includeIntronic, boolean includeIntergenic, boolean exonicOnly)
	{
		Map<String, Map<String, Hit>> hitmap = new HashMap<String, Map<String, Hit>>();
		
		List<GeneFunc> pelist = new LinkedList<GeneFunc>();
		pelist.add(GeneFunc.EXONIC);
		if (!exonicOnly)
		{
			pelist.add(GeneFunc.SPLICING);
			pelist.add(GeneFunc.NCRNA);
			pelist.add(GeneFunc.UTR5);
			pelist.add(GeneFunc.UTR3);
			pelist.add(GeneFunc.UPSTREAM);
			if (includeIntronic) pelist.add(GeneFunc.INTRONIC);
			if (includeIntergenic) {
				pelist.add(GeneFunc.DOWNSTREAM);
				pelist.add(GeneFunc.INTERGENIC);
			}	
		}
		
		Set<String> sampset = pathmap.keySet();
		for (String s : sampset)
		{
			String dir = pathmap.get(s);
			for(SVType t : types)
			{
				String tstr = t.getString();
				for (GeneFunc ef : pelist)
				{
					String estr = ef.toString();
					String segpath = dir + File.separator + tstr + File.separator + tstr + "_" + estr + "_genes.txt";
					String uphhpath = dir + File.separator + tstr + File.separator + tstr + "_" + estr + "_genes_uphh.txt";
					
					try 
					{
						Set<String> genes = new HashSet<String>();
						FileReader fr = new FileReader(segpath);
						BufferedReader br = new BufferedReader(fr);
						
						String line = null;
						while ((line = br.readLine()) != null)
						{
							String[] fields = line.split("\t");
							if (fields == null || fields.length < 1) continue;
							genes.add(fields[0]);
						}
						
						br.close();
						fr.close();
						
						for (String g : genes)
						{
							Map<String, Hit> gmap = hitmap.get(g);
							if(gmap == null)
							{
								gmap = new HashMap<String, Hit>();
								hitmap.put(g, gmap);
							}
							Hit h = gmap.get(s);
							if (h == null)
							{
								h = new Hit();
								h.poseff = ef;
								h.type = t;
								h.seg = true;
								gmap.put(s, h);
							}
							else
							{
								if (ef.getPriority() < h.poseff.getPriority())
								{
									h.poseff = ef;
									h.type = t;
									h.seg = true;
									continue;
								}
								if (!h.seg)
								{
									h.poseff = ef;
									h.type = t;
									h.seg = true;
									continue;
								}
								if (t == SVType.DEL && h.type != SVType.DEL)
								{
									h.poseff = ef;
									h.type = t;
									h.seg = true;
									continue;
								}
							}
						}
					} 
					catch (IOException e) {
						System.err.println("ERROR: File " + segpath + " could not be read!");
						System.err.println("No segregating genes (" + estr + " " + tstr + ") found for " + s);
						//e.printStackTrace();
					}
					
					
					try 
					{
						Set<String> genes = new HashSet<String>();
						FileReader fr = new FileReader(uphhpath);
						BufferedReader br = new BufferedReader(fr);
						
						String line = null;
						while ((line = br.readLine()) != null)
						{
							String[] fields = line.split("\t");
							if (fields == null || fields.length < 1) continue;
							genes.add(fields[0]);
						}
						
						br.close();
						fr.close();
						
						for (String g : genes)
						{
							Map<String, Hit> gmap = hitmap.get(g);
							if(gmap == null)
							{
								gmap = new HashMap<String, Hit>();
								hitmap.put(g, gmap);
							}
							Hit h = gmap.get(s);
							if (h == null)
							{
								h = new Hit();
								h.poseff = ef;
								h.type = t;
								h.seg = false;
								gmap.put(s, h);
							}
							else
							{
								if (ef.getPriority() < h.poseff.getPriority())
								{
									h.poseff = ef;
									h.type = t;
									h.seg = false;
									continue;
								}
								if (t == SVType.DEL && h.type != SVType.DEL)
								{
									h.poseff = ef;
									h.type = t;
									h.seg = false;
									continue;
								}
							}
						}
					} 
					catch (IOException e) {
						System.err.println("ERROR: File " + uphhpath + " could not be read!");
						System.err.println("No unpaired half-het genes (" + estr + " " + tstr + ") found for " + s);
						//e.printStackTrace();
					}
					
				}
			}
		}
		
		
		return hitmap;
	}

	public static void printHitMap(Map<String, Map<String, Hit>> hitmap, GeneSet genes)
	{
		List<GeneHit> list = new LinkedList<GeneHit>();
		Set<String> geneset = hitmap.keySet();
		for (String g : geneset)
		{
			GeneHit gh = new GeneHit();
			gh.gene = g;
			int h = 0;
			Map<String, Hit> hmap = hitmap.get(g);
			if (hmap != null) h = hmap.size();
			gh.hits = h;
			list.add(gh);
		}
		
		Collections.sort(list);
		
		//GeneName	#Hits	HitDetails	OMIM
		//Hit details: Sample,Type,PosEff,Seg;Sample,Type,PosEff,Seg;(etc)
		System.out.println("GeneName\t#Hits\tHitDetails\tOMIM");
		
		for (GeneHit gh : list)
		{
			String g = gh.gene;
			System.out.print(g + "\t" + gh.hits + "\t");
			Map<String, Hit> hmap = hitmap.get(g);
			Set<String> sset = hmap.keySet();
			boolean first = true;
			for (String s : sset)
			{
				if(!first)System.out.print(";");
				Hit h = hmap.get(s);
				System.out.print(s + "," + h.type.toString() + "," + h.poseff.toString() + "," + h.seg);
				first = false;
			}
			System.out.print("\t");
			//Try to get gene
			List<Gene> gmatches = genes.getGeneByName(g);
			if(gmatches == null || gmatches.isEmpty()) System.out.print("[N/A]");
			else
			{
				first = true;
				for(Gene go : gmatches)
				{
					String anno = go.getAnnotation(OMIMGeneMapImporter.ANNO_KEY);
					if(anno != null && !anno.isEmpty())
					{
						if(!first) System.out.print(";");
						System.out.print(anno);
						first = false;
						break;
					}
				}
				if(first) System.out.print("[N/A]"); //No genes had annotations
			}
			System.out.print("\n");
		}
		
		
	}
	
	public static void runGeneTally(String[] args, GenomeBuild gb)
	{
		String inPath = null;
		String tlistraw = null;
		boolean use_intronic = false;
		boolean use_intergenic = false;
		boolean exonic_only = false;
		String omimPath = null;
		
		for (int i = 0; i < args.length; i++)
		{
			String s = args[i];
			if (s.equals(OP_INFILE))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_INFILE + " flag MUST be followed by input table path!");
					printUsage();
					System.exit(1);
				}
				inPath = args[i+1];
			}
			if (s.equals(OP_TYPELIST))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_TYPELIST + " flag MUST be followed by comma delimited list of SV types!");
					printUsage();
					System.exit(1);
				}
				tlistraw = args[i+1];
			}
			if (s.equals(OP_OMIM_TABLE))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_OMIM_TABLE + " flag MUST be followed by path to OMIM genemap2.txt file!");
					printUsage();
					System.exit(1);
				}
				omimPath = args[i+1];
			}
			if (s.equals(OP_INTRONIC))
			{
				use_intronic = true;
			}
			if (s.equals(OP_INTERGENIC))
			{
				use_intergenic = true;
			}
			if (s.equals(OP_EXONIC_ONLY))
			{
				exonic_only = true;
			}
		}
		
		
		boolean pass = true;
		if (inPath == null || inPath.isEmpty())
		{
			System.err.println("ERROR: Input path is required!");
			pass = false;
		}
		List<SVType> tlist = new LinkedList<SVType>();
		if (tlistraw == null || tlistraw.isEmpty())
		{
			tlist.add(SVType.DEL);
			tlist.add(SVType.DUP);
		}
		else
		{
			String[] fields = tlistraw.split(",");
			for (String s : fields)
			{
				SVType t = SVType.getType(s);
				if (t == null)
				{
					System.err.println("ERROR: SV Type \"" + t + "\" not recognized!");
					pass = false;
					break;
				}
				tlist.add(t);
			}
		}
		
		if (!pass)
		{
			printUsage();
			System.exit(1);
		}
		
		//Load genes
		GeneSet gs = GeneSet.loadRefGene(gb);
		if (gb == null)
		{
			System.err.println("ERROR: Gene set could not be opened!");
			System.exit(1);
		}
		System.err.println("Gene info loaded!");
		
		//Load omim, if specified
		if(omimPath != null && !omimPath.isEmpty())
		{
			OMIMGeneMapImporter omimImport = new OMIMGeneMapImporter(omimPath);
			if (omimImport.importTable(gs)) System.err.println("OMIM Table read!");
			else System.err.println("OMIM Table read failed!");
		}
		
		
		//Parse input table
		Map<String,String> pathmap = new HashMap<String, String>();
		
		try
		{
			FileReader fr = new FileReader(inPath);
			BufferedReader br = new BufferedReader(fr);
			
			String line = null;
			while ((line = br.readLine()) != null)
			{
				String[] fields = line.split("\t");
				if (fields == null || fields.length < 2) continue;
				pathmap.put(fields[0], fields[1]);
			}
			
			br.close();
			fr.close();
		}
		catch (IOException e)
		{
			System.err.println("ERROR: File " + inPath + " could not be read!");
			e.printStackTrace();
			System.exit(1);
		}
		
		//Actually run it
		Map<String, Map<String, Hit>> hitmap = generateHitMap(pathmap, tlist, use_intronic, use_intergenic, exonic_only);
		
		//Print
		printHitMap(hitmap, gs);
		
	}

}
