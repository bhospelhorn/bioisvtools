package hospelhornbg_genomeBuild;


import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import hospelhornbg_bioinformatics.BreakendPair;
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.Translocation;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_bioinformatics.VariantPool.InfoDefinition;
import hospelhornbg_genomeBuild.Gene.Exon;
import waffleoRai_Utils.CompositeBuffer;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.Huffman;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

/*
 * UPDATE LOG
 * 	Version 1.0.0 created April 3, 2018
 * 
 * 1.0.0 -> 1.0.1 | April 5, 2018
 * 		Added a method to print a text table file.
 * 		Debugging.
 * 
 * 1.0.1 -> 1.1.0 | April 9, 2018
 * 		Debugging.
 * 		Added standard build loaders!
 * 		First pass at making annotation threadsafe
 * 
 * 1.1.0 -> 1.2.0 | April 19, 2018
 * 		More debugging
 * 
 * 1.2.0 -> 1.2.1 | July 18, 2018
 * 		Added an inputstream constructor
 * 
 * 1.2.1 -> 1.2.2 | July 29, 2018
 * 		Added ability to look up gene by name
 * 		Annotation methods also return list of genes for non-intergenic variants.
 * 
 * 1.2.2 -> 1.2.3 | August 30, 2018
 * 		Fixed bug for reading file from package- needed inputstream
 */


/**
 * A set of gene annotations for a genome build.
 * @author Blythe Hospelhorn
 * @version 1.2.3
 * @since August 30, 2018
 *
 */
public class GeneSet 
{
	
	/* --- Constants --- */
	
	/**
	 * Expected first four bytes (as ASCII string) in a compressed gbgd
	 * (genome build gene database) file.
	 */
	public static final String MAGIC_GBGD_COMPRESSED = "gbgc";
	
	/**
	 * Expected first four bytes (as ASCII string) in an uncompressed gbgd
	 * (genome build gene database) file.
	 */
	public static final String MAGIC_GBGD = "gbgd";
	
	/**
	 * Four byte ASCII marker for a chromosome/contig block in an uncompressed
	 * gbgd file.
	 */
	public static final String MAGIC_CHROMBLOCK = "CHRM";
	
	/**
	 * Distance between contig index points in basepairs.
	 * Should make for faster searching...
	 */
	public static final int INDEX_DIST = 1000000;
	
	/**
	 * Distance from the end of a gene (in basepairs) that is considered
	 * close enough to be called "upstream" or "downstream" rather than simply
	 * intergenic.
	 */
	public static final int NEARGENE_DIST = 500;
	
	/**
	 * ID = REFGENE_GENE
	 * <br>Number = .
	 * <br>Type = String
	 * <br>Description = "RefSeq genes associated with variant"
	 */
	public static final InfoDefinition INFODEF_INFO_GENES = new InfoDefinition("REFGENE_GENE", VariantPool.INFODEF_STRING, "RefSeq genes associated with variant", VariantPool.INFODEF_NARGS_VARIABLE);
	
	/**
	 * ID = REFGENE_GFUNC
	 * <br>Number = 1
	 * <br>Type = String
	 * <br>Description = "Highest priority location relative to associated gene(s)"
	 */
	public static final InfoDefinition INFODEF_INFO_GFUNC = new InfoDefinition("REFGENE_GFUNC", VariantPool.INFODEF_STRING, "Highest priority location relative to associated gene(s)", 1);
	
	/**
	 * ID = REFGENE_GFLANKR
	 * <br>Number = .
	 * <br>Type = String
	 * <br>Description = "Right flanking gene(s) for intergenic variants"
	 */
	public static final InfoDefinition INFODEF_INFO_RFLANK = new InfoDefinition("REFGENE_GFLANKR", VariantPool.INFODEF_STRING, "Right flanking gene(s) for intergenic variants", VariantPool.INFODEF_NARGS_VARIABLE);
	
	/**
	 * ID = REFGENE_GFLANKL
	 * <br>Number = .
	 * <br>Type = String
	 * <br>Description = "Left flanking gene(s) for intergenic variants"
	 */
	public static final InfoDefinition INFODEF_INFO_LFLANK = new InfoDefinition("REFGENE_GFLANKL", VariantPool.INFODEF_STRING, "Left flanking gene(s) for intergenic variants", VariantPool.INFODEF_NARGS_VARIABLE);
	
	/**
	 * ID = REFGENE_GDISTR
	 * <br>Number = .
	 * <br>Type = Integer
	 * <br>Description = "Number of base pairs from nearest right flanking gene(s) for intergenic variants"
	 */
	public static final InfoDefinition INFODEF_INFO_RDIST = new InfoDefinition("REFGENE_GDISTR", VariantPool.INFODEF_INT, "Number of base pairs from nearest right flanking gene(s) for intergenic variants", VariantPool.INFODEF_NARGS_VARIABLE);
	
	/**
	 * ID = REFGENE_GDISTL
	 * <br>Number = .
	 * <br>Type = Integer
	 * <br>Description = "Number of base pairs from nearest left flanking gene(s) for intergenic variants"
	 */
	public static final InfoDefinition INFODEF_INFO_LDIST = new InfoDefinition("REFGENE_GDISTL", VariantPool.INFODEF_INT, "Number of base pairs from nearest left flanking gene(s) for intergenic variants", VariantPool.INFODEF_NARGS_VARIABLE);

	/* --- Instance Variables --- */
	
	private String name;
	
	private GenomeBuild genome;
	
	private Map<Contig, ChromSet> genemap;
	
	/* --- Inner Structures --- */
	
	private static class ChromSet
	{
		private int[] indexTable;
		private boolean adjustSize;
		
		private GeneList genes;
		private OverlapMap overlapMap;
		
		private static class GeneList
		{
			private List<Gene> genes;
			
			public GeneList(int initSize)
			{
				genes = new ArrayList<Gene>(initSize);
			}
			
			public synchronized Gene get(int index)
			{
				return genes.get(index);
			}
			
			public synchronized int size()
			{
				return genes.size();
			}
		
			public synchronized void addAll(Collection<Gene> c)
			{
				genes.addAll(c);
			}
		
			public synchronized void sort()
			{
				Collections.sort(genes);
			}
			
			public synchronized List<Gene> subList(int fromIndex, int toIndex)
			{
				return genes.subList(fromIndex, toIndex);
			}
			
			public synchronized List<Gene> getGenes()
			{
				List<Gene> glist = new ArrayList<Gene>(genes.size());
				glist.addAll(genes);
				return glist;
			}
			
		}
		
		private static class OverlapMap
		{
			private Map<Integer, List<Integer>> overlapMap;
			
			public OverlapMap()
			{
				overlapMap = new HashMap<Integer, List<Integer>>();
			}
			
			public synchronized void put(Integer key, List<Integer> value)
			{
				overlapMap.put(key, value);
			}
			
			public synchronized List<Integer> get(Integer key)
			{
				return overlapMap.get(key);
			}
			
		}
		
		public ChromSet(long chrSize)
		{
			//System.err.println("GeneSet.ChromSet.<init> || Called! chrSize = " + chrSize);
			int megabases = (int)(chrSize / INDEX_DIST);
			indexTable = new int[megabases + 1];
			int geneCountEst = (int)(chrSize / 30000L);
			genes = new GeneList(geneCountEst);
			adjustSize = false;
			overlapMap = new OverlapMap();
		}
		
		public ChromSet(int nGenes)
		{
			genes = new GeneList(nGenes);
			long szEstimate = nGenes * 30000;
			int megabases = (int)(szEstimate / 1000000L);
			indexTable = new int[megabases + 1];
			adjustSize = true;
			overlapMap = new OverlapMap();
		}
		
		private void updateIndexTable()
		{
			int listi = 0;
			int pos = INDEX_DIST;
			
			if (adjustSize)
			{
				Gene lastg = genes.get(genes.size() - 1);
				long sz = lastg.getTranscriptEnd() + 1;
				int megabases = (int)(sz / INDEX_DIST);
				indexTable = new int[megabases + 1];
			}
			
			indexTable[0] = 0;
			int gnum = genes.size();
			
			for (int i = 1; i < indexTable.length; i++)
			{
				/*
				System.err.println("GeneSet.ChromSet.updateIndexTable || pos: " + pos);
				boolean match = false;
				while (!match)
				{
					if (listi >= gnum) break;
					Gene g = genes.get(listi);
					System.err.println("GeneSet.ChromSet.updateIndexTable || Checking pos against gene start " + g.getTranscriptStart());
					if (g.getTranscriptStart() <= pos) match = true;
					else listi++;
				}*/
				
				//Run up index until it finds a gene that starts after position
				boolean breaker = false;
				while (!breaker)
				{
					if (listi >= gnum) break;
					Gene g = genes.get(listi);
					//System.err.println("GeneSet.ChromSet.updateIndexTable || Checking pos against gene start " + g.getTranscriptStart());
					if (g.getTranscriptStart() > pos) breaker = true;
					else listi++;
				}
				
				//Back up index until it finds a gene that ends before position
				breaker = false;
				while (!breaker)
				{
					if (listi <= 0 || listi >= gnum) break;
					Gene g = genes.get(listi);
					//System.err.println("GeneSet.ChromSet.updateIndexTable || Checking pos against gene start " + g.getTranscriptStart());
					if (g.getTranscriptEnd() <= pos) breaker = true;
					else listi--;
				}
				if (listi < gnum) listi++;
				
				indexTable[i] = listi;
				//System.err.println("GeneSet.ChromSet.updateIndexTable || indexTable[" + i + "] set to " + listi);
				pos += INDEX_DIST;
			}
		}
		
		private void updateOverlapMap()
		{
			//System.err.println("GeneSet.ChromSet.updateOverlapMap || Called");
			overlapMap = new OverlapMap();
			//System.err.println("GeneSet.ChromSet.updateOverlapMap || Map instantiated");
			
			int gNum = genes.size();
			
			//DEBUG BLOCK
			//System.err.println("GeneSet.ChromSet.updateOverlapMap || gNum = " + gNum);
			//System.out.println("SUPER DUPER DEBUG!! GeneSet.ChromSet.updateOverlapMap || genes: ");
			//for (Gene g : genes) System.out.println(g.printInfo());
			//System.exit(0);
			
			//int hCount = 0;
			//int mCount = 0;
			//int tCount = 0;
			for (int i = 0; i < gNum - 1; i++)
			{
				//System.err.println("GeneSet.ChromSet.updateOverlapMap || -- i = " + i);
				Gene g = genes.get(i);
				for (int j = 1; j < (gNum - i); j++)
				{
					//tCount++;
					//System.err.println("GeneSet.ChromSet.updateOverlapMap || -- j = " + j);
					//System.err.println("GeneSet.ChromSet.updateOverlapMap || -- o = " + (i+j));
					Gene o = genes.get(i + j);
					//Check if o starts after g ends
					if (o.getTranscriptStart() > g.getTranscriptEnd()){
						//mCount++;
						continue;
					}
					//System.err.println("GeneSet.ChromSet.updateOverlapMap || \t o does not start after g ends...");
					//Check if g and o overlap
					if (g.overlaps(o))
					{
						//System.err.println("GeneSet.ChromSet.updateOverlapMap || \t g and o overlap...");
						List<Integer> glist = overlapMap.get(i);
						List<Integer> olist = overlapMap.get(i + j);
						if (glist == null) glist = new LinkedList<Integer>();
						if (olist == null) olist = new LinkedList<Integer>();
						glist.add(new Integer(i+j));
						olist.add(i);
						overlapMap.put(i, glist);
						overlapMap.put((i+j), olist);
						//hCount++;
						//System.err.println("GeneSet.ChromSet.updateOverlapMap || \t references saved...");
					}
					//else mCount++;
				}
			}
			//System.err.println("GeneSet.ChromSet.updateOverlapMap || Complete");
			//System.err.println("GeneSet.ChromSet.updateOverlapMap || Total Crosses: " + tCount);
			//System.err.println("GeneSet.ChromSet.updateOverlapMap || Hits: " + hCount);
			//System.err.println("GeneSet.ChromSet.updateOverlapMap || Misses: " + mCount);
		}
		
		public void addGenes(Collection<Gene> geneset)
		{
			//System.err.println("GeneSet.ChromSet.addGenes || Called");
			genes.addAll(geneset);
			//System.err.println("GeneSet.ChromSet.addGenes || Genes added");
			genes.sort();
			//System.err.println("GeneSet.ChromSet.addGenes || Genes sorted");
			updateIndexTable();
			//System.err.println("GeneSet.ChromSet.addGenes || Genes indexed");
			updateOverlapMap();
			//System.err.println("GeneSet.ChromSet.addGenes || Gene overlap mapped");
		}
		
		public FileBuffer serializeMe(int UID)
		{
			//System.err.println("GeneSet.ChromSet.serializeMe || UID " + Integer.toHexString(UID));
			int ngenes = genes.size();
			List<FileBuffer> glist = new ArrayList<FileBuffer>(ngenes);
			
			//Serialize genes
			for (int i = 0; i < ngenes; i++)
			{
				glist.add(genes.get(i).serializeMe());
			}
			//System.err.println("GeneSet.ChromSet.serializeMe || Genes Serialized: " + ngenes);
			
			//Generate the position table. 
			//Note - offsets are from START OF GENE RECORDS!!!
			int tentries = indexTable.length;
			//System.err.println("GeneSet.ChromSet.serializeMe || Index Entries: " + tentries);
			int pos = 0;
			int listi = 0;
			FileBuffer idx = new FileBuffer((tentries * 4) + 4, true);
			idx.addToFile(tentries);
			
			for (int i = 0; i < tentries; i++)
			{
				int ind = indexTable[i];
				//listi = 0;
				//System.err.println("GeneSet.ChromSet.serializeMe || Index entry " + i + " : " + ind);
				while (listi < ind)
				{
					pos += glist.get(listi).getFileSize();
					listi++;
				}
				idx.addToFile(pos);
				//System.err.println("GeneSet.ChromSet.serializeMe || Entry: 0x" + Integer.toHexString(pos));
			}
			
			//Calculate chunk size
			int csz = (int)idx.getFileSize();
			for (FileBuffer g : glist) csz += g.getFileSize();
			csz += 4 + 4;
			//System.err.println("GeneSet.ChromSet.serializeMe || Chunk Size: 0x" + Integer.toHexString(csz));
			
			//Serialize the header
			FileBuffer header = new FileBuffer((4*4), true);
			header.printASCIIToFile(MAGIC_CHROMBLOCK);
			header.addToFile(csz);
			header.addToFile(UID);
			header.addToFile(ngenes);
			
			//Make a composite
			FileBuffer CHRM = new CompositeBuffer(2 + ngenes);
			CHRM.addToFile(header);
			CHRM.addToFile(idx);
			for (FileBuffer g : glist) CHRM.addToFile(g);
			
			return CHRM;
		}
		
		/*
		public HitRecord search(int position)
		{
			return search(position, false, false);
		}*/
		
		public HitRecord search(int position, boolean getLeftFlanking, boolean getRightFlanking)
		{
			//If this contig has no genes, the main body of this method will crash hard!
			if (this.genes.size() <= 0)
			{
				HitRecord rec = new HitRecord();
				rec.index = -1;
				rec.lindex = -1;
				rec.location = GeneFunc.INTERGENIC;
				return rec;
			}
			
			//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Called - position = " + position + " | left flank: " + getLeftFlanking + " | right flank: " + getRightFlanking);
			int chunk = position/INDEX_DIST;
			if (chunk < 0 || chunk >= indexTable.length) return null;
			//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Chunk check passed! Chunk = " + chunk);
			int listi = indexTable[chunk];
			if(listi < 0) System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || WARNING!! listi < 0 | Chunk = " + chunk);
			if(listi < 0) System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || WARNING!! listi < 0 | listi = " + listi);
			//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Initial list index = " + listi);
			//The index lookup may put the position AFTER the index it should actually be at. Check and backpedal if needed...
			if (listi == genes.size()) listi--;
			while ((listi > 0) && (genes.get(listi).getTranscriptStart() > position)) {
					listi--;
			}
			if(listi < 0) System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || WARNING!! listi < 0 | After backpedal | listi = " + listi);
			//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || List index rewound to: " + listi);
				
			//First, we must check to see if it's BEFORE the first gene!!
			if (listi == 0 && genes.get(0).getTranscriptStart() > position)
			{
				//Process it on the spot.
				//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Position falls before the first gene.");
				HitRecord rec = new HitRecord();
				rec.index = -1;
				rec.lindex = 0;
				//Close enough to first gene?
				int rightDist = genes.get(0).getTranscriptStart() - position;
				//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || rightDist = " + rightDist);
				if (rightDist <= NEARGENE_DIST){
					//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Position is close to nearest gene.");
					if (genes.get(0).getStrand()) rec.location = GeneFunc.UPSTREAM;
					else
					{
						//Check to see if 0 overlaps with anything on the + strand...
						//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Nearest gene is on minus strand. Checking for overlapping genes that are on plus strand.");
						rec.location = GeneFunc.DOWNSTREAM;
						List<Integer> overlap = overlapMap.get(0);
						if (overlap != null && !overlap.isEmpty())
						{
							//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Overlapping genes found.");
							for (Integer n : overlap)
							{
								Gene g = genes.get(n);
								int dist = g.getTranscriptStart() - position;
								//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Overlapping gene: " + g.getName() + " | Dist = " + dist);
								if (dist <= NEARGENE_DIST)
								{
									if (g.getStrand()){
										//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Gene is on plus strand and within near distance of position.");
										rec.location = GeneFunc.UPSTREAM;
										rec.lindex = n;
										break;
									}
								}
							}
						}
					}
				}
				else rec.location = GeneFunc.INTERGENIC;
				//rec.dist = rightDist;
				
				//Get the right flanking gene(s), if requested.
				if (getRightFlanking)
				{
					//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Getting right flanking gene(s)...");
					Gene g0 = getGene(0);
					//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Gene 0 = " + g0.getName());
					List<Integer> overlap = overlapMap.get(0);
					//Gene 0 should always be the one with the lowest start position
					int dist = g0.getTranscriptStart() - position;
					//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Gene 0 dist = " + dist);
					rec.addRightFlankingGene(g0, dist);
					if (overlap != null && !overlap.isEmpty())
					{
						//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Overlapping genes found.");
						for (Integer i : overlap)
						{
							Gene ogene = getGene(i);
							//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Overlapping gene: " + ogene.getName());
							if (ogene.getTranscriptStart() <= g0.getTranscriptStart())
							{
								//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Overlapping gene " + ogene.getName() + " added as right flanking.");
								dist = ogene.getTranscriptStart() - position;
								rec.addRightFlankingGene(ogene, dist);
							}
						}	
					}
				}
			
				//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Returning... (Pre-gene 0)");
				return rec;
			}
			
			int gnum = genes.size();
			//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Gene number = " + gnum);
			
			//Advance listi until hit gene with start point above position... (List should be sorted!!)
			while ((listi + 1 < gnum) && genes.get(listi+1).getTranscriptStart() <= position) listi++;
			if(listi < 0) System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || WARNING!! listi < 0 | After forward | listi = " + listi);
			//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Index set = " + listi);
			
			//Initialize return object
			HitRecord rec = new HitRecord();
			rec.index = listi;
			if(listi < 0) System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || WARNING!! listi < 0 | After rec set | listi = " + listi);
			
			//See if the position is in the previous gene, and in what capacity.
			Gene g = genes.get(listi);
			//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Gene: " + g.getName() + " | " + g.getID());
			//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Gene Position: " + g.getTranscriptStart() + "-" + g.getTranscriptEnd());
			//boolean ingene = g.positionInGene(position);
			
			//See how position is related to its indexed gene
			rec.location = g.getRelativeLocationEffect(position);
			//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Initial location effect: " + rec.location.toString());
			rec.lindex = listi;
			
			//See if can get anything better for overlapping genes
			List<Integer> olgenes = overlapMap.get(listi);
			if (olgenes != null && !olgenes.isEmpty())
			{
				//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Overlapping gene(s) found!");
				for (Integer n : olgenes)
				{
					GeneFunc f = genes.get(n).getRelativeLocationEffect(position);
					//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Overlapping gene: " + genes.get(n).getName());
					//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Position effect: " + f.toString());
					if (f.getPriority() < rec.location.getPriority())
					{
						rec.location = f;
						rec.lindex = n;
						//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Position effect updated to: " + f.toString());
					}
				}
			}
			
			//If intergenic, see if can get anything better for right flanking and its overlapping genes
			if ((listi < gnum - 1) && rec.location.isIntergenic())
			{
				//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Position is intergenic. Checking overlapping right flanking genes...");
				Gene right = genes.get(listi + 1);
				//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Right gene: " + right.getName() + " | " + right.getID());
				GeneFunc rEff = right.getRelativeLocationEffect(position);
				//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Right gene relative effect: " + rEff.toString());
				if (rEff.getPriority() < rec.location.getPriority())
				{
					rec.location = rEff;
					rec.lindex = listi + 1;
					//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Position effect updated to: " + rEff.toString());
				}
				olgenes = overlapMap.get(listi + 1);
				if (olgenes != null && !olgenes.isEmpty())
				{
					//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Right flanking overlapping gene(s) found!");
					for (Integer n : olgenes)
					{
						GeneFunc f = genes.get(n).getRelativeLocationEffect(position);
						//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Overlapping gene: " + genes.get(n).getName());
						//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Position effect: " + f.toString());
						if (f.getPriority() < rec.location.getPriority())
						{
							rec.location = f;
							rec.lindex = n;
							//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Position effect updated to: " + f.toString());
						}
					}	
				}
			}
			//Nab flanking if intergenic and requested...
			if (rec.location.isIntergenic())
			{
				//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Position is intergenic. Looking for flanking genes...");
				if (getLeftFlanking)
				{
					//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Finding left flanking gene(s)...");
					Gene gi = getGene(listi);
					int d = position - gi.getTranscriptEnd();
					rec.addLeftFlankingGene(gi, d);
					//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Left flanking gene found: " + gi.getName() + " | " + gi.getID());
					olgenes = overlapMap.get(listi);
					if(olgenes != null && !olgenes.isEmpty())
					{
						//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Left flanking overlapping gene(s) found!");
						for (Integer n : olgenes)
						{
							Gene go = getGene(n);
							//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Overlapping gene: " + go.getName());
							if (go.getTranscriptEnd() >= gi.getTranscriptEnd())
							{
								d = position - go.getTranscriptEnd();
								rec.addLeftFlankingGene(go, d);
								//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Overlapping gene added: " + go.getName() + " | " + go.getID() + " | dist = " + d);
							}
						}		
					}
				}
				if (getRightFlanking && (listi < gnum - 1))
				{
					//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Finding right flanking gene(s)...");
					Gene gi = getGene(listi + 1);
					int d = gi.getTranscriptStart() - position;
					rec.addRightFlankingGene(gi, d);
					//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Right flanking gene found: " + gi.getName() + " | " + gi.getID());
					olgenes = overlapMap.get(listi + 1);
					if(olgenes!=null && !olgenes.isEmpty())
					{
						//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Right flanking overlapping gene(s) found!");
						for (Integer n : olgenes)
						{
							Gene go = getGene(n);
							//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Overlapping gene: " + go.getName());
							if (go.getTranscriptStart() <= gi.getTranscriptStart())
							{
								d = go.getTranscriptStart() - position;
								rec.addLeftFlankingGene(go, d);
								//System.err.println(Thread.currentThread().getName() + " || GeneSet.ChromSet.search || Overlapping gene added: " + go.getName() + " | " + go.getID() + " | dist = " + d);
							}
						}		
					}
				}
			}
			
			return rec;
		}
		
		public Gene getGene(int index)
		{
			return genes.get(index);
		}
		
		public List<Integer> getOverlappingIndicies(int gIndex)
		{
			return overlapMap.get(gIndex);
		}
		
		public List<Gene> getGenes(int stIndex, int edIndex)
		{
			List<Gene> sublist = new ArrayList<Gene>(edIndex - stIndex + 1);
			sublist.addAll(genes.subList(stIndex, edIndex + 1));
			return sublist;
		}
		
		public List<Gene> getAllGenes()
		{
			return genes.getGenes();
		}
		
		/*public int geneCount()
		{
			return genes.size();
		}*/
		
		public boolean hasGenes()
		{
			return (genes.size() > 0);
		}
	
	}

	private static class HitRecord
	{
		public int index;
		
		public int lindex; //index of gene that location is relative to
		public GeneFunc location;
		
		private List<FlankingGene> lflank;
		private List<FlankingGene> rflank;
		
		public static class FlankingGene
		{
			public Gene gene;
			public int dist;
			
			public FlankingGene(Gene g, int d)
			{
				gene = g;
				dist = d;
			}
		}
		
		public HitRecord()
		{
			index = -2; //-1 means before gene 0!!
			lindex = -2;
			location = null;
			//exon = -1;
			//dist = -1;
			lflank = new LinkedList<FlankingGene>();
			rflank = new LinkedList<FlankingGene>();
		}

		public void addLeftFlankingGene(Gene g, int d)
		{
			lflank.add(new FlankingGene(g, d));
		}
		
		public void addRightFlankingGene(Gene g, int d)
		{
			rflank.add(new FlankingGene(g, d));
		}
		
		public String[] getLFLANK_anno()
		{
			if (lflank == null) return null;
			if (lflank.isEmpty()) return null;
			
			int fcount = lflank.size();
			String[] sarr = new String[fcount];
			for (int i = 0; i < fcount; i++)
			{
				//sarr[i] = lflank.get(i).gene.getName();
				sarr[i] = lflank.get(i).gene.getID();
			}
			
			return sarr;
		}
		
		public String[] getRFLANK_anno()
		{
			if (rflank == null) return null;
			if (rflank.isEmpty()) return null;
			
			int fcount = rflank.size();
			String[] sarr = new String[fcount];
			for (int i = 0; i < fcount; i++)
			{
				//sarr[i] = rflank.get(i).gene.getName();
				sarr[i] = rflank.get(i).gene.getID();
			}
			
			return sarr;
		}
		
		public int[] getLDIST_anno()
		{
			if (lflank == null) return null;
			if (lflank.isEmpty()) return null;
			
			int fcount = lflank.size();
			int[] iarr = new int[fcount];
			for (int i = 0; i < fcount; i++)
			{
				iarr[i] = lflank.get(i).dist;
			}
			
			return iarr;
		}
		
		public int[] getRDIST_anno()
		{
			if (rflank == null) return null;
			if (rflank.isEmpty()) return null;
			
			int fcount = rflank.size();
			int[] iarr = new int[fcount];
			for (int i = 0; i < fcount; i++)
			{
				iarr[i] = rflank.get(i).dist;
			}
			
			return iarr;
		}
		
	}
	
	private static class AnnoRecord
	{
		public String[] genes;
		public GeneFunc effect;
		
		public String[] lflank;
		public int[] ldist;
		public String[] rflank;
		public int[] rdist;
		
		public AnnoRecord()
		{
			genes = null;
			effect = null;
			lflank = null;
			ldist = null;
			rflank = null;
			rdist = null;
		}
	
		public void loadGeneList(Collection<Gene> glist)
		{
			if (glist == null) return;
			if (glist.isEmpty()) return;
			Set<String> strset = new HashSet<String>();
			//To remove redundant names eg. same gene, different transcript.
			//While knowing the different transcripts is important, this is mostly for efficiency
			for (Gene g : glist) strset.add(g.getName());
			int sz = strset.size();
			genes = new String[sz];
			int i = 0;
			for (String s : strset)
			{
				genes[i] = s;
				i++;
			}
		}
		
	}
	
	/* --- Construction --- */
	
	/**
	 * Construct an empty GeneSet referencing existing genome build gb. 
	 * @param setname Name to set as the set name.
	 * @param gb GenomeBuild this gene set is based off of. This is required
	 * so that contigs/chromosomes can be referenced.
	 * This parameter may be null, but if it is never set explicitly, then 
	 * a new empty genome build will be generated.
	 * <br>It is strongly recommended that you provide an existing genome build.
	 */
	public GeneSet(String setname, GenomeBuild gb)
	{
		name = setname;
		genome = gb;
		genemap = new HashMap<Contig, ChromSet>();
		
		List<Contig> buildcontigs = gb.getChromosomes();
		for (Contig c : buildcontigs)
		{
			genemap.put(c, new ChromSet(c.getLength()));
		}
	}
	
	/**
	 * Construct a GeneSet object from a gbgd file on disk, and an existing
	 * genome build.
	 * @param gbgdPath Path to the file to read.
	 * @param gb Genome build the file uses. If this is not provided and the match
	 * boolean is false, an empty build will be created and populated.
	 * @param strictBuildMatch Whether to reject the file if the genome build provided
	 * as an argument doesn't match the one requested in the file. If this flag is on
	 * and there is a mismatch, the parser will throw an exception.
	 * @throws UnsupportedFileTypeException If the file could not be opened, but is not formatted
	 * in such a way that it can be read.
	 * <br>This exception is also thrown if there is a mismatch between the build provided and the build
	 * listed in the file.
	 * @throws IOException If there is an error accessing the file on disk.
	 */
	public GeneSet(String gbgdPath, GenomeBuild gb, boolean strictBuildMatch) throws UnsupportedFileTypeException, IOException
	{
		genemap = new HashMap<Contig, ChromSet>();
		parseGBGD(gbgdPath, gb, strictBuildMatch);
		
	}
	
	/**
	 * Construct a GeneSet object from a gbgd file as an inputstream and an existing genome build.
	 * @param gbgdStream The gbgd file as an input stream.
	 * @param gb Genome build the file uses. If this is not provided and the match
	 * boolean is false, an empty build will be created and populated.
	 * @param strictBuildMatch Whether to reject the file if the genome build provided
	 * as an argument doesn't match the one requested in the file. If this flag is on
	 * and there is a mismatch, the parser will throw an exception.
	 * @throws UnsupportedFileTypeException If the file could not be opened, but is not formatted
	 * in such a way that it can be read.
	 * <br>This exception is also thrown if there is a mismatch between the build provided and the build
	 * listed in the file.
	 * @throws IOException If there is an error with the input stream or decompressing the file.
	 */
	public GeneSet(InputStream gbgdStream, GenomeBuild gb, boolean strictBuildMatch) throws UnsupportedFileTypeException, IOException
	{
		genemap = new HashMap<Contig, ChromSet>();
		
		FileBuffer myFile = new FileBuffer(1024 * 500); //500KB
	//	int sz = 0;
		int b = gbgdStream.read();
		while (b != -1)
		{
			//sz++;
			//System.err.println("GenomeBuild.<init> || DEBUG: sz = " + sz + " b = " + b);
			byte y = (byte)b;
			myFile.addToFile(y);
			b = gbgdStream.read();
		}
		
		parseGBGD(myFile, gb, strictBuildMatch);
	}
	
	/* --- I/O --- */
	
	private void parseGBGD(String gbgdPath, GenomeBuild gb, boolean strictBuildMatch) throws UnsupportedFileTypeException, IOException
	{
		if (gbgdPath == null || gbgdPath.isEmpty()) throw new FileBuffer.UnsupportedFileTypeException();
		if (gb == null) throw new UnsupportedFileTypeException();
		
		//Open file
		FileBuffer myFile = FileBuffer.createBuffer(gbgdPath, true);
		
		parseGBGD(myFile, gb, strictBuildMatch);
	}
	
	private void parseGBGD(FileBuffer myFile, GenomeBuild gb, boolean strictBuildMatch) throws UnsupportedFileTypeException, IOException
	{
		// --- Header
				//Check magic, decompress, and check for inner magic
				long cPos = myFile.findString(0, 0x10, MAGIC_GBGD_COMPRESSED);
				if (cPos < 0) throw new UnsupportedFileTypeException();
				FileBuffer compressed = myFile;
				myFile = Huffman.HuffDecodeFile(compressed, cPos + 4);
				cPos = myFile.findString(0, 0x10, MAGIC_GBGD);
				if (cPos < 0) throw new UnsupportedFileTypeException();
				//DEBUG
				//myFile.writeFile("C:\\Users\\Blythe\\Documents\\Work\\Bioinformatics\\References\\grch37_refSeq_redecompressed.gbgd");
				
				//Get genome build name and confirm it with desired build
				cPos = 4;
				String gbName = myFile.getASCII_string(cPos, 12); cPos += 12;
				//System.err.println("GeneSet.parseGBDB || Build name: " + gbName);
				boolean namematch = gbName.equals(gb.getBuildName());
				//System.err.println("GeneSet.parseGBDB || Namematch = " + namematch);
				if (strictBuildMatch && !namematch) throw new UnsupportedFileTypeException();
				if (namematch) genome = gb;
				else genome = new GenomeBuild("unknown", gbName);
				
				//Get database name
				name = myFile.getASCII_string(cPos, 16); cPos += 16;
				//System.err.println("GeneSet.parseGBDB || Database Name: " + name);
				
				//--- Contig Table
				//System.err.println("GeneSet.parseGBDB || Parsing Contig Table ---");
				Map<Integer, Contig> contigUIDmap = new HashMap<Integer, Contig>();
				int cCount = myFile.intFromFile(cPos); cPos += 4;
				//System.err.println("GeneSet.parseGBDB || \tContig Count: " + cCount);
				for (int i = 0; i < cCount; i++)
				{
					int UID = myFile.intFromFile(cPos); cPos += 4;
					String cName = myFile.getASCII_string(cPos, 40); cPos += 40;
					cPos += 4; //This int is the offset of the contig's CHRM block from this value. 
								//Used if you don't want to/ can't read file into memory in one go (ie. indexing)
					//System.err.println("GeneSet.parseGBDB || \tContig 0x" + Integer.toHexString(UID));
					//System.err.println("GeneSet.parseGBDB || \t\tName: " + cName);
					Contig c = genome.getContig(cName);
					//if (c!= null) System.err.println("GeneSet.parseGBDB || \t\tMatched to build contig " + c.getUDPName());
					if ((c == null) && !namematch){
						//System.err.println("GeneSet.parseGBDB || \t\tContig not found in build. Creating...");
						c = new Contig();
						c.setUCSCName(cName);
						c.setUDPName(cName);
						c.setType(Contig.SORTCLASS_UNKNOWN);
						//Length estimated from block readings...
						genome.addContig(c);
					}
					else if (c == null && namematch) throw new UnsupportedFileTypeException();
					
					contigUIDmap.put(UID, c);
				}
				
				//--- Chrom blocks
				//System.err.println("GeneSet.parseGBDB || Parsing CHRM Blocks ---");
				for (int i = 0; i < cCount; i++)
				{
					List<Gene> glist = new LinkedList<Gene>();
					cPos = myFile.findString(cPos, cPos + 0x10, MAGIC_CHROMBLOCK);
					if (cPos < 0) throw new UnsupportedFileTypeException();
					else cPos += 4;
					//System.err.println("GeneSet.parseGBDB || \tCHRM block found at 0x" + Long.toHexString(cPos));
					//int chunkSize = myFile.intFromFile(cPos); cPos += 4;
					cPos += 4; //Skip chunk size
					int cUID = myFile.intFromFile(cPos); cPos += 4; //Chrom UID
					//System.err.println("GeneSet.parseGBDB || \tContig 0x" + Integer.toHexString(cUID));
					Contig c = contigUIDmap.get(cUID);
					//if (c!= null) System.err.println("GeneSet.parseGBDB || \t\tMatched to build contig " + c.getUDPName());
					if (c == null) throw new UnsupportedFileTypeException();
					int geneCount = myFile.intFromFile(cPos); cPos += 4;
					//System.err.println("GeneSet.parseGBDB || \t\tGene Count: " + geneCount);
					int idxEntries = myFile.intFromFile(cPos); cPos += 4;
					//System.err.println("GeneSet.parseGBDB || \t\tIndex Entries: " + idxEntries);
					//We'll actually skip the index table, since that is for reading off disk...
					cPos += (idxEntries * 4);
					
					for (int j = 0; j < geneCount; j++)
					{
						int gSt = myFile.intFromFile(cPos); cPos += 4;
						int gEd = myFile.intFromFile(cPos); cPos += 4;
						int tSt = myFile.intFromFile(cPos); cPos += 4;
						int tEd = myFile.intFromFile(cPos); cPos += 4;
						int flags = Byte.toUnsignedInt(myFile.getByte(cPos)); cPos += 2;
						//There is one byte of flags, and one byte "0" padding, so skip ahead 2.
						int nExons = Short.toUnsignedInt(myFile.shortFromFile(cPos)); cPos += 2;
						List<Exon> eList = new LinkedList<Exon>();
						for (int k = 0; k < nExons; k++)
						{
							int eSt = myFile.intFromFile(cPos); cPos += 4;
							int eEd = myFile.intFromFile(cPos); cPos += 4;
							eList.add(new Exon(eSt, eEd));
						}
						String gID = myFile.getASCII_string(cPos, '\0');
						cPos += gID.length() + 1;
						//Skip padding, if needed...
						if (gID.length() % 2 == 0) cPos++;
						String gName = myFile.getASCII_string(cPos, '\0');
						cPos += gName.length() + 1;
						if (gName.length() % 2 == 0) cPos++;
						
						Gene g = new Gene(nExons);
						g.setChromosome(c);
						g.setStart(gSt);
						g.setEnd(gEd);
						g.setTranslationStart(tSt);
						g.setTranslationEnd(tEd);
						g.setID(gID);
						g.setName(gName);
						g.setStrand((flags & 0x1) != 0);
						g.setNCRNA((flags & 0x2) != 0);
						g.addExons(eList);
						
						glist.add(g);
					}
					
					ChromSet chr = null;
					if (namematch) chr = new ChromSet(c.getLength());
					else chr = new ChromSet(geneCount);
					
					chr.addGenes(glist);
					
					genemap.put(c, chr);
				}
	}
	
	/**
	 * Serialize this gene set and write it out to disk in a compressed gbgd file.
	 * @param gbdbPath Path on local file system to write file to.
	 * @throws UnsupportedFileTypeException If there is a problem serializing the set.
	 * @throws IOException If there is an error writing to disk.
	 */
	public void serializeGBGD(String gbdbPath) throws UnsupportedFileTypeException, IOException
	{
		//Calculate some sizes...
		int headersize = 4 + 12 + 16;
		//System.err.println("GeneSet.serializeGBGD || headersize = " + headersize);
		int cttblsize = 4 + (genemap.size() * (4 + 40 + 4));
		//System.err.println("GeneSet.serializeGBGD || cttblsize = " + cttblsize);
		//System.err.println("GeneSet.serializeGBGD || Sizes calculated");
		
		//Header
		FileBuffer header = new FileBuffer(headersize, true);
		header.printASCIIToFile(MAGIC_GBGD);
		String gname = genome.getBuildName();
		if (gname.length() > 12) gname = gname.substring(0, 12);
		header.printASCIIToFile(gname);
		for (int i = gname.length(); i < 12; i++) header.addToFile((byte)0x00);
		
		if (name.length() > 16) name = name.substring(0, 16);
		header.printASCIIToFile(name);
		for (int i = name.length(); i < 16; i++) header.addToFile((byte)0x00);
		//int hsz = (int)header.getFileSize();
		//System.err.println("GeneSet.serializeGBGD || Header generated - size = 0x" + Integer.toHexString(hsz));
		//System.exit(1);
		
		//Chrom initialize
		int nChrom = genemap.size();
		List<Contig> clist = new ArrayList<Contig>(nChrom);
		clist.addAll(genemap.keySet());
		Collections.sort(clist);
		//System.err.println("GeneSet.serializeGBGD || Chrom initialize");
		
		//Contig Table/Blocks
		FileBuffer contigmap = new FileBuffer(cttblsize, true);
		FileBuffer chrmchunk = new CompositeBuffer(nChrom);
		contigmap.addToFile(nChrom);
		int pos = 0;
		//System.err.println("GeneSet.serializeGBGD || Contig table/blocks initialize");
		for (Contig c : clist)
		{
			int UID = c.hashCode() ^ gname.hashCode();
			
			//Map
			contigmap.addToFile(UID);
			String cname = c.getUDPName();
			if (cname.length() > 40) throw new UnsupportedFileTypeException();
			contigmap.printASCIIToFile(cname);
			for (int i = cname.length(); i < 40; i++) contigmap.addToFile((byte)0x00);
			contigmap.addToFile(pos);
			
			//Block
			FileBuffer CHRM = genemap.get(c).serializeMe(UID);
			chrmchunk.addToFile(CHRM);
			
			//Track position...
			pos += CHRM.getFileSize();
			//System.err.println("GeneSet.serializeGBGD || Chrom " + c.getUDPName() + " block generated");
		}
		//int msz = (int)contigmap.getFileSize();
		//int csz = (int)chrmchunk.getFileSize();
		//System.err.println("GeneSet.serializeGBGD || contigmap - size = 0x" + Integer.toHexString(msz));
		//System.err.println("GeneSet.serializeGBGD || chrmchunk - size = 0x" + Integer.toHexString(csz));
		//System.err.println("GeneSet.serializeGBGD || total expected size = 0x" + Integer.toHexString(csz + msz + hsz));
		
		//Composite
		FileBuffer fullfile = new CompositeBuffer(3);
		fullfile.addToFile(header);
		fullfile.addToFile(contigmap);
		fullfile.addToFile(chrmchunk);
		//System.err.println("GeneSet.serializeGBGD || Composite complete");
		//System.err.println("GeneSet.serializeGBGD || composite size = 0x" + Long.toHexString(fullfile.getFileSize()));
		
		//fullfile.writeFile(gbdbPath + "_decompressed.gbgd");
		//System.err.println("GeneSet.serializeGBGD || Decompressed file write complete");
		
		//Compress
		FileBuffer compressed = Huffman.HuffEncodeFile(fullfile, 8, MAGIC_GBGD_COMPRESSED);
		//System.err.println("GeneSet.serializeGBGD || Compression complete");
		
		//Write
		compressed.writeFile(gbdbPath);
		//System.err.println("GeneSet.serializeGBGD || Write complete");
		
	}
	
	/**
	 * Output a text table representation of this gene set to a file on disk.
	 * @param tablePath Path to write table to.
	 * @throws IOException If there is an error writing to disk.
	 */
	public void outputTable(String tablePath) throws IOException
	{
		FileWriter fw = new FileWriter(tablePath);
		BufferedWriter bw = new BufferedWriter(fw);
		
		List<Contig> contigList = new LinkedList<Contig>();
		contigList.addAll(genemap.keySet());
		Collections.sort(contigList);
		
		for (Contig c : contigList)
		{
			ChromSet cGenes = genemap.get(c);
			List<Gene> glist = cGenes.getAllGenes();
			for (Gene g : glist)
			{
				bw.write(g.printInfo());
				bw.newLine();
			}
		}
		
		bw.close();
		fw.close();
	}
	
	/* --- Annotation --- */
	
	/**
	 * Get all info definitions used by this class when annotating variants.
	 * @return LinkedList of InfoDefinitions used by this class.
	 */
	public static Collection<InfoDefinition> getInfoDefinitions()
	{
		List<InfoDefinition> coll = new LinkedList<InfoDefinition>();
		coll.add(INFODEF_INFO_GENES);
		coll.add(INFODEF_INFO_GFUNC);
		coll.add(INFODEF_INFO_LFLANK);
		coll.add(INFODEF_INFO_RFLANK);
		coll.add(INFODEF_INFO_LDIST);
		coll.add(INFODEF_INFO_RDIST);
		return coll;
	}
	
	private AnnoRecord annotateRegion(Contig c, int st, int ed, List<Gene> genelist)
	{
		if (c == null) return null;
		ChromSet cGenes = genemap.get(c);
		if (cGenes == null) return null;
		
		//int gnum = cGenes.geneCount();
		
		if (!cGenes.hasGenes())
		{
			AnnoRecord rec = new AnnoRecord();
			rec.effect = GeneFunc.INTERGENIC;
			return rec;
		}
		
		HitRecord posHit = cGenes.search(st, true, false);
		HitRecord endHit = cGenes.search(ed, false, true);
		
		AnnoRecord rec = new AnnoRecord();
		
		if (posHit.location.isIntergenic() && endHit.location.isIntergenic())
		{
			//If both breakends are intergenic (including upstream and downstream)...
			
			//If at this point, the breakends are both considered intergenic, then
			//	that means that neither lie inside of a gene, even one overlapping another.
			
			if (posHit.index == endHit.index)
			{
				//They must be between the same two genes.
				
				//Get effect
				GeneFunc eff = posHit.location;
				if (endHit.location.getPriority() < eff.getPriority()) eff = endHit.location;
				
				//Get left flanking genes and distances (EDIT: Don't need to anymore)
				//Get right flanking genes and distances (EDIT: Don't need to anymore)

				//Instead, let's just load the record.
				rec.effect = eff;
				rec.lflank = posHit.getLFLANK_anno();
				rec.ldist = posHit.getLDIST_anno();
				rec.rflank = endHit.getRFLANK_anno();
				rec.rdist = endHit.getRDIST_anno();
				
				return rec;
			}
			else
			{
				//They are straddling at least one gene.
				List<Gene> inGenes = cGenes.getGenes(posHit.index + 1, endHit.index); //hopefully, that gets them all!
				
				//See if there are any coding genes in there. If so, set effect as exonic.
				//If not, set effect ncRNA
				GeneFunc eff = GeneFunc.NCRNA;
				for (Gene g : inGenes)
				{
					if (!g.is_ncRNA())
					{
						eff = GeneFunc.EXONIC;
						break;
					}
				}
				//sv.addInfoField(eff.toString(), INFODEF_INFO_GFUNC);
				int igSize = inGenes.size();
				String[] genes = new String[igSize];
				for (int i = 0; i < igSize; i++){
					genelist.add(inGenes.get(i));
					genes[i] = inGenes.get(i).getName();
				}
				//sv.addInfoField(genes, INFODEF_INFO_GENES);
				
				rec.effect = eff;
				rec.lflank = posHit.getLFLANK_anno();
				rec.ldist = posHit.getLDIST_anno();
				rec.rflank = endHit.getRFLANK_anno();
				rec.rdist = endHit.getRDIST_anno();
				
				return rec;
			}	
		}
		
		//If it passes the intergenic check, then at least one breakend is inside a gene
		
		//Get list of all genes in between breakends.
		List<Gene> inGenes = null;
		if(posHit.location.isIntergenic()) inGenes = cGenes.getGenes(posHit.index + 1, endHit.index);
		else inGenes = cGenes.getGenes(posHit.index, endHit.index);
		if (inGenes == null) inGenes = new LinkedList<Gene>();
		
		//Get any genes overlapping the pos breakend gene that might have been missed.
		List<Integer> stOverlap = cGenes.getOverlappingIndicies(posHit.index);
		if (stOverlap != null && !stOverlap.isEmpty())
		{
			for (Integer n : stOverlap)
			{
				if (n >= posHit.index) continue;
				Gene g = cGenes.getGene(n);
				if (g.getTranscriptEnd() > st) inGenes.add(g);
			}	
		}
		
		//Determine the most detrimental effect.
		GeneFunc eff = posHit.location;
		if (endHit.location.getPriority() < eff.getPriority()) eff = endHit.location;
		if (eff != GeneFunc.getHighestPriority())
		{
			//One or neither is intergenic, and neither is exonic.
			//Here, we need to see if we can upgrade its effect urgency by 
				//seeing if the region spans anything higher priority than either breakend
			
			if (!inGenes.isEmpty())
			{
				//See if whole genes or gene parts are affected
				
				for (Gene g : inGenes)
				{
					boolean posIn = g.positionInGene(st);
					boolean endIn = g.positionInGene(ed);
					
					if (!posIn && !endIn)
					{
						//Assumed to span this gene.
						if (g.is_ncRNA() && eff.getPriority() > GeneFunc.NCRNA.getPriority()) eff = GeneFunc.NCRNA;
						else eff = GeneFunc.EXONIC;
					}
					else if (posIn && endIn)
					{
						//Assumed both in gene, and neither is in an exon (effect would be EXONIC, then)
						if (posHit.location == GeneFunc.UTR3) continue;
						if (endHit.location == GeneFunc.UTR5) continue;
						if (g.is_ncRNA() && eff.getPriority() > GeneFunc.NCRNA.getPriority()) eff = GeneFunc.NCRNA;
						else
						{
							int pnext = g.nextExon(st);
							int enext = g.nextExon(ed);
							if (pnext != enext) eff = GeneFunc.EXONIC;
							//else no upgrade possible - breakend would have to fall IN the region
						}
						
					}
					else if (posIn && !endIn)
					{
						//Is pos before the last exon?
						if (g.is_ncRNA() && eff.getPriority() > GeneFunc.NCRNA.getPriority()) eff = GeneFunc.NCRNA;
						else if (posHit.location != GeneFunc.UTR3) eff = GeneFunc.EXONIC;
					}
					else if (!posIn && endIn)
					{
						//Is end after first exon?
						if (g.is_ncRNA() && eff.getPriority() > GeneFunc.NCRNA.getPriority()) eff = GeneFunc.NCRNA;
						else if (endHit.location != GeneFunc.UTR5) eff = GeneFunc.EXONIC;
					}
					if (eff == GeneFunc.getHighestPriority()) break;
				}
			}
			else
			{
				//This shouldn't happen... getGenes should return SOMETHING as long as stInd <= edInd
				System.err.println("GeneSet.annotateRegion || Error: Invalid non-intergenic region?");
				System.err.println("GeneSet.annotateRegion || POS: " + st + " (" + posHit.index + ") || END: " + ed + " (" + endHit.index + ")");
				return null;
			}	
		}
		
		
		//Add gene list and effect to the record...
		rec.loadGeneList(inGenes);
		genelist.addAll(inGenes);
		rec.effect = eff;
		
		//Note any flanking if needed...
		if (posHit.location.isIntergenic())
		{
			rec.lflank = posHit.getLFLANK_anno();
			rec.ldist = posHit.getLDIST_anno();
		}
		if (endHit.location.isIntergenic())
		{
			rec.rflank = endHit.getRFLANK_anno();
			rec.rdist = endHit.getRDIST_anno();
		}
		
		return rec;
	}
	
	private AnnoRecord annotatePosition(Contig c, int pos, List<Gene> genelist)
	{
		if (c == null) return null;
		ChromSet cGenes = genemap.get(c);
		if (cGenes == null) return null;
		if (!cGenes.hasGenes())
		{
			AnnoRecord rec = new AnnoRecord();
			rec.effect = GeneFunc.INTERGENIC;
			return rec;
		}
		
		HitRecord hit = cGenes.search(pos, true, true);
		
		AnnoRecord rec = new AnnoRecord();
		rec.effect = hit.location;
		
		if (hit.location.isIntergenic())
		{
			rec.lflank = hit.getLFLANK_anno();
			rec.ldist = hit.getLDIST_anno();
			rec.rflank = hit.getRFLANK_anno();
			rec.rdist = hit.getRDIST_anno();
		}
		else
		{
			List<Integer> ogenes = cGenes.getOverlappingIndicies(hit.lindex);
			List<Gene> myGenes = new LinkedList<Gene>();
			myGenes.add(cGenes.getGene(hit.lindex));
			if (ogenes != null)
			{
				for (Integer n : ogenes) myGenes.add(cGenes.getGene(n));
			}
			genelist.addAll(myGenes);
			rec.loadGeneList(myGenes);
		}
		
		return rec;
	}
	
	private AnnoRecord annotatePointPair(Contig c1, int pos1, Contig c2, int pos2, List<Gene> genelist)
	{
		if (c1 == null) return null;
		if (c2 == null) return null;
		ChromSet chrom1 = genemap.get(c1);
		if (chrom1 == null) return null;
		ChromSet chrom2 = genemap.get(c2);
		if (chrom2 == null) return null;
		
		boolean c1genes = chrom1.hasGenes();
		boolean c2genes = chrom2.hasGenes();
		
		HitRecord hit1 = null;
		HitRecord hit2 = null;
		
		if (!c1genes)
		{
			//AnnoRecord rec = new AnnoRecord();
			//rec.effect = GeneFunc.INTERGENIC;
			//return rec;
			hit1 = new HitRecord();
			hit1.location = GeneFunc.INTERGENIC;
		}
		else hit1 = chrom1.search(pos1, true, true);
		if (!c2genes)
		{
			//AnnoRecord rec = new AnnoRecord();
			//rec.effect = GeneFunc.INTERGENIC;
			//return rec;
			hit2 = new HitRecord();
			hit2.location = GeneFunc.INTERGENIC;
		}
		else hit2 = chrom2.search(pos2, true, true);
		
		AnnoRecord rec = new AnnoRecord();
		GeneFunc eff = hit1.location;
		if (hit2.location.getPriority() < eff.getPriority()) eff = hit2.location;
		rec.effect = eff;
		
		List<Gene> myGenes = new LinkedList<Gene>();
		
		if(c1genes)
		{
			if (hit1.location.isIntergenic())
			{
				rec.lflank = hit1.getLFLANK_anno();
				rec.ldist = hit1.getLDIST_anno();
			}
			else
			{
				List<Integer> ogenes = chrom1.getOverlappingIndicies(hit1.lindex);
				myGenes.add(chrom1.getGene(hit1.lindex));
				if(ogenes != null)
				{
					for (Integer n : ogenes) myGenes.add(chrom1.getGene(n));
				}
				//rec.loadGeneList(myGenes);
			}	
		}
		
		if(c2genes)
		{
			if (hit2.location.isIntergenic())
			{
				rec.rflank = hit2.getRFLANK_anno();
				rec.rdist = hit2.getRDIST_anno();
			}
			else
			{
				List<Integer> ogenes = chrom2.getOverlappingIndicies(hit2.lindex);
				myGenes.add(chrom2.getGene(hit2.lindex));
				if(ogenes != null)
				{
					for (Integer n : ogenes) myGenes.add(chrom2.getGene(n));
				}
				//rec.loadGeneList(myGenes);
			}	
		}
		genelist.addAll(myGenes);
		rec.loadGeneList(myGenes);
		
		return rec;
	}
	
	private void annotateSV_CNV(StructuralVariant sv, List<Gene> genelist)
	{
		//Assumed to be on the same chromosome!
		
		AnnoRecord rec = annotateRegion(sv.getChromosome(), sv.getCIPosition(false, false, false), sv.getCIPosition(true, false, true), genelist);
		
		annotateVariant(sv, rec);
	}
	
	private void annotateSV_BND_Single(StructuralVariant sv, List<Gene> genelist)
	{
		if (sv == null) return;
		Contig c = sv.getChromosome();
		int pos = sv.getPosition();
		int st = sv.getCIPosition(false, false, false);
		int ed = sv.getCIPosition(false, false, true);
		
		AnnoRecord rec = null;
		
		int CIsz = ed - st;
		if (CIsz == 0) rec = annotatePosition(c, pos, genelist);
		else rec = annotateRegion(c, st, ed, genelist);
		
		annotateVariant(sv, rec);
	}
	
	private void annotateSV_BND_Pair(BreakendPair bp, List<Gene> genelist)
	{
		StructuralVariant v1 = bp.getBNDVariant(false);
		StructuralVariant v2 = bp.getBNDVariant(true);
		
		annotateSV_BND_Single(v1, genelist);
		annotateSV_BND_Single(v2, genelist);
	}
	
	private void annotateSV_INS(StructuralVariant sv, List<Gene> genelist)
	{
		if (sv == null) return;
		Contig c = sv.getChromosome();
		int pos = sv.getPosition();
		int st = sv.getCIPosition(false, false, false);
		int ed = sv.getCIPosition(false, false, true);
		
		AnnoRecord rec = null;
		
		int CIsz = ed - st;
		if (CIsz == 0) rec = annotatePosition(c, pos, genelist);
		else rec = annotateRegion(c, st, ed, genelist);
		
		annotateVariant(sv, rec);
	}
	
	private void annotateSV_INV(StructuralVariant sv, List<Gene> genelist)
	{
		if (sv == null) return;
		Contig c = sv.getChromosome();
		int pos = sv.getPosition();
		int end = sv.getEndPosition();
		AnnoRecord rec = null;
		
		if (sv.isImprecise())
		{
			int pSt = sv.getCIPosition(false, false, false);
			int pEd = sv.getCIPosition(false, false, true);
			int eSt = sv.getCIPosition(true, false, false);
			int eEd = sv.getCIPosition(true, false, true);
			
			AnnoRecord rec1 = null;
			AnnoRecord rec2 = null;
			if (pEd != pSt) rec1 = annotateRegion(c, pSt, pEd, genelist);
			else rec1 = annotatePosition(c, pos, genelist);
			if (eEd != eSt) rec2 = annotateRegion(c, eSt, eEd, genelist);
			else rec2 = annotatePosition(c, end, genelist);
			
			rec = new AnnoRecord();
			rec.effect = rec1.effect;
			if (rec2.effect.getPriority() < rec.effect.getPriority()) rec.effect = rec2.effect;
			
			int ngenes = 0;
			if(rec1.genes != null) ngenes += rec1.genes.length;
			if(rec2.genes != null) ngenes += rec2.genes.length;
			
			if (ngenes > 0)
			{
				List<String> geneSet = new ArrayList<String>(ngenes);
				if(rec1.genes != null)
				{
					for (int i = 0; i < rec1.genes.length; i++){
						if (!geneSet.contains(rec1.genes[i])) geneSet.add(rec1.genes[i]);
					}	
				}
				if(rec2.genes != null)
				{
					for (int i = 0; i < rec2.genes.length; i++){
						if (!geneSet.contains(rec2.genes[i])) geneSet.add(rec2.genes[i]);
					}
				}
				rec.genes = new String[geneSet.size()];
				int i = 0;
				for (String g : geneSet)
				{
					rec.genes[i] = g;
					i++;
				}
			}
			
			if (rec1.lflank != null) rec.lflank = rec1.lflank;
			if (rec1.ldist != null) rec.ldist = rec1.ldist;
			if (rec2.rflank != null) rec.rflank = rec2.rflank;
			if (rec2.rdist != null) rec.rdist = rec2.rdist;
		}
		else
		{
			rec = annotatePointPair(c, pos, c, end, genelist);
		}
		
		annotateVariant(sv, rec);
	}
	
	private void annotateSV_TRA(Translocation tra, List<Gene> genelist)
	{
		if (tra == null) return;
		Contig c1 = tra.getChromosome1();
		Contig c2 = tra.getChromosome2();
		int pos = tra.getPosition();
		int end = tra.getEndPosition();
		AnnoRecord rec = null;
		
		if (tra.isImprecise())
		{
			int pSt = tra.getCIPosition(false, false, false);
			int pEd = tra.getCIPosition(false, false, true);
			int eSt = tra.getCIPosition(true, false, false);
			int eEd = tra.getCIPosition(true, false, true);
			
			AnnoRecord rec1 = null;
			AnnoRecord rec2 = null;
			if (pEd != pSt) rec1 = annotateRegion(c1, pSt, pEd, genelist);
			else rec1 = annotatePosition(c1, pos, genelist);
			if (eEd != eSt) rec2 = annotateRegion(c2, eSt, eEd, genelist);
			else rec2 = annotatePosition(c2, end, genelist);
			
			rec = new AnnoRecord();
			rec.effect = rec1.effect;
			if (rec2.effect.getPriority() < rec.effect.getPriority()) rec.effect = rec2.effect;
			
			List<String> geneSet = new ArrayList<String>(rec1.genes.length + rec2.genes.length);
			for (int i = 0; i < rec1.genes.length; i++){
				if (!geneSet.contains(rec1.genes[i])) geneSet.add(rec1.genes[i]);
			}
			for (int i = 0; i < rec2.genes.length; i++){
				if (!geneSet.contains(rec2.genes[i])) geneSet.add(rec2.genes[i]);
			}
			
			rec.genes = new String[geneSet.size()];
			int i = 0;
			for (String g : geneSet)
			{
				rec.genes[i] = g;
				i++;
			}
			
			if (rec1.lflank != null) rec.lflank = rec1.lflank;
			if (rec1.ldist != null) rec.ldist = rec1.ldist;
			if (rec2.rflank != null) rec.rflank = rec2.rflank;
			if (rec2.rdist != null) rec.rdist = rec2.rdist;
		}
		else
		{
			rec = annotatePointPair(c1, pos, c2, end, genelist);
		}
		
		annotateVariant(tra, rec);
	}
	
	private void annotateVariant(Variant v, AnnoRecord rec)
	{
		if (v == null){
			System.err.println("GeneSet.annotateVariant || ERROR: cannot have null variant!");
			return;
		}
		if (rec == null){
			System.err.println("GeneSet.annotateVariant || ERROR: Variant " + v.getVarID() + " could not be annotated. Null record...");
			return;
		}
		if (rec.genes != null) v.addInfoField(rec.genes, INFODEF_INFO_GENES);
		if (rec.effect != null) v.addInfoField(rec.effect.toString(), INFODEF_INFO_GFUNC);
		if (rec.lflank != null) v.addInfoField(rec.lflank, INFODEF_INFO_LFLANK);
		if (rec.ldist != null) v.addInfoField(INFODEF_INFO_LDIST.getKey(), rec.ldist);
		if (rec.rflank != null) v.addInfoField(rec.rflank, INFODEF_INFO_RFLANK);
		if (rec.rdist != null) v.addInfoField(INFODEF_INFO_RDIST.getKey(), rec.rdist);
	}
	
	/**
	 * Compare the location of a variant with genes in this set, and annotate
	 * the variant (via INFO fields) with information about any overlapping or
	 * nearby genes.
	 * <br>Variant will remain un-annotated if lookup fails.
	 * @param v Variant to annotate.
	 */
	public List<Gene> annotateVariant(Variant v)
	{
		if (v == null) return null;
		if (v instanceof StructuralVariant)
		{
			StructuralVariant sv = (StructuralVariant)v;
			return annotateStructuralVariant(sv, true);
		}
		List<Gene> genelist = new LinkedList<Gene>();
		Contig c = v.getChromosome();
		int pos = v.getPosition();
		int len = v.getLargestAbsoluteLength();
		
		AnnoRecord rec = null;
		if (len < 2) rec = annotatePosition(c, pos, genelist);
		else rec = annotateRegion(c, pos, pos + (len - 1), genelist);
		
		annotateVariant(v, rec);
		Collections.sort(genelist);
		return genelist;
	}
	
	/**
	 * Compare the location of a structural variant with genes in this set, and annotate
	 * the variant (via INFO fields) with information about any overlapping or
	 * nearby genes. Approach to annotation will vary based upon SV type and precision.
	 * <br>Variant will remain un-annotated if lookup fails.
	 * @param sv Structural variant to annotate.
	 * @param inversionRegion If true, treat inversions as a region (ie. note all genes between breakends).
	 * If false, treat inversions as a pair of breakends (ie. ignore all genes outside of breakend CIs).
	 */
	public List<Gene> annotateStructuralVariant(StructuralVariant sv, boolean inversionRegion)
	{
		if (sv == null) return null;
		//System.err.println(Thread.currentThread().getName() + " || GeneSet.annotateStructuralVariant || Called - " + sv.getVarID() + " [" + sv.getChromosomeName() + ":" + sv.getPosition() + "-" + sv.getEndPosition() + "]");
		List<Gene> genelist = new LinkedList<Gene>();
		try {
			switch(sv.getType())
			{
			case BED_REGION: annotateSV_CNV(sv, genelist); break;
			case BND:
				if (sv instanceof BreakendPair) annotateSV_BND_Pair((BreakendPair)sv, genelist);
				else annotateSV_BND_Single(sv, genelist);
				break;
			case CNV: annotateSV_CNV(sv, genelist); break;
			case DEL: annotateSV_CNV(sv, genelist); break;
			case DELME: annotateSV_CNV(sv, genelist); break;
			case DUP: annotateSV_CNV(sv, genelist); break;
			case INS: annotateSV_INS(sv, genelist); break;
			case INSME: annotateSV_INS(sv, genelist); break;
			case INV: 
				if (inversionRegion) annotateSV_CNV(sv, genelist);
				else annotateSV_INV(sv, genelist); 
				break;
			case OTHER: annotateSV_CNV(sv, genelist); break;
			case TANDEM: annotateSV_CNV(sv, genelist); break;
			case TRA:
				if (sv instanceof Translocation) annotateSV_TRA((Translocation)sv, genelist);
				else annotateSV_BND_Single(sv, genelist);
				break;
			default: return null;
			}	
		}
		catch (Exception e)
		{
			System.err.println("ERROR: Exception caught when annotating the following variant:");
			System.err.println(sv.getVarID() + " | " + sv.toString());
			e.printStackTrace();
		}
		Collections.sort(genelist);
		return genelist;
	}
	
	/* --- Access --- */
	
	/**
	 * Add a collection of genes to this set.
	 * <br>Genes are requested in multiples because internal lists are resorted
	 * and reindexed every time there is an addition.
	 * @param genes Collection of genes to add to set.
	 */
	public void addGenes(Collection<Gene> genes)
	{
		Map<Contig, List<Gene>> addcoll = new HashMap<Contig, List<Gene>>();
		Set<Contig> cset = genemap.keySet();
		for (Contig c : cset){
			//System.out.println("|DEBUG| GeneSet.addGenes || Contig found: " + c.getUDPName());
			addcoll.put(c, new LinkedList<Gene>());
		}
		
		for (Gene g : genes)
		{
			List<Gene> glist = addcoll.get(g.getChromosome());
			if (glist == null) System.err.println("GeneSet.addGenes || WARNING: Contig requested, but not found: " + g.getChromosome());
			else glist.add(g);
			
			//System.out.println("|DEBUG| GeneSet.addGenes || Genes compiled for contig ");
		}
		
		//System.err.println("|DEBUG| GeneSet.addGenes || Genes sorted by contig.");
		for (Contig c : cset)
		{
			ChromSet cGenes = genemap.get(c);
			if (cGenes == null) continue;
			//int n = addcoll.get(c).size();
			cGenes.addGenes(addcoll.get(c));
			//System.err.println("|DEBUG| GeneSet.addGenes || " + n + " genes added to contig " + c.getUDPName());
			addcoll.remove(c);
		}
		
		//System.err.println("|DEBUG| GeneSet.addGenes || Add complete.");
	}
	
	/**
	 * Get a sorted list of all genes in this set.
	 * @return Sorted list of all genes.
	 */
	public List<Gene> getAllGenes()
	{
		Collection<ChromSet> csets = genemap.values();
		List<Gene> glist = new LinkedList<Gene>();
		for (ChromSet cs : csets)
		{
			glist.addAll(cs.getAllGenes());
		}
		Collections.sort(glist);
		
		return glist;
	}
	
	/**
	 * Get a random gene from this set.
	 * @return A random gene.
	 */
	public Gene getRandomGene()
	{
		List<Gene> genes = getAllGenes();
		Collections.shuffle(genes);
		return genes.get(0);
	}
	
	/**
	 * Get a list of all genes that have names matching the query (not case-sensitive).
	 * @param gene_name Name of gene to query.
	 * @return List of matching genes - empty list if none found.
	 */
	public List<Gene> getGeneByName(String gene_name)
	{
		List<Gene> allgenes = getAllGenes();
		List<Gene> results = new LinkedList<Gene>();
		for (Gene g : allgenes)
		{
			if (g.getName().equalsIgnoreCase(gene_name)) results.add(g);
		}
		Collections.sort(results);
		return results;
	}
	
	/* --- Load Standard --- */
	
	public static final String PACKAGEPATH_36 = "resources/ncbi36_refSeq.gbgd";
	public static final String PACKAGEPATH_37 = "resources/grch37_refSeq.gbgd";
	public static final String PACKAGEPATH_38 = "resources/grch38_refSeq.gbgd";
	
	public static final String DEBUGPATH_37 = "C:\\Users\\Blythe\\grch37_refSeq.gbgd";
	public static final String DEBUGPATH_38 = "C:\\Users\\Blythe\\grch38_refSeq.gbgd";
	
	private static Map<String, String> buildDict;
	private static Map<String, String> packPathMap;
	private static Map<String, String> debugPathMap;
	
	private static Map<String, GeneSet> loadedMap;
	
	private static void populateLoadMaps()
	{
		buildDict = new HashMap<String, String>();
		packPathMap = new HashMap<String, String>();
		debugPathMap = new HashMap<String, String>();
		loadedMap = new HashMap<String, GeneSet>();
		
		buildDict.put("hg18", "hg18");
		buildDict.put("ncbi36", "hg18");
		buildDict.put("hg19", "hg19");
		buildDict.put("grch37", "hg19");
		buildDict.put("hg38", "hg38");
		buildDict.put("grch38", "hg38");
		
		packPathMap.put("hg18", PACKAGEPATH_36);
		packPathMap.put("hg19", PACKAGEPATH_37);
		packPathMap.put("hg38", PACKAGEPATH_38);
		
		debugPathMap.put("hg19", DEBUGPATH_37);
		debugPathMap.put("hg38", DEBUGPATH_38);
	}
	
	private static void loadRefGeneIntoMap(String buildname)
	{	
		//Try to load genome build
		GenomeBuild gb = GenomeBuild.loadStandardBuild(buildname);
		if (gb == null) throw new IllegalArgumentException();
		
		//Try package local path
		String packPath = packPathMap.get(buildname);
		if (packPath == null) throw new IllegalArgumentException();
		GeneSet gs = null;
		InputStream is = GenomeBuild.class.getResourceAsStream(packPath);
		try 
		{
			gs = new GeneSet(is, gb, true);
			is.close();
			//gs = new GeneSet(packPath, gb, true);
		} 
		catch (UnsupportedFileTypeException e) 
		{
			System.err.println("GeneSet.loadRefGeneIntoMap || ERROR! " + buildname + " (" + packPath + ") could not be parsed!!");
			e.printStackTrace();
		}
		catch (IOException e) 
		{
			System.err.println("GeneSet.loadRefGeneIntoMap || ERROR! " + buildname + " (" + packPath + ") could not be read!!");
			e.printStackTrace();
		}
		if (gs != null)
		{
			loadedMap.put(buildname, gs);
			return;
		}
		
		//Try debug path
		System.err.println("GeneSet.loadRefGeneIntoMap || Standard path for " + buildname + " (" + packPath + ") did not work. Trying debug path...");
		String debugPath = debugPathMap.get(buildname);
		if (debugPath == null) throw new IllegalArgumentException();
		try 
		{
			gs = new GeneSet(debugPath, gb, true);
		} 
		catch (UnsupportedFileTypeException e) 
		{
			System.err.println("GeneSet.loadRefGeneIntoMap || ERROR! " + buildname + " (" + debugPath + ") could not be parsed!!");
			e.printStackTrace();
		}
		catch (IOException e) 
		{
			System.err.println("GeneSet.loadRefGeneIntoMap || ERROR! " + buildname + " (" + debugPath + ") could not be read!!");
			e.printStackTrace();
		}
		if (gs != null)
		{
			loadedMap.put(buildname, gs);
			return;
		}
		
	}
	
	private static void loadRefGeneIntoMap(GenomeBuild gb)
	{	
		//Try to load genome build
		if (gb == null) throw new IllegalArgumentException();
		String buildname = gb.getBuildName();
		buildname = buildDict.get(buildname);
		if (buildname == null) return;
		
		//Try package local path
		String packPath = packPathMap.get(buildname);
		if (packPath == null) throw new IllegalArgumentException();
		GeneSet gs = null;
		InputStream is = GenomeBuild.class.getResourceAsStream(packPath);
		try 
		{
			gs = new GeneSet(is, gb, true);
			is.close();
		} 
		catch (UnsupportedFileTypeException e) 
		{
			System.err.println("GeneSet.loadRefGeneIntoMap || ERROR! " + buildname + " (" + packPath + ") could not be parsed!!");
			e.printStackTrace();
		}
		catch (IOException e) 
		{
			System.err.println("GeneSet.loadRefGeneIntoMap || ERROR! " + buildname + " (" + packPath + ") could not be read!!");
			e.printStackTrace();
		}
		if (gs != null)
		{
			loadedMap.put(buildname, gs);
			return;
		}
		
		//Try debug path
		System.err.println("GeneSet.loadRefGeneIntoMap || Standard path for " + buildname + " (" + packPath + ") did not work. Trying debug path...");
		String debugPath = debugPathMap.get(buildname);
		if (debugPath == null) throw new IllegalArgumentException();
		try 
		{
			gs = new GeneSet(debugPath, gb, true);
		} 
		catch (UnsupportedFileTypeException e) 
		{
			System.err.println("GeneSet.loadRefGeneIntoMap || ERROR! " + buildname + " (" + debugPath + ") could not be parsed!!");
			e.printStackTrace();
		}
		catch (IOException e) 
		{
			System.err.println("GeneSet.loadRefGeneIntoMap || ERROR! " + buildname + " (" + debugPath + ") could not be read!!");
			e.printStackTrace();
		}
		if (gs != null)
		{
			loadedMap.put(buildname, gs);
			return;
		}
		
	}
	
	/**
	 * Get the path to the database file relative to the hospelhornbg_genomeBuild
	 * package. Use this for reading file from the JAR.
	 * <br>Obviously, problems arise if the database is not kept in hospelhornbg_genomeBuild.resources
	 * with its standard name!
	 * @param buildname Valid name of genome build to load. Case insensitive.
	 * @return JAR path of refSeq database for requested genome build, if present.
	 */
	public static String getStandardDB_packagePath(String buildname)
	{
		if (loadedMap == null) populateLoadMaps();
		buildname = buildname.toLowerCase();
		
		String bName = buildDict.get(buildname);
		if (bName == null) return null;
		
		return packPathMap.get(bName);
	}
	
	/**
	 * Load the refSeq gene database for the specified genome build, if
	 * present in the jar. The loaded build will remain in memory until all threads
	 * terminate or is explicitly unloaded - so keep this in mind!
	 * @param buildname Valid name of genome build to load. Case insensitive.
	 * @return refSeq database as a GeneSet object if the build name is recognized,
	 * null if the build or a refSeq database for the build could not be found.
	 */
	public static GeneSet loadRefGene(String buildname)
	{
		if (loadedMap == null) populateLoadMaps();
		buildname = buildname.toLowerCase();
		
		String bName = buildDict.get(buildname);
		if (bName == null) return null;
		
		//Look to see if table is already loaded...
		GeneSet gs = loadedMap.get(bName);
		if (gs != null) return gs;
		
		loadRefGeneIntoMap(bName);
		gs = loadedMap.get(bName);
		if (gs != null) return gs;
		
		return null;
	}
	
	/**
	 * Load the refSeq gene database for the specified genome build, if
	 * present in the jar. The loaded build will remain in memory until all threads
	 * terminate or is explicitly unloaded - so keep this in mind!
	 * @param build GenomeBuild to load accompanying refGene db of.
	 * @return refSeq database as a GeneSet object if the build name is recognized,
	 * null if the build or a refSeq database for the build could not be found.
	 */
	public static GeneSet loadRefGene(GenomeBuild build)
	{
		if (loadedMap == null) populateLoadMaps();
		String buildname = build.getBuildName();
		
		String bName = buildDict.get(buildname);
		if (bName == null) return null;
		
		//Look to see if table is already loaded...
		GeneSet gs = loadedMap.get(bName);
		if (gs != null) return gs;
		
		loadRefGeneIntoMap(build);
		gs = loadedMap.get(bName);
		if (gs != null) return gs;
		
		return null;
	}
	
	/**
	 * Unload a standard refSeq database from memory. If the build name
	 * is not recognized or the database was never loaded, this function
	 * will do nothing.
	 * @param buildname Valid name of genome build to load. Case insensitive.
	 */
	public static void unloadRefGene(String buildname)
	{
		if (loadedMap == null) return;
		buildname = buildname.toLowerCase();
		
		String bName = buildDict.get(buildname);
		if (bName == null) return;
		
		loadedMap.remove(bName);
	}
	
	/**
	 * Unload all standard refSeq databases that have been loaded from memory.
	 * This will only set the GeneSet class mappings to null; if any other
	 * part of the program is still referencing any given GeneSet, the object
	 * will NOT be freed for the Java garbage collector to have at.
	 */
	public static void unloadAllRefGene()
	{
		if (loadedMap == null) return;
		loadedMap.clear();
	}
	
	
}