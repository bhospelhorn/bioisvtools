package hospelhornbg_svdb;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ConcurrentMap;

import hospelhornbg_bioinformatics.Genotype;
import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.VCFReadStreamer;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.Gene;
import hospelhornbg_genomeBuild.GeneFunc;
import hospelhornbg_genomeBuild.GeneSet;
import hospelhornbg_genomeBuild.GenomeBuild;
import hospelhornbg_segregation.FamilyMember;
import hospelhornbg_segregation.Population;
import hospelhornbg_svdb.DBVariantTable.GeneHitCounter;
import waffleoRai_Utils.Arunnable;

public class SQLVTVarAdder {
	
	//Read and parse variants from VCF
	//Check DB for matching variants
	//Perform merge or create new dbvar
	//Update variant table
	
	//Update sample geno table (maybe after all added?)
	
	/*--- Constants ---*/
	
	public static final int MAX_Q_SIZE_PARSED = 256;
	
	public static final String INFO_MATCHED_UID_KEY = "svdb_uid";
	
	/*--- Instance Variables ---*/
	
	private boolean verbose;
	
	private GenomeBuild genome;
	private GeneSet genes;
	private SQLVariantTable table;
	
	private ConcurrentMap<String, FamilyMember> sampleMap;
	
	private ConcurrentLinkedQueue<StructuralVariant> parsedQ;
	private ConcurrentLinkedQueue<MergePair> mergeQ;
	private ConcurrentLinkedQueue<StructuralVariant> newQ;
	
	private Object regGetLock;
	private volatile StructuralVariant nextSV;
	private volatile Collection<DBVariant> nextRegion;
	
	private volatile boolean error_kill;
	private volatile boolean reader_done;
	private volatile boolean getter_done;
	private volatile boolean merger_done;
	
	private VCFReader parseRunner;
	private RegGetter regGetter;
	private VarMerger varMerger;
	private NewVarAdder varAdder;
	private VarUpdater varUpdater;
	
	private volatile long c_parsed;
	private volatile long c_reg;
	private volatile long c_merged;
	private volatile long c_added;
	private volatile long c_updated;
	
	/*--- Other Objects ---*/
	
	private static class MergePair
	{
		public StructuralVariant newSV;
		public DBVariant oldSV;
		
		public MergePair(StructuralVariant sv, DBVariant dbv)
		{
			newSV = sv;
			oldSV = dbv;
		}
	}
	
	/*--- Runner Objects ---*/
	
	public class VCFReader extends Arunnable
	{
		//Blocks when queue is full.
		//Can be interrupted by the region checker
		
		private Contig lastChrom;
		
		private boolean ignore_translocations;
		private VCFReadStreamer reader;
		private Iterator<StructuralVariant> itr;
		
		public VCFReader(String vcfPath, boolean no_tra)
		{
			super.setName("vcf_parser");
			super.sleeps = true;
			super.delay = 0;
			super.sleeptime = 10;
			ignore_translocations = no_tra;
			
			reader = new VCFReadStreamer(vcfPath, genome);
		}
		
		public void open() throws IOException
		{
			reader.open();
		}
		
		public void close() throws IOException
		{
			reader.close();
		}

		@Override
		public void doSomething() 
		{
			try{
			if(itr == null)
			{
				//Get the iterator
				itr = reader.getSVIterator();
			}
			
			if(!itr.hasNext())
			{
				//Close and kys
				try {close();}
				catch(Exception e)
				{
					e.printStackTrace();
				}
				
				reader_done = true;
				this.requestTermination();
			}
			
			//Read loop
			while(itr.hasNext() && (parsedQ.size() <= MAX_Q_SIZE_PARSED))
			{
				StructuralVariant sv = itr.next();
				if(verbose)
				{
					Contig c = sv.getChromosome();
					if(lastChrom == null)
					{
						lastChrom = c;
						System.err.println("Now parsing variants on chrom " + c.getUDPName() + "...");
					}
					else
					{
						if(!lastChrom.equals(c))
						{
							lastChrom = c;
							System.err.println("Now parsing variants on chrom " + c.getUDPName() + "...");
						}
					}
				}
				if(ignore_translocations && ((sv.getType() == SVType.BND) || (sv.getType() == SVType.TRA))) continue;
				parsedQ.add(itr.next());
				c_parsed++;
			}
			
			if(regGetter.anyThreadsAlive()){
				//System.err.println(Thread.currentThread().getName() + " || Parser - Interrupting RegGetter");
				regGetter.interruptThreads();
			}
			}
			catch(Exception e){killAll(e);}
	}
		
	}

	public class RegGetter extends Arunnable
	{
		
		//private boolean ignore_translocations;
		private Contig lastChrom;
		
		public RegGetter()
		{
			super.setName("region_checker");
			super.sleeps = true;
			super.delay = 0;
			super.sleeptime = 10;
			//ignore_translocations = no_tra;
		}

		@Override
		public void doSomething() 
		{
			try{
			//System.err.println(Thread.currentThread().getName() + " || RegGetter - doSomething() called!");
			super.disableInterrupts();
			while(nextRegion == null && !parsedQ.isEmpty())
			{
				StructuralVariant sv = parsedQ.poll();
				Contig c = sv.getChromosome();
				int st = sv.getCIPosition(false, false, false);
				int ed = sv.getCIPosition(true, false, true);
				SVType type = sv.getType();
				//System.err.println(Thread.currentThread().getName() + " || RegGetter - looking for variants in region!");
				Collection<DBVariant> results = table.getVariantsInRegion(c, st, ed, type);
				//System.err.println(Thread.currentThread().getName() + " || RegGetter - variant search done!");
				
				//Note
				if(verbose)
				{
					if(lastChrom == null)
					{
						lastChrom = c;
						System.err.println("Now processing variants on chrom " + c.getUDPName() + "...");
					}
					else
					{
						if(!lastChrom.equals(c))
						{
							lastChrom = c;
							System.err.println("Now processing variants on chrom " + c.getUDPName() + "...");
						}
					}
				}
				
				//Set
				synchronized(regGetLock)
				{
					nextSV = sv;
					nextRegion = results;
				}
				
				if(parseRunner.anyThreadsAlive()){
					//System.err.println(Thread.currentThread().getName() + " || RegGetter - Interrupting Parser");
					parseRunner.interruptThreads();
				}
				if(varMerger.anyThreadsAlive()){
					//System.err.println(Thread.currentThread().getName() + " || RegGetter - Interrupting VarMerger");
					varMerger.interruptThreads();
				}
				c_reg++;
			}
			if(reader_done) {getter_done = true; this.requestTermination();}
			super.enableInterrupts();
			//System.err.println(Thread.currentThread().getName() + " || RegGetter - doSomething() returning!");
			}
			catch(Exception e){killAll(e);}
			}
	}
	
	public class VarMerger extends Arunnable
	{
		
		private double leeway;

		public VarMerger(double percentLeeway)
		{
			leeway = percentLeeway;
			
			super.setName("VariantRegionScanner");
			super.sleeps = true;
			super.sleeptime = 100;
		}
		
		@Override
		public void doSomething() 
		{
			try{
			//Looks at the last region search and SV in queue to 
			//determine if new SV can be merged to existing variant
			
			while(nextRegion != null)
			{
				StructuralVariant sv = nextSV;
				Collection<DBVariant> regvars = nextRegion;
				
				synchronized(regGetLock)
				{
					nextSV = null;
					nextRegion = null;
				}
				//Let the reg getter have at the next variant
				if(regGetter.anyThreadsAlive()){
					//System.err.println(Thread.currentThread().getName() + " || VarMerger - Interrupting RegGetter");
					regGetter.interruptThreads();
				}
				
				DBVariant match = null;
				
				for(DBVariant dbv : regvars)
				{
					if(dbv.svIsEquivalent(sv, leeway))
					{
						match = dbv;
						break; //Hopefully won't match to more than one!
					}
				}
				
				if(match != null)
				{
					mergeQ.add(new MergePair(sv, match));
					if(varUpdater.anyThreadsAlive()){
						//System.err.println(Thread.currentThread().getName() + " || VarMerger - Interrupting VarUpdater");
						varUpdater.interruptThreads();
					}
				}
				else
				{
					newQ.add(sv);
					if(varAdder.anyThreadsAlive()){
						//System.err.println(Thread.currentThread().getName() + " || VarMerger - Interrupting NewVarAdder");
						varAdder.interruptThreads();
					}
				}
				c_merged++;
			}
			
			if(getter_done) {merger_done = true; this.requestTermination();}
			}
			catch(Exception e){killAll(e);}
		}
		
	}
	
	public class NewVarAdder extends Arunnable
	{
		
		private Contig lastChrom;
		
		public NewVarAdder()
		{
			super.setName("NewVarAdder");
			super.sleeps = true;
			super.sleeptime = 1000;
		}

		@Override
		public void doSomething() 
		{
			try{
			while(!newQ.isEmpty())
			{
				//Pop!
				StructuralVariant sv = newQ.poll();
				Contig c = sv.getChromosome();
				int stpos = sv.getCIPosition(false, false, false);
				
				//Message
				if(verbose)
				{
					if(lastChrom == null)
					{
						lastChrom = c;
						System.err.println("Now adding variants on chrom " + c.getUDPName() + "...");
					}
					else
					{
						if(!lastChrom.equals(c))
						{
							lastChrom = c;
							System.err.println("Now adding variants on chrom " + c.getUDPName() + "...");
						}
					}
				}
				
				//Generate a new UID
				super.disableInterrupts();
				//System.err.println(Thread.currentThread().getName() + " || NewVarAdder - getting new UID");
				long vuid = table.generateUID(c, stpos);
				//System.err.println(Thread.currentThread().getName() + " || NewVarAdder - getting new UID -- Done!");
				super.enableInterrupts();
				
				//Convert to DBV
				DBVariant dbv = DBVariant.getFromVariant(sv, sv.getVarID());
				dbv.setLongUID(vuid);
				dbv.noteGenes(genes);
				
				//Generate genotype & update pop counts
				VariantGenotype vg = new VariantGenotype(dbv.getLongID());
				List<String> samplelist = new ArrayList<String>(sampleMap.size()+1);
				samplelist.addAll(sampleMap.keySet());
				for(String s : samplelist)
				{
					Genotype g = sv.getSampleGenotype(s);
					if(g == null || (g.isHomozygous() && g.hasAllele(0))) continue;
					//Generate genotype
					FamilyMember mem = sampleMap.get(s);
					SVDBGenotype gt = SVDBGenotype.generateGenotype(mem.getUID(), g, sv);
					vg.addGenotype(gt);
					
					//Population counts...
					boolean hom = gt.isHomozygous();
					dbv.incrementTotalCount();
					if(hom) dbv.incrementHomozygoteCount();
					
					Collection<Population> plist = mem.getPopulationTags();
					for(Population p : plist)
					{
						dbv.incrementTotalCount(p);
						if(hom) dbv.incrementHomozygoteCount(p);
					}
				}

				//Note gene hits
				List<Gene> glist = dbv.getGeneListReference();
				super.disableInterrupts();
				if(glist != null)
				{
					for(Gene g : glist)
					{
						//System.err.println(Thread.currentThread().getName() + " || NewVarAdder - Getting gene hit record: " + g.getGUID());
						GeneHitCounter ghc = table.getGeneHitRecord(g.getGUID());
						//System.err.println(Thread.currentThread().getName() + " || NewVarAdder - Getting gene hit record -- Done! ");
						boolean exonic = (g.getRelativeRegionLocationEffect(dbv.getStartPosition().getStart(), dbv.getEndPosition().getEnd()) == GeneFunc.EXONIC);
						//ghc.total_hits_var++;
						//ghc.total_hits_indiv.addAll(vg.getAllIndividuals());
						ghc.incrementTotalHits_sync();
						ghc.addIndividuals_Total_sync(vg.getAllIndividuals());
						
						if(exonic)
						{
							//ghc.exon_hits_var++;
							//ghc.exon_hits_indiv.addAll(vg.getAllIndividuals());
							ghc.incrementExonHits_sync();
							ghc.addIndividuals_Exon_sync(vg.getAllIndividuals());
						}
					}
					table.tickGHCCDirty();
				}
				
				//Send to table to add	
				//System.err.println(Thread.currentThread().getName() + " || NewVarAdder - Adding variant... ");
				if (!table.addOrUpdateVariant(dbv, vg, true)) {killAll(null); return;}
				//System.err.println(Thread.currentThread().getName() + " || NewVarAdder - Adding variant -- Done!");
				super.enableInterrupts();
				c_added++;
			}
			if(merger_done) this.requestTermination();
			}
			catch(Exception e){killAll(e);}
		}
		
	}
	
	public class VarUpdater extends Arunnable
	{

		private Contig lastChrom;
		
		public VarUpdater()
		{
			super.setName("VarUpdater");
			super.sleeps = true;
			super.sleeptime = 1000;
		}
		
		@Override
		public void doSomething() 
		{
			try{
			while(!mergeQ.isEmpty())
			{
				//Pop
				MergePair mp = mergeQ.poll();
				
				//Progress message
				if(verbose)
				{
					Contig c = mp.newSV.getChromosome();
					if(lastChrom == null)
					{
						lastChrom = c;
						System.err.println("Now updating variants on chrom " + c.getUDPName() + "...");
					}
					else
					{
						if(!lastChrom.equals(c))
						{
							lastChrom = c;
							System.err.println("Now updating variants on chrom " + c.getUDPName() + "...");
						}
					}
				}
				
				//Update variant ends (and genes covered)
				int in_st = mp.newSV.getCIPosition(false, false, false);
				int in_ed = mp.newSV.getCIPosition(true, false, true);
				boolean alteredEnds = false;
				int oldStart = mp.oldSV.getStartPosition().getStart();
				int oldEnd = mp.oldSV.getEndPosition().getEnd();
				if(in_st < mp.oldSV.getStartPosition().getStart())
				{
					mp.oldSV.getStartPosition().setStart(in_st);
					alteredEnds = true;
				}
				if(in_st > mp.oldSV.getStartPosition().getEnd())
				{
					mp.oldSV.getStartPosition().setEnd(in_st);
				}
				if(in_ed < mp.oldSV.getEndPosition().getStart())
				{
					mp.oldSV.getStartPosition().setStart(in_ed);
				}
				if(in_ed > mp.oldSV.getEndPosition().getEnd())
				{
					mp.oldSV.getStartPosition().setEnd(in_ed);
					alteredEnds = true;
				}
				
				Set<Integer> newGenes = new TreeSet<Integer>();
				Map<Integer, GeneFunc> oldposeff = new HashMap<Integer, GeneFunc>();
				if(alteredEnds)
				{
					List<Gene> glist = mp.oldSV.getGeneListReference();
					Set<Integer> oldgenes = new TreeSet<Integer>();
					for(Gene g : glist) {
						GeneFunc poseff = g.getRelativeRegionLocationEffect(oldStart, oldEnd);
						oldgenes.add(g.getGUID());
						oldposeff.put(g.getGUID(), poseff);
					}
					
					mp.oldSV.noteGenes(genes);
					glist = mp.oldSV.getGeneListReference();
					for(Gene g : glist)
					{
						if(!oldgenes.contains(g.getGUID())) newGenes.add(g.getGUID());
					}
				}
				
				//Update genotypes & population counts
				super.disableInterrupts();
				//System.err.println(Thread.currentThread().getName() + " || VarUpdater - getting genotype");
				VariantGenotype vg = table.getGenotype(mp.oldSV.getLongID());
				//System.err.println(Thread.currentThread().getName() + " || VarUpdater - getting genotype -- Done!");
				super.enableInterrupts();
				List<String> samplelist = new ArrayList<String>(sampleMap.size()+1);
				samplelist.addAll(sampleMap.keySet());
				Set<Integer> removed = new TreeSet<Integer>();
				for(String s : samplelist)
				{
					Genotype g = mp.newSV.getSampleGenotype(s);
					if(g == null) continue;
					
					FamilyMember mem = sampleMap.get(s);
					
					//See if already has genotype
					boolean has = vg.hasGenotypeFor(mem.getUID());
					
					//See if needs to add
					g = mp.newSV.getSampleGenotype(s);
					if(g.isHomozygous() && g.hasAllele(0))
					{
						if(!has) continue; //No problems, just continue
						//Otherwise we are presumably switching to homref and we need to remove
						SVDBGenotype gt = vg.getGenotype(mem.getUID());
						boolean hom = gt.isHomozygous();
						mp.oldSV.decrementTotalCount();
						if(hom)mp.oldSV.decrementHomozygoteCount();
						Collection<Population> plist = mem.getPopulationTags();
						for(Population p : plist)
						{
							mp.oldSV.decrementTotalCount(p);
							if(hom) mp.oldSV.decrementHomozygoteCount(p);
						}
						removed.add(mem.getUID());
						continue;
					}
					
					//Add new genotype
					SVDBGenotype gt = SVDBGenotype.generateGenotype(mem.getUID(), g, mp.newSV);
					vg.addGenotype(gt);
					
					//Population counts...
					boolean hom = gt.isHomozygous();
					mp.oldSV.incrementTotalCount();
					if(hom) mp.oldSV.incrementHomozygoteCount();
					
					Collection<Population> plist = mem.getPopulationTags();
					for(Population p : plist)
					{
						mp.oldSV.incrementTotalCount(p);
						if(hom) mp.oldSV.incrementHomozygoteCount(p);
					}
					
				}
				
				//Update gene hit counts
				List<Gene> glist = mp.oldSV.getGeneListReference();
				super.disableInterrupts();
				for(Gene g : glist)
				{
					boolean isnew = newGenes.contains(g.getGUID());
					//System.err.println(Thread.currentThread().getName() + " || VarUpdater - getting gene hit info");
					GeneHitCounter ghc = table.getGeneHitRecord(g.getGUID());
					//System.err.println(Thread.currentThread().getName() + " || VarUpdater - getting gene hit info -- Done!");
					boolean exonic = (g.getRelativeRegionLocationEffect(mp.oldSV.getStartPosition().getStart(), mp.oldSV.getEndPosition().getEnd()) == GeneFunc.EXONIC);
					//if(isnew) ghc.total_hits_var++;
					//ghc.total_hits_indiv.addAll(vg.getAllIndividuals());
					//ghc.total_hits_indiv.removeAll(removed);
					if(isnew) ghc.incrementTotalHits_sync();
					ghc.addIndividuals_Total_sync(vg.getAllIndividuals());
					ghc.removeIndividuals_Total_sync(removed);
					
					if(exonic)
					{
						//if(isnew || oldposeff.get(g.getGUID()) != GeneFunc.EXONIC) ghc.exon_hits_var++;
						//ghc.exon_hits_indiv.addAll(vg.getAllIndividuals());
						//ghc.exon_hits_indiv.removeAll(removed);
						if(isnew || oldposeff.get(g.getGUID()) != GeneFunc.EXONIC) ghc.incrementExonHits_sync();
						ghc.addIndividuals_Exon_sync(vg.getAllIndividuals());
						ghc.removeIndividuals_Exon_sync(removed);
					}
				}
				table.tickGHCCDirty();
				
				//Update in database table
				//System.err.println(Thread.currentThread().getName() + " || VarUpdater - updating variant in db");
				if (!table.addOrUpdateVariant(mp.oldSV, vg, false)){killAll(null); return;}
				//System.err.println(Thread.currentThread().getName() + " || VarUpdater - updating variant in db -- Done!");
				super.enableInterrupts();
				c_updated++;
			}
			if(merger_done) this.requestTermination();
			}
			catch(Exception e){killAll(e);}
		}
		
	}
	
	/*--- Construction ---*/
	
	public SQLVTVarAdder(SQLVariantTable mytable, Map<String, FamilyMember> samples, boolean v)
	{
		verbose = v;
		table = mytable;
		genome = table.getGenomeBuild();
		genes = table.getTranscriptSet();
		
		sampleMap = new ConcurrentHashMap<String, FamilyMember>();
		Set<String> keys = samples.keySet();
		for(String s : keys)
		{
			sampleMap.put(s, samples.get(s));
		}
		
		error_kill = false;
	}
	
	/*--- Running ---*/
	
	public synchronized void start(String vcfpath, boolean noTRA) throws IOException
	{
		table.setThreadlock(true);	
		reader_done = false;
		getter_done = false;
		merger_done = false;
		
		error_kill = false;
		
		c_parsed = 0;
		c_merged = 0;
		c_added = 0;
		c_reg = 0;
		c_updated = 0;
		
		parsedQ = new ConcurrentLinkedQueue<StructuralVariant>();
		mergeQ = new ConcurrentLinkedQueue<MergePair>();
		newQ = new ConcurrentLinkedQueue<StructuralVariant>();
		
		regGetLock = new Boolean(false);
		
		parseRunner = new VCFReader(vcfpath, noTRA);
		regGetter = new RegGetter();
		varMerger = new VarMerger(table.getPercentLeeway());
		varAdder = new NewVarAdder();
		varUpdater = new VarUpdater();
		
		Thread t_parser = new Thread(parseRunner);
		parseRunner.open();
		t_parser.setName("VCFParserThread");
		t_parser.setDaemon(true);
		
		Thread t_regGetter = new Thread(regGetter);
		t_regGetter.setName("RegionVariantGetterThread");
		t_regGetter.setDaemon(true);
		
		Thread t_varMerger = new Thread(varMerger);
		t_varMerger.setName("VariantMergerThread");
		t_varMerger.setDaemon(true);
		
		Thread t_varAdder = new Thread(varAdder);
		t_varAdder.setName("NewVariantThread");
		t_varAdder.setDaemon(true);
		
		Thread t_varUpdater = new Thread(varUpdater);
		t_varUpdater.setName("UpdateVariantThread");
		t_varUpdater.setDaemon(true);
		
		t_varUpdater.start();
		t_varAdder.start();
		t_varMerger.start();
		t_regGetter.start();
		t_parser.start();
		
	}
	
	public synchronized boolean isDone()
	{
		if(error_kill) return true;
		if(!reader_done || !this.getter_done || !this.merger_done) return false;
		if(parseRunner.anyThreadsAlive()) return false;
		if(regGetter.anyThreadsAlive()) return false;
		if(varMerger.anyThreadsAlive()) return false;
		if(varAdder.anyThreadsAlive()) return false;
		if(varUpdater.anyThreadsAlive()) return false;
		
		return true;
	}
	
	private synchronized void killAll(Exception e)
	{
		System.err.println("An unhandled exception has been detected! Terminating variant addition...");
		System.err.println("Parsed: " + c_parsed);
		System.err.println("Region Searched: " + c_reg);
		System.err.println("Merged: " + c_merged);
		System.err.println("Added: " + c_added);
		System.err.println("Updated: " + c_updated);
		error_kill = true;
		parseRunner.requestTermination();
		regGetter.requestTermination();
		varMerger.requestTermination();
		varAdder.requestTermination();
		varUpdater.requestTermination();
		if(e != null) e.printStackTrace();
	}

}
