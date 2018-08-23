package hospelhornbg_sviewerGUI;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.DefaultListModel;
import javax.swing.ListModel;

import hospelhornbg_bioinformatics.BreakendPair;
import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_bioinformatics.VariantFilter;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_bioinformatics.VariantPool.InfoDefinition;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.Gene;
import hospelhornbg_genomeBuild.GeneFunc;
import hospelhornbg_genomeBuild.GeneSet;
import hospelhornbg_genomeBuild.GenomeBuild;
import hospelhornbg_segregation.Candidate;
import hospelhornbg_segregation.Inheritance;
import hospelhornbg_segregation.Inheritor;
import hospelhornbg_segregation.Pedigree;

public class ViewManager {
	
	/* --- Constants --- */
	
	/* --- Instance Variables --- */
	
	private VariantPool mainpool;
	private VariantPool truthset;
	
	private Pedigree family;
	
	private List<Candidate> candidateList;
	private Map<Integer, Candidate> visibleCandidateLinks;
	
		// Genome Build
	
	private GenomeBuild build;
	private GeneSet genes;
	
		// Basic filters
	
	private boolean confirmedOnly;
	
	private List<Contig> allChromosomes;
	private Set<Contig> includedChromosomes;
	
	private boolean includeSV;
	private boolean includeNonSV;
	private Set<SVType> includedSVTypes;
	
	private int minimumSize;
	private int maximumSize;
	
	private boolean filterQual;
	private double minimumQual;
	private boolean filterNoQual;
	
	private boolean harshBND; // BND chromsome inclusion mode
	
	private Set<VariantFilter> infoFilters;
	
		// Advanced filters
	
	private Set<GeneFunc> passedEffects;
	private Set<Inheritance> passedInheritance;
	private boolean passUnpairedHalfHets;
	
		// Columns
	
	private List<ViewColumn> includedColumns;
	private List<String> infoColumns;
	private List<String> sampleColumns;
	
	/* --- Inner Classes --- */
	
	public static enum ViewColumn
	{
		CHROM("Chromosome"),
		POSITION("Position"),
		SVEND("End Position"),
		REFALLELE("Ref. Allele"),
		ALTALLELES("Alt. Allele(s)"),
		SVTYPE("SV Type"),
		SVLEN("SV Length"),
		CONFIRMED("Confirmed"),
		VARIANTID("Variant ID"),
		FILTERS("Filters"),
		QUALITY("Quality"),
		SVCI("SV Confidence Intervals"),
		SVCI95("SV Confidence Intervals (95%)"),
		OTHERINFO("Info Field [GENERAL]"), //Lumpy stuff will just have to go here.
		SAMPLEGENO("Sample Genotype [GENERAL]"),
		SEGREGATION("Inheritance Pattern"),
		POSITION_EFF("Position Effect"),
		GENE("Gene");
		
		private String label;
		
		private ViewColumn(String enName)
		{
			label = enName;
		}
		
		public static List<ViewColumn> getStandardSet()
		{
			List<ViewColumn> all = new ArrayList<ViewColumn>(13);
			all.add(CHROM);
			all.add(POSITION);
			all.add(SVEND);
			all.add(REFALLELE);
			all.add(ALTALLELES);
			all.add(SVTYPE);
			all.add(SVLEN);
			all.add(CONFIRMED);
			all.add(VARIANTID);
			all.add(FILTERS);
			all.add(QUALITY);
			all.add(SVCI);
			all.add(SVCI95);
			return all;
		}
		
		public String toString()
		{
			return label;
		}
	}
	
	/* --- Construction --- */
	
	public ViewManager(String genomebuild)
	{
		mainpool = null;
		truthset = null;
		confirmedOnly = false;
		allChromosomes = new ArrayList<Contig>(64);
		includedChromosomes = new HashSet<Contig>();
		includeSV = true;
		includeNonSV = true;
		includedSVTypes = new HashSet<SVType>();
		includedSVTypes.addAll(SVType.allTypes());
		minimumSize = 0;
		maximumSize = Integer.MAX_VALUE;
		includedColumns = new LinkedList<ViewColumn>();
		infoColumns = new LinkedList<String>();
		sampleColumns = new LinkedList<String>();
		infoFilters = new HashSet<VariantFilter>();
		filterQual = false;
		filterNoQual = false;
		minimumQual = -2.0;
		harshBND = false;
		family = null;
		candidateList = null;
		visibleCandidateLinks = new HashMap<Integer, Candidate>();
		
		build = GenomeBuild.loadStandardBuild(genomebuild);
		if (build == null)
		{
			throw new IllegalArgumentException();
		}
		genes = GeneSet.loadRefGene(build);
		
		passedEffects = new HashSet<GeneFunc>();
		passedInheritance = new HashSet<Inheritance>();
		passUnpairedHalfHets = true;
	}
	
	/* --- Build --- */
	
	public void setGenome(String buildname)
	{
		try
		{
			GenomeBuild gb = GenomeBuild.loadStandardBuild(buildname);
			if (gb != null)
			{
				build = gb;
				genes = GeneSet.loadRefGene(build);
			}
			else
			{
				System.err.println("ViewManager.setGenome || ERROR: Build " + buildname + " could not be found!");
				throw new IllegalArgumentException();
			}
		}
		catch (Exception e)
		{
			System.err.println("ViewManager.setGenome || ERROR: Build " + buildname + " could not be found!");
			throw new IllegalArgumentException();
		}
		
		//Refresh if needed - unload variant pool if loaded
		mainpool = null;
		truthset = null;
		candidateList = null;
		readAllChromosomes();
	}
	
	/* --- Pools --- */
	
	public void loadVariantPool(VariantPool pool, boolean defoSV)
	{
		if (pool == null) return;
		mainpool = pool;
		truthset = null;
		mainpool.unconfirmAll();
		readAllChromosomes();
		if (defoSV) {
			setColumnsSVDefault();
			minimumSize = 50;
			maximumSize = 5000000;
			excludeSVType(SVType.BND);
		}
		else {
			setColumnsDefault();
			minimumSize = 0;
			maximumSize = Integer.MAX_VALUE;
		}
		
		//See if ped is loaded. If so, convert to candidates...
		if (family != null)
		{
			candidateList = Inheritor.getCandidates(mainpool, family, genes);	
		}
		
	}
	
	public void loadTruthSet(VariantPool truth, boolean stringent_bothEnds, boolean stringent_CI, boolean stringent_overlap)
	{
		mainpool.unconfirmAll();
		if (truth == null) return;
		mainpool.confirm(truth, stringent_bothEnds, stringent_CI, stringent_overlap);
		truthset = truth;
	}
	
	public void reconfirmTruthSet(boolean stringent_bothEnds, boolean stringent_CI, boolean stringent_overlap)
	{
		mainpool.unconfirmAll();
		if (truthset == null) return;
		mainpool.confirm(truthset, stringent_bothEnds, stringent_CI, stringent_overlap);
	}
	
	private void readAllChromosomes()
	{
		if (mainpool == null) return;
		//allChromosomes.addAll(mainpool.getAllChromosomes());
		allChromosomes = build.getChromosomes();
		Collections.sort(allChromosomes);
		includedChromosomes.addAll(allChromosomes);
	}
	
	public boolean poolLoaded()
	{
		return (mainpool != null);
	}
	
	public boolean truthSetLoaded()
	{
		return (truthset != null);
	}
	
	/* --- Pedigree --- */
	
	public boolean pedigreeLoaded()
	{
		return (family != null);
	}

	public void loadPedigree(Pedigree ped)
	{
		family = ped;
		if (family != null)
		{
			candidateList = Inheritor.getCandidates(mainpool, family, genes);	
			if(!includedColumns.contains(ViewColumn.SEGREGATION)) includedColumns.add(ViewColumn.SEGREGATION);
			if(!includedColumns.contains(ViewColumn.GENE))includedColumns.add(ViewColumn.GENE);
			if(!includedColumns.contains(ViewColumn.POSITION_EFF))includedColumns.add(ViewColumn.POSITION_EFF);
		}
	}
	
	public Pedigree getPedigree()
	{
		return family;
	}
	
	/* --- Columns --- */
	
	public Set<GeneFunc> getPassedCodingEffects()
	{
		Set<GeneFunc> copy = new HashSet<GeneFunc>();
		copy.addAll(this.passedEffects);
		return copy;
	}
	
	public void addCodingEffect(GeneFunc eff)
	{
		passedEffects.add(eff);
	}
	
	public void addAllCodingEffects()
	{
		GeneFunc[] alleff = GeneFunc.values();
		for (GeneFunc e: alleff)
		{
			passedEffects.add(e);
		}
	}
	
	public void clearCodingEffects()
	{
		passedEffects.clear();
	}
	
	public Set<Inheritance> getPassedInheritancePatterns()
	{
		Set<Inheritance> copy = new HashSet<Inheritance>();
		copy.addAll(this.passedInheritance);
		return copy;
	}
	
	public void addInheritancePattern(Inheritance ip)
	{
		passedInheritance.add(ip);
	}
	
	public void addAllInheritancePatterns()
	{
		Inheritance[] allip = Inheritance.values();
		for (Inheritance i: allip)
		{
			passedInheritance.add(i);
		}
	}
	
	public boolean includeUnpairedHalfHets()
	{
		return this.passUnpairedHalfHets;
	}
	
	public void setIncludeUnpairedHalfHets(boolean b)
	{
		passUnpairedHalfHets = b;
	}

	public void clearInheritancePatterns()
	{
		passedInheritance.clear();
		passUnpairedHalfHets = false;
	}
	
	public String[] getColumnHeader()
	{
		String[] columns = new String[getNumberColumns()];
		int lSz = includedColumns.size();
		int i = 0;
		for (int j = 0; j < lSz; j++)
		{
			ViewColumn c = includedColumns.get(j);
			if (c == ViewColumn.OTHERINFO)
			{
				//System.out.println("ViewManager.getColumnHeader || Other Info flag found ");
				int infofields = infoColumns.size();
				for (int k = 0; k < infofields; k++)
				{
					columns[i] = infoColumns.get(k);
					i++;
				}
			}
			else if (c == ViewColumn.SAMPLEGENO)
			{
				//System.out.println("ViewManager.getColumnHeader || Sample fields flag found ");
				int genofields = sampleColumns.size();
				for (int k = 0; k < genofields; k++)
				{
					columns[i] = sampleColumns.get(k);
					i++;
				}
			}
			else 
			{
				columns[i] = getStandardColumnHeader(c);
				i++;
			}
		}
		//System.out.println("ViewManager.getColumnHeader || Number Columns: " + columns.length);
		//for (String s : columns) System.out.println("ViewManager.getColumnHeader || Header: " + s);
		//System.out.println();
		return columns;
	}
	
	private String getStandardColumnHeader(ViewColumn c)
	{
		switch (c)
		{
		case ALTALLELES: return "ALT";
		case CHROM: return "CHROM";
		case CONFIRMED: return "CONFIRMED";
		case FILTERS: return "FILTERS";
		case OTHERINFO: return "INFO";
		case POSITION: return "POS";
		case QUALITY: return "QUAL";
		case REFALLELE: return "REF";
		case SAMPLEGENO: return "GENOTYPE";
		case SVCI: return "SVCI";
		case SVCI95: return "SVCI95";
		case SVEND: return "ENDPOS";
		case SVLEN: return "SVLEN";
		case SVTYPE: return "SVTYPE";
		case VARIANTID: return "ID";
		default:
			break;
		}
		return "";
	}
	
	public int getNumberColumns()
	{
		int t = 0;
		for (ViewColumn c : includedColumns)
		{
			if (c == ViewColumn.OTHERINFO)
			{
				t += infoColumns.size();
			}
			else if (c == ViewColumn.SAMPLEGENO)
			{
				t += sampleColumns.size();
			}
			else t++;
		}
		return t;
	}
	
	public void setColumnsDefault()
	{
		includedColumns.clear();
		
		includedColumns.add(ViewColumn.CHROM);
		includedColumns.add(ViewColumn.POSITION);
		includedColumns.add(ViewColumn.REFALLELE);
		includedColumns.add(ViewColumn.ALTALLELES);
		includedColumns.add(ViewColumn.CONFIRMED);
		includedColumns.add(ViewColumn.SAMPLEGENO);
		setAllSamples();
	}
	
	public void setColumnsSVDefault()
	{
		includedColumns.clear();
		
		includedColumns.add(ViewColumn.CHROM);
		includedColumns.add(ViewColumn.POSITION);
		includedColumns.add(ViewColumn.SVEND);
		includedColumns.add(ViewColumn.SVTYPE);
		includedColumns.add(ViewColumn.SVLEN);
		includedColumns.add(ViewColumn.CONFIRMED);
		includedColumns.add(ViewColumn.REFALLELE);
		includedColumns.add(ViewColumn.ALTALLELES);
		includedColumns.add(ViewColumn.SAMPLEGENO);
		setAllSamples();
	}
	
	public void setAllSamples()
	{
		sampleColumns.clear();
		if (mainpool != null)
		{
			List<String> samps = mainpool.getAllSamples();
			sampleColumns.addAll(samps);
		}
	}
	
	public void setAllInfoFields()
	{
		infoColumns.clear();
		if (mainpool != null)
		{
			infoColumns.addAll(mainpool.getOrderedInfoKeys());
		}
	}
	
	public void clearColumns()
	{
		includedColumns.clear();
		infoColumns.clear();
		sampleColumns.clear();
	}
	
	public List<ViewColumn> getAllPossibleStandardFields()
	{
		return ViewColumn.getStandardSet();
	}
	
	public final List<ViewColumn> getIncludedStandardFields()
	{
		List<ViewColumn> copy = new ArrayList<ViewColumn>(includedColumns.size());
		for (ViewColumn c : includedColumns)
		{
			copy.add(c);
		}
		return copy;
	}
	
	public List<String> getAllPossibleInfoFields()
	{
		if (mainpool == null) return null;
		return mainpool.getOrderedInfoKeys();
	}
	
	public final List<String> getIncludedInfoFields()
	{
		return infoColumns;
	}
	
	public List<String> getAllSamples()
	{
		if (mainpool == null) return null;
		return mainpool.getAllSamples();
	}
	
	public final List<String> getIncludedSamples()
	{
		return sampleColumns;
	}
	
	public boolean columnIncluded(ViewColumn c)
	{
		return includedColumns.contains(c);
	}
	
	public boolean infoColumnIncluded(String key)
	{
		return infoColumns.contains(key);
	}
	
	public boolean sampleColumnIncluded(String sample)
	{
		return sampleColumns.contains(sample);
	}
	
	protected void includeInfoColumn(String key)
	{
		if (infoColumns.contains(key)) return;
		//System.out.println("ViewManager.includeInfoColumn || Including " + key);
		infoColumns.add(key);
	}
	
	protected void includeSampleColumn(String sampleName)
	{
		if (sampleColumns.contains(sampleName)) return;
		sampleColumns.add(sampleName);
	}
	
	public void includeStandardColumn(ViewColumn column)
	{
		if (includedColumns.contains(column)) return;
		//System.out.println("ViewManager.includeStandardColumn || Including " + column.toString());
		includedColumns.add(column);
	}
	
	/* --- Chrom Sorting --- */
	
	private static class ChromWrap implements Comparable<ChromWrap>
	{
		private Contig chrom;
		
		public ChromWrap(Contig c)
		{
			chrom = c;
		}

		public int compareTo(ChromWrap o) 
		{
			return chrom.compareTo(o.chrom);
		}
	}
	
	/* --- Filter Management --- */
	
	public boolean confirmedOnly()
	{
		return confirmedOnly;
	}
	
	public void setConfirmedOnly(boolean b)
	{
		confirmedOnly = b;
	}
	
	public List<Contig> getAllChromosomes()
	{
		return allChromosomes;
	}
	
	public List<Contig> getChromosomeList()
	{
		List<ChromWrap> sorter = new ArrayList<ChromWrap>(includedChromosomes.size());
		List<Contig> chroms = new ArrayList<Contig>(sorter.size());
		for (Contig c : includedChromosomes) sorter.add(new ChromWrap(c));
		Collections.sort(sorter);
		for (ChromWrap c : sorter) chroms.add(c.chrom);
		return chroms;
	}
	
	public void includeChromosome(Contig chrom)
	{
		if (allChromosomes.contains(chrom))
		{
			includedChromosomes.add(chrom);
		}
	}
	
	public void excludeChromosome(Contig chrom)
	{
		includedChromosomes.remove(chrom);
	}

	public void clearIncludedChromosomes()
	{
		includedChromosomes.clear();
	}
	
	public boolean chromosomeIncluded(Contig chrom)
	{
		return includedChromosomes.contains(chrom);
	}
	
	public boolean variantChromosomeIncluded(Variant v)
	{
		if (chromosomeIncluded(v.getChromosome())) return true;
		for (Contig c : includedChromosomes)
		{
			if (v.isOnChromosome(c)) return true;
		}
		return false;
	}
	
	public boolean includeSV()
	{
		return includeSV;
	}
	
	public void setIncludeSV(boolean b)
	{
		includeSV = b;
	}
	
	public boolean includeNonSV()
	{
		return includeNonSV;
	}
	
	public void setIncludeNonSV(boolean b)
	{
		includeNonSV = b;
	}
	
	public boolean SVTypeIncluded(SVType t)
	{
		return includedSVTypes.contains(t);
	}
	
	public void includeSVType(SVType t)
	{
		includedSVTypes.add(t);
	}
	
	public void excludeSVType(SVType t)
	{
		includedSVTypes.remove(t);
	}
	
	public void clearIncludedSVTypes()
	{
		includedSVTypes.clear();
	}
	
	public int getMinimumSize()
	{
		return minimumSize;
	}
	
	public void setMinimumSize(int min)
	{
		if (min < 0) return;
		if (min > maximumSize) return;
		minimumSize = min;
	}
	
	public int getMaximumSize()
	{
		return maximumSize;
	}
	
	public void setMaximumSize(int max)
	{
		if (max < 1 || max < minimumSize) return;
		maximumSize = max;
	}
	
	public boolean passesAllCustomFilters(Variant v)
	{
		for (VariantFilter f : infoFilters)
		{
			if (!f.passes(v)) return false;
		}
		return true;
	}
	
	public boolean passesAllFilters(Variant v)
	{
		if (confirmedOnly && !v.isConfirmed()) return false;
		if (filterQual && (v.getQuality() < minimumQual)) {
			if (filterNoQual) return false;
			else if (v.getQuality() >= 0) return false;
		}
		if (!variantChromosomeIncluded(v)) return false;
		if (!includeSV && (v instanceof StructuralVariant)) return false;
		if (!includeNonSV && !(v instanceof StructuralVariant)) return false;
		if ((v instanceof StructuralVariant))
		{
			StructuralVariant sv = (StructuralVariant)v;
			if (!SVTypeIncluded(sv.getType())) return false;
			if (sv.getType() != SVType.BND)
			{
				if (sv.getAbsoluteSVLength() < minimumSize) return false;
				if (sv.getAbsoluteSVLength() > maximumSize) return false;	
			}
			if (harshBND && sv instanceof BreakendPair)
			{
				Collection<Contig> chroms = sv.getAllChromosomes();
				boolean badchrom = false;
				for (Contig c : chroms)
				{
					if (!this.chromosomeIncluded(c)) {
						badchrom = true;
						break;
					}
				}
				if (badchrom) return false;
			}
		}
		else
		{
			if (v.getLargestAbsoluteLength() < minimumSize) return false;
			if (v.getSmallestAbsoluteLength() > maximumSize) return false;
		}
		if (!passesAllCustomFilters(v)) return false;
		
		return true;
	}
	
	public List<Variant> getFilteredSet(Collection<Variant> sourceset)
	{
		List<Variant> filtered = new LinkedList<Variant>();
		
		for (Variant v : sourceset)
		{
			boolean passes = passesAllFilters(v);
			if (passes) filtered.add(v);
		}
		
		return filtered;
	}
	
	public List<Variant> getFilteredSet()
	{
		if (pedigreeLoaded())
		{
			List<Candidate> cset = getFilteredCandidateSet();
			int sz = cset.size();
			List<Variant> vset = new ArrayList<Variant>(sz + 1);
			for (Candidate c : cset) vset.add(c.getVariant());
			return vset;
		}
		List<Variant> sourceset = mainpool.getVariants();
		return getFilteredSet(sourceset);
	}
	
	public boolean candidatePassesFilters(Candidate c)
	{
		if (c == null) return false;
		//See if internal variant passes
		Variant v = c.getVariant();
		if (!passesAllFilters(v)) return false;
		//See if variant effect passes
		boolean epass = false;
		for (GeneFunc eff : passedEffects)
		{
			if (v.getGeneFunction() == eff)
			{
				epass = true;
				break;
			}
		}
		if (!epass) return false;
		//See if segregation passes
		Inheritance ip = c.getInheritancePattern();
		boolean spass = false;
		for (Inheritance i : passedInheritance)
		{
			if (ip == i)
			{
				spass = true;
				break;
			}
		}
		if (!spass) return false;
		//See if unpaired halfhet
		if(!passUnpairedHalfHets)
		{
			if (ip == Inheritance.HALF_HET || ip == Inheritance.HALF_HET_SV)
			{
				List<Candidate> partners = c.getAllPartners();
				if (partners == null) return false;
				if (partners.isEmpty()) return false;
			}
		}
		return true;
	}
	
	public List<Candidate> getFilteredCandidateSet()
	{
		List<Candidate> passed = new LinkedList<Candidate>();
		if (candidateList == null) return passed;
		
		//List<Candidate> rescuecand = new LinkedList<Candidate>();
		
		for (Candidate c : candidateList)
		{
			boolean cpassed = candidatePassesFilters(c);
			if (cpassed) passed.add(c);
			//else rescuecand.add(c);
		}
		
		return passed;
	}

	public boolean filterByQuality()
	{
		return filterQual;
	}
	
	public void setFilterByQuality(boolean b)
	{
		filterQual = b;
	}
	
	public double getMinimumQuality()
	{
		return this.minimumQual;
	}
	
	public void setMinimumQuality(double q)
	{
		minimumQual = q;
	}
	
	public boolean excludeNoQuality()
	{
		return this.filterNoQual;
	}
	
	public void setExcludeNoQuality(boolean b)
	{
		filterNoQual = b;
	}
	
	public boolean compositeChromosomeExclusive()
	{
		return harshBND;
	}
	
	public void setCompositeChrosomsomeExclusionMode(boolean harsh)
	{
		harshBND = harsh;
	}
	
	/* --- Custom Filters --- */
	
	public ListModel<VariantFilter> getCustomFilters_Swing()
	{
		List<VariantFilter> fList = new ArrayList<VariantFilter>(infoFilters.size());
		fList.addAll(infoFilters);
		Collections.sort(fList);
		DefaultListModel<VariantFilter> model = new DefaultListModel<VariantFilter>();
		for (VariantFilter f : fList) model.addElement(f);
		return model;
	}
	
	public void addCustomFilter(VariantFilter filter)
	{
		infoFilters.add(filter);
	}
	
	public void removeCustomFilter(VariantFilter filter)
	{
		infoFilters.remove(filter);
	}
	
	public void clearCustomFilters()
	{
		infoFilters.clear();
	}
	
	public int getInfoFieldType(String fieldKey)
	{
		InfoDefinition id = mainpool.getInfoDef(fieldKey);
		if (id == null) return VariantPool.INFODEF_UNK;
		return id.getType();
	}
	
	/* --- Variant View --- */
	
	private String[][] getTable_var()
	{
		List<Variant> filtered = getFilteredSet();
		int cols = getNumberColumns();
		int sCols = includedColumns.size();
		String[][] data = new String[filtered.size()][cols];
	
		int i = 0;
		for (Variant v : filtered)
		{
			int l = 0;
			for (int j = 0; j < sCols; j++)
			{
				ViewColumn c = includedColumns.get(j);
				if (c == ViewColumn.OTHERINFO)
				{
					for (String k : infoColumns)
					{
						data[i][l] = getInfoString(k, v);
						l++;
					}
				}
				else if (c == ViewColumn.SAMPLEGENO)
				{
					for (String s : sampleColumns)
					{
						data[i][l] = getSampleGenoString(s, v);
						l++;
					}
				}
				else 
				{
					data[i][l] = getPropertyString(c, v);
					l++;
				}
			}
			i++;
		}
		
		return data;
	}
	
	private String[][] getTable_cand()
	{
		if (candidateList == null) return null;
		visibleCandidateLinks.clear();
		List<Candidate> filtered = getFilteredCandidateSet();
		int cols = getNumberColumns();
		int sCols = includedColumns.size();
		String[][] data = new String[filtered.size()][cols];
		
		int i = 0;
		for (Candidate n : filtered)
		{
			Variant v = n.getVariant();
			int l = 0;
			for (int j = 0; j < sCols; j++)
			{
				ViewColumn c = includedColumns.get(j);
				if (c == ViewColumn.OTHERINFO)
				{
					for (String k : infoColumns)
					{
						data[i][l] = getInfoString(k, v);
						l++;
					}
				}
				else if (c == ViewColumn.SAMPLEGENO)
				{
					for (String s : sampleColumns)
					{
						data[i][l] = getSampleGenoString(s, v);
						l++;
					}
				}
				else 
				{
					data[i][l] = getPropertyString(c, n);
					l++;
				}
			}
			i++;
			visibleCandidateLinks.put(i, n);
		}
		
		return data;
	}
	
	public String[][] getTable()
	{
		if (candidateList == null) return getTable_var();
		return getTable_cand();
	}
	
	public Candidate getLinkedCandidate(int tblIndex)
	{
		if (candidateList == null) return null;
		if (tblIndex < 0) return null;
		return this.visibleCandidateLinks.get(tblIndex);
	}
	
	public String getPropertyString(ViewColumn c, Variant v)
	{
		if (v == null) return "";
		switch(c)
		{
		case ALTALLELES:
			int nAlts = v.countAltAlleles();
			if (nAlts < 1) return "";
			if (nAlts == 1)
			{
				return v.getAltAllele(0);
			}
			String[] alts = v.getAllAltAlleles();
			if (alts != null && alts.length > 0)
			{
				String s = "";
				for (String a : alts)
				{
					s += a + "\n";
				}
				s = s.substring(0, s.length() - 1);
				return s;
			}
			else return "";
		case CHROM: return v.getChromosome().getUDPName();
		case CONFIRMED:
			if (v.isConfirmed()) return "#";
			else return "";
		case FILTERS:
			if (v.passedAllFilters()) return "PASS";
			String[] filters = v.getFiltersFailed();
			if (filters != null && filters.length > 0)
			{
				String s = "";
				for (String f : filters)
				{
					s += f + "\n";
				}
				s = s.substring(0, s.length() - 1);
				return s;
			}
			else return "";
		case OTHERINFO:
			return ""; //This should be accessed another way.
		case POSITION:
			return formatInteger(v.getPosition());
		case QUALITY:
			return Double.toString(v.getQuality());
		case REFALLELE:
			return v.getRefAllele();
		case SAMPLEGENO:
			return ""; //This should be accessed another way.
		case SVCI:
			if (v instanceof StructuralVariant)
			{
				StructuralVariant sv = (StructuralVariant)v;
				int sl = sv.getCIDiff(false, false, false);
				int sh = sv.getCIDiff(false, false, true);
				int el = sv.getCIDiff(true, false, false);
				int eh = sv.getCIDiff(true, false, true);
				return sl + " : " + sh + " || " + el + " : " + eh;
			}
			else return "";
		case SVCI95:
			if (v instanceof StructuralVariant)
			{
				StructuralVariant sv = (StructuralVariant)v;
				int sl = sv.getCIDiff(false, true, false);
				int sh = sv.getCIDiff(false, true, true);
				int el = sv.getCIDiff(true, true, false);
				int eh = sv.getCIDiff(true, true, true);
				return sl + " : " + sh + " || " + el + " : " + eh;
			}
			else return "";
		case SVEND:
			if (v instanceof StructuralVariant)
			{
				StructuralVariant sv = (StructuralVariant)v;
				return formatInteger(sv.getEndPosition());
			}
			else return "";
		case SVLEN:
			if (v instanceof StructuralVariant)
			{
				StructuralVariant sv = (StructuralVariant)v;
				return formatInteger(sv.getSVLength());
			}
			else return "";
		case SVTYPE:
			if (v instanceof StructuralVariant)
			{
				StructuralVariant sv = (StructuralVariant)v;
				return sv.getType().toString();
			}
			else return "";
		case VARIANTID:
			return v.getVarID();
		case GENE:
			return ""; //Candidate only
		case POSITION_EFF:
			return ""; //Candidate only
		case SEGREGATION:
			return ""; //Candidate only
		default:
			break;
		}
		return "";
	}

	public String getPropertyString(ViewColumn c, Candidate n)
	{
		if (n == null) return "";
		Variant v = n.getVariant();
		if (v == null) return "";
		switch(c)
		{
		case ALTALLELES:
			int nAlts = v.countAltAlleles();
			if (nAlts < 1) return "";
			if (nAlts == 1)
			{
				return v.getAltAllele(0);
			}
			String[] alts = v.getAllAltAlleles();
			if (alts != null && alts.length > 0)
			{
				String s = "";
				for (String a : alts)
				{
					s += a + "\n";
				}
				s = s.substring(0, s.length() - 1);
				return s;
			}
			else return "";
		case CHROM: return v.getChromosome().getUDPName();
		case CONFIRMED:
			if (v.isConfirmed()) return "#";
			else return "";
		case FILTERS:
			if (v.passedAllFilters()) return "PASS";
			String[] filters = v.getFiltersFailed();
			if (filters != null && filters.length > 0)
			{
				String s = "";
				for (String f : filters)
				{
					s += f + "\n";
				}
				s = s.substring(0, s.length() - 1);
				return s;
			}
			else return "";
		case OTHERINFO:
			return ""; //This should be accessed another way.
		case POSITION:
			return formatInteger(v.getPosition());
		case QUALITY:
			return Double.toString(v.getQuality());
		case REFALLELE:
			return v.getRefAllele();
		case SAMPLEGENO:
			return ""; //This should be accessed another way.
		case SVCI:
			if (v instanceof StructuralVariant)
			{
				StructuralVariant sv = (StructuralVariant)v;
				int sl = sv.getCIDiff(false, false, false);
				int sh = sv.getCIDiff(false, false, true);
				int el = sv.getCIDiff(true, false, false);
				int eh = sv.getCIDiff(true, false, true);
				return sl + " : " + sh + " || " + el + " : " + eh;
			}
			else return "";
		case SVCI95:
			if (v instanceof StructuralVariant)
			{
				StructuralVariant sv = (StructuralVariant)v;
				int sl = sv.getCIDiff(false, true, false);
				int sh = sv.getCIDiff(false, true, true);
				int el = sv.getCIDiff(true, true, false);
				int eh = sv.getCIDiff(true, true, true);
				return sl + " : " + sh + " || " + el + " : " + eh;
			}
			else return "";
		case SVEND:
			if (v instanceof StructuralVariant)
			{
				StructuralVariant sv = (StructuralVariant)v;
				return formatInteger(sv.getEndPosition());
			}
			else return "";
		case SVLEN:
			if (v instanceof StructuralVariant)
			{
				StructuralVariant sv = (StructuralVariant)v;
				return formatInteger(sv.getSVLength());
			}
			else return "";
		case SVTYPE:
			if (v instanceof StructuralVariant)
			{
				StructuralVariant sv = (StructuralVariant)v;
				return sv.getType().toString();
			}
			else return "";
		case VARIANTID:
			return v.getVarID();
		case GENE:
			Gene g = n.getGene();
			if (g == null) return "[None]";
			return g.getName();
		case POSITION_EFF:
			GeneFunc e = v.getGeneFunction();
			if (e == null) return "[Unknown]";
			return e.toString();
		case SEGREGATION:
			Inheritance ip = n.getInheritancePattern();
			if (ip == null) return "[Unknown]";
			return ip.toString();
		default:
			break;
		}
		return "";
	}
	
	public String getInfoString(String fieldKey, Variant v)
	{
		if (v == null) return "";
		if (fieldKey == null) return "";
		if (fieldKey.isEmpty()) return "";
		if (v.getInfoFlag(fieldKey))
		{
			return "#";
		}
		String[] vals = v.getInfoEntry(fieldKey);
		if (vals == null) return "";
		String s = "";
		int nfields = vals.length;
		for (int i = 0; i < nfields; i++)
		{
			String e = vals[i];
			s += e;
			if (i < nfields - 1)
			{
				s += ", ";
			}
		}
		s = s.substring(0, s.length() - 1);
		return s;
	}
	
	public String getSampleGenoString(String sample, Variant v)
	{
		if (v == null) return "";
		if (sample == null) return "";
		if (sample.isEmpty()) return "";
		return v.getSampleGenotypeString(sample);
	}
	
	private String formatInteger(int i)
	{
		String s = String.format("%,d", i);
		return s;
	}
	

}
