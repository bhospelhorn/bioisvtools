package hospelhornbg_svproject;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import hospelhornbg_bioinformatics.LiteSV;
import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.VCF;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.GeneSet;
import hospelhornbg_genomeBuild.GenomeBuild;
import hospelhornbg_genomeBuild.TwoSexChromSegModel;
import hospelhornbg_segregation.Candidate;
import hospelhornbg_segregation.Family;
import hospelhornbg_segregation.Individual;
import hospelhornbg_segregation.Inheritor;
import hospelhornbg_svproject.VartableFile.GenomeBuildMismatchException;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class SVProject {
	
	/* --- Constants --- */
	
	public static final String MAGIC = "SV_PROJ_";
	public static final int CURRENT_VERSION = 1;
	public static final String EXT = "svproj";
	
	/* --- Static Variables --- */
	
	private static GenomeBuild genomeBuild;
	private static GeneSet geneSet;
	
	/* --- Instance Variables --- */
	
	private String sProjectName;
	
	private Family iFamily;
	private List<Integer> lCohorts;
	
	private List<Candidate> lCandidates;
	
	/* --- Construction --- */
	
	public SVProject(String name)
	{
		sProjectName = name;
		iFamily = null;
		lCohorts = new LinkedList<Integer>();
		lCandidates = null;
	}
	
	public static SVProject loadProject(String projectDirectory, boolean verbose) throws UnsupportedFileTypeException, IOException, GenomeBuildMismatchException
	{
		if (projectDirectory == null || projectDirectory.isEmpty()) return null;
		int lastslash = projectDirectory.lastIndexOf(File.separator);
		String pname = projectDirectory;
		if (lastslash >= 0) pname = projectDirectory.substring(lastslash);
		
		SVProject proj = new SVProject(pname);
		proj.parseSVProject(projectDirectory, verbose);
		
		return proj;
	}
	
	/* --- Parsing --- */
	
	private void parseSVProject(String pDir, boolean verbose) throws UnsupportedFileTypeException, IOException, GenomeBuildMismatchException
	{
		//Magic [8]
		//Version [4]
		//Genome Build UID [4]
		//#Cohort UIDs [4]
			// Cohort UID [4]
		
		//Family & VarTable files should be in the same dir
		
		//First, open the project info file...
		String projectFile = pDir + File.separator + sProjectName + "." + EXT;
		if (!FileBuffer.fileExists(projectFile)) throw new FileBuffer.UnsupportedFileTypeException();
		FileBuffer pfile = FileBuffer.createBuffer(projectFile, true);
		
		long cpos = 0;
		cpos = pfile.findString(0, 0x10, MAGIC);
		if (cpos != 0) throw new FileBuffer.UnsupportedFileTypeException();
		
		//Skip magic and version
		cpos += 12;
		
		//Get Genome UID
		int gbuid = pfile.intFromFile(cpos); cpos += 4;
		genomeBuild = CommonLoader.loadGenomeBuild(gbuid);
		if (genomeBuild == null)
		{
			if (verbose) System.err.println("SVProject.parseSVProject || WARNING: Genome with UID 0x" + Integer.toHexString(gbuid) + " could not be found!");
		}
		else geneSet = CommonLoader.loadGeneSet(genomeBuild, verbose);
		
		//Get cohort UIDs...
		int cohortCount = pfile.intFromFile(cpos); cpos += 4;
		for (int i = 0; i < cohortCount; i++)
		{
			int cid = pfile.intFromFile(cpos); cpos += 4;
			lCohorts.add(cid);
		}
		
		//Find and load family file
		String fampath = pDir + File.separator + sProjectName + "." + Family.FAM_EXT;
		if (!FileBuffer.fileExists(fampath))
		{
			if (verbose) System.err.println("SVProject.parseSVProject || WARNING: Family file could not be found!");
		}
		else
		{
			iFamily = Family.readFromFAMI(fampath);
		}
		
		//If genome build and gene set are non-null, find and load vartable
		if (genomeBuild != null && geneSet != null && iFamily != null)
		{
			String varpath = pDir + File.separator + sProjectName + "." + VartableFile.EXT;
			List<String> suppStrings = CommonLoader.getSupportStrings_ordered();
			VariantPool pool = VartableFile.readVartable(genomeBuild, iFamily, suppStrings, varpath);
			lCandidates = VartableFile.processCandidates(pool, iFamily, geneSet);
		}
		else
		{
			if (verbose) System.err.println("SVProject.parseSVProject || WARNING: VarTable could not be read because genome build, gene set, or family is null!");
		}
	}
	
	/* --- Serialization --- */
	
	public void saveSVProject(String pDir, boolean verbose) throws IOException
	{
		//Write project file
		String projectFile = pDir + File.separator + sProjectName + "." + EXT;
		int ccount = lCohorts.size();
		int pfsize = 8 + 4 + 4 + 4 + (ccount * 4);
		FileBuffer pfile = new FileBuffer(pfsize, true);
		pfile.printASCIIToFile(MAGIC);
		pfile.addToFile(CURRENT_VERSION);
		int gbuid = -1;
		if (genomeBuild != null)
		{
			if (genomeBuild.getUIDEnum() != null) gbuid = genomeBuild.getUIDEnum().getUID();
			else gbuid = genomeBuild.getBuildName().hashCode();
		}
		pfile.addToFile(gbuid);
		pfile.addToFile(ccount);
		for(Integer i : lCohorts) pfile.addToFile(i);
	
		pfile.writeFile(projectFile);
		
		//Write family file
		if(iFamily != null)
		{
			String famFile = pDir + File.separator + sProjectName + "." + Family.FAM_EXT;
			Family.writeToFAMI(iFamily, famFile, true);
		}
		
		//Write vartable file
		if (genomeBuild != null && lCandidates != null && !lCandidates.isEmpty())
		{
			String varpath = pDir + File.separator + sProjectName + "." + VartableFile.EXT;
			List<String> supportStrings = CommonLoader.getSupportStrings_ordered();
			VartableFile.writeVartable(genomeBuild, lCandidates, iFamily, supportStrings, varpath);
			VartableFile.indexVarTable(varpath);
		}
	}
	
	/* --- Getters --- */
	
	public List<StructuralVariant> getSVs()
	{
		LinkedList<StructuralVariant> svlist = new LinkedList<StructuralVariant>();
		for (Candidate c : lCandidates)
		{
			Variant v = c.getVariant();
			//Cast as structural variant
			if (v instanceof StructuralVariant)
			{
				StructuralVariant sv = (StructuralVariant)v;
				//Check to see if already in list
				boolean has = false;
				for (StructuralVariant i : svlist)
				{
					if (i == sv)
					{
						has = true;
						break;
					}
				}
				if (!has)
				{
					//Push.
					//Candidate list is by default sorted by variant
					//So, candidates with the same variant are more likely to be
					//	next to each other.
					//Iterating through the list normally starts at the head.
					//That way, if we check the head first, we check what was
					//	added most recently.
					svlist.push(sv);
				}
			}
		}
		Collections.sort(svlist);
		return svlist;
	}
	
	/* --- Setters --- */
	
	/* --- Import --- */
	
	public void importFamilyFromPED(String pedPath, String fname) throws IOException
	{
		iFamily = Family.readFromPED(pedPath).get(fname);
	}
	
	public void importVariantsFromVCF(String vcfPath) throws IOException, UnsupportedFileTypeException
	{
		VariantPool pool = VCF.readVCF(vcfPath, genomeBuild, true);
		iFamily.adjustSexChromGenotypes(pool.getVariants(), new TwoSexChromSegModel(genomeBuild.getContig("X"), genomeBuild.getContig("Y"), genomeBuild));
		lCandidates = Inheritor.getCandidates(pool, iFamily, geneSet);
		//Dump all allele 0 candidates...
		List<Candidate> nlist = new LinkedList<Candidate>();
		for (Candidate c : lCandidates)
		{
			if (c.getAllele() != 0) nlist.add(c);
		}
		lCandidates = nlist;
	}
	
	public void readEvidenceVCF(String vcfPath, String evidenceKey, String indivSampleName, int leeway, boolean typeMatch, String reportPath) throws IOException
	{
		//Reads in a call set for an individual from a caller like lumpy or manta
		//Scans through project variants to see if there is a match in input callset.
		//If so, flags that var as having evidence from that set
		//The leeway arg is the number of bases beyond the CI90 that should be considered for merging
		
		List<LiteSV> evlist = readEvidenceVCF(vcfPath);
		List<StructuralVariant> svlist = getSVs();
		
		//Index evlist by chrom1 (evlist is an ArrayList)
		int evcount = evlist.size();
		Map<Contig, Integer> evIndex = new HashMap<Contig, Integer>();
		List<Contig> clist = genomeBuild.getChromosomes();
		for (Contig c : clist)
		{
			for (int i = 0; i < evcount; i++)
			{
				LiteSV esv = evlist.get(i);
				if (esv.getChrom1().equals(c))
				{
					evIndex.put(c, i);
					break;
				}
			}
		}
		
		//Open report stream if path is provided
		BufferedWriter bw = null;
		if (reportPath != null)
		{
			bw = new BufferedWriter(new FileWriter(reportPath));
			bw.write("#Leeway=" + leeway + "bp\n");
			bw.write("#TypeMatch=" + typeMatch + "\n");
			bw.write("#ProjectSV\tEvidenceSV\n");
		}
		
		//Scan project svs
		for(StructuralVariant sv : svlist)
		{
			Contig c1 = sv.getChromosome();
			int stIndex = evIndex.get(c1);
			if (stIndex < 0) stIndex = 0;
			for(int i = stIndex; i < evcount; i++)
			{
				//Shouldn't hit list end; breaks when it finds a mismatch...
				LiteSV comp = evlist.get(i);
				if (comp.getChrom1() != c1) break; //Change to .equals if causing issues
				boolean match = comp.svIsEquivalent(sv, leeway, typeMatch);
				if(match)
				{
					//Note in variant
					sv.addSupportMark(indivSampleName, evidenceKey);
					//Note in report (if stream open)
					if (bw != null)
					{
						bw.write(c1.getUDPName() + ":");
						bw.write(sv.getPosition() + "-");
						Contig c2 = sv.getEndChromosome();
						if (c2 == null || c2 == c1)
						{
							bw.write(sv.getEndPosition() + ";");
						}
						else
						{
							bw.write(c2.getUDPName() + ":");
							bw.write(sv.getEndPosition() + ";");
						}
						bw.write(sv.getVarID() + "\t");
						
						Contig e1 = comp.getChrom1();
						bw.write(e1.getUDPName() + ":");
						bw.write(comp.getPosition() + "-");
						Contig e2 = comp.getChrom2();
						if (e2 == null || e2 == c1)
						{
							bw.write(comp.getEndPosition() + ";");
						}
						else
						{
							bw.write(e2.getUDPName() + ":");
							bw.write(comp.getEndPosition() + ";");
						}
						bw.write(comp.getVariantID() + "\t");
						
					}
					break; //Only record first match
				}
			}
		}
		
		if(bw != null) bw.close();
		
	}
	
	private static List<LiteSV> readEvidenceVCF(String path) throws IOException
	{
		Map<String, LiteSV> svNameMap = new HashMap<String, LiteSV>(); //For quick BND pairing
		List<LiteSV> bnds = new LinkedList<LiteSV>();
		
		BufferedReader br = new BufferedReader(new FileReader(path));
		String line = null;
		while ((line = br.readLine()) != null)
		{
			if (line.isEmpty() || line.charAt(0) == '#') continue;
			LiteSV lsv = LiteSV.readFromVCFLine(line, genomeBuild);
			if (lsv != null)
			{
				svNameMap.put(lsv.getVariantID(), lsv);
				if (lsv.getSVType() == SVType.BND) bnds.add(lsv);
			}
		}
		br.close();
		
		//Pair BNDs
		for(LiteSV bnd : bnds)
		{
			//Check if still in map...
			if (svNameMap.containsKey(bnd.getVariantID()))
			{
				//Remove mate from map and dump in this one
				LiteSV mate = svNameMap.remove(bnd.getMateID());
				if (mate != null) bnd.setAsEnd(mate);
			}
		}
		
		//Dump in list and sort
		Collection<LiteSV> all = svNameMap.values();
		int count = all.size() + 1;
		List<LiteSV> svlist = new ArrayList<LiteSV>(count);
		svlist.addAll(all);
		Collections.sort(svlist);
		
		return svlist;
	}
	
	/* --- Export --- */
	
	public VariantPool getAsVariantPool()
	{
		List<StructuralVariant> svlist = getSVs();
		VariantPool pool = new VariantPool(iFamily.countMembers());
		
		//Add infodefs...
		StructuralVariant.addStandardDefs(pool, true);
		//Add samples
		List<Individual> ilist = iFamily.getAllMembers();
		for (Individual i : ilist) pool.addSample(i.getName());
		
		//Add variants
		//(Have to loop manually - won't take list of SVs, picky bastard
		for (StructuralVariant sv : svlist) pool.addVariant(sv);
		
		return pool;
	}
	
	public void exportFamilyToPED(String pedPath) throws IOException
	{
		Family.writeToPED(iFamily, pedPath);
	}
	
	public void exportVariantsToVCF(String vcfPath) throws IOException
	{
		VariantPool pool = getAsVariantPool();
		VCF.writeVCF(pool, "bioi_svAnalyze", vcfPath);
	}

}
