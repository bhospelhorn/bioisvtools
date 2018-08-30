package hospelhornbg_segregation;

import java.io.IOException;
import java.util.Collections;
import java.util.List;

import hospelhornbg_bioinformatics.BreakendPair;
import hospelhornbg_bioinformatics.Genotype;
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.Translocation;
import hospelhornbg_bioinformatics.VCF;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_genomeBuild.Gene;
import hospelhornbg_genomeBuild.GeneSet;
import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class FamTest {

	public static void main(String[] args) 
	{
		if (args.length < 4)
		{
			System.err.println("ERROR: Insufficient arguments detected...");
			System.exit(1);
		}
		//Parse args
		// 1 - Pedigree file
		// 2 - Family VCF file
		// 3 - Genome build
		
		System.err.println("Reading arguments...");
		
		String pedpath = args[1];
		String vcfpath = args[2];
		String gbuild = args[3];
		
		if (pedpath == null || pedpath.isEmpty())
		{
			System.err.println("ERROR: Invalid PED path!");
			System.exit(1);
		}
		else System.err.println("PED Path: " + pedpath);
		if (vcfpath == null || vcfpath.isEmpty())
		{
			System.err.println("ERROR: Invalid VCF path!");
			System.exit(1);
		}
		else System.err.println("VCF Path: " + vcfpath);
		if (gbuild == null || gbuild.isEmpty())
		{
			System.err.println("ERROR: Invalid Genome Build!");
			System.exit(1);
		}
		else System.err.println("Genome Build: " + gbuild);
		
		//Parse Family
		System.err.println("Reading pedigree...");
		Pedigree fam = null;
				
		try 
		{
			fam = new Pedigree(pedpath);
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: IO Error - Input PED file could not be opened!");
			e.printStackTrace();
			System.exit(1);
		} 
		catch (UnsupportedFileTypeException e) 
		{
			System.err.println("ERROR: Parsing Error - Input PED file could not be read!");
			e.printStackTrace();
			System.exit(1);
		}
		System.err.println("Pedigree read!");
		
		//Output family table
		
		//Family Table
		//(As in info form)
		//Fam table lines are preceded by #
		//Name	relation	affected	sex
		List<Individual> famlist = fam.getAllMembers();
		System.out.println("##Family Name: " + fam.getFamilyName());
		System.out.println("##Name\tRelation\tAffected\tSex");
		//Individual proband = fam.getProband();
		for (Individual i : famlist)
		{
			String iname = i.getName();
			String relation = fam.getRelationshipString_ENG(i);
			String affected = i.getAffectedStatus().toString();
			String isex = i.getSex().toString();
			System.out.println("#" + iname + "\t" + relation + "\t" + affected + "\t" + isex);
		}
		
		//Read in genome build
		System.err.println("Loading genome build " + gbuild + " ...");
		GenomeBuild gb = GenomeBuild.loadStandardBuild(gbuild);
		if (gb == null)
		{
			System.err.println("ERROR: Genome build by name\"" + gbuild + "\" not recognized!");
			System.exit(1);
		}
		System.err.println("Build loaded!");
		
		//Load geneset
		System.err.println("Loading gene information...");
		GeneSet gs = GeneSet.loadRefGene(gb);
		if (gb == null)
		{
			System.err.println("ERROR: Gene set for build by name\"" + gbuild + "\" could not be opened!");
			System.exit(1);
		}
		System.err.println("Gene info loaded!");
		
		//Load VCF
		System.err.println("Loading input VCF...");
		VariantPool pool = null;
		try 
		{
			pool = VCF.readVCF(vcfpath, gb, true);
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: IO Error - Input VCF file could not be opened!");
			e.printStackTrace();
			System.exit(1);
		} 
		catch (UnsupportedFileTypeException e) 
		{
			System.err.println("ERROR: Parsing Error - Input VCF file could not be read!");
			e.printStackTrace();
			System.exit(1);
		}
		System.err.println("Input VCF read!");
		
		//Get candidates (and cry so much...)
		System.err.println("Performing inheritance analysis...");
		List<Candidate> clist = Inheritor.getCandidates(pool, fam, gs);
		Collections.sort(clist);
		System.err.println("Inheritance analysis and annotation complete!");
		
		//Table
		//[Normal]
		//VarID	Chr:Pos-End	Type	Size	PosEff	Gene(s)	Transcript(s)	Allele	[SegPattern	Partners]...	Genotypes...
		//[BND/TRA]
		//VarID	Chr:Pos|Chr2:Pos2	Type	Size	PosEff	Gene(s)	Transcript(s)	Allele	[SegPattern	Partners]...	Genotypes...
		int ncand = clist.size();
		System.out.println("-----------------------------");
		System.out.println("##Candidates found: " + ncand);
		String header = "VARID\tPOS\tTYPE\tSIZE\tPOSEFF\tGENE\tTID\tALLELE";
		List<Individual> affected = fam.getAllAffected();
		for (Individual aff : affected){
			header += "\t" + "SEG_" + aff.getName();
			header += "\t" + "SEGPARTNERS" + aff.getName();
		}
		for (Individual i : famlist) header += "\t" + i.getName();
		System.out.println("##" + header);
		
		for (Candidate c : clist)
		{
			Variant v = c.getVariant();
			String rec = "";
			
			//VARID
			rec += v.getVarID() + "\t";
			
			//CHROM/POS
			//TYPE
			//SIZE
			//This field depends upon what v is...
			if (v instanceof StructuralVariant)
			{
				StructuralVariant sv = (StructuralVariant) v;
				if (sv instanceof BreakendPair)
				{
					rec += sv.getChromosome().getUDPName() + ":";
					rec += sv.getPosition() + "|";
					rec += sv.getEndChromosome().getUDPName() + ":";
					rec += sv.getEndPosition() + "\t";
				}
				else if (sv instanceof Translocation)
				{
					rec += sv.getChromosome().getUDPName() + ":";
					rec += sv.getPosition() + "|";
					rec += sv.getEndChromosome().getUDPName() + ":";
					rec += sv.getEndPosition() + "\t";
				}
				else
				{
					rec += sv.getChromosome().getUDPName() + ":";
					rec += sv.getPosition() + "-" + sv.getEndPosition() + "\t";
				}
				rec += sv.getType().getString() + "\t";
				rec += sv.getAbsoluteSVLength() + "bp\t";
			}
			else
			{
				rec += v.getChromosome().getUDPName() + ":" + v.getPosition() + "\t";
				rec += "[None]\t";
				rec += "1bp\t";
			}
			
			//POSEFF
			rec += v.getGeneFuncINFO() + "\t";
			
			//GENE, TID
			Gene g = c.getGene();
			rec += g.getName() + "\t";
			rec += g.getID() + "\t";
			
			//ALLELE
			rec += c.getAllele() + "\t";
			
			//SEG & SEGPARTNERS
			for (Individual aff : affected)
			{
				rec += c.getInheritancePattern(aff).toString() + "\t";
				
				List<Candidate> plist = c.getAllPartners(aff);
				boolean first = true;
				for (Candidate p : plist)
				{
					if(!first) rec += ";";
					rec += p.getVariant().getVarID() + "," + c.getAllele();
					first = false;
				}
				rec += "\t";
			}
			
			//GENOTYPES
			boolean first = true;
			for (Individual i : famlist)
			{
				Genotype gt = v.getSampleGenotype(i.getName());
				if (!first) rec += "\t";
				if (gt == null)
				{
					rec += "null";
				}
				else
				{
					rec += gt.getField(Genotype.INFODEF_GT.getKey());
				}
				first = false;
			}
			
			//Print
			System.out.println(rec);
		}
	
		//GIIIIIIT
	}

}
