package hospelhornbg_genomeBuild;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import hospelhornbg_genomeBuild.Gene.Exon;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class GBGDMaker {

	public static void main(String[] args) 
	{
		/*
		String inpath = "C:\\Users\\Blythe\\Documents\\Work\\Bioinformatics\\References\\refGene38.txt"; //Hardcode - I'm lazy
		String outpath = "C:\\Users\\Blythe\\Documents\\Work\\Bioinformatics\\References\\grch38_refSeq.gbgd";
		String tblpath = "C:\\Users\\Blythe\\Documents\\Work\\Bioinformatics\\References\\grch38_refSeq.txt";
		String checkpath = "C:\\Users\\Blythe\\Documents\\Work\\Bioinformatics\\References\\grch38_refSeq_check.txt";*/
		
		//Update Sep 2018
		String inpath = "X:\\usr\\hospelhornbg\\Java\\db\\refGene_hg38_Aug2018.txt";
		String outpath = "X:\\usr\\hospelhornbg\\Java\\db\\grch38_refSeq.gbgd";
		String tblpath = "X:\\usr\\hospelhornbg\\Java\\db\\grch38_refSeq.txt";
		String checkpath = "X:\\usr\\hospelhornbg\\Java\\db\\grch38_refSeq_check.txt";
		
		//GenomeBuild ncbi36 = GenomeBuild.loadStandardBuild("hg18");
		//GenomeBuild grch37 = GenomeBuild.loadStandardBuild("grch37");
		GenomeBuild grch38 = GenomeBuild.loadStandardBuild("grch38");

		try 
		{
			FileReader fr = new FileReader(inpath);
			BufferedReader br = new BufferedReader(fr);
			
			GeneSet set38 = new GeneSet("refGene", grch38);
			List<Gene> glist = new LinkedList<Gene>();
			
			String line = null;
			while ((line = br.readLine()) != null)
			{
				/*
				 * 0	bin
				 * 1	name
				 * 2	chrom
				 * 3	strand
				 * 4	txStart
				 * 5	txEnd
				 * 6	cdsStart
				 * 7	cdsEnd
				 * 8	exonCount
				 * 9	exonStarts
				 * 10	exonEnds
				 * 11	score
				 * 12	name2
				 * 13	cdsStartStat
				 * 14	cdsEndStat
				 * 15	exonFrames
				 */
				
				String[] fields = line.split("\t");
				
				String name = fields[12];
				Contig chrom = grch38.getContig(fields[2]);
				boolean strand = true;
				if (fields[3].equals("-")) strand = false;
				int txStart = Integer.parseInt(fields[4]);
				int txEnd = Integer.parseInt(fields[5]);
				int cdsStart = Integer.parseInt(fields[6]);
				int cdsEnd = Integer.parseInt(fields[7]);
				int exonCount = Integer.parseInt(fields[8]);
				String[] exonStarts = fields[9].split(",");
				String[] exonEnds = fields[10].split(",");
				String ID = fields[1];
				boolean ncrna = false;
				if (ID.startsWith("NR")) ncrna = true;
				
				Gene gene = new Gene(exonCount);
				gene.setChromosome(chrom);
				gene.setID(ID);
				gene.setName(name);
				gene.setStrand(strand);
				gene.setStart(txStart);
				gene.setEnd(txEnd);
				gene.setTranslationStart(cdsStart);
				gene.setTranslationEnd(cdsEnd);
				gene.setNCRNA(ncrna);
				
				//int tSt = 0;
				//int tEd = 0;
				
				for (int i = 0; i < exonCount; i++)
				{
					int eSt = Integer.parseInt(exonStarts[i]);
					int eEd = Integer.parseInt(exonEnds[i]);
					//if (i == 0) eSt = cdsStart;
					//if (i == exonCount - 1) eEd = cdsEnd;
					Exon e = new Exon(eSt, eEd);
					gene.addExon(e);
				}
				
				//gene.setTranslationStart(tSt);
				//gene.setTranslationEnd(tEd);
				
				glist.add(gene);
				
			}
			
			br.close();
			fr.close();
			System.out.println("Parse completed.");
			
			set38.addGenes(glist);
			
			//Now, let's output a table.
			set38.outputTable(tblpath);
			System.out.println("Table output completed.");
			//System.exit(0);
			
			set38.serializeGBGD(outpath);
			System.out.println("Serialization completed.");
			
			//Check
			GeneSet check36 = new GeneSet(outpath, grch38, true);
			System.out.println("Reparse completed.");
			check36.outputTable(checkpath);
			System.out.println("Reparse table output completed.");
			
		} 
		catch (IOException e) 
		{
			System.out.println("GBGDMaker.main || IOError with hg38...");
			e.printStackTrace();
		}
		catch (NumberFormatException e) 
		{
			System.out.println("GBGDMaker.main || Integer parsing error with hg38...");
			e.printStackTrace();
		} 
		catch (UnsupportedFileTypeException e) 
		{
			System.out.println("GBGDMaker.main || Serialization error with hg38...");
			e.printStackTrace();
		}
		
		
		
		
	}

}
