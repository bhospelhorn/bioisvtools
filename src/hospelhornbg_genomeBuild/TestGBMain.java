package hospelhornbg_genomeBuild;

import java.io.IOException;

import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class TestGBMain {

	public static void main(String[] args) 
	{
		
		//String tbl37 = "C:\\Users\\Blythe\\GRCh37.txt";
		//String tbl38 = "C:\\Users\\Blythe\\GRCh38.txt";
		
		//GenomeBuild gb37 = GenomeBuild.loadStandardBuild("grch37");
		//GenomeBuild gb38 = GenomeBuild.loadStandardBuild("grch38");
		
		//gb37.printMe();
		
		//For converting files
		//String ingb = "C:\\Users\\Blythe\\eclipse-workspace\\bioisvtools\\src\\hospelhornbg_genomeBuild\\resources\\GRCh37.gbdh";
		String ingb = "C:\\Users\\Blythe\\eclipse-workspace\\bioisvtools\\src\\hospelhornbg_genomeBuild\\resources\\GRCh37.gbdh";
		String outgb = "C:\\Users\\Blythe\\Desktop\\GRCh37.gbdh";
		
		String ings = "C:\\Users\\Blythe\\eclipse-workspace\\bioisvtools\\src\\hospelhornbg_genomeBuild\\resources\\grch37_refSeq.gbgd";
		String outgs = "C:\\Users\\Blythe\\Desktop\\grch37_refSeq.gbgd";
		String outtbli = "C:\\Users\\Blythe\\Desktop\\grch37_refSeq_in.out";
		String outtblo = "C:\\Users\\Blythe\\Desktop\\grch37_refSeq_out.out";
		
		/*String ingb = "C:\\Users\\Blythe\\eclipse-workspace\\bioisvtools\\src\\hospelhornbg_genomeBuild\\resources\\NCBI36.gbdh";
		String outgb = "C:\\Users\\Blythe\\Desktop\\NCBI36.gbdh";
		
		String ings = "C:\\Users\\Blythe\\eclipse-workspace\\bioisvtools\\src\\hospelhornbg_genomeBuild\\resources\\ncbi36_refSeq.gbgd";
		String outgs = "C:\\Users\\Blythe\\Desktop\\ncbi36_refSeq.gbgd";
		String outtbli = "C:\\Users\\Blythe\\Desktop\\ncbi36_refSeq_in.out";
		String outtblo = "C:\\Users\\Blythe\\Desktop\\ncbi36_refSeq_out.out";*/
		
		/*String ingb = "C:\\Users\\Blythe\\eclipse-workspace\\bioisvtools\\src\\hospelhornbg_genomeBuild\\resources\\GRCh38.gbdh";
		String outgb = "C:\\Users\\Blythe\\Desktop\\GRCh38.gbdh";
		
		String ings = "C:\\Users\\Blythe\\eclipse-workspace\\bioisvtools\\src\\hospelhornbg_genomeBuild\\resources\\grch38_refSeq.gbgd";
		String outgs = "C:\\Users\\Blythe\\Desktop\\grch38_refSeq.gbgd";
		String outtbli = "C:\\Users\\Blythe\\Desktop\\grch38_refSeq_in.out";
		String outtblo = "C:\\Users\\Blythe\\Desktop\\grch38_refSeq_out.out";*/
		
		try 
		{
			System.err.println("Read input genome build...");
			GenomeBuild gb = new GenomeBuild(ingb);
			System.err.println("Writing genome build version 4...");
			gb.saveGLBD(outgb, true);
			
			System.err.println("Read input gene set...");
			GeneSet gs = new GeneSet(ings, gb, true);
			System.err.println("Writing gene set table...");
			gs.outputTable(outtbli);
			System.err.println("Writing gene set version 3...");
			gs.serializeGBGD(outgs);
			
			//Read back in
			System.err.println("Reading back in genome build...");
			gb = new GenomeBuild(outgb);
			System.err.println("Genome build read:");
			gb.printMe();
			System.err.println("Reading back input gene set...");
			gs = new GeneSet(outgs, gb, true);
			System.err.println("Writing gene set table...");
			gs.outputTable(outtblo);
			System.err.println("Done!");
		} 
		catch (IOException e) 
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		catch (UnsupportedFileTypeException e) 
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

}
