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
		String ingb = "C:\\Users\\Blythe\\eclipse-workspace\\bioisvtools\\src\\hospelhornbg_genomeBuild\\resources\\NCBI36.gbdh";
		String outgb = "C:\\Users\\Blythe\\Desktop\\NCBI36.gbdh";
		
		String ings = "C:\\Users\\Blythe\\eclipse-workspace\\bioisvtools\\src\\hospelhornbg_genomeBuild\\resources\\ncbi36_refSeq.gbgd";
		String outgs = "C:\\Users\\Blythe\\Desktop\\ncbi36_refSeq.gbgd";
		//String outtbli = "C:\\Users\\Blythe\\Desktop\\grch37_refSeq_in.out";
		//String outtblo = "C:\\Users\\Blythe\\Desktop\\grch37_refSeq_out.out";
		
		try 
		{
			GenomeBuild gb = new GenomeBuild(ingb);
			gb.printMe();
			gb.setUID(GenomeBuildUID.NCBI36);
			gb.saveGLBD(outgb, true);
			
			GeneSet gs = new GeneSet(ings, gb, true);
			//gs.outputTable(outtbli);
			gs.serializeGBGD(outgs);
			
			//Read back in
			gb = new GenomeBuild(outgb);
			System.err.println("\nReading back in...");
			gb.printMe();
			
			gs = new GeneSet(ings, gb, true);
			//gs.outputTable(outtblo);
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
