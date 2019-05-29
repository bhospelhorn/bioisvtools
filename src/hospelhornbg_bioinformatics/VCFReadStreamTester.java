package hospelhornbg_bioinformatics;

import hospelhornbg_genomeBuild.GenomeBuild;

public class VCFReadStreamTester {

	public static void main(String[] args) 
	{
		String vcfpath = "C:\\Users\\Blythe\\Documents\\Bioinformatics\\GIAB\\NA24385\\GNA002PB_SVmerged.vcf";
		String genomepath = "C:\\Users\\Blythe\\Documents\\Bioinformatics\\bioisvtools\\GRCh37.gbdh";
		
		try
		{
			GenomeBuild gb = new GenomeBuild(genomepath);
			
			VCFReadStreamer vcfReader = new VCFReadStreamer(vcfpath, gb);
			vcfReader.open();
			
		}
		catch(Exception e)
		{
			e.printStackTrace();
			System.exit(1);
		}
		
		
	}

}
