package hospelhornbg_genomeBuild;

public class TestGBMain {

	public static void main(String[] args) 
	{
		//String tbl37 = "C:\\Users\\Blythe\\GRCh37.txt";
		//String tbl38 = "C:\\Users\\Blythe\\GRCh38.txt";
		
		GenomeBuild gb37 = GenomeBuild.loadStandardBuild("grch37");
		//GenomeBuild gb38 = GenomeBuild.loadStandardBuild("grch38");
		
		gb37.printMe();
		
	}

}
