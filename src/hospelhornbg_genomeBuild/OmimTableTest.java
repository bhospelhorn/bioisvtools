package hospelhornbg_genomeBuild;

import java.io.IOException;

import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class OmimTableTest {

	public static void main(String[] args) 
	{
		String tablepath = "C:\\Users\\hospelhornbg\\Desktop\\genemap2.txt";
		String gbpath = "Z:\\svref\\bioi\\GRCh37.gbdh";
		String gspath = "Z:\\svref\\bioi\\grch37_refSeq.gbgd";
		
		try 
		{
			GenomeBuild gb = new GenomeBuild(gbpath);
			GeneSet gs = new GeneSet(gspath, gb, true);
			OMIMGeneMapImporter im = new OMIMGeneMapImporter(tablepath);
			im.importTable(gs);
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
		} 
		catch (UnsupportedFileTypeException e) 
		{
			e.printStackTrace();
		}
		
		
		
	}

}
