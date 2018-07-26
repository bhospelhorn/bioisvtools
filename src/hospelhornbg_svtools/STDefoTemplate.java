package hospelhornbg_svtools;

import java.io.IOException;

import hospelhornbg_bioinformatics.VCF;
import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

class STDefoTemplate extends SVDefoTemplate{
	
	public STDefoTemplate(GenomeBuild g)
	{
		super(g);
	}
	
	protected VCF readVCF(String inPath) throws UnsupportedFileTypeException, IOException
	{
		return new VCF(inPath, false, genome);
	}
	

}
