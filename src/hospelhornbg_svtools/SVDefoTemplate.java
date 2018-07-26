package hospelhornbg_svtools;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Random;

import hospelhornbg_bioinformatics.UCSCGVBED;
import hospelhornbg_bioinformatics.VCF;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

class SVDefoTemplate implements TrackTemplate{

	private String trName;
	private String trDesc;
	
	protected GenomeBuild genome;
	
	public SVDefoTemplate(GenomeBuild g)
	{
		Random rand = new Random();
		trName = "Track_svdefo_" + Integer.toHexString(rand.nextInt());
		trDesc = "[No description]";
		genome = g;
	}
	
	public void setTrackName(String trackName) 
	{
		if (trackName == null) return;
		if (trackName.isEmpty()) return;
		trName = trackName;
	}

	public void setTrackDescription(String trackDesc) 
	{
		if (trackDesc == null) return;
		if (trackDesc.isEmpty()) return;
		trDesc = trackDesc;
	}

	protected VCF readVCF(String inPath) throws UnsupportedFileTypeException, IOException
	{
		return new VCF(inPath, true, genome);
	}
	
	protected UCSCGVBED generateContainer(String inPath) throws UnsupportedFileTypeException, IOException
	{
		VCF reader = readVCF(inPath);
		VariantPool var = reader.getVariants();
		
		UCSCGVBED myTrack = new UCSCGVBED(trName, var.getVariants());
		myTrack.setDescription(trDesc);
		myTrack.setSampleList(var.getAllSamples());
		
		return myTrack;
	}
	
	public void generateSingle(String inPath, String outPath, String sampleName) throws UnsupportedFileTypeException, IOException 
	{
		UCSCGVBED myTrack = generateContainer(inPath);
		myTrack.write(outPath, sampleName);
	}

	public void generateGroup(String inPath, String outDir) throws UnsupportedFileTypeException, IOException 
	{
		UCSCGVBED myTrack = generateContainer(inPath);
		Collection<String> sampList = myTrack.getSampleList();

		for (String s : sampList)
		{
			String tname = myTrack.getName();
			tname = tname.replaceAll(" ", "");
			tname = tname.replaceAll("?", "");
			tname = tname.replaceAll("!", "");
			tname = tname.replaceAll("/", "");
			tname = tname.replaceAll("\\", "");
			tname = tname.replaceAll("+", "");
			tname = tname.replaceAll("*", "");
			tname = tname.replaceAll("&", "");
			tname = tname.replaceAll("\"", "");
			tname = tname.replaceAll("|", "");
			tname = tname.replaceAll("<", "");
			tname = tname.replaceAll(">", "");
			tname = tname.replaceAll(";", "");
			tname = tname.replaceAll(":", "");
			String outpath = outDir + File.separator + s + "_" + tname + ".bed";
			myTrack.setName(trName + " [" + s + "]");
			myTrack.write(outpath, s);
		}
		
	}

	
	
}
