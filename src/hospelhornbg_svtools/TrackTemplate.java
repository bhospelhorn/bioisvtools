package hospelhornbg_svtools;

import java.io.IOException;

import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public interface TrackTemplate {
	
	public void setTrackName(String trackName);
	public void setTrackDescription(String trackDesc);
	public void generateSingle(String inPath, String outPath, String sampleName) throws UnsupportedFileTypeException, IOException;
	public void generateGroup(String inPath, String outDir) throws UnsupportedFileTypeException, IOException;

}
