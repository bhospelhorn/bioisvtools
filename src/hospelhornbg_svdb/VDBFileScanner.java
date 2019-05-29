package hospelhornbg_svdb;

import java.io.IOException;
import java.util.Iterator;

import hospelhornbg_genomeBuild.GeneSet;
import hospelhornbg_genomeBuild.GenomeBuild;
import hospelhornbg_svdb.DBVariant.ParsedVariant;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;
import waffleoRai_Utils.StreamBuffer;

public class VDBFileScanner implements Iterable<DBVariant>{
	
	public static final String VARTBL_EXT = "vdb";
	public static final int VARTBL_VERSION = 1;
	public static final String VARTBL_MAGIC = "VDBt";
	
	private String path;
	private long stOff;
	private long edOff;
	
	private FileBuffer file; 
	private long cpos;
	private long fsz;
	
	private GenomeBuild genome;
	private GeneSet genes;
	
	public VDBFileScanner(GenomeBuild gb, GeneSet gs, String filepath)
	{
		path = filepath;
		stOff = 8; //Skips the header
		edOff = FileBuffer.fileSize(filepath);
		fsz = edOff - 8;
		cpos = -1;
		genome = gb;
		genes = gs;
	}
	
	public VDBFileScanner(GenomeBuild gb, GeneSet gs, String filepath, long stpos, long edpos)
	{
		path = filepath;
		stOff = stpos;
		if(stOff < 8) stOff = 8; //Skips the header
		edOff = edpos;
		cpos = -1;
		fsz = edpos - stpos;
		genome = gb;
		genes = gs;
	}
	
	public void open() throws IOException, UnsupportedFileTypeException
	{
		file = new StreamBuffer(path, stOff, edOff);
		//Check header
		cpos = file.findString(0, 0x10, VARTBL_MAGIC);
		if (cpos != 0) throw new FileBuffer.UnsupportedFileTypeException("VDBFileScanner - File is not valid VDB file! (Magic number not found)");
		//Check version
		int version = file.intFromFile(4L);
		if(version > VARTBL_VERSION) throw new FileBuffer.UnsupportedFileTypeException("VDBFileScanner - Version of VDB file is not recognized!");
		cpos = 8;
	}
	
	public void close()
	{
		file = null;
		cpos = -1;
	}
	
	public DBVariant parseNextVariant()
	{
		if(file == null) return null;
		if(cpos >= fsz) return null;
		
		System.err.println("VDBFileScanner.parseNextVariant | -DEBUG- cpos (before): 0x" + Long.toHexString(cpos));
		ParsedVariant pv = DBVariant.getFromVDBRecord(file, genome, genes, cpos);
		if(pv == null) return null;
		cpos += pv.getSize();
		System.err.println("VDBFileScanner.parseNextVariant | -DEBUG- File Path: " + path);
		System.err.println("VDBFileScanner.parseNextVariant | -DEBUG- Variant Read: " + pv.getVariant().getName());
		System.err.println("VDBFileScanner.parseNextVariant | -DEBUG- Variant Size: 0x" + Long.toHexString(pv.getSize()));
		System.err.println("VDBFileScanner.parseNextVariant | -DEBUG- cpos: 0x" + Long.toHexString(cpos));
		
		return pv.getVariant();
	}
	
	public class VDBIterator implements Iterator<DBVariant>
	{
		
		@Override
		public boolean hasNext() 
		{
			if(cpos < 0) return false;
			return (cpos < fsz);
		}

		@Override
		public DBVariant next() 
		{
			return parseNextVariant();
		}
		
	}
	
	@Override
	public Iterator<DBVariant> iterator() 
	{
		return new VDBIterator();
	}
	
}
