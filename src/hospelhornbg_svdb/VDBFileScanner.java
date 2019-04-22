package hospelhornbg_svdb;

import java.io.IOException;
import java.util.Iterator;

import hospelhornbg_genomeBuild.GeneSet;
import hospelhornbg_genomeBuild.GenomeBuild;
import hospelhornbg_svdb.DBVariant.ParsedVariant;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.StreamBuffer;

public class VDBFileScanner implements Iterable<DBVariant>{

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
	}
	
	public VDBFileScanner(GenomeBuild gb, GeneSet gs, String filepath, long stpos, long edpos)
	{
		path = filepath;
		stOff = stpos;
		if(stOff < 8) stOff = 8; //Skips the header
		edOff = edpos;
		cpos = -1;
		fsz = edpos - stpos;
	}
	
	public void open() throws IOException
	{
		cpos = 0;
		file = new StreamBuffer(path, stOff, edOff);
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
		
		ParsedVariant pv = DBVariant.getFromVDBRecord(file, genome, genes, cpos);
		if(pv == null) return null;
		cpos += pv.getSize();
		
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
