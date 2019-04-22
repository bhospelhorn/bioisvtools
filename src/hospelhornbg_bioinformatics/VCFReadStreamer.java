package hospelhornbg_bioinformatics;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ConcurrentMap;

import hospelhornbg_bioinformatics.VariantPool.InfoDefinition;
import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class VCFReadStreamer {
	
private String vcf_path;
	
	private GenomeBuild genome;
	
	private List<String> genoSampleList;
	
	private ConcurrentMap<String, String> shortHeaderFields;
	private ConcurrentMap<String, InfoDefinition> infoFields;
	private ConcurrentMap<String, String> customAlts;
	private ConcurrentMap<String, String> filters;
	private ConcurrentMap<String, InfoDefinition> genotypeFields;
	
	private ConcurrentLinkedQueue<String> lineBuffer;
	private BufferedReader readBuffer;
	
	public class VIterator implements Iterator<Variant>
	{

		@Override
		public boolean hasNext() 
		{
			return bufferHasNext();
		}

		@Override
		public Variant next() 
		{
			String record = nextLine();
			if(record == null) return null;
			try 
			{
				return VCF.parseVCFLine(record, genoSampleList, genome);
			} 
			catch (UnsupportedFileTypeException e) 
			{
				e.printStackTrace();
				return null;
			}
		}
		
	}
	
	public class SVIterator implements Iterator<StructuralVariant>
	{

		@Override
		public boolean hasNext() 
		{
			return bufferHasNext();
		}

		@Override
		public StructuralVariant next() 
		{
			String record = nextLine();
			if(record == null) return null;
			try 
			{
				return VCF.parseVCFLineAsSV(record, genoSampleList, genome);
			} 
			catch (UnsupportedFileTypeException e) 
			{
				e.printStackTrace();
				return null;
			}
		}
		
	}
	
	public VCFReadStreamer(String vcfPath, GenomeBuild gb)
	{
		vcf_path = vcfPath;
		genome = gb;
		
		shortHeaderFields = new ConcurrentHashMap<String, String>();
		infoFields = new ConcurrentHashMap<String, InfoDefinition>();
		customAlts = new ConcurrentHashMap<String, String>();
		filters = new ConcurrentHashMap<String, String>();
		genotypeFields = new ConcurrentHashMap<String, InfoDefinition>();
		
		genoSampleList = new LinkedList<String>();
		lineBuffer = new ConcurrentLinkedQueue<String>();
	}
	
	private boolean bufferHasNext()
	{
		if(lineBuffer.isEmpty())
		{
			String line = null;
			try 
			{
				line = readBuffer.readLine();
			} 
			catch (IOException e) 
			{
				return false;
			}
			if (line == null) return false;
			lineBuffer.add(line);
			return true;
		}
		return true;
	}
	
	private String nextLine()
	{
		String preread = lineBuffer.poll();
		if(preread != null) return preread;
		try {
			return readBuffer.readLine();
		} catch (IOException e) {
			e.printStackTrace();
			return null;
		}
	}
	
	public Iterator<Variant> getIterator()
	{
		return new VIterator();
	}
	
	public Iterator<StructuralVariant> getSVIterator()
	{
		return new SVIterator();
	}
	
	private void readHeader() throws IOException
	{
		String line = null;
		while((line = readBuffer.readLine()) != null)
		{
			if(line.isEmpty()) continue;
			if(!line.startsWith("##")) {
				String[] fields = line.split("\t");
				if(fields.length > 9)
				{
					for(int i = 9; i < fields.length; i++)
					{
						genoSampleList.add(fields[i]);
					}
				}
				return; //We hit the end of the header
			}
			
			//Split into key and value
			String[] kv = line.split("=");
			if(kv.length < 1) continue;
			String key = kv[0].substring(2); //Chop off the ##
			String value = null;
			if(kv.length >= 2) value = kv[1];
			
			//Figure out what type of header line it is
			if(key.equals("INFO"))
			{
				InfoDefinition def = parseInfoDef(value);
				if(def != null) infoFields.put(def.getKey(), def);
			}
			else if (key.equals("ALT"))
			{
				InfoDefinition def = parseInfoDef(value);
				if(def != null) customAlts.put(def.getKey(), def.getDescription());
			}
			else if(key.equals("FORMAT"))
			{
				InfoDefinition def = parseInfoDef(value);
				if(def != null) genotypeFields.put(def.getKey(), def);
			}
			else if(key.equals("FILTER"))
			{
				InfoDefinition def = parseInfoDef(value);
				if(def != null) filters.put(def.getKey(), def.getDescription());
			}
			else shortHeaderFields.put(key, value);
		}
	}
	
	private static InfoDefinition parseInfoDef(String hlineValue)
	{
		//copypasted from VCF class because I am lazy
		hlineValue = hlineValue.substring(hlineValue.indexOf("<") + 1, hlineValue.lastIndexOf(">"));
		String[] fields = hlineValue.split(",");
		String key = "";
		String type = "";
		String desc = "";
		String num = "";
		int nArgs = 0;
		for (String s : fields)
		{
			String[] fsplit = s.split("=");
			if (fsplit[0].equals("ID")) key = fsplit[1];
			if (fsplit[0].equals("Number")) num = fsplit[1];
			if (fsplit[0].equals("Type")) type = fsplit[1];
			if (fsplit[0].equals("Description")) desc = fsplit[1];
			//System.err.println("VCF.parseHeaderLine || Description found: " + desc);
			//String quote = "\"";
			desc = desc.replaceAll("\"", "");
			//System.err.println("VCF.parseHeaderLine || Description after quote strip: " + desc);
		}
		if (!num.isEmpty())
		{
			nArgs = VariantPool.InfoDefinition.parseNumberString(num);
			if (nArgs < -4) return null;
			if (!key.isEmpty() && !type.isEmpty() && !desc.isEmpty())
			{
				int ifieldtype = VariantPool.getInfoDefType(type);
				InfoDefinition def = new InfoDefinition(key, ifieldtype, desc, nArgs);
				return def;
			}
		}
		return null;
	}
	
	public void open() throws IOException
	{
		readBuffer = new BufferedReader(new FileReader(vcf_path));
		readHeader();
	}
	
	public void close() throws IOException
	{
		readBuffer.close();
		readBuffer = null;
	}

}
