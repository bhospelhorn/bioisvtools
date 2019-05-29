package hospelhornbg_bioinformatics;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collection;
import java.util.Deque;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import hospelhornbg_bioinformatics.VariantPool.InfoDefinition;
import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

/*
 * UPDATES
 * 
 * Initial Version: 1.0.0 (December 14, 2017)
 * 
 * 1.0.0 -> 1.1.0 | January 17, 2018
 * 	Javadoc
 * 	Update for 1.1.0 compatibility with the rest of the package
 * 	Added more setters for improved accessibility
 * 
 * 1.1.0 -> 1.2.0 | January 23, 2018
 * 	Multithreaded the parser.
 * 
 * 1.2.0 -> 1.3.0 | January 24, 2018
 * 	Altered the parser to just use the native Java library for line reading instead of
 * 	my underdeveloped LineBuffer class.
 * 
 * 1.3.0 -> 1.4.0 | January 31, 2018
 * 	Overhauled parser multithreading framework to hopefully eliminate memory access conflicts.
 * 
 * 1.4.0 -> 1.4.1 | February 1, 2018
 * 	Commented out simultaneous parsing/reading framework due to undiagnosed threading issues - made
 * 	reading and parsing sequential for the time being.
 *  There don't appear to be any obvious race conditions, so why the simultaneous approach isn't working
 *  has yet to be determined. (All testing performed on multi-core machines)
 * 
 * 1.4.1 -> 1.5.0 | February 20, 2018
 * 	GenomeBuild/ Contig update
 *
 * 1.5.0 -> 1.5.1 | March 1, 2018
 * 	Commented out some of the parser error messages
 * 
 * 1.5.1 -> 1.5.2 | April 18, 2018
 * 	Added readVCF overload that takes genome build
 *
 * 1.5.2 -> 1.6.0 | April 19, 2018
 * 	Rewrote parser. Again.
 * 
 * 1.6.0 -> 1.6.1 | January 17, 2019
 * 	GenomeBuild format version 2.0
 *
 * 1.6.1 -> 1.6.2 | February 28, 2019
 * 	Static method for parsing INFO field into a key/value map
 * 
 * 1.6.2 -> 1.7.0 | April 22, 2019
 * 	Added static single variant parser methods
 */


/**
 * A container for reading and writing a collection of variants and annotation metadata to
 * VCF format.
 * @author Blythe Hospelhorn
 * @version 1.7.0
 * @since April 22, 2019
 *
 */
public class VCF {
	
	/* --- Constants --- */
	
	/* --- Instance Variables --- */
	
	private VariantPool variants;
	
	private String fileformat;
	private String dateString;
	private String source;
	
	private List<String[]> strayHeaderLines; 
	
	private GenomeBuild gBuild;
	
	/* --- Construction/Parsing --- */
	
	/**
	 * Construct a VCF container from an existing VariantPool and note
	 * the name of the application the VCF is originating from.
	 * @param pool An existing collection of variants and VCF style metadata to wrap.
	 * @param sourceApp A string denoting the name of the application creating
	 * this VCF container.
	 */
	public VCF(VariantPool pool, String sourceApp)
	{
		variants = pool;
		source = sourceApp;
		fileformat = "fileformat=VCFv4.2";
		dateString = getHeaderDate();
		strayHeaderLines = new LinkedList<String[]>();
		gBuild = pool.getGenomeBuild();
	}
	
	/**
	 * Construct a VCF container and internal VariantPool by reading in an existing VCF file.
	 * @param filePath Path on local file system of file to read in.
	 * @param readSVs Whether to additionally attempt to parse out structural variants
	 * into StructuralVariant objects, or to just leave them as annotated Variants.
	 * @throws IOException If there is an error accessing the file on disk.
	 * @throws UnsupportedFileTypeException If there is an error parsing the file as a VCF.
	 */
	public VCF(String filePath, boolean readSVs) throws IOException, UnsupportedFileTypeException
	{
		//FileBuffer vcf = FileBuffer.createBuffer(filePath);
		variants = new VariantPool(256);
		strayHeaderLines = new LinkedList<String[]>();
		gBuild = null;
		parseVCF(filePath, readSVs);
	}
	
	/**
	 * Constructor added with version 1.5.0
	 * <br>Default parse constructor will attempt to determine reference build from VCF header,
	 * and create a new build made from what it finds in the VCF if it can't figure it out.
	 * <br>This is an alternative constructor that allows one to explicitly pass a GenomeBuild
	 * to the parser.
	 * @param filePath Path on local file system of file to read in.
	 * @param readSVs Whether to additionally attempt to parse out structural variants
	 * into StructuralVariant objects, or to just leave them as annotated Variants.
	 * @param genome GenomeBuild to reference when parsing variant coordinates.
	 * @throws IOException If there is an error accessing the file on disk.
	 * @throws UnsupportedFileTypeException If there is an error parsing the file as a VCF.
	 */
	public VCF(String filePath, boolean readSVs, GenomeBuild genome) throws UnsupportedFileTypeException, IOException
	{
		variants = new VariantPool(256);
		strayHeaderLines = new LinkedList<String[]>();
		gBuild = genome;
		parseVCF(filePath, readSVs);
	}
	
	private class ParserThread extends Thread
	{
		private static final String killsignal = "#KILL";
		private Deque<String> linequeue;
		private List<Variant> variants;
		private List<String> slist;
		
		private List<String> rejected;
		
		private boolean running;
		
		public ParserThread(int i)
		{
			linequeue = new LinkedList<String>();
			variants = new LinkedList<Variant>();
			running = false;
			slist = new LinkedList<String>();
			rejected = new LinkedList<String>();
			//slist = sampleList;
			super.setName("VCFParserThread_" + i);
			super.setDaemon(true);
		}
		
		public void run()
		{
			running = true;
			while(running)
			{
				while(!queueEmpty())
				{
					String line = popQueue();
					if (line.equals(killsignal)) {
						running = false;
						break;
					}
					try 
					{
						//if(gBuild == null) System.err.println("VCF.ParserThread.run || gBuild is null!!");
						//else System.err.println("VCF.ParserThread.run || gBuild is not null!!");
						Variant v = new Variant(line, slist, gBuild, !presetGenome);
						addParsedVariant(v);
					} 
					catch (UnsupportedFileTypeException e) 
					{
						parsingErr = true;
						running = false;
						System.err.println("VCF.ParserThread.run || Parsing error! Variant Line: ");
						System.err.println(line);
						e.printStackTrace();
						break;
					}
					catch (NullPointerException e)
					{
						//Rejected contig
						addRejectedRecord(line);
					}
				}
				try {
					Thread.sleep(10);
				} catch (InterruptedException e) {
					Thread.interrupted();
					//e.printStackTrace();
				}
			}
		}
		
		private synchronized void addRejectedRecord(String record)
		{
			rejected.add(record);
		}
		
		public synchronized boolean queueEmpty()
		{
			return linequeue.isEmpty();
		}
		
		public synchronized void addLine(String line)
		{
			linequeue.add(line);
		}
		
		private synchronized String popQueue()
		{
			return linequeue.pop();
		}
		
		public synchronized void killWhenDone()
		{
			linequeue.add(killsignal);
		}
		
		public Collection<Variant> getVariants()
		{
			if(running) return null;
			return variants;
		}
		
		public synchronized void setSampleList(Collection<String> sampleList)
		{
			slist.addAll(sampleList);
		}
		
		private synchronized void addParsedVariant(Variant v)
		{
			variants.add(v);
		}
		
		public synchronized int countParsed()
		{
			return variants.size();
		}
		
		public synchronized int countRejected()
		{
			return rejected.size();
		}
		
	}
	
	private boolean parsingErr = false;
	private boolean presetGenome = false;
	
	/**
	 * VCF core parser.
	 * @param vcfpath Path of file to open and parse.
	 * @param readSVs Whether to attempt additional structural variant parsing.
	 * @throws UnsupportedFileTypeException If there is an error parsing the provided file as a VCF.
	 * @throws IOException If there is an error accessing the file on disk or retrieving data from the file.
	 */
	private void parseVCF(String vcfpath, boolean readSVs) throws UnsupportedFileTypeException, IOException
	{
		//Credit to http://www.avajava.com/tutorials/lessons/how-do-i-read-a-string-from-a-file-line-by-line.html
				// where I got the FileReader to BufferedReader tip from
		//System.err.println(Thread.currentThread().getName() + " || VCF.parseVCF || Called!");
		if (vcfpath == null || vcfpath.isEmpty()) throw new UnsupportedFileTypeException();
		FileReader freader = new FileReader(vcfpath);
		BufferedReader reader = new BufferedReader(freader);
		variants.setGenomeBuild(gBuild);
		if(gBuild != null) System.err.println(Thread.currentThread().getName() + " || VCF.parseVCF || Genome Build set by user: " + gBuild.getBuildName());
		else System.err.println(Thread.currentThread().getName() + " || VCF.parseVCF || Genome Build not set by user!");
		
		presetGenome = (gBuild != null);
		//if(gBuild != null) gBuild.printMe();
		
		int c = 0;
		
		boolean flipper = false;
		ParserThread t1 = new ParserThread(1);
		ParserThread t2 = new ParserThread(2);
		
		//Read header
		String line = null;
		try {line = reader.readLine();}
		catch(IOException e) {
			reader.close();
			freader.close();
			throw e;
		}
		
		while (line != null && (!line.isEmpty()))
		{
			if (line.charAt(0) == '#')
			{
				//It's a header line
				//Look for an = sign and split if find one. 
					//If not, then header line is unrecognized, or is the CHROM line
				if (line.length() > 2 && line.charAt(1) == '#')
				{
					if (line.indexOf('=') >= 0)
					{
						line = line.substring(line.lastIndexOf('#') + 1);
						String[] fields = line.split("=");
						if (fields == null) continue;
						if (fields.length < 2) {
							String k = fields[0];
							fields = new String[2];
							fields[0] = k;
							fields[1] = "";
						}
						
						if (fields[0].equals("INFO"))
						{
							InfoDefinition def = parseHeaderLine(line);
							String key = def.getKey();
							variants.addInfoFieldDefinition(key, def);
							variants.addInfoKeyToActiveList(key);
							//System.err.println(Thread.currentThread().getName() + " || VCF.parseVCF || INFO field found: " + key);
						}
						else if (fields[0].equals("FORMAT"))
						{
							InfoDefinition def = parseHeaderLine(line);
							String key = def.getKey();
							variants.addFormatFieldDefinition(key, def);
							variants.addFormatKeyToActiveList(key);
							//System.err.println(Thread.currentThread().getName() + " || VCF.parseVCF || FORMAT field found: " + key);
						}
						else if (fields[0].equals("ALT"))
						{
							String[] pieces = parseShortLine(line);
							if (pieces.length == 2)
							{
								variants.addCustomAlt(pieces[0], pieces[1]);
							}
							//System.err.println(Thread.currentThread().getName() + " || VCF.parseVCF || ALT field found: " + pieces[0]);
						}
						else if (fields[0].equals("FILTER"))
						{
							String[] pieces = parseShortLine(line);
							if (pieces.length == 2)
							{
								variants.addFilter(pieces[0], pieces[1]);
							}
							//System.err.println(Thread.currentThread().getName() + " || VCF.parseVCF || FILTER field found: " + pieces[0]);
						}
						else if (fields[0].equals("fileformat"))
						{
							fileformat = fields[1];
							//System.err.println(Thread.currentThread().getName() + " || VCF.parseVCF || Format: " + fields[1]);
						}
						else if (fields[0].equals("fileDate"))
						{
							dateString = fields[1];
							//System.err.println(Thread.currentThread().getName() + " || VCF.parseVCF || Date: " + fields[1]);
						}
						else if (fields[0].equals("source"))
						{
							source = fields[1];
							//System.err.println(Thread.currentThread().getName() + " || VCF.parseVCF || Source: " + fields[1]);
						}
						else if (fields[0].equals("reference") && (gBuild == null))
						{
							String ref = fields[1];
							System.err.println(Thread.currentThread().getName() + " || VCF.parseVCF || Reference Genome Build Detected: " + ref);
							if (!ref.isEmpty()) gBuild = GenomeBuild.loadStandardBuild(ref);
						}
						else
						{
							strayHeaderLines.add(fields);
						}
					}
					else
					{
						String[] arr = new String[1];
						arr[0] = line;
						strayHeaderLines.add(arr);
					}
				}
				else
				{
					//Assume it's the #CHROM etc. column header line
						//Extract sample names, if present
					String[] cols = line.split("\t");
					//System.err.println(Thread.currentThread().getName() + " || VCF.parseVCF || Header line found: " + line);
					//System.err.println(Thread.currentThread().getName() + " || VCF.parseVCF || Columns: " + cols.length);
					if (cols.length > 8)
					{
						//System.err.println(Thread.currentThread().getName() + " || VCF.parseVCF || This VCF contains genotypes calls.");
						for (int i = 9; i < cols.length; i++) variants.addSample(cols[i]);
					}
					List<String> slist1 = variants.getAllSamples();
					List<String> slist2 = variants.getAllSamples();
					//System.err.println(Thread.currentThread().getName() + " || VCF.parseVCF || Sample List: ");
					for (String s : slist1) System.err.println("\t" + s);
					//System.err.println(Thread.currentThread().getName() + " || VCF.parseVCF || Sample List Size: " + slist1.size());
					t1.setSampleList(slist1);
					t2.setSampleList(slist2);
					//Make a fillable genome build if one wasn't provided or detected.
					if (gBuild == null)
					{
						Random r = new Random();
						gBuild = new GenomeBuild("Unknown", "VCF_" + String.format("%08x", r.nextInt()), null);	
					}
					t1.start();
					t2.start();
				}
			}
			else
			{
				//Assume it's a variant record.
				if(flipper) t2.addLine(line);
				else t1.addLine(line);
				flipper = !flipper;
				if (parsingErr) {
					reader.close();
					freader.close();
					throw new UnsupportedFileTypeException();
				}
				c++;
			}
			/*if (c % 5000 == 0 && c != 0)
			{
				System.err.println("VCF.parseVCF || Variants read: " + c);
				int tot = 0;
				tot += t1.countParsed();
				tot += t2.countParsed();
				System.err.println("VCF.parseVCF || Variants parsed: " + tot);
				tot = 0;
				tot += t1.countRejected();
				tot += t2.countRejected();
				System.err.println("VCF.parseVCF || Variants rejected: " + tot);
			}*/
			try {line = reader.readLine();}
			catch(IOException e) {
				reader.close();
				freader.close();
				throw e;
			}
		}
		
		System.err.println("VCF.parseVCF || VCF Reading complete. Variants read: " + c);
		reader.close();
		freader.close();
		
		t1.killWhenDone();
		t2.killWhenDone();
		
		while(t1.isAlive() || t2.isAlive())
		{
			//int tot = 0;
			//tot += t1.countParsed();
			//tot += t2.countParsed();
			//System.err.println("VCF.parseVCF || Variants parsed: " + tot);
			//tot = 0;
			//tot += t1.countRejected();
			//tot += t2.countRejected();
			//System.err.println("VCF.parseVCF || Variants rejected: " + tot);
			try {
				Thread.sleep(5000);
			} catch (InterruptedException e) {
				Thread.interrupted();
				System.err.println("VCF.parseVCF || Main thread sleep interrupted. Checking worker thread status... ");
				//e.printStackTrace();
			}
		}
		
		System.err.println("VCF.parseVCF || Variant parsing complete!");
		int tot = 0;
		tot += t1.countParsed();
		tot += t2.countParsed();
		System.err.println("VCF.parseVCF || Variants parsed: " + tot);
		tot = 0;
		tot += t1.countRejected();
		tot += t2.countRejected();
		System.err.println("VCF.parseVCF || Variants rejected: " + tot);
		
		//Get the parsed variants.
		variants.addVariants(t1.getVariants());
		variants.addVariants(t2.getVariants());
		
		variants.sortVariants();
		System.err.println(Thread.currentThread().getName() + " || VCF.parseVCF || Variants found: " + variants.getVariants().size());
		
		if (readSVs)
		{
			variants.castStructuralVariants();
			System.err.println(Thread.currentThread().getName() + " VCF.parseVCF || SVs parsed. New variant count: " + variants.getVariants().size());
		}
		
	}
	
	private InfoDefinition parseHeaderLine(String line) throws UnsupportedFileTypeException
	{
		//System.err.println("VCF.parseHeaderLine || Line: " + line);
		//System.err.println("VCF.parseHeaderLine || < " + line.indexOf("<") + " > " + line.lastIndexOf(">"));
		line = line.substring(line.indexOf("<") + 1, line.lastIndexOf(">"));
		String[] fields = line.split(",");
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
			if (nArgs < -4) throw new UnsupportedFileTypeException();
			if (!key.isEmpty() && !type.isEmpty() && !desc.isEmpty())
			{
				int ifieldtype = VariantPool.getInfoDefType(type);
				InfoDefinition def = new InfoDefinition(key, ifieldtype, desc, nArgs);
				return def;
			}
		}
		return null;
	}
	
	private String[] parseShortLine(String line)
	{
		line = line.substring(line.indexOf("<") + 1, line.lastIndexOf(">"));
		String[] fields = line.split(",");
		String key = "";
		String desc = "";
		for (String s : fields)
		{
			String[] fsplit = s.split("=");
			if (fsplit[0].equals("ID")) key = fsplit[1];
			if (fsplit[0].equals("Description")) desc = fsplit[1];
			desc = desc.replaceAll("\"", "");
		}
		fields = new String[2];
		fields[0] = key;
		fields[1] = desc;
		return fields;
	}
	
	public static Map<String, String[]> mapINFOValues(String INFO_field)
	{
		Map<String, String[]> infoMap = new HashMap<String, String[]>();
		String[] fields = INFO_field.split(";");
		
		for(String s : fields)
		{
			String[] kv = s.split("=");
			if (kv.length != 2) continue;
			String key = kv[0];
			String val = kv[1];
			String[] vals = val.split(",");
			infoMap.put(key, vals);
		}
		
		return infoMap;
	}
	
	/* --- Inner Classes --- */
	
	/* --- Getters --- */
	
	/**
	 * Check the fileformat header line that was grabbed from the beginning of the file.
	 * If it doesn't indicate VCF 4.2, or there wasn't one, then this method returns false.
	 * @return True - If fileformat line read from the file matches the VCF 4.2 standard.
	 * <br>False - Otherwise.
	 */
	public boolean formatValid()
	{
		return (fileformat.equals("VCFv4.2"));
	}
	
	/**
	 * Get the string representing the date this VCF was originally written.
	 * It should be in YYYYMMDD format.
	 * @return 8 character unparsed string representing the date (YYYYMMDD), null
	 * if none was read from this file.
	 */
	public String getDateString()
	{
		return dateString;
	}
	
	/**
	 * Get the name of the application that wrote the VCF, if a header line saying 
	 * what it is is present.
	 * @return Source application name as found in the VCF file, null if none found.
	 */
	public String getSourceAppName()
	{
		return source;
	}
	
	/**
	 * Get the VariantPool contained within this object. Variant Pool should contain
	 * all variants, sample names, and info field definitions.
	 * @return Contained VariantPool, or null if there is none.
	 */
	public VariantPool getVariants()
	{
		return variants;
	}
	
	/**
	 * Get an ordered list of all of the field sets found as unparsed header lines in 
	 * this VCF.
	 * @return List (copy of internal collection) of unparsed header lines as field sets.
	 */
	public List<String[]> getStrayHeaderLines()
	{
		List<String[]> copy = new ArrayList<String[]>(strayHeaderLines.size());
		copy.addAll(strayHeaderLines);
		return copy;
	}
	
	/* --- Setters --- */
	
	/**
	 * Set the fileformat string for this VCF container to that used for the current
	 * version this class best recognizes.
	 * <br>Currently, this is version 4.2
	 */
	public void setFormatToCurrent()
	{
		fileformat = "VCFv4.2";
	}
	
	/**
	 * Update the fileDate string to match the date this method is called on.
	 */
	public void stampDate()
	{
		dateString = getHeaderDate();
	}
	
	/**
	 * Set the application name tag to the specified string. If null, then
	 * source application tag will be treated as unset.
	 * @param appName Name of application writing the VCF.
	 */
	public void setSourceApplicationName(String appName)
	{
		source = appName;
	}
	
	/**
	 * Add a one or two field header line as a String array. One field header
	 * lines will be copied to the header as is, two field header lines will
	 * have a '=' in between the key and value. Arrays longer than 2 will be not
	 * be added.
	 * @param fields Array of length 1 or 2 containing strings to use in a header field
	 * that is not explicitly specified otherwise.
	 */
	public void addUnparsedHeaderLine(String[] fields)
	{
		if (fields == null) return;
		if (fields.length > 2) return;
		strayHeaderLines.add(fields);
	}
	
	/**
	 * Clear the list of additional header lines for this VCF.
	 */
	public void clearUnparsedHeaderLines()
	{
		strayHeaderLines.clear();
	}
	
	/* --- Serialization --- */
	
	/**
	 * Generate a VCF fileDate header line style date string (YYYYMMDD) representing
	 * the date this method is called.
	 * @return YYYYMMDD formatted date string for the current date.
	 */
 	public String getHeaderDate()
	{
		Calendar date = new GregorianCalendar();
		int year = date.get(Calendar.YEAR);
		int month = date.get(Calendar.MONTH) + 1;
		int day = date.get(Calendar.DAY_OF_MONTH);
		String d = String.format("%d%02d%02d", year, month, day);
		return d;
	}
	
 	private String writeStrayHeaderLines()
 	{
 		//Includes newline at the end
 		String s = "";
 		for (String[] arr : strayHeaderLines)
 		{
 			if (arr.length == 1)
 			{
 				s += "##" + arr[0];
 			}
 			else if (arr.length == 2)
 			{
 				s += "##" + arr[0] + "=" + arr[1];
 			}
 			s += "\n";
 		}
 		return s;
 	}
 	
 	/**
 	 * Write this VCF to a file on the local disk system.
 	 * @param path The path to write the file to.
 	 * @throws IOException If there is an error accessing disk or writing to the provided
 	 * path.
 	 */
	public void writeToDisk(String path) throws IOException
	{
		FileWriter writer = new FileWriter(path);
		writer.write("##fileformat=" + fileformat + "\n");
		writer.write("##fileDate=" + dateString + "\n");
		if (source != null && !source.isEmpty())
		{
			source = source.replaceAll(" ", "_");
			writer.write("##source=" + source + "\n");
		}
		if (gBuild != null)
		{
			writer.write("##reference=" + gBuild.getBuildName() + "\n");
		}
		//Other header lines
		writer.write(writeStrayHeaderLines());
		//Info fields
		List<String> infolist = variants.getOrderedInfoKeys();
		for (String k : infolist)
		{
			InfoDefinition def = variants.getInfoDef(k);
			writer.write("##INFO=<" + def.formatForVCFHeader() + ">\n");
		}
		//Filter fields
		List<String> filterlist = variants.getOrderedFilterKeys();
		for (String k : filterlist)
		{
			String desc = variants.getFilterDescription(k);
			writer.write("##FILTER=<ID=" + k + ",Description=\"" + desc + "\">\n");
		}
		//Format fields
		List<String> genolist = variants.getOrderedGenotypeKeys();
		for (String k : genolist)
		{
			InfoDefinition def = variants.getGenoDef(k);
			writer.write("##FORMAT=<" + def.formatForVCFHeader() + ">\n");
		}
		//Alt fields
		List<String> altlist = variants.getAllAltKeys();
		for (String k : altlist)
		{
			String desc = variants.getAltDescription(k);
			writer.write("##ALT=<ID=" + k + ",Description=\"" + desc + "\">\n");
		}
		//Header line (#CHROM)
		List<String> samples = variants.getAllSamples();
		String cline = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
		if (samples != null && !samples.isEmpty())
		{
			cline += "\tFORMAT";
			for (String smpl : samples)
			{
				cline += "\t" + smpl;
			}
		}
		cline += "\n";
		writer.write(cline);
		//Variants
		List<Variant> varList = variants.getVariants();
		for (Variant v : varList)
		{
			List<InfoDefinition> deflist = variants.getOrderedInfoDefs();
			writer.write(v.toVCFLine(samples, deflist) + "\n");
		}
		//Close
		writer.close();
	}
	
	/* --- Static --- */
	
	/**
	 * Read a VCF file from disk and return a VariantPool representing the information
	 * in that VCF file.
	 * <br>WARNING: The VariantPool class does not contain fields for certain VCF specific metadata!
	 * Using this static method may result in the loss of information encoded in the header!
	 * @param path Path on local file system of file to read.
	 * @param readSVs Whether to additionally attempt additional parsing of structural variants
	 * to a more specialized structure.
	 * @return VariantPool representation of the VCF file contents, if parsing is successful.
	 * @throws IOException If the file could not be accessed on disk.
	 * @throws UnsupportedFileTypeException If the file could not be parsed as a VCF.
	 */
	public static VariantPool readVCF(String path, boolean readSVs) throws IOException, UnsupportedFileTypeException
	{
		VCF myVCF = new VCF(path, readSVs);
		return myVCF.variants;
	}
	
	/**
	 * Read a VCF file from disk and return a VariantPool representing the information
	 * in that VCF file.
	 * <br>WARNING: The VariantPool class does not contain fields for certain VCF specific metadata!
	 * Using this static method may result in the loss of information encoded in the header!
	 * @param path Path on local file system of file to read.
	 * @param genome Genome build to use for VCF parsing.
	 * @param readSVs Whether to additionally attempt additional parsing of structural variants
	 * to a more specialized structure.
	 * @return VariantPool representation of the VCF file contents, if parsing is successful.
	 * @throws IOException If the file could not be accessed on disk.
	 * @throws UnsupportedFileTypeException If the file could not be parsed as a VCF.
	 */
	public static VariantPool readVCF(String path, GenomeBuild genome, boolean readSVs) throws IOException, UnsupportedFileTypeException
	{
		VCF myVCF = new VCF(path, readSVs, genome);
		return myVCF.variants;
	}
	
	/**
	 * Write a variant pool out as a VCF file. Source application can be specified, but
	 * date is set to date the method is called and format is set to the current default 4.2.
	 * The only other VCF header metadata included is that which is in the VariantPool.
	 * <br>If additional header metadata is desired, then this static method should not be used.
	 * @param pool VariantPool to write.
	 * @param sourceApp Name of the application writing the VCF. This is optional.
	 * @param outpath Path on local file system to write VCF to.
	 * @throws IOException If the path is invalid or the disk cannot be written to.
	 */
	public static void writeVCF(VariantPool pool, String sourceApp, String outpath) throws IOException
	{
		VCF myVCF = new VCF(pool, sourceApp);
		myVCF.setFormatToCurrent();
		myVCF.stampDate();
		myVCF.writeToDisk(outpath);
	}
	
	public static Variant parseVCFLine(String line, List<String> genoSamples, GenomeBuild gb) throws UnsupportedFileTypeException
	{
		if(gb == null) return null;
		Variant v = new Variant();
		
		String[] fields = line.split("\t");
		try 
		{
			v.setChromosome(gb.getContig(fields[0]));
			v.setPosition(Integer.parseInt(fields[1]));
			v.setVariantName(fields[2]);
			v.setRefAllele(fields[3]);
			//Alt allele(s)
			String[] alts = fields[4].split(",");
			for(String a : alts) v.addAltAllele(a);
			//QUAL
			v.setQuality(Double.parseDouble(fields[5]));
			//FILTER
			if(fields[6].equals("PASS")) v.setFilterPass(true);
			else
			{
				v.setFilterPass(false);
				String[] filters = fields[6].split(";");
				if(filters != null && filters.length > 0)
				{
					for(String f : filters) v.addFailedFilter(f);
				}	
			}
			//INFO
			String[] infoFields = fields[7].split(";");
			if(infoFields != null)
			{
				for (String i : infoFields)
				{
					String[] kv = i.split("=");
					if(kv.length < 2)
					{
						v.addInfoFlag(i);
						continue;
					}
					String key = kv[0];
					String[] values = kv[1].split(",");
					v.addInfoField(key, values);
				}
			}
			//FORMAT
			if(fields.length < 9) return v;
			String formatString = fields[8];
			//Genotypes
			int i = 9;
			for(String s : genoSamples)
			{
				if (i >= fields.length) break;
				String rawGeno = fields[i];
				
				Genotype g = new Genotype(formatString, rawGeno);
				v.addGenotype(s, g);
				
				i++;
			}
		}
		catch(NullPointerException e)
		{
			e.printStackTrace();
			throw new FileBuffer.UnsupportedFileTypeException("VCF.parseVCFLine || NPE - Likely VCF record has insufficient tab separated fields");
		}
		catch(NumberFormatException e)
		{
			e.printStackTrace();
			throw new FileBuffer.UnsupportedFileTypeException("VCF.parseVCFLine || One or more integers could not be read as such!");
		}
		
		return v;
	}
	
	public static StructuralVariant parseVCFLineAsSV(String line, List<String> genoSamples, GenomeBuild gb) throws UnsupportedFileTypeException
	{
		if(gb == null) return null;
		StructuralVariant sv = null;
		
		String[] fields = line.split("\t");
		
		try 
		{
			//This time, we're gonna parse INFO first in case we need to
			// switch to TRA!
			
			Set<String> flags = new HashSet<String>();
			Map<String, String[]> infos = new HashMap<String, String[]>();
			
			String[] infoFields = fields[7].split(";");
			if(infoFields != null)
			{
				for (String i : infoFields)
				{
					String[] kv = i.split("=");
					if(kv.length < 2)
					{
						flags.add(i);
						continue;
					}
					String key = kv[0];
					String[] values = kv[1].split(",");
					infos.put(key, values);
				}
			}
			
			//Get the SV related INFO fields out of the way!
			String typeraw = "BND";
			String[] val = infos.remove(StructuralVariant.INFODEF_INFO_SVTYPE.getKey());
			if (val.length >= 1) typeraw = val[0];
			SVType t = SVType.getType(typeraw);
			
			//Initialize the SV
			if(t == SVType.TRA) sv = new Translocation();
			else sv = new StructuralVariant();
			sv.setType(t);
			
			//Parse the other interesting SV fields
				//End
				//Chr2
				//SVLEN
				//CI
				//MateID
				//Imprecise
				//Secondary
			
			sv.setChromosome(gb.getContig(fields[0]));
			sv.setPosition(Integer.parseInt(fields[1]));
			
			try
			{
				val = infos.remove(StructuralVariant.INFODEF_INFO_END.getKey());
				if(val != null && val.length >= 1)
				{
					int end = Integer.parseInt(val[0]);
					sv.setEndPosition(end);
				}
			}
			catch(NumberFormatException e)
			{
				e.printStackTrace();
				throw new FileBuffer.UnsupportedFileTypeException("VCF.parseVCFLineAsSV || SV end position could not be read! (Number parsing error)");
			}
			
			
			if(t == SVType.TRA)
			{
				val = infos.remove(Translocation.INFODEF_INFO_CHR2.getKey());
				if(val == null) sv.setEndChromosome(sv.getChromosome());
				else sv.setEndChromosome(gb.getContig(val[0]));
			}
			
			
			val = infos.remove(StructuralVariant.INFODEF_INFO_SVLEN.getKey());
			if(val != null)
			{
				try
				{
					for(int i = 0; i < val.length; i++)
					{
						int len = Integer.parseInt(val[i]);
						sv.setSVLength(i, len);
					}
				}
				catch(NumberFormatException e)
				{
					e.printStackTrace();
					throw new FileBuffer.UnsupportedFileTypeException("VCF.parseVCFLineAsSV || SVLEN could not be read! (Number parsing error)");
				}
			}
			
			try
			{
				val = infos.remove(StructuralVariant.INFODEF_INFO_CIPOS.getKey());
				if(val != null)
				{
					String[] rng = val[0].split(",");
					if(rng.length >= 2)
					{
						sv.setCIDiff(Integer.parseInt(rng[0]), false, false, false);
						sv.setCIDiff(Integer.parseInt(rng[1]), false, false, true);
					}
				}
				
				val = infos.remove(StructuralVariant.INFODEF_INFO_CIEND.getKey());
				if(val != null)
				{
					String[] rng = val[0].split(",");
					if(rng.length >= 2)
					{
						sv.setCIDiff(Integer.parseInt(rng[0]), true, false, false);
						sv.setCIDiff(Integer.parseInt(rng[1]), true, false, true);
					}
				}
				
				val = infos.remove(StructuralVariant.INFODEF_INFO_CIPOS95.getKey());
				if(val != null)
				{
					String[] rng = val[0].split(",");
					if(rng.length >= 2)
					{
						sv.setCIDiff(Integer.parseInt(rng[0]), false, true, false);
						sv.setCIDiff(Integer.parseInt(rng[1]), false, true, true);
					}
				}
				
				val = infos.remove(StructuralVariant.INFODEF_INFO_CIEND95.getKey());
				if(val != null)
				{
					String[] rng = val[0].split(",");
					if(rng.length >= 2)
					{
						sv.setCIDiff(Integer.parseInt(rng[0]), true, true, false);
						sv.setCIDiff(Integer.parseInt(rng[1]), true, true, true);
					}
				}
			}
			catch(NumberFormatException e)
			{
				e.printStackTrace();
				throw new FileBuffer.UnsupportedFileTypeException("VCF.parseVCFLineAsSV || CIs could not be read! (Number parsing error)");
			}
			
			val = infos.remove(StructuralVariant.INFODEF_INFO_MATEID.getKey());
			if(val != null) sv.addMate(val[0]);
			
			String fkey = StructuralVariant.INFODEF_INFO_IMPRECISE.getKey();
			if(flags.remove(fkey)) sv.setImprecise(true);
			else sv.setImprecise(false);
			
			fkey = StructuralVariant.INFODEF_INFO_SECONDARY.getKey();
			if(flags.remove(fkey)) sv.setSecondary(true);
			else sv.setSecondary(false);
			
			//Copy the remainder of the INFO fields
			
			if(!flags.isEmpty())
			{
				for(String f : flags) sv.addInfoFlag(f);
			}
			
			if(!infos.isEmpty())
			{
				for(String k : infos.keySet()) sv.addInfoField(k, infos.get(k));
			}
			
			
			/*--------------------*/
			
			sv.setVariantName(fields[2]);
			sv.setRefAllele(fields[3]);
			//Alt allele(s)
			String[] alts = fields[4].split(",");
			for(String a : alts) sv.addAltAllele(a);
			//QUAL
			if(fields[5].equals(".")) sv.setQuality(-1);
			else
			{
				try 
				{
					sv.setQuality(Double.parseDouble(fields[5]));
				}
				catch(NumberFormatException e)
				{
					e.printStackTrace();
					throw new FileBuffer.UnsupportedFileTypeException();
				}
			}
			//FILTER
			if(fields[6].equals("PASS")) sv.setFilterPass(true);
			else
			{
				sv.setFilterPass(false);
				String[] filters = fields[6].split(";");
				if(filters != null && filters.length > 0)
				{
					for(String f : filters) sv.addFailedFilter(f);
				}	
			}
			//INFO (Did before...)
			
			//FORMAT
			if(fields.length < 9) return sv;
			String formatString = fields[8];
			//Genotypes
			int i = 9;
			for(String s : genoSamples)
			{
				if (i >= fields.length) break;
				String rawGeno = fields[i];
				
				Genotype g = new Genotype(formatString, rawGeno);
				sv.addGenotype(s, g);
				
				i++;
			}
		}
		catch(NullPointerException e)
		{
			e.printStackTrace();
			throw new FileBuffer.UnsupportedFileTypeException("VCF.parseVCFLineAsSV || NPE - Likely VCF record has insufficient tab separated fields");
		}
		catch(NumberFormatException e)
		{
			e.printStackTrace();
			throw new FileBuffer.UnsupportedFileTypeException("VCF.parseVCFLineAsSV || One or more integers could not be read as such!");
		}
		
		return sv;
	}

}
