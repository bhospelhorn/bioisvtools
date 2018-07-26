package hospelhornbg_bioinformatics;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Collection;
import java.util.Random;

import hospelhornbg_bioinformatics.VariantPool.InfoDefinition;
import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

/*
 * UPDATE NOTES
 * 
 * Initial version date: December 14, 2017 (1.0.0)
 * 
 * 1.0.0 -> 1.1.0 | January 16, 2018
 * 		Javadoc annotation
 * 		Tweaked BEDPE parser
 * 		Allowed for parsing of all fields, not just default/required
 * 		Added more methods for conversion and access
 * 		Overhauled InfoDef adding - should now be handled by StructuralVariant class
 * 
 * 1.1.0 -> 1.2.0 | February 6, 2018
 * 		Made parsing to StructuralVariant cleaner - set more fields.
 * 		Made input stream safer (using Java classes instead)
 * 		
 * 1.2.0 -> 1.2.1 | February 8, 2018
 * 		Fixed a dumb in the BEDPE parser.
 * 
 * 1.2.1 -> 1.2.2 | February 9, 2018
 * 		Added a handler in the parser for stupidly formatted headers.
 * 
 * 1.2.2 -> 1.3.0 | February 20, 2018
 * 		Update to be compatible with Contig object chromosomes instead of strings.
 */

/*
 * Possible future updates:
 * 	- Write LUMPY BEDPE parser/ writer ?
 */

/**
 * Parsing/ serialization manager for the BED variant format.
 * @author Blythe Hospelhorn
 * @version 1.3.0
 * @since February 20, 2018
 */
public class BED {
	
	/* --- Constants --- */
	
	public enum BEDType
	{
		Standard,
		BEDPE_AllBNDPairs,
		BEDPE_Default,
		BEDPE_Lumpy;
	}
	
	/* --- Instance Variables --- */
	
	/**
	 * Collection of parsed variants
	 */
	private VariantPool variants;
	
	/* --- Construction/Parsing --- */
	
	/**
	 * Constructor for creating a new BED container from a variant pool.
	 * Pool is referenced by this new BED, not copied.
	 * @param varPool Variant pool to construct a BED container around.
	 */
	public BED(VariantPool varPool)
	{
		variants = varPool;
	}
	
	/**
	 * Constructor used for opening a file at the given path as a BED file and
	 * parsing out the variants.
	 * <br>Variants are by default parsed as Structural Variants of the "CNV"
	 * type with confidence intervals of size 0 (precise).
	 * <br>This parser only handles the required BED fields.
	 * @param path String representing the path on the file system to the file to open as BED.
	 * @throws IOException If file cannot be found or opened.
	 * @throws UnsupportedFileTypeException If file cannot be parsed as a BED file.
	 */
	public BED(String path, GenomeBuild genome) throws IOException, UnsupportedFileTypeException
	{
		this(path, SVType.CNV, genome);
	}
	
	/**
	 * Constructor used for opening a file at the given path as a BED file and
	 * parsing out the variants.
	 * <br>Variants are by default parsed as Structural Variants of the given type
	 * (defaulting to BED Region if null) with confidence intervals of size 0 (precise).
	 * <br>This parser only handles the required BED fields. 
	 * @param path String representing the path on the file system to the file to open as BED.
	 * @param varType Structural Variant type to parse variants as.
	 * @throws IOException If file cannot be found or opened.
	 * @throws UnsupportedFileTypeException If file cannot be parsed as a BED file.
	 */
	public BED(String path, SVType varType, GenomeBuild genome) throws IOException, UnsupportedFileTypeException
	{
		this(path, varType, BEDType.Standard, genome);
	}
	
	/**
	 * Constructor used for opening a file at the given path as a BED or BEDPE file and
	 * parsing out the variants.
	 * <br>Variants are by default parsed as Structural Variants of the given type
	 * (defaulting to BED Region if null) with confidence intervals of size 0 (precise).
	 * <br>This parser only handles the required BED or BEDPE fields.  
	 * @param path String representing the path on the file system to the file to open as BED or BEDPE.
	 * @param varType Structural Variant type to parse variants as.
	 * @param bedpe Whether input file is in BEDPE format (true), or standard BED format (false).
	 * @throws IOException If file cannot be found or opened.
	 * @throws UnsupportedFileTypeException If file cannot be parsed as the requested file type.
	 */
	public BED(String path, SVType varType, BEDType parser, GenomeBuild genome) throws IOException, UnsupportedFileTypeException
	{
		variants = new VariantPool(32);
		//FileBuffer bed = FileBuffer.createBuffer(path);
		if (varType == null) varType = SVType.BED_REGION;
		parseBED(path, varType, parser, genome);
	}
	
	/**
	 * Core parser. Reads as tab delimited text file line by line (one line per variant).
	 * @param bed BED file as a wrapped byte array.
	 * @param varType Structural Variant type to label variants as.
	 * @param bedpe Whether input file is in BEDPE format (true), or standard BED format (false).
	 * @throws UnsupportedFileTypeException If file cannot be parsed as the requested file type.
	 * @throws IOException If the file could not be found, or there was an error with the stream.
	 */
	private void parseBED(String path, SVType varType, BEDType parser, GenomeBuild genome) throws UnsupportedFileTypeException, IOException
	{
		populateINFODefs();
		Random rand = new Random();
		String bedname = "BED" + Integer.toHexString(rand.nextInt());
		int c = 1;
		if (parser == BEDType.BEDPE_AllBNDPairs) defaultParseBEDPE(path, varType, true, genome);
		else if (parser == BEDType.BEDPE_Default) defaultParseBEDPE(path, varType, false, genome);
		else if (parser == BEDType.BEDPE_Lumpy) throw new UnsupportedFileTypeException();
		else
		{
			FileReader fr = new FileReader(path);
			BufferedReader br = new BufferedReader(fr);
			String line = br.readLine();
			while (line != null)
			{
				if (line.charAt(0) == '#')
				{
					line = br.readLine();
					continue;
				}
				String[] fields = line.split("\t");
				try
				{
					String chrom = fields[0];
					if (chrom.indexOf("chr") >= 0)
					{
						chrom = chrom.replaceFirst("chr", "");
					}
					String name = bedname + "_" + c;
					int start = -1;
					int end = -1;
					try
					{
						start = Integer.parseInt(fields[1]);
						end = Integer.parseInt(fields[2]);
					}
					catch (NumberFormatException e)
					{
						//It could be a stupid header line that is improperly formatted
						line = br.readLine();
						continue;
					}
					
					start++;
					if (fields.length > 3) name = fields[3];	
					StructuralVariant sv = new StructuralVariant();
					sv.setChromosome(genome.getContig(chrom));
					sv.setPosition(start);
					sv.setEndPosition(end);
					sv.setType(varType);
					sv.setRefAllele("N");
					sv.addAltAllele("<" + varType.toString() + ">");
					sv.setVariantName(name);
					sv.setSVLength(end - start);
					if (sv.getType() == SVType.DEL) sv.setSVLength(start - end);
					variants.addVariant(sv);
				}
				catch(ArrayIndexOutOfBoundsException | NumberFormatException e)
				{
					e.printStackTrace();
					fr.close();
					br.close();
					throw new UnsupportedFileTypeException();
				}
				line = br.readLine();
				c++;
			}
			fr.close();
			br.close();
		}
	}
	
	private void defaultParseBEDPE(String path, SVType varType, boolean pairAll, GenomeBuild genome) throws UnsupportedFileTypeException, IOException
	{
		Random rand = new Random();
		String bedname = "BEDPE" + Integer.toHexString(rand.nextInt());
		int c = 1;
		InfoDefinition def = StructuralVariant.INFODEF_INFO_CIPOS;
		variants.addInfoFieldDefinition(def.getKey(), def);
		variants.addInfoKeyToActiveList(def.getKey());
		def = StructuralVariant.INFODEF_INFO_CIEND;
		variants.addInfoFieldDefinition(def.getKey(), def);
		variants.addInfoKeyToActiveList(def.getKey());
		def = StructuralVariant.INFODEF_INFO_SECONDARY;
		variants.addInfoFieldDefinition(def.getKey(), def);
		variants.addInfoKeyToActiveList(def.getKey());
		def = StructuralVariant.INFODEF_INFO_EVENT;
		variants.addInfoFieldDefinition(def.getKey(), def);
		variants.addInfoKeyToActiveList(def.getKey());
		def = StructuralVariant.INFODEF_INFO_MATEID;
		variants.addInfoFieldDefinition(def.getKey(), def);
		variants.addInfoKeyToActiveList(def.getKey());
		
		FileReader fr = new FileReader(path);
		BufferedReader br = new BufferedReader(fr);
		String line = br.readLine();
		//int c = 1;
		while (line != null)
		{
			if (line.charAt(0) == '#')
			{
				line = br.readLine();
				continue;
			}
			String[] fields = line.split("\t");
			try
			{
				String chrom1 = fields[0];
				if (chrom1.indexOf("chr") >= 0)
				{
					chrom1 = chrom1.replaceFirst("chr", "");
				}
				int start1 = -1;
				int end1 = -1;
				try
				{
					start1 = Integer.parseInt(fields[1]);
					end1 = Integer.parseInt(fields[2]);
				}
				catch (NumberFormatException e)
				{
					//It could be a stupid header line that is improperly formatted
					line = br.readLine();
					continue;
				}
				String chrom2 = fields[3];
				if (chrom2.indexOf("chr") >= 0)
				{
					chrom2 = chrom2.replaceFirst("chr", "");
				}
				int start2 = Integer.parseInt(fields[4]);
				int end2 = Integer.parseInt(fields[5]);
				String name = bedname + "_"  + c;
				int qual = -1;
				if (fields.length >= 8)
				{
					try
					{
						qual = Integer.parseInt(fields[7]);
					}
					catch (NumberFormatException e)
					{
						qual = -1;
					}
				}
				start1++; //BED coordinates are 0
				start2++; //BED coordinates are 0
				if (fields.length > 6) name = fields[6];
				if (pairAll || !chrom1.equals(chrom2))
				{
					StructuralVariant var1 = new StructuralVariant();
					StructuralVariant var2 = new StructuralVariant();
					var1.setRefAllele("N");
					var2.setRefAllele("N");
					var1.setType(SVType.BND);
					var2.setType(SVType.BND);
					int pos1 = (end1 + start1) / 2;
					int pos2 = (end2 + start2) / 2;
					var1.addAltAllele("N]" + chrom2 + ":" + pos2 + "]");
					var2.addAltAllele("[" + chrom1 + ":" + pos1 + "[N");
					var1.setChromosome(genome.getContig(chrom1));
					var2.setChromosome(genome.getContig(chrom2));
					var1.setPosition(pos1);
					var2.setPosition(pos2);
					if (qual >= 0) var1.setQuality(qual);
					if (qual >= 0) var2.setQuality(qual);
					var1.setCIDiff(start1-pos1, false, false, false);
					var1.setCIDiff(end1-pos1, false, false, true);
					var1.setCIDiff(start1-pos1, false, true, false);
					var1.setCIDiff(end1-pos1, false, true, true);
					var2.setCIDiff(start2-pos2, false, false, false);
					var2.setCIDiff(end2-pos2, false, false, true);
					var2.setCIDiff(start2-pos2, false, true, false);
					var2.setCIDiff(end2-pos2, false, true, true);
					var1.setImprecise(true);
					var2.setImprecise(true);
					var2.setSecondary(true);
					if (name.isEmpty() || name == null || name.equals(".")) name = Integer.toString(c);
					var1.setVariantName(name + "_1");
					var2.setVariantName(name + "_2");
					var1.setEventID(Integer.toString(c));
					var2.setEventID(Integer.toString(c));
					var1.addMate(name + "_2");
					var2.addMate(name + "_1");
					//Are BND end positions taken care of?
					BreakendPair bp = new BreakendPair(var1, var2);
					bp.setVariantName(name);
					variants.addVariant(bp);
				}
				else
				{
					StructuralVariant sv = new StructuralVariant();
					sv.setChromosome(genome.getContig(chrom1));
					sv.setRefAllele("N");
					sv.addAltAllele("<" + varType.toString() + ">");
					sv.setVariantName(name);
					sv.setType(varType);
					int st = (start1 + end1)/2;
					int ed = (start2 + end2)/2;
					sv.setSVLength(ed - st);
					if (sv.getType() == SVType.DEL) sv.setSVLength(st - ed);
					sv.setPosition(st);
					sv.setEndPosition(ed);
					sv.setCIDiff(start1 - st, false, false, false);
					sv.setCIDiff(end1 - st, false, false, true);
					sv.setCIDiff(start2 - ed, true, false, false);
					sv.setCIDiff(end2 - ed, true, false, true);
					sv.setCIDiff(start1 - st, false, true, false);
					sv.setCIDiff(end1 - st, false, true, true);
					sv.setCIDiff(start2 - ed, true, true, false);
					sv.setCIDiff(end2 - ed, true, true, true);
					sv.setImprecise(true);
					variants.addVariant(sv);
				}
			}
			catch(ArrayIndexOutOfBoundsException | NumberFormatException e)
			{
				e.printStackTrace();
				fr.close();
				br.close();
				throw new UnsupportedFileTypeException();
			}
			line = br.readLine();
			c++;
		}
		fr.close();
		br.close();
	}

	private void populateINFODefs()
	{
		variants.addCustomAlt(StructuralVariant.INFODEF_ALT_DEL, StructuralVariant.getAltDefDescription(StructuralVariant.INFODEF_ALT_DEL));
		variants.addCustomAlt(StructuralVariant.INFODEF_ALT_DUP, StructuralVariant.getAltDefDescription(StructuralVariant.INFODEF_ALT_DUP));
		variants.addCustomAlt(StructuralVariant.INFODEF_ALT_INS, StructuralVariant.getAltDefDescription(StructuralVariant.INFODEF_ALT_INS));
		variants.addCustomAlt(StructuralVariant.INFODEF_ALT_INV, StructuralVariant.getAltDefDescription(StructuralVariant.INFODEF_ALT_INV));
		variants.addCustomAlt(StructuralVariant.INFODEF_ALT_CNV, StructuralVariant.getAltDefDescription(StructuralVariant.INFODEF_ALT_CNV));
		variants.addCustomAlt(StructuralVariant.INFODEF_ALT_TANDEM, StructuralVariant.getAltDefDescription(StructuralVariant.INFODEF_ALT_TANDEM));
		variants.addCustomAlt(StructuralVariant.INFODEF_ALT_DELME, StructuralVariant.getAltDefDescription(StructuralVariant.INFODEF_ALT_DELME));
		variants.addCustomAlt(StructuralVariant.INFODEF_ALT_INSME, StructuralVariant.getAltDefDescription(StructuralVariant.INFODEF_ALT_INSME));
		
		InfoDefinition def = StructuralVariant.INFODEF_INFO_SVTYPE;
		variants.addInfoFieldDefinition(def.getKey(), def);
		variants.addInfoKeyToActiveList(def.getKey());
		def = StructuralVariant.INFODEF_INFO_SVLEN;
		variants.addInfoFieldDefinition(def.getKey(), def);
		variants.addInfoKeyToActiveList(def.getKey());
		def = StructuralVariant.INFODEF_INFO_END;
		variants.addInfoFieldDefinition(def.getKey(), def);
		variants.addInfoKeyToActiveList(def.getKey());
		
		def = StructuralVariant.INFODEF_INFO_IMPRECISE;
		variants.addInfoFieldDefinition(def.getKey(), def);
		variants.addInfoKeyToActiveList(def.getKey());
		

	}
	
	/* --- Inner Classes --- */
	
	/* --- Getters --- */
	
	/**
	 * Get a reference to the variant pool used by this BED container.
	 * @return VariantPool object reference that container is using. In other words,
	 * the pool of variants from this container.
	 */
	public VariantPool getVariantPool()
	{
		return variants;
	}
	
	/* --- Setters --- */
	
	/* --- Serialization --- */
	
	/**
	 * Write the variant pool contained by this BED object out to disk at
	 * the specified path in BED format.
	 * <br>IMPORTANT: This method only writes structural variants to the BED
	 * file. SNV's and smaller variants are ignored.
	 * @param path Path to write file to.
	 * @throws FileNotFoundException If path is invalid.
	 */
	public void writeBED(String path) throws FileNotFoundException
	{
		if (variants == null) return;
		writeBED(variants.getVariants(), path);
	}
	
	/**
	 * Write the variant pool contained by this BED object out to disk at
	 * the specified path in BEDPE format.
	 * <br>IMPORTANT: This method only writes structural variants to the new
	 * file. SNV's and smaller variants are ignored.
	 * @param path Path to write file to.
	 * @throws FileNotFoundException If path is invalid.
	 */
	public void writeBEDPE(String path) throws FileNotFoundException
	{
		if (variants == null) return;
		writeBEDPE(variants.getVariants(), path);
	}
	
	/* --- Static --- */

	/**
	 * Read a file on disk and parse it as a BED file. Return the contents as a VariantPool
	 * of structural variants.
	 * @param path Path on local file system to file to read.
	 * @param type Type of structural variant to set all regions read in the BED as. Although
	 * this field must be set when calling this method, it can be ignored by the user when looking
	 * through variants.
	 * @param parser The BED parser to call. There is the standard BED parser as well as different
	 * ways to read BEDPE files.
	 * @return A VariantPool object containing a collection of the regions/variants found in the
	 * BED file as StructuralVariant objects.
	 * @throws IOException If the file could not be accessed, either because there is a problem
	 * with the disk or because the path string is invalid.
	 * @throws UnsupportedFileTypeException If the file could not be read be the specified parser.
	 */
	public static VariantPool readBED(String path, SVType type, BEDType parser, GenomeBuild genome) throws IOException, UnsupportedFileTypeException
	{
		BED myfile = new BED(path, type, parser, genome);
		return myfile.getVariantPool();
	}
	
	/**
	 * Write all structural variants in a collection of variants (other variants are ignored)
	 * to a file on disk in BED format.
	 * @param variants Collection of variants to include in the BED output.
	 * @param path Path of file to write.
	 * @throws FileNotFoundException If the output file path is invalid.
	 */
	public static void writeBED(Collection<Variant> variants, String path) throws FileNotFoundException
	{
		if (variants == null) return;
		PrintStream ps = new PrintStream(path);
		for (Variant v : variants)
		{
			if (v instanceof StructuralVariant)
			{
				StructuralVariant sv = (StructuralVariant)v;
				ps.println(sv.toBEDLine());
			}
		}
		ps.close();
	}
	
	/**
	 * Write all structural variants in a collection of variants (other variants are ignored)
	 * to a file on disk in BEDPE format.
	 * @param variants Collection of variants to include in the BED output.
	 * @param path Path of file to write.
	 * @throws FileNotFoundException If the output file path is invalid.
	 */
	public static void writeBEDPE(Collection<Variant> variants, String path) throws FileNotFoundException
	{
		if (variants == null) return;
		PrintStream ps = new PrintStream(path);
		for (Variant v : variants)
		{
			if (v instanceof StructuralVariant)
			{
				StructuralVariant sv = (StructuralVariant)v;
				ps.println(sv.toBEDPELine());
			}
		}
		ps.close();
	}
	
}
