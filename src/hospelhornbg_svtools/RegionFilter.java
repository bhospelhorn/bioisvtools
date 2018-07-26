package hospelhornbg_svtools;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import hospelhornbg_bioinformatics.VCF;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class RegionFilter {

	//Need an input BED file (to know which regions)
	//Input VCF to filter
	//Output path
	//Inclusive/exclusive - for SVs, should filter out if one or both end points falls in a bad region?
	
	public static final String OP_VCFIN = "-i"; 
	public static final String OP_VCFOUT = "-o"; 
	public static final String OP_BEDIN = "-r"; 
	public static final String OP_EXCLUSIVE = "-x";
	
	private static class Region implements Comparable<Region>
	{
		public Contig chrom;
		public int start;
		public int end;
		
		public boolean equals(Object o)
		{
			if (o == null) return false;
			if (o == this) return true;
			if (!(o instanceof Region)) return false;
			Region r = (Region)o;
			if (this.start != r.start) return false;
			if (this.end != r.end) return false;
			if (!this.chrom.equals(r.chrom)) return false;
			return true;
		}
		
		public int hashCode()
		{
			return chrom.hashCode() ^ start ^ end;
		}
		
		@Override
		public int compareTo(Region o) {
			if (o == null) return 1;
			if (o == this) return 0;
			if (!this.chrom.equals(o.chrom))
			{
				return this.chrom.compareTo(o.chrom);
			}
			if (this.start != o.start) return this.start - o.start;
			if (this.end != o.end) return this.end - o.end;
			
			return 0;
		}
	}
	
	public static void printUsage()
	{
		
	}
	
	public static List<Region> readBED(String bedpath, GenomeBuild gb) throws IOException, UnsupportedFileTypeException
	{
		List<Region> rlist = new LinkedList<Region>();
		
		FileReader fr = new FileReader(bedpath);
		BufferedReader br = new BufferedReader(fr);
		
		int i = 0;
		String line = null;
		while ((line = br.readLine()) != null)
		{
			i++;
			if (line.isEmpty()) continue;
			if (line.charAt(0) == '#') continue;
			String[] fields = line.split("\t");
			if (fields.length < 3){
				br.close();
				fr.close();
				throw new FileBuffer.UnsupportedFileTypeException();
			}
			try
			{
				Region r = new Region();
				r.chrom = gb.getContig(fields[0]);
				r.start = Integer.parseInt(fields[1]);
				r.end = Integer.parseInt(fields[2]);
				rlist.add(r);
			}
			catch (NumberFormatException e)
			{
				System.err.println("RegionFilter.readBED || ERROR: Line " + i + " could not be read!");
				br.close();
				fr.close();
				throw new FileBuffer.UnsupportedFileTypeException();
			}
		}
		
		br.close();
		fr.close();
		
		Collections.sort(rlist);
		return rlist;
	}
	
	public static void runFilter(VariantPool variants, List<Region> outregions, boolean exclusive)
	{
		List<Variant> vlist = variants.getVariants();
		List<Variant> passed = new LinkedList<Variant>();
		boolean good = false;
		
		for (Variant v : vlist)
		{
			good = true;
			for (Region r : outregions)
			{
				boolean inregion = v.inRegion(r.chrom, r.start, r.end, exclusive);	
				if (inregion)
				{
					good = false;
					break;
				}
			}
			if (good) passed.add(v);
		}
		
		variants.clearVariants();
		variants.addVariants(passed);
		variants.sortVariants();
		
	}
	
	public static void filterRegions(String[] args, GenomeBuild gb)
	{
		String inFile = null;
		String outFile = null;
		String bedFile = null;
		boolean ex = false;
		
		for (int i = 0; i < args.length; i++)
		{
			String s = args[i];
			if (s.equals(OP_VCFIN))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_VCFIN + " flag MUST be followed by input VCF path!");
					printUsage();
					System.exit(1);
				}
				inFile = args[i+1];
			}
			else if (s.equals(OP_VCFOUT))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_VCFOUT + " flag MUST be followed by output VCF path!");
					printUsage();
					System.exit(1);
				}
				outFile = args[i+1];
			}
			else if (s.equals(OP_BEDIN))
			{
				if (i+1 >= args.length)
				{
					System.err.println("ERROR: " + OP_BEDIN + " flag MUST be followed by input BED file!");
					printUsage();
					System.exit(1);
				}
				bedFile = args[i+1];
			}
			else if (s.equals(OP_EXCLUSIVE))
			{
				ex = true;
			}
		}
		
		
		boolean pass = true;
		if (inFile == null || inFile.isEmpty())
		{
			System.err.println("ERROR: Input path is required!");
			pass = false;
		}
		if (outFile == null || outFile.isEmpty())
		{
			System.err.println("ERROR: Output path is required!");
			pass = false;
		}
		if (bedFile == null || bedFile.isEmpty())
		{
			System.err.println("ERROR: BED file path is required!");
			pass = false;
		}
		if (!pass)
		{
			printUsage();
			System.exit(1);
		}
		
		List<Region> rlist = null;
		VariantPool pool = null;
		
		try 
		{
			pool = VCF.readVCF(inFile, gb, true);
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: Input VCF file " + inFile + " could not be opened!");
			e.printStackTrace();
			System.exit(1);
		} 
		catch (UnsupportedFileTypeException e) {
			System.err.println("ERROR: Input VCF file " + inFile + " could not be read!");
			e.printStackTrace();
			System.exit(1);
		}
		
		try 
		{
			rlist = readBED(bedFile, gb);
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: Input BED file " + bedFile + " could not be opened!");
			e.printStackTrace();
			System.exit(1);
		} 
		catch (UnsupportedFileTypeException e) 
		{
			System.err.println("ERROR: Input BED file " + bedFile + " could not be read!");
			e.printStackTrace();
			System.exit(1);
		}
		
		runFilter(pool, rlist, ex);
		
		try 
		{
			VCF.writeVCF(pool, "bioisvtools", outFile);
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: Output VCF file " + outFile + " could not be written!");
			e.printStackTrace();
			System.exit(1);
		}
		
	}
	
}
