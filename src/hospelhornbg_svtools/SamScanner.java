package hospelhornbg_svtools;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import hospelhornbg_bioinformatics.SAMHeaderLine;
import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.FileBuffer;

public class SamScanner {
	
	public static final String OP_INPUT = "-i"; 
	public static final String OP_TEMPDIR = "-T"; 
	public static final String OP_THREADS = "-t"; 
	
	private static class Counts
	{
		private volatile int err_syntax_general;
		private volatile int err_nullseq_nnqual;
		private volatile int err_bad_customFieldType;
		private Set<String> bad_customFieldTypes;
		
		private volatile int err_invalid_RNAME;
		private volatile int err_invalid_POS;
		private volatile int err_invalid_RNEXT;
		private volatile int err_invalid_PNEXT;
		private Set<String> bad_RNAME;
		private Set<String> bad_RNEXT;
		
		private volatile int err_qualstr_len_bad;
		
		public Counts()
		{
			err_syntax_general = 0;
			err_nullseq_nnqual = 0;
			err_bad_customFieldType = 0;
			bad_customFieldTypes = new HashSet<String>();
			
			err_invalid_RNAME = 0;
			err_invalid_POS = 0;
			err_invalid_RNEXT = 0;
			err_invalid_PNEXT = 0;
			bad_RNAME = new HashSet<String>();
			bad_RNEXT = new HashSet<String>();
			
			err_qualstr_len_bad = 0;
		}
		
		public void incrementSyntaxErrCount()
		{
			err_syntax_general++;
		}
	}
	
	private static Contig generateContig(SAMHeaderLine sq)
	{
		if (sq.getFieldCount("SN") > 1)
		{
			//Duplication! Bad!
			System.out.println("ERROR   || @SQ line has multiple sequence name records!");
			System.out.println("\t" + sq);
			return null;
		}
		String cname = sq.getValue("SN");
		if (sq.getFieldCount("LN") > 1)
		{
			//Duplication! Bad!
			System.out.println("ERROR   || @SQ line has multiple sequence length records!");
			System.out.println("\t" + sq);
			return null;
		}
		String clenstr = sq.getValue("LN");
		long clen = -1;
		if(clenstr != null)
		{
			try
			{
				clen = Long.parseUnsignedLong(clenstr);
			}
			catch(NumberFormatException e)
			{
				System.out.println("ERROR   || @SQ LN value invalid (" + clenstr + ")");
				System.out.println("\t" + sq);
				return null;
			}
		}
		if (cname == null)
		{
			System.out.println("ERROR   || @SQ line has no SN field!");
			System.out.println("\t" + sq);
			return null;
		}
		else
		{
			if (clen < 0)
			{
				System.out.println("ERROR   || @SQ line has no LN field!");
				System.out.println("\t" + sq);
				return null;
			}
			else
			{
					Contig c = new Contig();	
					c.setUDPName(cname);
					c.setUCSCName(cname);
					c.setLength(clen);
					return c;
			}
		}
	}
	
	private static GenomeBuild scanSamHeader(BufferedReader reader, GenomeBuild gb) throws IOException
	{
		//Read header
		//Check for HD line (which is mandatory)
		//If there is a provided GenomeBuild, report any SQ lines that have a contig not in genome build
		//	and any genome build contigs not in SQ (remove these from local gb copy so that warnings will be thrown if these contigs show up
		//	in records
		//If there is no provided GenomeBuild, build one from the SQ lines and regurgitate for user
		//Check for RG line
		System.out.println("===================== HEADER =====================");
		String line = reader.readLine();
		//This line should be the @HD. If not, report.
		if (!(line.startsWith("@HD")))
		{
			System.out.println("ERROR   || @HD line was NOT found!");
			if (!line.startsWith("@")) {
				System.out.println("ERROR   || SAM lacks header!");
				reader.reset();
			}
		}
		else
		{
			System.out.println("INFO   || @HD line found!");
			String[] fields = line.split("\t");
			int fc = fields.length;
			for (int i = 1; i < fc; i++) System.out.println("\t" + fields[i]);
		}
		
		Set<SAMHeaderLine> sqlines = new HashSet<SAMHeaderLine>();
		boolean rgfound = false;
		boolean rgsmfound = false;
		line = reader.readLine();
		while(line.startsWith("@"))
		{
			SAMHeaderLine hl = null;
			try {hl = new SAMHeaderLine(line);}
			catch(FileBuffer.UnsupportedFileTypeException e) {
				System.out.println("ERROR   || Syntax Error - Improperly formatted header line:");
				System.out.println("\t" + line);
				continue;
			};
			
			if (hl.getKey().equals("SQ")) sqlines.add(hl);
			else if (hl.getKey().equals("RG"))
			{
				rgfound = true;
				String sm = hl.getValue("SM");
				if (sm != null)
				{
					System.out.println("INFO   || @RG line found! Sample Name: " + sm);
					rgsmfound = true;
				}
			}
			line = reader.readLine();
		}
		
		//Handle RG
		if (!rgfound) System.out.println("WARNING || No @RG line was found! Some packages require this line to read a SAM/BAM file.");
		if (rgfound && !rgsmfound) System.out.println("WARNING || @RG line(s) found, but there were no SM (sample) fields. Some packages require this to obtain the sample name!");
		
		//Handle SQ set
		
		if (gb != null)
		{
			Set<Contig> matched = new HashSet<Contig>();
			for (SAMHeaderLine sq : sqlines)
			{
				Contig c = generateContig(sq);
				if (c != null)
				{
					String cname = c.getUDPName();
					Contig o = gb.getContig(cname);
					if (o != null)
					{
						//Lengths must match!
						long mylen = c.getLength();
						long olen = o.getLength();
						if (mylen != olen)
						{
							System.out.println("ERROR   || Contig found in build, but length is mismatched: " + cname + " [Build:" + olen + "][Local:" + mylen + "]");
						}
						else matched.add(o);
					}
					else
					{
						//No match
						System.out.println("ERROR   || Contig not found in build: " + cname);
					}
				}
			}
			
			//Find contigs in build that are NOT in SAM
			List<Contig> allc = gb.getChromosomes();
			for (Contig c : allc)
			{
				if(!matched.contains(c))
				{
					System.out.println("WARNING || Build contig " + c.getUDPName() + " was not found in SAM file header!");
					gb.removeContig(c.getUDPName());
				}
			}
			
			reader.reset();
			return gb;
		}
		else
		{
			GenomeBuild ngb = new GenomeBuild("UNKNOWN", "UNKNOWN");
			for (SAMHeaderLine sq : sqlines)
			{
				Contig c = generateContig(sq);
				if (c != null)
				{
					String cname = c.getUDPName();
					Contig o = ngb.getContig(cname);
					if (o != null)
					{
						System.out.println("ERROR   || Duplicate contig (" + cname + ")");
					}
					else ngb.addContig(c);
				}
			}
			
			reader.reset();
			return ngb;
		}
	}
	
	private static Counts getCounts(BufferedReader reader, GenomeBuild gb, boolean verbose)
	{
		return null;
	}
	
	private static Counts getCounts(BufferedReader reader, GenomeBuild gb, boolean verbose, int threads)
	{
		return null;
	}
	
	public static void scanSam(BufferedReader reader, GenomeBuild gb, boolean verbose, int threads) throws IOException
	{
		//Read header (and reset stream)
		gb = scanSamHeader(reader, gb);
		
		//Dump contig list to stdout before looking at reads
		System.out.println("--------- Generated Genome Build --------- ");
		gb.printMe();
		
		//Scan reads to look for major syntactical errors
		Counts c = null;
		if (threads < 2) c = getCounts(reader, gb, verbose);
		else c = getCounts(reader, gb, verbose, threads);
		
		//Print stats on errors found
		
	}
	
	public static void runSamScanner(String[] args, GenomeBuild gb, boolean verbose)
	{
		//Genome build may be null. If it is null, genome build is generated from SAM header
		
	}

}
