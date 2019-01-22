package hospelhornbg_svtools;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.List;

import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.GenomeBuild;

public class GenomeIntervals {

	public static void generateIntervals(GenomeBuild gb, int intervalSize, String outpath) throws IOException
	{
		if (gb == null) return;
		if (intervalSize < 1) return;
		
		List<Contig> clist = gb.getChromosomes();
		Collections.sort(clist);
		
		//Open stream
		FileWriter fw = new FileWriter(outpath);
		BufferedWriter bw = new BufferedWriter(fw);
		
		for (Contig c : clist)
		{
			if (c.getUDPName().equals("M")) continue;
			if (c.getUDPName().equals("GL000250.1")) break;
			System.err.println("Breaking down contig " + c.getUDPName());
			long len = c.getLength();
			int start = 0;
			int end = intervalSize;
			while(start < len)
			{
				bw.write(c.getUDPName() + "\t" + start + "\t" + end + "\n");
				start = end;
				end += intervalSize;
				if (end > len) end = (int)len;
			}
		}
		
		bw.close();
	}
	
	public static void main(String[] args) 
	{
		String gbpath = "Z:\\svref\\bioi\\GRCh37.gbdh";
		String outpath = "Z:\\svref\\gatk\\grch37_50bp_intervals.bed";
		int interval = 50;
		
		try 
		{
		GenomeBuild gb = new GenomeBuild(gbpath);
		generateIntervals(gb, interval, outpath);
		}
		catch(Exception e)
		{
			System.err.println("Exception caught. Terminating...");
		}

	}

}
