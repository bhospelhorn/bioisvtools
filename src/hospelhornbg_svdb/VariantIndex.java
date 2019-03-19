package hospelhornbg_svdb;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.GenomeBuild;

public class VariantIndex {

	//Gives the number (1 based index) of the first line
	//of entries for each chrom.
	//For TRA variants, sorted by start chr.
	
	//Saved on disk as a tsv (chrom	line)
	
	public static final int CHR1_COLUMN_INDEX = 3;
	
	private ConcurrentHashMap<Contig, Integer> lineMap;
	
	private VariantIndex()
	{
		lineMap = new ConcurrentHashMap<Contig, Integer>();
	}
	
	public int getFirstLineOfContig(Contig c)
	{
		Integer i = lineMap.get(c);
		if (i == null) return -1;
		return i;
	}
	
	public static VariantIndex buildIndexFromTable(String tablePath, GenomeBuild gb) throws IOException
	{
		if(gb == null) return null;
		VariantIndex idx = new VariantIndex();
		BufferedReader br = new BufferedReader(new FileReader(tablePath));
		int c = 0;
		String line = null;
		Contig lastc = null;
		while((line = br.readLine()) != null)
		{
			c++;
			if (line.isEmpty()) continue;
			if (line.charAt(0) == '#') continue;
			//Get contig, always in the third column
			String[] fields = line.split("\t");
			if (fields.length < CHR1_COLUMN_INDEX + 1) continue;
			String ctg = fields[CHR1_COLUMN_INDEX];
			Contig mychr = gb.getContig(ctg);
			if (lastc == null || !lastc.equals(mychr))
			{
				lastc = mychr;
				idx.lineMap.put(lastc, c);
			}
		}
		br.close();
		return idx;
	}
	
	public static VariantIndex readIndexFromDisk(String indexPath, GenomeBuild gb) throws IOException
	{
		if(gb == null) return null;
		VariantIndex idx = new VariantIndex();
		BufferedReader br = new BufferedReader(new FileReader(indexPath));
		String line = null;
		while((line = br.readLine()) != null)
		{
			if (line.isEmpty()) continue;
			String[] fields = line.split("\t");
			if (fields.length < 2) continue;
			String cname = fields[0];
			String i = fields[1];
			try
			{
				Contig c = gb.getContig(cname);
				if (c == null) continue;
				int n = Integer.parseInt(i);
				idx.lineMap.put(c, n);
			}
			catch(NumberFormatException e)
			{
				System.err.println("VariantIndex.readIndexFromDisk || Could not read line. Line number \"" + i + "\" invalid integer.");
				continue;
			}
		}
		br.close();
		return idx;
	}
	
	public void writeIndexToDisk(String path) throws IOException
	{
		Set<Contig> keyset = lineMap.keySet();
		List<Contig> clist = new LinkedList<Contig>();
		clist.addAll(keyset);
		Collections.sort(clist);
		BufferedWriter bw = new BufferedWriter(new FileWriter(path));
		for(Contig c : clist)
		{
			int i = lineMap.get(c);
			bw.write(c.getUDPName() + "\t" + i + "\n");
		}
		bw.close();
	}
	
}
