package hospelhornbg_svproject;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.Gene;
import hospelhornbg_genomeBuild.GeneSet;
import hospelhornbg_genomeBuild.GenomeBuild;

public class GenomeUIDTable {
	
	private Map<Integer, Contig> contigMap;
	private Map<Integer, Gene> transcriptMap;
	
	public GenomeUIDTable()
	{
		contigMap = new HashMap<Integer, Contig>();
		transcriptMap = new HashMap<Integer, Gene>();
	}
	
	public void loadContigMap(GenomeBuild gb, String tablepath) throws IOException
	{
		if (gb == null) return;
		
		FileReader fr = new FileReader(tablepath);
		BufferedReader br = new BufferedReader(fr);
		
		String line = null;
		while((line = br.readLine()) != null)
		{
			if (line.startsWith("#")) continue;
			String[] fields = line.split("\t");
			if (fields.length != 2) continue;
			int id = -1;
			String cname = null;
			try
			{
				//id = Integer.parseInt(fields[0]);
				id = Integer.parseUnsignedInt(fields[0], 16);
			}
			catch (NumberFormatException e)
			{
				continue;
			}
			cname = fields[1];
			if (id != -1)
			{
				Contig c = gb.getContig(cname);
				if (c != null)
				{
					contigMap.put(id, c);
				}
			}
		}
		
		br.close();
		
	}
	
	public void loadTranscriptMap(GeneSet gs, String tablepath) throws IOException
	{
		if (gs == null) return;
		
		//Remap genes by transcript name...
		List<Gene> allgenes = gs.getAllGenes();
		Map<String, Gene> tmap = new HashMap<String, Gene>();
		for (Gene g : allgenes)
		{
			tmap.put(g.getID(), g);
		}
		
		FileReader fr = new FileReader(tablepath);
		BufferedReader br = new BufferedReader(fr);
		
		String line = null;
		while((line = br.readLine()) != null)
		{
			if (line.startsWith("#")) continue;
			String[] fields = line.split("\t");
			if (fields.length != 2) continue;
			int id = -1;
			String tname = null;
			try
			{
				//id = Integer.parseInt(fields[0]);
				id = Integer.parseUnsignedInt(fields[0], 16);
			}
			catch (NumberFormatException e)
			{
				continue;
			}
			tname = fields[1];
			if (id != -1)
			{
				//Look up transcript...
				Gene g = tmap.get(tname);
				if (g != null)
				{
					transcriptMap.put(id, g);
				}
			}
		}
		
		br.close();
		
	}

	public Contig getContigByUID(int id)
	{
		return contigMap.get(id);
	}
	
	public Gene getTranscriptByUID(int id)
	{
		return transcriptMap.get(id);
	}

	public static void generateContigUIDTable(GenomeBuild gb, String outpath) throws IOException
	{
		List<Contig> clist = gb.getChromosomes();
		
		FileWriter fw = new FileWriter(outpath);
		BufferedWriter bw = new BufferedWriter(fw);
		
		bw.write("#" + gb.getBuildName() + "\n");
		
		for(Contig c : clist)
		{
			String name = c.getUCSCName();
			int id = name.hashCode();
			bw.write(Integer.toHexString(id) + "\t" + name + "\n");
		}
		
		bw.close();
	}
	
	public static void generateTranscriptUIDTable(GeneSet gs, String outpath) throws IOException
	{
		List<Gene> tlist = gs.getAllGenes();
		
		FileWriter fw = new FileWriter(outpath);
		BufferedWriter bw = new BufferedWriter(fw);
		
		for(Gene t : tlist)
		{
			String name = t.getID();
			int id = name.hashCode();
			bw.write(Integer.toHexString(id) + "\t" + name + "\n");
		}
		
		bw.close();
	}
	
}
