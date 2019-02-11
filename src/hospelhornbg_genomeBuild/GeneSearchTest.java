package hospelhornbg_genomeBuild;

import java.util.List;

public class GeneSearchTest {

	public static void main(String[] args) {
		
		String gbPath = "C:\\Users\\Blythe\\Documents\\GitHub\\bioisvtools\\src\\hospelhornbg_genomeBuild\\resources\\GRCh37.gbdh";
		String gsPath = "C:\\Users\\Blythe\\Documents\\GitHub\\bioisvtools\\src\\hospelhornbg_genomeBuild\\resources\\grch37_refSeq.gbgd";
		
		String ctg = "6";
		int pos = 162329287;
		int end = 162396659;
		
		try
		{
			System.out.println("Loading Genome Build...");
			GenomeBuild gb = new GenomeBuild(gbPath);
			//gb.printMe();
			System.out.println("Loading Gene Set...");
			GeneSet gs = new GeneSet(gsPath, gb, true);	
			
			System.out.println("Getting Genes...");
			List<Gene> glist = gs.getGenesInRegion(gb.getContig(ctg), pos, end);
			if (glist == null)
			{
				System.out.println("No genes found :(");
			}
			else
			{
				System.out.println("Genes found: " + glist.size());
				int i = 0;
				for (Gene g : glist)
				{
					System.out.println("Gene " + i + ": " + g.getName());
					System.out.println(g.printInfo());
					i++;
				}
			}
			
		}
		catch (Exception e)
		{
			System.err.println("Exception!");
			e.printStackTrace();
		}
		

	}

}
