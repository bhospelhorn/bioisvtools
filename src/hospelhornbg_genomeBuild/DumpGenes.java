package hospelhornbg_genomeBuild;

import java.io.IOException;

public class DumpGenes {

	public static void main(String[] args) {
		
		if (args.length < 1)
		{
			System.err.println("Genome build name required.");
			System.exit(1);
		}
		if (args.length < 2)
		{
			System.err.println("Output path required.");
			System.exit(1);
		}
		
		
		String bname = args[0];
		String outpath = args[1];
		
		System.err.println("Loading genes for build " + bname);
		GeneSet gs = GeneSet.loadRefGene(bname);
		if (gs == null)
		{
			System.err.println("No gene set found for build " + bname + ". Now exiting...");
			System.exit(1);
		}
		
		System.err.println("Writing table to " + outpath);
		try 
		{
			gs.outputTable(outpath);
		} 
		catch (IOException e) 
		{
			System.err.println("ERROR: Table could not be written to " + outpath);
			System.exit(1);
		}
		System.err.println("Table successfully dumped to " + outpath);

	}

}
