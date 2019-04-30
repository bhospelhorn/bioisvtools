package hospelhornbg_genomeBuild;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;

public class OMIMGeneMapImporter {
	
	public static final String ANNO_KEY = "OMIM_PHENO";
	
	private String table_path;
	
	public OMIMGeneMapImporter(String tablePath)
	{
		table_path = tablePath;
	}
	
	public boolean importTable(GeneSet genes)
	{
		if(genes == null) return false;
		if(table_path == null || table_path.isEmpty()) return false;
		
		final int GENES_COLUMN = 6;
		final int PHENO_COLUMN = 12;
		
		BufferedReader br;
		try 
		{
			br = new BufferedReader(new FileReader(table_path));
			String line = null;
			while((line = br.readLine()) != null)
			{
				if(line.isEmpty()) continue;
				if(line.startsWith("#")) continue;
				String[] fields = line.split("\t");
				//System.err.println("Line: " + line);
				//System.err.println("Field count: " + fields.length);
				if(fields.length <= GENES_COLUMN) continue;
				if(fields.length <= PHENO_COLUMN) continue;
				//System.err.println("Pheno pass!: " + fields.length);
				String pheno = fields[PHENO_COLUMN];
				if(pheno.isEmpty()) continue;
				//System.err.println("Pheno: " + pheno);
				String genelist = fields[GENES_COLUMN];
				//System.err.println("Genelist: " + genelist);
				String gclean = genelist.replace(" ", ""); //Delete spaces
				String[] glsplit = gclean.split(",");
				//System.err.println("Genelist (list): ");
				for(String gname : glsplit)
				{
					//System.err.println("\t" + gname);
					List<Gene> cand = genes.getGeneByName(gname);
					if (cand == null || cand.isEmpty()) continue;
					for(Gene g : cand) {
						//System.err.println("Gene Match Found: " + g.getName());
						g.addAnnotation(ANNO_KEY, pheno);
					}
				}
			}
			br.close();
		}
		catch(IOException e) {e.printStackTrace(); return false;}
		
		
		
		return true;
	}

}
