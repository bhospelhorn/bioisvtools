package hospelhornbg_genomeBuild;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public enum GeneFunc {
	
	EXONIC("exonic", 0, false),
	SPLICING("splicing", 1, false),
	NCRNA("ncRNA", 2, false),
	UTR5("UTR5", 3, false),
	UTR3("UTR3", 4, false),
	INTRONIC("intronic", 5, false),
	UPSTREAM("upstream", 6, true),
	DOWNSTREAM("downstream", 7, true),
	INTERGENIC("intergenic", 8, true);
	
	private String annotation;
	private int p;
	private boolean intergenic;
	
	private GeneFunc(String anno, int priority, boolean ig)
	{
		annotation = anno;
		p = priority;
		intergenic = ig;
	}
	
	public String toString()
	{
		return annotation;
	}
	
	public int getPriority()
	{
		return p;
	}
	
	public boolean isIntergenic()
	{
		return intergenic;
	}

	public static GeneFunc getHighestPriority()
	{
		return EXONIC;
	}
	
	public static GeneFunc getLowestPriority()
	{
		return INTERGENIC;
	}

	public static GeneFunc getFunction(String infoValue)
	{
		if (infoValue == null) return null;
		if (infoValue.isEmpty()) return null;
		
		if (infoValue.equals(NCRNA.toString())) return NCRNA;
		if (infoValue.equals(EXONIC.toString())) return EXONIC;
		if (infoValue.equals(INTRONIC.toString())) return INTRONIC;
		if (infoValue.equals(INTERGENIC.toString())) return INTERGENIC;
		if (infoValue.equals(UTR3.toString())) return UTR3;
		if (infoValue.equals(UTR5.toString())) return UTR5;
		if (infoValue.equals(DOWNSTREAM.toString())) return DOWNSTREAM;
		if (infoValue.equals(UPSTREAM.toString())) return UPSTREAM;
		if (infoValue.equals(SPLICING.toString())) return SPLICING;
		
		return null;
	}

	public static Collection<GeneFunc> getAll()
	{
		List<GeneFunc> all = new ArrayList<GeneFunc>(9);
		all.add(EXONIC);
		all.add(SPLICING);
		all.add(NCRNA);
		all.add(UTR5);
		all.add(UTR3);
		all.add(INTRONIC);
		all.add(UPSTREAM);
		all.add(DOWNSTREAM);
		all.add(INTERGENIC);
		
		return all;
	}
	
}
