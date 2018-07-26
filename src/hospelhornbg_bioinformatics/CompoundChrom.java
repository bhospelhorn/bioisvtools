package hospelhornbg_bioinformatics;

import java.util.ArrayList;
import java.util.List;

import hospelhornbg_genomeBuild.Contig;

public class CompoundChrom extends Contig{

	//The only difference is that it contains multiple contigs and compares accordingly.
	
	private List<Contig> contigs;
	
	public CompoundChrom(int nChrom)
	{
		contigs = new ArrayList<Contig>(nChrom);
		super.setType(SORTCLASS_COMPOUND);
	}
	
	public void addChrom(Contig chrom)
	{
		contigs.add(chrom);
	}
	
	public void clearChrom()
	{
		contigs.clear();
	}
	
	public Contig getChrom(int index)
	{
		if (index < 0) return null;
		if (index >= contigs.size()) return null;
		return contigs.get(index);
	}
	
	public boolean equals(Object o)
	{
		if (o == null) return false;
		if (this == o) return true;
		if (!(o instanceof CompoundChrom)) return false;
		
		CompoundChrom c = (CompoundChrom)o;
		if (this.contigs.size() != c.contigs.size()) return false;
		if (!this.contigs.containsAll(c.contigs)) return false;
		if (!c.contigs.containsAll(this.contigs)) return false;
		
		return true;
	}
	
	public int hashCode()
	{
		int hash = 0;
		for (Contig c : contigs) hash ^= c.getUDPName().hashCode();
		return hash;
	}

	public int compareTo(Contig o) 
	{
		if (o == null) return 1;
		//First, is it a compound?
		if (this.getType() != o.getType()) return this.getType() - o.getType();
		
		//Else we assume that it is a compound
		if (!(o instanceof CompoundChrom)) return -1; //Something is typed wrong. Shove it to the bottom.
		CompoundChrom c = (CompoundChrom)o;
		
		//Compare sizes
		if (this.contigs.size() != c.contigs.size()) return this.contigs.size() - c.contigs.size();
		
		//Compare actual contigs
		for (int i = 0; i < this.contigs.size(); i++)
		{
			int comp = contigs.get(i).compareTo(c.getChrom(i));
			if (comp != 0) return comp;
		}
		
		return 0;
	}
	
	public String toString()
	{
		String s = "";
		for (int i = 0; i < contigs.size(); i++)
		{
			s += contigs.get(i).toString();
			if (i < contigs.size() - 1) s += ":";
		}
		return s;
	}
	
	
}
