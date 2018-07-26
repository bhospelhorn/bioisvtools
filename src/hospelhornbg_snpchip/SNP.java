package hospelhornbg_snpchip;

import hospelhornbg_genomeBuild.Contig;

public class SNP implements Comparable<SNP>{
	
	private String ID;
	
	private Contig chrom;
	private int position;

	private int alleleCount;
	
	public int hashCode()
	{
		return 0;
	}
	
	public boolean equals(Object o)
	{
		return false;
	}

	public int compareTo(SNP o) 
	{
	
		return 0;
	}

}
