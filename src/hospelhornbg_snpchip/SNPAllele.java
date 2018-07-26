package hospelhornbg_snpchip;

public class SNPAllele {
	
	private int allele;
	
	private boolean parent1;
	private boolean parent2;
	private boolean denovo;
	
	//private boolean isMaternal;
	//private boolean isPaternal;
	
	public SNPAllele(int a)
	{
		allele = a;
		parent1 = false;
		parent2 = false;
		denovo = false;
	}
	
	public int getAllele()
	{
		return allele;
	}
	
	public boolean phased_parent1()
	{
		return parent1;
	}
	
	public boolean phased_parent2()
	{
		return parent2;
	}
	
	public boolean flaggedDeNovo()
	{
		return denovo;
	}
	
	public boolean isPhased()
	{
		//Mm, don't think it should ever be flagged as both...
		return (parent1 ^ parent2);
	}
	
	public void setAllele(int a)
	{
		allele = a;
	}
	
	public void flagParent1(boolean b)
	{
		parent1 = b;
		if(b) parent2 = false;
	}
	
	public void flagParent2(boolean b)
	{
		parent2 = b;
		if(b) parent1 = false;
	}

	public void flagDeNovo(boolean b)
	{
		denovo = b;
		if (b)
		{
			parent1 = false;
			parent2 = false;
		}
	}
	
	public void clearPhasingFlags()
	{
		parent1 = false;
		parent2 = false;
		denovo = false;
	}
	
}
