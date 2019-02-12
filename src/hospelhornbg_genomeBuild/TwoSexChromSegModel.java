package hospelhornbg_genomeBuild;

import java.awt.Point;
import java.util.ArrayList;
import java.util.List;

public class TwoSexChromSegModel {
	
	private Contig homogameticChrom;
	private Contig heterogameticChrom;
	
	private ArrayList<Point> homPARs;
	private ArrayList<Point> hetPARs;

	public TwoSexChromSegModel(Contig homgam, Contig hetgam, GenomeBuild gb)
	{
		homogameticChrom = homgam;
		heterogameticChrom = hetgam;
		
		//Get PAR definitions
		List<Point> pars1 = gb.getPARsForContig(homgam);
		List<Point> pars2 = gb.getPARsForContig(hetgam);
		
		homPARs = new ArrayList<Point>(pars1.size());
		homPARs.addAll(pars1);
		hetPARs = new ArrayList<Point>(pars2.size());
		hetPARs.addAll(pars2);
	}
	
	public Contig getHomogameticChrom()
	{
		return this.homogameticChrom;
	}
	
	public Contig getHeterogameticChrom()
	{
		return this.heterogameticChrom;
	}
	
	public boolean inHomChromPAR(int pos)
	{
		for(Point p : homPARs)
		{
			if (pos >= p.x && pos < p.y) return true;
		}
		return false;
	}
	
	public boolean inHetChromPAR(int pos)
	{
		for(Point p : hetPARs)
		{
			if (pos >= p.x && pos < p.y) return true;
		}
		return false;
	}
	
	public int mapHetPosToHom(int pos)
	{
		int pari = -1;
		int offset = -1;
		int i = 0;
		for(Point p : hetPARs)
		{
			if (pos >= p.x && pos < p.y)
			{
				pari = i;
				offset = pos - p.x;
				break;
			}
			i++;
		}
		if(pari < 0) return -1;
		
		Point other = homPARs.get(pari);
		
		return offset + other.x;
	}
	
	public int mapHomPosToHet(int pos)
	{
		int pari = -1;
		int offset = -1;
		int i = 0;
		for(Point p : homPARs)
		{
			if (pos >= p.x && pos < p.y)
			{
				pari = i;
				offset = pos - p.x;
				break;
			}
			i++;
		}
		if(pari < 0) return -1;
		
		Point other = hetPARs.get(pari);
		
		return offset + other.x;
	}

}
