package hospelhornbg_bioinformatics;

public class Utilities {

	public static double phredToProb(int phredScore)
	{
		double p;
		double fac = ((double)phredScore / 10.0) * -1.0;
		p = Math.pow(10.0, fac);
		return p;
	}
	
	public static int probToPhred(double probability)
	{
		double phred = (-10.0) * Math.log10(probability);
		int i = (int)Math.round(phred);
		return i;
	}
	
	public static double floatScaleToProb(double scaled)
	{
		double val = Math.pow(10.0, scaled); // May lose precision - can't really go smaller than 10^ -300 or so.
		return val;
	}
	
	public static double probToFloatScale(double probability)
	{
		double val = Math.log10(probability);
		return val;
	}
	
	public static double phredToFloatScale(int phredScore)
	{
		double p = (double)phredScore / (-10.0);
		return p;
	}
	
	public static int floatScaleToPhred(double scaled)
	{
		double p = scaled * -10.0;
		int i = (int)Math.round(p);
		return i;
	}
	
}
