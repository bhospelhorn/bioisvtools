package waffleoRai_Utils;


/**
 * Method collection for simple statistics calculations.
 * @author Blythe Hospelhorn
 * @version 1.0.0
 * @since January 4, 2018
 * 
 */
public class Statistics {
	
	public static double average(int[] values)
	{
		if (values == null) return Double.NaN;
		if (values.length < 1) return Double.NaN;
		int n = values.length;
		int sum = 0;
		for (int v : values) sum += v;
		return (double)sum / (double)n;
	}
	
	public static double average(double[] values)
	{
		if (values == null) return Double.NaN;
		if (values.length < 1) return Double.NaN;
		double n = (double)values.length;
		double sum = 0.0;
		for (double v : values) sum += v;
		return sum / n;
	}
	
	public static double stdev(int[] values)
	{
		if (values == null) return Double.NaN;
		if (values.length < 1) return Double.NaN;
		double avg = average(values);
		double n = (double)values.length;
		double sum = 0.0;
		for (int v : values)
		{
			double diff = (double)v - avg;
			sum += Math.pow(diff, 2.0);
		}
		return Math.sqrt(sum / (n - 1.0));
	}
	
	public static double stdev(double[] values)
	{
		if (values == null) return Double.NaN;
		if (values.length < 1) return Double.NaN;
		double avg = average(values);
		double n = (double)values.length;
		double sum = 0.0;
		for (double v : values)
		{
			double diff = v - avg;
			sum += Math.pow(diff, 2.0);
		}
		return Math.sqrt(sum / (n - 1.0));
	}

	private static double[] getCI(int[] values, final double crit)
	{
		if (values == null) return null;
		if (values.length < 1) return null;
		//final double crit = 1.645; //wikipedia
		int n = values.length;
		
		double mean = average(values);
		double stdev = stdev(values);
		
		double sqrtn = Math.sqrt((double)n);
		
		double bot = mean - (crit * (stdev/sqrtn));
		double top = mean + (crit * (stdev/sqrtn));
		
		double[] rng = new double[2];
		rng[0] = bot;
		rng[1] = top;
		
		return rng;
	}
	
	public static double[] getCI90(int[] values)
	{
		final double crit = 1.645; //wikipedia
		return getCI(values, crit);
	}
	
	public static double[] getCI95(int[] values)
	{
		final double crit = 1.96; //wikipedia
		return getCI(values, crit);
	}
	
	
}
