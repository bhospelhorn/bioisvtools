package hospelhornbg_bioinformatics;

import java.util.ArrayList;
import java.util.List;

import waffleoRai_Utils.Statistics;

public class LumpyEvidence {
	
	private int highestTotal;
	private double SU_avg;
	private double SU_stdev;
	
	private boolean PE_include;
	private boolean SR_include;
	private boolean BD_include;
	
	private double PE_avg;
	private double PE_stdev;
	private int PE_max;
	private double SR_avg;
	private double SR_stdev;
	private int SR_max;
	private double BD_avg;
	private double BD_stdev;
	private int BD_max;
	
	public LumpyEvidence(VariantPool initialPool, boolean includePE, boolean includeSR, boolean includeBD)
	{
		calibrate(initialPool);
		PE_include = includePE;
		SR_include = includeSR;
		BD_include = includeBD;
	}
	
	private void calibrate(VariantPool pool)
	{
		//Because of filtering, the total score will have to be
			//the difference from "max" total, as some may be over
		setDefaultValues();
		if (pool == null) return;
		final List<StructuralVariant> svList = filterPool(pool);
		int[] allPE = new int[svList.size()];
		int[] allSR = new int[svList.size()];
		int[] allBD = new int[svList.size()];
		int[] allSU = new int[svList.size()];
		int i = 0;
		for (StructuralVariant sv : svList)
		{
			int PE = sv.getPEEvidenceCount();
			int SR = sv.getSREvidenceCount();
			int BD = sv.getBDEvidenceCount();
			int SU = sv.getTotalEvidenceCount();
			allPE[i] = PE;
			allSR[i] = SR;
			allBD[i] = BD;
			allSU[i] = SU;
			if (SU > highestTotal) highestTotal = SU;
			if (PE > PE_max) PE_max = PE;
			if (SR > SR_max) SR_max = SR;
			if (BD > BD_max) BD_max = BD;
			i++;
		}
		PE_avg = Statistics.average(allPE);
		PE_stdev = Statistics.stdev(allPE);
		SR_avg = Statistics.average(allSR);
		SR_stdev = Statistics.average(allSR);
		BD_avg = Statistics.average(allBD);
		BD_stdev = Statistics.average(allBD);
		SU_avg = Statistics.average(allSU);
		SU_stdev = Statistics.stdev(allSU);
	}
	
	private void setDefaultValues()
	{
		highestTotal = -1;
		PE_avg = -1;
		PE_stdev = -1;
		SR_avg = -1;
		SR_stdev = -1;
		BD_avg = -1;
		BD_stdev = -1;
		PE_max = -1;
		SR_max = -1;
		BD_max = -1;
	}
	
	private List<StructuralVariant> filterPool(VariantPool pool)
	{
		final List<Variant> vList = pool.getVariants();
		List<StructuralVariant> svList = new ArrayList<StructuralVariant>(vList.size());
		for (Variant v : vList)
		{
			if (v instanceof StructuralVariant)
			{
				StructuralVariant sv = (StructuralVariant)v;
				if (sv.getTotalEvidenceCount() >= 0)
				{
					if (sv.getType() != SVType.BND)
					{
						if (sv.getSVLength() >= 50 && sv.getSVLength() <= 5000000)
						{
							svList.add(sv);
						}
					}
				}
			}
		}
		return svList;
	}
	
	private static double singleTypeScore(int eCount, double avg, double stdev, int max)
	{
		//TODO: May need to normalize
		if (stdev == 0) return 1.0;
		final double TSTDEV_FACTOR = stdev;
		final double TDIST_FACTOR = 2.0;
		final double ADIST_FACTOR = 2.0;
		final double TBASE = 2.0;
		final double ABASE = 1.1;
		final double AWEIGHT = 0.35;
		
		double d = (double)eCount;
		
		double aDist = Math.pow((avg - d) / stdev, ADIST_FACTOR);
		double tDist = Math.pow(((double)max - d)/(stdev * TSTDEV_FACTOR), TDIST_FACTOR);
		
		double aScore = Math.pow(ABASE, aDist * -1.0);
		double tScore = Math.pow(TBASE, tDist * -1.0);
		double score = (AWEIGHT * aScore) + ((1 - AWEIGHT) * tScore);
		
		return score;
	}
	
	public double calc_PE_score(StructuralVariant sv)
	{
		//A measure of the PE evidence contribution
			//Values close to the max and average are scored the highest
		return singleTypeScore(sv.getPEEvidenceCount(), PE_avg, PE_stdev, PE_max);
	}
	
	public double calc_SR_score(StructuralVariant sv)
	{
		//A measure of the SR evidence contribution
			//Values close to the max and average are scored the highest
		return singleTypeScore(sv.getSREvidenceCount(), SR_avg, SR_stdev, SR_max);
	}
	
	public double calc_BD_score(StructuralVariant sv)
	{
		//A measure of the BD evidence contribution
			//Values close to the max and average are scored the highest
		return singleTypeScore(sv.getBDEvidenceCount(), BD_avg, BD_stdev, BD_max);
	}
	
	public double calc_PSR_score(StructuralVariant sv)
	{
		//A measure of the lopsidedness of PE-SR evidence counts
			//Variants with a near equal amount of PE and SR evidence score the highest
		int PE = sv.getPEEvidenceCount();
		int SR = sv.getSREvidenceCount();
		/*int tot = PE + SR;
		if (tot == 0) return 1.0; //Technically, it's not lopsided
		double fPE = (double)PE / (double)tot;
		double fSR = (double)SR / (double)tot;
		double score = Math.pow((fPE - fSR), 2.0);
		score = 1.0 - score;	
		if (score < 0.0) score = 0.0;*/
		double avgRatio = PE_avg/ SR_avg;
		double ratio = (double)PE / (double)SR;
		double score = Math.pow((avgRatio - ratio), 2.0);
		score = 1.0 - score;
		if (score < 0.0) score = 0.0;
		
		return score;
	}
	
	public double calc_SU_score(StructuralVariant sv)
	{
		//A measure of the all evidence relative to evidence amounts in pool
		//Values close to the max and average are scored the highest
	return singleTypeScore(sv.getTotalEvidenceCount(), SU_avg, SU_stdev, highestTotal);
	}
	
	public double calc_RawLumpyEvidenceScore(StructuralVariant sv)
	{
		double initScore = 0.0;
		if (PE_include && SR_include)
		{
			double cAvg = (calc_PE_score(sv) + calc_SR_score(sv)) / 2.0;
			double pairScore = calc_PSR_score(sv);
			double tScore = calc_SU_score(sv);
			
			final double CWEIGHT = 0.2;
			final double PWEIGHT = 0.7;
			final double TWEIGHT = 0.1;
			
			initScore = (cAvg * CWEIGHT) + (tScore * TWEIGHT) + (pairScore * PWEIGHT);
		}
		else if (PE_include || SR_include)
		{
			double cAvg = 0.5;
			if (PE_include) cAvg = calc_PE_score(sv);
			else cAvg = calc_SR_score(sv);
			double tScore = calc_SU_score(sv);
			
			final double CWEIGHT = 0.6;
			final double TWEIGHT = 0.4;
			
			initScore = (cAvg * CWEIGHT) + (tScore * TWEIGHT);
		}
		else
		{
			initScore = calc_SU_score(sv);
		}
		
		//Figure BED bonus
		if (BD_include)
		{
			double bScore = calc_BD_score(sv);
			final double BWEIGHT = 0.1;
			
			double bweight = BWEIGHT * bScore;
			double score = ((1.0 - bweight) * initScore) + (bweight * bScore);
			
			return score;
		}
		else return initScore;
		
	}

	public int calc_LumpyEvidenceScore(StructuralVariant sv)
	{
		return (int)Math.round((this.calc_RawLumpyEvidenceScore(sv) * 100.0));
	}

}
