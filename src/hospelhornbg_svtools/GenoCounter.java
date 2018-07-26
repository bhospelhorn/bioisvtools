package hospelhornbg_svtools;

import java.util.ArrayList;
import java.util.List;

import hospelhornbg_bioinformatics.Genotype;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_bioinformatics.VariantPool;

public class GenoCounter {
	
	private List<String> sampleList;
	private int sampleCount;
	private int[] tally;
	private int[] homTally;
	private int[] hetTally;
	private int[] otherTally;
	
	public GenoCounter(VariantPool pool)
	{
		sampleCount = -1;
		if (pool == null) return;
		runCount(pool);
	}
	
	private void runCount(VariantPool pool)
	{
		List<String> samplelist = pool.getAllSamples();
		if (samplelist == null) return;
		if (samplelist.isEmpty()) return;
		List<Variant> varList = pool.getVariants();
		
		sampleCount = samplelist.size();
		if (sampleCount > 16) throw new IllegalArgumentException();
		sampleList = new ArrayList<String>(sampleCount);
		sampleList.addAll(samplelist);
		
		int comboCount = (int)Math.pow(2, sampleCount);
		tally = new int[comboCount];
		homTally = new int[comboCount];
		hetTally = new int[comboCount];
		otherTally = new int[comboCount];
		
		
		for (int i = 0; i < comboCount; i++)
		{
			tally[i] = 0;
			homTally[i] = 0;
			hetTally[i] = 0;
			otherTally[i] = 0;
		}
		
		for (Variant v: varList)
		{
			int vecT = 0;
			int vecHM = 0;
			int vecHT = 0;
			int vecO = 0;
			
			for (String s : samplelist)
			{
				int zy = v.getSampleGenotype(s).getZygosity();
				switch(zy)
				{
				case Genotype.ZYGOSITY_HOMOREF: break; //Nothing
				case Genotype.ZYGOSITY_HOMOALT:
					vecT |= 0x1;
					vecHM |= 0x1;
					break;
				case Genotype.ZYGOSITY_HETERORA:
					vecT |= 0x1;
					vecHT |= 0x1;
					break;
				case Genotype.ZYGOSITY_HETEROAA:
					vecT |= 0x1;
					vecO |= 0x1;
					break;
				case Genotype.ZYGOSITY_UNKNOWN: break; //Nothing
				case Genotype.ZYGOSITY_HOMOREF_CNV: break; //Nothing
				case Genotype.ZYGOSITY_HOMOALT_CNV:
					vecT |= 0x1;
					vecHM |= 0x1;
					vecO |= 0x1;
					break;
				case Genotype.ZYGOSITY_HETERORA_CNV:
					vecT |= 0x1;
					vecHT |= 0x1;
					vecO |= 0x1;
					break;
				case Genotype.ZYGOSITY_HETEROAA_CNV:
					vecT |= 0x1;
					vecO |= 0x1;
					break;
				case Genotype.ZYGOSITY_HETERORAA_CNV:
					vecT |= 0x1;
					vecO |= 0x1;
					break;
				default: break; //Nothing
				}
				
				vecT = vecT << 1;
				vecHM = vecHM << 1;
				vecHT = vecHT << 1;
				vecO = vecO << 1;
			}
			
			//Tally record for each vector
			tally[vecT]++;
			homTally[vecHM]++;
			hetTally[vecHT]++;
			otherTally[vecO]++;
			
		}
	}
	
	public int getSampleNumber()
	{
		return sampleCount;
	}
	
	public int getPresentTally(int vector)
	{
		if (vector < 0) return 0;
		if (vector > tally.length) return 0;
		return tally[vector];
	}
	
	public int getHomTally(int vector)
	{
		if (vector < 0) return 0;
		if (vector > homTally.length) return 0;
		return homTally[vector];
	}
	
	public int getHetTally(int vector)
	{
		if (vector < 0) return 0;
		if (vector > hetTally.length) return 0;
		return hetTally[vector];
	}
	
	public int getOutlierTally(int vector)
	{
		if (vector < 0) return 0;
		if (vector > otherTally.length) return 0;
		return otherTally[vector];
	}

	public void printReport()
	{
		System.out.println("Samples (In Order)");
		System.out.println("------------------");
		for (String s : sampleList)
		{
			System.out.println(s);
		}
		System.out.println();
		
		System.out.println("General Tally");
		System.out.println("-------------");
		for (int i = 0; i < tally.length; i++)
		{
			String binString = Integer.toBinaryString(i);
			while (binString.length() < sampleCount)
			{
				binString = "0" + binString;
			}
			System.out.println(binString + "\t" + tally[i]);
		}
		System.out.println();
		
		System.out.println("Homozygous Alt Tally");
		System.out.println("--------------------");
		for (int i = 0; i < homTally.length; i++)
		{
			String binString = Integer.toBinaryString(i);
			while (binString.length() < sampleCount)
			{
				binString = "0" + binString;
			}
			System.out.println(binString + "\t" + homTally[i]);
		}
		System.out.println();
		
		System.out.println("Heterozygous Alt Tally");
		System.out.println("----------------------");
		for (int i = 0; i < hetTally.length; i++)
		{
			String binString = Integer.toBinaryString(i);
			while (binString.length() < sampleCount)
			{
				binString = "0" + binString;
			}
			System.out.println(binString + "\t" + hetTally[i]);
		}
		System.out.println();
		
		System.out.println("Outlier Tally");
		System.out.println("----------------------");
		for (int i = 0; i < otherTally.length; i++)
		{
			String binString = Integer.toBinaryString(i);
			while (binString.length() < sampleCount)
			{
				binString = "0" + binString;
			}
			System.out.println(binString + "\t" + otherTally[i]);
		}
		System.out.println();
		
		
	}
	
}
