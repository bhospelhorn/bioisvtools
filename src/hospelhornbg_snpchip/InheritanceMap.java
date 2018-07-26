package hospelhornbg_snpchip;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class InheritanceMap {
	
	/* ----- Instance Variables ----- */
	
	private Map <Integer, Allele> imap;
	
	/* ----- Inner Classes ----- */
	
	public static class Allele
	{
		private int allele;
		
		private int count;
		private int p1Count;
		private int p2Count;
		
		private int flagP1; //Notes for number to flag later
		private int flagP2;
		private int flagError; //So as to not count as "unphased"
		
		private Set<SNPAllele> aobjects;
		
		public Allele(int a)
		{
			allele = a;
			count = 0;
			p1Count = -1;
			p2Count = -2;
			flagP1 = 0;
			flagP2 = 0;
			flagError = 0;
			aobjects = new HashSet<SNPAllele>();
		}
		
		public synchronized void incrementCount()
		{
			count++;
		}
		
		public synchronized void tallyParent1(SNPGeno p1Geno)
		{
			p1Count = p1Geno.countAllele(allele);
		}
		
		public synchronized void tallyParent2(SNPGeno p2Geno)
		{
			p2Count = p2Geno.countAllele(allele);
		}		
		
		public synchronized int getAlleleCount()
		{
			return count;
		}
		
		public synchronized int getUnflaggedCount()
		{
			return count - flagP1 - flagP2 - flagError;
		}
		
		public synchronized int getParent1Count()
		{
			return p1Count;
		}
		
		public synchronized int getParent2Count()
		{
			return p2Count;
		}
		
		public synchronized int getAllele()
		{
			return allele;
		}
	
		public synchronized int getFlag1Count()
		{
			return flagP1;
		}
		
		public synchronized void setFlag1Count(int i)
		{
			flagP1 = i;
		}
		
		public synchronized int getFlag2Count()
		{
			return flagP2;
		}
		
		public synchronized void setFlag2Count(int i)
		{
			flagP2 = i;
		}
	
		public synchronized int getUnphasedCount()
		{
			return count - flagP1 - flagP2 - flagError;
		}

		public synchronized int getErrorFlagCount()
		{
			return flagError;
		}

		public synchronized void setErrorFlagCount(int i)
		{
			flagError = i;
		}
		
		public synchronized boolean linkAlleleObject(SNPAllele a)
		{
			if (a.getAllele() != allele) return false;
			aobjects.add(a);
			return true;
		}
		
		public synchronized boolean flagLinkedAlleles()
		{
			if (aobjects.size() < (flagP1 + flagP2)) return false;
			int i = 0;
			int j = 0;
			for (SNPAllele a : aobjects)
			{
				a.clearPhasingFlags();
				if (i < flagP1){
					a.flagParent1(true);
					i++;
					continue;
				}
				if (j < flagP2)
				{
					a.flagParent2(true);
					j++;
				}
			}
			return true;
		}
		
	}

	/* ----- Construction ----- */
	
	public InheritanceMap()
	{
		imap = new HashMap<Integer, Allele>();
	}
	
	public InheritanceMap(SNPGeno geno, SNPGeno parent1, SNPGeno parent2)
	{
		imap = new HashMap<Integer, Allele>();
		List<SNPAllele> alist = geno.getAlleles();
		for (SNPAllele a : alist)
		{
			Allele as = imap.get(a.getAllele());
			if (as == null)
			{
				as = new Allele(a.getAllele());
				as.incrementCount();
				as.linkAlleleObject(a);
				if(parent1 != null) as.setFlag1Count(parent1.countAllele(a.getAllele()));
				if(parent2 != null) as.setFlag2Count(parent2.countAllele(a.getAllele()));
			}
			else{
				as.incrementCount();
				as.linkAlleleObject(a);
			}
		}
	}
	
	/* ----- Getters ----- */
	
	public synchronized Allele getAllele(int a)
	{
		return imap.get(a);
	}
	
	public synchronized Set<Integer> getKeySet()
	{
		return imap.keySet();
	}
	
	/* ----- Setters ----- */
	
	public synchronized void putAllele(int a, Allele as)
	{
		imap.put(a, as);
	}
	
	/* ----- Utility ----- */
	
	public static boolean possible_idup(List<Integer> possibleIDUP, List<Allele> pflagged, List<Allele> unflagged, boolean p1)
	{
		if (possibleIDUP == null) return false;
		//Between pflagged and unphased, see if there are at least two of any alleles in the dup list
		for (Integer i : possibleIDUP)
		{
			int act = 0;
			for (Allele a : pflagged)
			{
				if (a.getAllele() == i)
				{
					if (p1) act += a.getFlag1Count();
					else a.getFlag2Count();	
				}
			}
			for (Allele a : unflagged)
			{
				if (a.getAllele() == i)
				{
					act += a.getUnflaggedCount();	
				}
			}
			if (act > 1) return true;
		}
		
		return false;
	}
	
	public static boolean possible_idup(List<Integer> possibleIDUP, List<Allele> unflagged)
	{
		if (possibleIDUP == null) return false;
		//Between pflagged and unphased, see if there are at least two of any alleles in the dup list
		for (Integer i : possibleIDUP)
		{
			int act = 0;
			for (Allele a : unflagged)
			{
				if (a.getAllele() == i)
				{
					act += a.getUnflaggedCount();	
				}
			}
			if (act > 1) return true;
		}
		
		return false;
	}
	
	public static boolean possible_upid(List<Integer> possibleUPID, List<Allele> pflagged, List<Allele> unflagged, boolean p1)
	{
		if (possibleUPID == null) return false;
		for (Integer i : possibleUPID)
		{
			int act = 0;
			for (Allele a : pflagged)
			{
				if (a.getAllele() == i)
				{
					if (p1) act += a.getFlag1Count();
					else a.getFlag2Count();	
				}
			}
			for (Allele a : unflagged)
			{
				if (a.getAllele() == i)
				{
					act += a.getUnflaggedCount();	
				}
			}
			if (act > 1) return true;
		}
		
		return false;
	}
	
	public static boolean possible_upid(List<Integer> possibleUPID, List<Allele> unflagged)
	{
		if (possibleUPID == null) return false;
		for (Integer i : possibleUPID)
		{
			int act = 0;
			for (Allele a : unflagged)
			{
				if (a.getAllele() == i)
				{
					act += a.getUnflaggedCount();	
				}
			}
			if (act > 1) return true;
		}
		
		return false;
	}
	
	public static boolean possible_uphd(List<Integer> possibleUPHD, List<Allele> pflagged, List<Allele> unflagged, boolean p1)
	{
		if (possibleUPHD == null) return false;
		int act = 0;
		for (Integer i : possibleUPHD)
		{
			//Need to find at least two alleles of any combination.
			for (Allele a : pflagged)
			{
				if (a.getAllele() == i)
				{
					if (p1) act += a.getFlag1Count();
					else a.getFlag2Count();	
				}
			}
			for (Allele a : unflagged)
			{
				if (a.getAllele() == i)
				{
					act += a.getUnflaggedCount();	
				}
			}
		}
		if (act > 1) return true;
		
		return false;
	}
	
	public static boolean possible_uphd(List<Integer> possibleUPHD, List<Allele> unflagged)
	{
		if (possibleUPHD == null) return false;
		int act = 0;
		for (Integer i : possibleUPHD)
		{
			//Need to find at least two alleles of any combination.
			for (Allele a : unflagged)
			{
				if (a.getAllele() == i)
				{
					act += a.getUnflaggedCount();	
				}
			}
		}
		if (act > 1) return true;
		
		return false;
	}
	
	/* ----- Queries ----- */
	
	public boolean anyPhasedP1()
	{
		Set<Integer> keys = imap.keySet();
		for (Integer i : keys)
		{
			Allele as = imap.get(i);
			if (as.getFlag1Count() > 0) return true;
		}
		return false;
	}
	
	public boolean anyPhasedP2()
	{
		Set<Integer> keys = imap.keySet();
		for (Integer i : keys)
		{
			Allele as = imap.get(i);
			if (as.getFlag2Count() > 0) return true;
		}
		return false;
	}
	
	public int countPhasedP1()
	{
		int c = 0;
		Set<Integer> keys = imap.keySet();
		for (Integer i : keys)
		{
			Allele as = imap.get(i);
			c += as.getFlag1Count();
		}
		return c;
	}
	
	public int countPhasedP2()
	{
		int c = 0;
		Set<Integer> keys = imap.keySet();
		for (Integer i : keys)
		{
			Allele as = imap.get(i);
			c += as.getFlag2Count();
		}
		return c;
	}
	
	public int countUnphased()
	{
		int c = 0;
		Set<Integer> keys = imap.keySet();
		for (Integer i : keys)
		{
			Allele as = imap.get(i);
			c += as.getUnflaggedCount();
		}
		return c;
	}

	public int countErrorAlleles()
	{
		int c = 0;
		Set<Integer> keys = imap.keySet();
		for (Integer i : keys)
		{
			Allele as = imap.get(i);
			c += as.getErrorFlagCount();
		}
		return c;
	}

	public void flagAlleleAsError(int a)
	{
		Allele as = imap.get(a);
		if (as != null)
		{
			as.setErrorFlagCount(as.getAlleleCount());
		}
	}
	
	public List<Allele> getUniquePhasedP1()
	{
		Set<Integer> keys = imap.keySet();
		List<Allele> alist = new LinkedList<Allele>(); 
		for (Integer i : keys)
		{
			Allele as = imap.get(i);
			if (as.getFlag1Count() > 0) alist.add(as);
		}
		if (alist.size() > 0) return alist;
		return null;
	}
	
	public List<Allele> getUniquePhasedP2()
	{
		Set<Integer> keys = imap.keySet();
		List<Allele> alist = new LinkedList<Allele>(); 
		for (Integer i : keys)
		{
			Allele as = imap.get(i);
			if (as.getFlag2Count() > 0) alist.add(as);
		}
		if (alist.size() > 0) return alist;
		return null;
	}
	
	public List<Allele> getUniqueUnphased()
	{
		Set<Integer> keys = imap.keySet();
		List<Allele> alist = new LinkedList<Allele>(); 
		for (Integer i : keys)
		{
			Allele as = imap.get(i);
			if (as.getUnflaggedCount() > 0) alist.add(as);
		}
		if (alist.size() > 0) return alist;
		return null;
	}

	public void flagLinkedAlleles()
	{
		Set<Integer> keys = imap.keySet();
		for (Integer i : keys)
		{
			Allele as = imap.get(i);
			as.flagLinkedAlleles();
		}
	}

}
