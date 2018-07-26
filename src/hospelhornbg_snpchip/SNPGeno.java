package hospelhornbg_snpchip;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import hospelhornbg_bioinformatics.Sex;
import hospelhornbg_snpchip.InheritanceMap.Allele;
import waffleoRai_Utils.TallyMap;

//TODO: 1. Add mosaicism functionality

/*
 * UPDATES
 */

/**
 * A class that for holding allele objects and special phasing flags.
 * Also contains static methods for phasing single positions.
 * <br>Parent class does not handle mosaicism.
 * <br>Needs re-organizing, but will have to do for now.
 * @author Blythe Hospelhorn
 * @version 1.2.0
 * @since June 21, 2018
 *
 */
public class SNPGeno {
	
	/* ----- Constants ----- */
	
	/**
	 * Default value for the phasepath variable when phasing has not been attempted.
	 */
	public static final int PHASEPATH_NONE = -1;
	
	/* ----- Instance Variables ----- */
	
	private Set<SNPAllele> alleles;
	
	private int copyNumber; //The reason this is included is to distinguish CNV3s from mosaicism
	
	// Inheritance unresolved flags
	private boolean phasingAttempted;
	private boolean iu_allAllelesUnknown;
	private boolean iu_anyAllelesUnknown;
	private boolean iu_unexpectedSegregation;
	private boolean iu_unknownGenoAssumed;
	private boolean iu_denovoCandidate;
	private boolean iu_UPID_Candidate_p1;
	private boolean iu_UPID_Candidate_p2;
	private boolean iu_UPHD_Candidate_p1;
	private boolean iu_UPHD_Candidate_p2;
	private boolean iu_inheritedDEL_candidate_p1;
	private boolean iu_inheritedDEL_candidate_p2;
	private boolean iu_inheritedDUP_candidate_p1;
	private boolean iu_inheritedDUP_candidate_p2;
	private boolean iu_mosaic_candidate;
	private boolean iu_unlikelyCall;
	private boolean iu_PC_err;
	private boolean iu_PPC_err;
	
	private int phasePath; //Enum to track phasing path
	
	/* ----- Construction ----- */

	/**
	 * Construct an empty SNPGeno object. The copy number is set to 2 by default.
	 */
	public SNPGeno()
	{
		alleles = new HashSet<SNPAllele>();
		copyNumber = 2;
		this.clearIUFlags();
	}
	
	/* ----- Inner Classes ----- */
	
	/* ----- Getters ----- */
	
	public int getCopyNumber()
	{
		return copyNumber;
	}
	
	public int getAlleleCount()
	{
		return alleles.size();
	}
	
	public List<SNPAllele> getAlleles()
	{
		List<SNPAllele> alist = new ArrayList<SNPAllele>(alleles.size());
		alist.addAll(alleles);
		return alist;
	}
	
	public boolean hasAllele(SNPAllele a)
	{
		for(SNPAllele al : alleles)
		{
			if (a.getAllele() == al.getAllele()) return true;
		}
		return false;
	}
	
	public int countAllele(SNPAllele a)
	{
		int c = 0;
		for(SNPAllele al : alleles)
		{
			if (a.getAllele() == al.getAllele()) c++;
		}
		return c;
	}
	
	public int countAllele(int a)
	{
		int c = 0;
		for(SNPAllele al : alleles)
		{
			if (a == al.getAllele()) c++;
		}
		return c;
	}
	
	public int countUniqueAlleles()
	{
		Set<Integer> aset = new HashSet<Integer>();
		for(SNPAllele a : alleles)
		{
			aset.add(a.getAllele());
		}
		return aset.size();
	}
	
	public Set<Integer> getUniqueAlleles()
	{
		Set<Integer> aset = new HashSet<Integer>();
		for(SNPAllele a : alleles)
		{
			aset.add(a.getAllele());
		}
		return aset;
	}
	
	public TallyMap getAlleleCountBreakdown()
	{
		TallyMap tm = new TallyMap();
		for (SNPAllele a : alleles) tm.increment(a.getAllele());
		return tm;
	}
	
	public InheritanceMap buildInheritanceMap(SNPGeno p1Geno, SNPGeno p2Geno)
	{
		InheritanceMap imap = new InheritanceMap(this, p1Geno, p2Geno);
		return imap;
	}
	
		// -- Inheritance Unresolved Flags
	
	public boolean phasingAttempted()
	{
		return this.phasingAttempted;
	}
	
	public boolean allAllelesUnknown()
	{
		return iu_allAllelesUnknown;
	}
	
	public boolean anyAllelesUnknown()
	{
		return iu_anyAllelesUnknown;
	}
	
	public boolean unexpectedSegregation()
	{
		return iu_unexpectedSegregation;
	}
	
	public boolean unknownGenotypeAssumed()
	{
		return iu_unknownGenoAssumed;
	}
	
	public boolean denovoCandidate()
	{
		return iu_denovoCandidate;
	}
	
	public boolean isUniparentalIsodisomyCandidate_p1()
	{
		return iu_UPID_Candidate_p1;
	}
	
	public boolean isUniparentalIsodisomyCandidate_p2()
	{
		return iu_UPID_Candidate_p2;
	}
	
	public boolean isUniparentalHeterodisomyCandidate_p1()
	{
		return iu_UPHD_Candidate_p1;
	}
	
	public boolean isUniparentalHeterodisomyCandidate_p2()
	{
		return iu_UPHD_Candidate_p2;
	}
	
	public boolean isInheritedDeletionCandidate_p1()
	{
		return iu_inheritedDEL_candidate_p1;
	}
	
	public boolean isInheritedDeletionCandidate_p2()
	{
		return iu_inheritedDEL_candidate_p2;
	}
	
	public boolean isInheritedDuplicationCandidate_p1()
	{
		return iu_inheritedDUP_candidate_p1;
	}
	
	public boolean isInheritedDuplicationCandidate_p2()
	{
		return iu_inheritedDUP_candidate_p2;
	}
	
	public boolean isMosaicCandidate()
	{
		return iu_mosaic_candidate;
	}
	
	public boolean unlikelyPhasingCall()
	{
		return iu_unlikelyCall;
	}
	
	public boolean PC_error()
	{
		return iu_PC_err;
	}
	
	public boolean PPC_error()
	{
		return iu_PPC_err;
	}
	
	public int getPhasePath()
	{
		return this.phasePath;
	}
	
	/* ----- Setters ----- */
	
	private void clearIUFlags()
	{
		phasingAttempted = false;
		iu_allAllelesUnknown = false;
		iu_anyAllelesUnknown = false;
		iu_unexpectedSegregation = false;
		iu_unknownGenoAssumed = false;
		iu_denovoCandidate = false;
		//iu_UPD_Candidate = false;
		iu_UPID_Candidate_p1 = false;
		iu_UPID_Candidate_p2 = false;
		iu_UPHD_Candidate_p1 = false;
		iu_UPHD_Candidate_p2 = false;
		iu_inheritedDEL_candidate_p1 = false;
		iu_inheritedDUP_candidate_p1 = false;
		iu_inheritedDEL_candidate_p2 = false;
		iu_inheritedDUP_candidate_p2 = false;
		iu_unlikelyCall = false;
		iu_PC_err = false;
		iu_PPC_err = false;
		iu_mosaic_candidate = false;
		phasePath = PHASEPATH_NONE;
	}
	
	public void clearPhasingMarks()
	{
		for (SNPAllele a : alleles)
		{
			a.clearPhasingFlags();
		}
		clearIUFlags();
	}

	/* ----- Phasing ----- */
	
		// -- Helper Methods
	
	public boolean anyPhasedP1()
	{
		for (SNPAllele a : alleles)
		{
			if (a.phased_parent1()) return true;
		}
		return false;
	}
	
	public boolean anyPhasedP2()
	{
		for (SNPAllele a : alleles)
		{
			if (a.phased_parent2()) return true;
		}
		return false;
	}
	
	public int countPhasedP1()
	{
		int tot = 0;
		for (SNPAllele a : alleles)
		{
			if (a.phased_parent1()) tot++;
		}
		return tot;
	}
	
	public int countPhasedP2()
	{
		int tot = 0;
		for (SNPAllele a : alleles)
		{
			if (a.phased_parent2()) tot++;
		}
		return tot;
	}
	
	public int countUnphased()
	{
		int tot = 0;
		for (SNPAllele a : alleles)
		{
			if (!a.isPhased()) tot++;
		}
		return tot;
	}
	
	public boolean isDiploid()
	{
		return (copyNumber == 2);
	}
	
	public boolean isHomozygous()
	{
		if (alleles.isEmpty()){
			copyNumber = 0;
			return true; //Technically!
		}
		Set<Integer> aset = new HashSet<Integer>();
		for (SNPAllele a : alleles) aset.add(a.getAllele());
		return (aset.size() == 1);
	}
	
	public static Sample reconstructMissingParent(Collection<Sample> children, Sample knownParent)
	{
		//Need at least 3 children
		return null;
	}
	
	public static boolean canInherit_PC(SNPGeno parentGeno, SNPGeno childGeno)
	{
		//Returns false if there is a PC error
		//Considers copy number
		if (parentGeno == null) return false;
		if (childGeno == null) return false;
		
		List<SNPAllele> pall = parentGeno.getAlleles();
		List<SNPAllele> call = childGeno.getAlleles();
		
		//There must be at least one allele both share
		for (SNPAllele ca: call)
		{
			int cai = ca.getAllele();
			for (SNPAllele pa: pall)
			{
				int pai = pa.getAllele();
				if (cai == pai) return true;
			}
		}
		
		return false;
	}
	
	public static boolean canInherit_PPC(SNPGeno p1Geno, SNPGeno p2Geno, SNPGeno childGeno)
	{
		if(!canInherit_PC(p1Geno, childGeno)) return false;
		if(!canInherit_PC(p2Geno, childGeno)) return false;
		
		return check_for_uninheritableAlleles_PPC(p1Geno, p2Geno, childGeno);
	}
	
	public static boolean check_for_uninheritableAlleles_PPC(SNPGeno p1Geno, SNPGeno p2Geno, SNPGeno childGeno)
	{
		List<SNPAllele> call = childGeno.getAlleles();
		for (SNPAllele ca: call)
		{
			int p1c = p1Geno.countAllele(ca);
			int p2c = p2Geno.countAllele(ca);
			if (p1c < 1) return false;
			if (p2c < 1) return false;
		}
		return true;
	}
	
	public static List<Integer> possibleInheritedDuplications(SNPGeno parentGeno, SNPGeno childGeno)
	{
		if (childGeno.getCopyNumber() < 3) return null;
		if (parentGeno.getCopyNumber() < 3) return null;
		TallyMap cAlleles = childGeno.getAlleleCountBreakdown();
		TallyMap pAlleles = parentGeno.getAlleleCountBreakdown();
		
		List<Integer> alist = new LinkedList<Integer>();
		//If any have >= 2 at the same value, then it passes
		List<Integer> ckeys = cAlleles.getAllValues();
		for (Integer a : ckeys)
		{
			int ccount = cAlleles.getCount(a);
			int pcount = pAlleles.getCount(a);
			if ((ccount > 1) && (pcount > 1)) alist.add(a);
		}
		
		if(alist.size() > 0) return alist;
		return null;
	}
	
	public static List<Integer> UPID_possible(SNPGeno parentGeno, SNPGeno childGeno)
	{
		//Must be two copies of same allele
		TallyMap cAlleles = childGeno.getAlleleCountBreakdown();
		
		List<Integer> alist = new LinkedList<Integer>();
		List<SNPAllele> pall = parentGeno.getAlleles();
		for (SNPAllele a : pall)
		{
			int ccount = cAlleles.getCount(a.getAllele());
			if (ccount > 1) alist.add(a.getAllele());
		}
		
		if(alist.size() > 0) return alist;
		return null;
	}
	
	public static boolean UPHD_possible(SNPGeno parentGeno, SNPGeno childGeno)
	{
		//Must have at least one copy of at least two parental alleles
		Set<Integer> caset = childGeno.getUniqueAlleles();
		Set<Integer> paset = parentGeno.getUniqueAlleles();
		
		int matchcount = 0;
		for (Integer i : caset)
		{
			if (paset.contains(i)) matchcount++;
			if (matchcount >= 2) return true;
		}
		
		return false;
	}
	
	public static List<Integer> canInherit_PC(InheritanceMap childimap, boolean p1)
	{
		Set<Integer> call = childimap.getKeySet();
		List<Integer> PCpass = new LinkedList<Integer>();
		
		for (Integer i : call)
		{
			Allele as = childimap.getAllele(i);
			if(p1 && (as.getParent1Count() > 0)) PCpass.add(i);
			if(!p1 && (as.getParent2Count() > 0)) PCpass.add(i);
		}
		
		if (PCpass.size() == 0) return null;
		return PCpass;
	}
	
	public static List<Integer> get_uninheritableAlleles_PPC(InheritanceMap childimap)
	{
		Set<Integer> call = childimap.getKeySet();
		List<Integer> PPCfail = new LinkedList<Integer>();
		
		for (Integer i : call)
		{
			Allele as = childimap.getAllele(i);
			int p1c = as.getParent1Count();
			int p2c = as.getParent2Count();
			if ((p1c < 1) && (p2c < 1)) PPCfail.add(i);
		}
		
		if(PPCfail.size() == 0) return null;
		return PPCfail;
	}
	
	public static List<Integer> get_uninheritableAlleles_PC(InheritanceMap childimap)
	{
		Set<Integer> call = childimap.getKeySet();
		List<Integer> PCfail = new LinkedList<Integer>();
		
		for (Integer i : call)
		{
			Allele as = childimap.getAllele(i);
			int p1c = as.getParent1Count();
			if (p1c < 1) PCfail.add(i);
		}
		
		if(PCfail.size() == 0) return null;
		return PCfail;
	}
	
	public static List<Integer> possibleInheritedDuplications(InheritanceMap childimap, boolean p1, int ccnv, int pcnv)
	{
		if (ccnv < 3) return null;
		if (pcnv < 3) return null;
		
		Set<Integer> call = childimap.getKeySet();
		List<Integer> dupCand = new LinkedList<Integer>();
		
		for (Integer i : call)
		{
			Allele as = childimap.getAllele(i);
			if (as.getAlleleCount() > 1)
			{
				if (p1 && (as.getParent1Count() > 1)) dupCand.add(i);
				if (!p1 && (as.getParent2Count() > 1)) dupCand.add(i);
			}
		}
		
		if (dupCand.size() > 0) return dupCand;
		return null;
	}
	
	public static List<Integer> UPID_possible(InheritanceMap childimap, boolean p1)
	{
		Set<Integer> call = childimap.getKeySet();
		List<Integer> updCand = new LinkedList<Integer>();
		
		for (Integer i : call)
		{
			Allele as = childimap.getAllele(i);
			if (as.getAlleleCount() > 1)
			{
				if (p1 && (as.getParent1Count() > 0)) updCand.add(i);
				if (!p1 && (as.getParent2Count() > 0)) updCand.add(i);
			}
		}
		
		if (updCand.size() > 0) return updCand;
		return null;
	}
	
	public static List<Integer> UPHD_possible(InheritanceMap childimap, boolean p1)
	{
		Set<Integer> call = childimap.getKeySet();
		List<Integer> updCand = new LinkedList<Integer>();
		
		for (Integer i : call)
		{
			Allele as = childimap.getAllele(i);
			//If parent has, add to list
			if (p1 && (as.getParent1Count() > 0)) updCand.add(i);
			if (!p1 && (as.getParent2Count() > 0)) updCand.add(i);
		}
		
		//Assess list size
		if (updCand.size() < 2) return null;
		
		return updCand;
	}
	
	public void flagAllOneParent(boolean p1)
	{
		for (SNPAllele a : alleles)
		{
			if (p1) a.flagParent1(true);
			else a.flagParent2(true);
		}
	}
	
	public void flagRemainingOneParent(boolean p1)
	{
		for (SNPAllele a : alleles)
		{
			if (a.isPhased()) continue;
			if (p1) a.flagParent1(true);
			else a.flagParent2(true);
		}
	}
	
	public void flagAllOneParent(List<Integer> pall, boolean p1)
	{
		if (pall == null) return;
		if (pall.isEmpty()) return;
		for (SNPAllele a : alleles)
		{
			Integer ai = a.getAllele();
			if (pall.contains(ai))
			{
				if(p1) a.flagParent1(true);
				else a.flagParent2(true);
			}
		}
	}
	
	public void flagOneAllele(int allele, boolean p1)
	{
		for (SNPAllele a : alleles)
		{
			if (!a.isPhased() && (a.getAllele() == allele))
			{
				if (p1) a.flagParent1(true);
				if (!p1) a.flagParent2(true);
				return;
			}
		}
	}
	
	public void flagOneEachParent(int allele)
	{
		flagOneAllele(allele, true);
		flagOneAllele(allele, false);
	}
	
	public void flagAllNotOtherParent(int allele, boolean p1)
	{
		for (SNPAllele a : alleles)
		{
			if (a.getAllele() != allele)
			{
				if(p1) a.flagParent2(true);
				else a.flagParent1(true);
			}
		}
	}
	
	public void flagAllNotOtherParent(List<Integer> pall, boolean p1)
	{
		for (SNPAllele a : alleles)
		{
			Integer ai = a.getAllele();
			if (!pall.contains(ai))
			{
				if(p1) a.flagParent2(true);
				else a.flagParent1(true);
			}
		}
	}
	
		// -- Phasing Methods (Specific)
	
	private static void phase_PC_CNV0_auto(SNPGeno parentGeno, SNPGeno childGeno)
	{
		//Child is CNV0, so no alleles to flag.
		//However, can try to flag genotype based on whether deletions could be inherited
		
		childGeno.phasingAttempted = true;
		
		//Does known parent have a deletion?
		if (parentGeno.getCopyNumber() < 2)
		{
			//Possibilities (Technically any):
				//P1[Del] / DN[Del] *
				//P1[Del] / P1[Del] *
				//P1[Del] / P2[Del] *
				//P2[Del] / P2[Del]
				//P2[Del] / DN[Del]
				//DN[Del] / DN[Del]
			
			//Most likely P1 del was inherited
			//Others are POSSIBLE, but we'll flag as if it was P1 inherited.
			childGeno.iu_denovoCandidate = true;
			childGeno.iu_inheritedDEL_candidate_p1 = true;
			childGeno.iu_UPID_Candidate_p1 = true;
			childGeno.phasePath = 0;
		}
		else
		{
			//Possibilities (Technically any):
				//P2[Del] / P2[Del]
				//P2[Del] / DN[Del]
				//DN[Del] / DN[Del]
			childGeno.iu_unknownGenoAssumed = true;
			childGeno.iu_denovoCandidate = true;
			childGeno.iu_inheritedDEL_candidate_p2 = true;
			childGeno.iu_UPID_Candidate_p2 = true;
			childGeno.phasePath = 1;
		}
	}
	
	private static void phase_PPC_CNV0_auto(SNPGeno p1Geno, SNPGeno p2Geno, SNPGeno childGeno)
	{
		//Child is CNV0, so no alleles to flag.
		//However, can try to flag genotype based on whether deletions could be inherited
		boolean p1Del = p1Geno.getCopyNumber() < 2;		
		boolean p2Del = p2Geno.getCopyNumber() < 2;
		
		childGeno.phasingAttempted = true;
		
		if (p1Del)
		{
			if (p2Del)
			{
				//Most likely:
				//P1[Del] / P2[Del] *
				//P1[Del] / P1[Del]
				//P2[Del] / P2[Del]
				
				//Denovo's always possible of course.
				childGeno.iu_inheritedDEL_candidate_p1 = true;
				childGeno.iu_inheritedDEL_candidate_p2 = true;
				childGeno.phasePath = 2;
			}
			else
			{
				//Most likely:
				//P1[Del] / DN[Del] *
				//P1[Del] / P1[Del] *
				//DN[Del] / DN[Del]
				
				childGeno.iu_inheritedDEL_candidate_p1 = true;
				childGeno.iu_unexpectedSegregation = true;
				childGeno.iu_denovoCandidate = true;
				childGeno.iu_UPID_Candidate_p1 = true;
				childGeno.phasePath = 3;
			}
		}
		else
		{
			if (p2Del)
			{
				//Most likely:
				//P2[Del] / DN[Del] *
				//P2[Del] / P2[Del] *
				//DN[Del] / DN[Del]
				
				childGeno.iu_inheritedDEL_candidate_p2 = true;
				childGeno.iu_unexpectedSegregation = true;
				childGeno.iu_denovoCandidate = true;
				childGeno.iu_UPID_Candidate_p2 = true;
				childGeno.phasePath = 4;
			}
			else
			{
				//Must be denovo
				//DN[Del] / DN[Del]
				
				childGeno.iu_denovoCandidate = true;
				childGeno.iu_unexpectedSegregation = true;
				childGeno.phasePath = 5;
			}
		}
		
		
	}
	
	private static void phase_PC_CNV1_auto(SNPGeno parentGeno, SNPGeno childGeno)
	{
		//See if known parent has a deletion
		//If not, see if child could have inherited an allele from known parent
		
		childGeno.phasingAttempted = true;
		int pcnv = parentGeno.getCopyNumber();
		SNPAllele call = childGeno.getAlleles().get(0);
		
		if (pcnv < 2)
		{
			if (pcnv == 0)
			{
				//Possibilities:
					//P1[del]/P2 *
				call.flagParent2(true);
				childGeno.iu_inheritedDEL_candidate_p1 = true;
				childGeno.iu_unknownGenoAssumed = true;
				childGeno.phasePath = 6;
			}
			else
			{
				if (canInherit_PC(parentGeno, childGeno))
				{
					//Possibilities:
						//P1[del]/P2 *
						//P1[del]/P1 *
						//P2[del]/P2
						//DN[del]/P1 *
						//DN[del]/P2
					//Do not phase known allele
					childGeno.iu_allAllelesUnknown = true;
					childGeno.iu_inheritedDEL_candidate_p1 = true;
					childGeno.iu_UPHD_Candidate_p1 = true;
					childGeno.phasePath = 7;
				}
				else
				{
					//Possibilities:
						//P1[del]/P2 *Most likely
						//P2[del]/P2
						//DN[del]/P2
					call.flagParent2(true);
					childGeno.iu_unknownGenoAssumed = true;
					childGeno.iu_inheritedDEL_candidate_p1 = true;
					childGeno.phasePath = 8;
				}
			}				
		}
		else
		{
			//Known parent is at least diploid
			if (canInherit_PC(parentGeno, childGeno))
			{
				//Possibilities:
					//P1/P2[del] *Most likely
					//P1/DN[del] *
					//DN[del]/P2 *
					//P2/P2[del]
				childGeno.iu_allAllelesUnknown = true;
				childGeno.iu_unknownGenoAssumed = true;
				childGeno.iu_denovoCandidate = true;
				childGeno.iu_inheritedDEL_candidate_p2 = true;
				childGeno.iu_UPHD_Candidate_p2 = true;
				childGeno.phasePath = 9;
			}
			else
			{
				//Possibilities:
					//DN[del]/P2
					//P2/P2[del]
				call.flagParent2(true);
				childGeno.iu_unexpectedSegregation = true;
				childGeno.iu_unknownGenoAssumed = true;
				childGeno.iu_denovoCandidate = true;
				childGeno.iu_inheritedDEL_candidate_p2 = true;
				childGeno.iu_UPHD_Candidate_p2 = true;
				childGeno.phasePath = 10;
			}
		}
		
	}

	private static void phase_PPC_CNV1_auto(SNPGeno p1Geno, SNPGeno p2Geno, SNPGeno childGeno)
	{
		childGeno.phasingAttempted = true;
		int pcnv1 = p1Geno.getCopyNumber();
		int pcnv2 = p2Geno.getCopyNumber();
		SNPAllele ca = childGeno.getAlleles().get(0);
		
		boolean p1_hasDel = (pcnv1 < 2);
		boolean p2_hasDel = (pcnv2 < 2);
		
		boolean p1i = canInherit_PC(p1Geno, childGeno);
		boolean p2i = canInherit_PC(p2Geno, childGeno);
		
		if (p1_hasDel)
		{
			if(p2_hasDel)
			{
				//Both parents have a deletion.
				if(p1i)
				{
					if (p2i)
					{
						//Remaining allele could have been inherited from either parent.
						//Do not flag the allele. Flag both inherited del
						childGeno.iu_allAllelesUnknown = true;
						childGeno.iu_inheritedDEL_candidate_p1 = true;
						childGeno.iu_inheritedDEL_candidate_p2 = true;
						childGeno.phasePath = 11;
						return;
					}
					else
					{
						//Only could have inherited from parent 1
						ca.flagParent1(true);
						childGeno.iu_inheritedDEL_candidate_p2 = true;
						childGeno.phasePath = 12;
						return;
					}
				}
				else
				{
					//Can't have inherited from p1
					if (p2i)
					{
						//Only could have inherited from parent 2
						ca.flagParent2(true);
						childGeno.iu_inheritedDEL_candidate_p1 = true;
						childGeno.phasePath = 13;
						return;
					}
					else
					{
						//Can't have inherited from either...
						childGeno.iu_PPC_err = true;
						childGeno.iu_allAllelesUnknown = true;
						childGeno.iu_denovoCandidate = true;
						childGeno.iu_unlikelyCall = true;
						childGeno.iu_unexpectedSegregation = true;
						childGeno.phasePath = 14;
						return;
					}
				}
			}
			else
			{
				//Only p1 has a deletion
				if(p1i)
				{
					if (p2i)
					{
						//Remaining allele could have been inherited from either parent
						//Assumed it's P2 and deletion is P1
						ca.flagParent2(true);
						childGeno.iu_inheritedDEL_candidate_p1 = true;
						childGeno.phasePath = 15;
						return;
					}
					else
					{
						//Remaining allele could have only been inherited from P1
						ca.flagParent1(true);
						childGeno.iu_denovoCandidate = true;
						childGeno.iu_inheritedDEL_candidate_p1 = true;
						childGeno.iu_unexpectedSegregation = true;
						childGeno.iu_UPHD_Candidate_p1 = true;
						childGeno.phasePath = 16;
						return;
					}
				}
				else
				{
					if (p2i)
					{
						//Can only have inherited from p2
						ca.flagParent2(true);
						childGeno.iu_inheritedDEL_candidate_p1 = true;
						childGeno.phasePath = 17;
						return;
					}
					else
					{
						//Can't have inherited from either
						childGeno.iu_PPC_err = true;
						childGeno.iu_allAllelesUnknown = true;
						childGeno.iu_denovoCandidate = true;
						childGeno.iu_unlikelyCall = true;
						childGeno.iu_unexpectedSegregation = true;
						childGeno.phasePath = 18;
						return;
					}
				}
			}
		}
		else
		{
			if(p2_hasDel)
			{
				//Only p2 has a deletion
				if(p1i)
				{
					if (p2i)
					{
						//Could have inherited from either
						//Assumed deletion inherited from p2
						ca.flagParent1(true);
						childGeno.iu_inheritedDEL_candidate_p2 = true;
						childGeno.phasePath = 19;
						return;
					}
					else
					{
						//Could have only inherited from p1
						//Assumed deletion inherited from p2
						ca.flagParent1(true);
						childGeno.iu_inheritedDEL_candidate_p2 = true;
						childGeno.phasePath = 20;
						return;
					}
				}
				else
				{
					//Can't have inherited from p1
					if (p2i)
					{
						//Could have only inherited from p2
						ca.flagParent2(true);
						childGeno.iu_unexpectedSegregation = true;
						childGeno.iu_denovoCandidate = true;
						childGeno.iu_inheritedDEL_candidate_p2 = true;
						childGeno.iu_UPHD_Candidate_p2 = true;
						childGeno.phasePath = 21;
						return;
					}
					else
					{
						//Can't have inherited from either
						childGeno.iu_PPC_err = true;
						childGeno.iu_allAllelesUnknown = true;
						childGeno.iu_denovoCandidate = true;
						childGeno.iu_unlikelyCall = true;
						childGeno.iu_unexpectedSegregation = true;
						childGeno.phasePath = 22;
						return;
					}
				}
			}
			else
			{
				//Neither parent has a deletion
				//Deletion must be denovo?
				childGeno.iu_unexpectedSegregation = true;
				childGeno.iu_denovoCandidate = true;
				if(p1i)
				{
					if (p2i)
					{
						//Unknown
						childGeno.iu_allAllelesUnknown = true;
						childGeno.phasePath = 23;
					}
					else
					{
						//Only could be inherited from p1
						ca.flagParent1(true);
						childGeno.phasePath = 24;
					}
				}
				else
				{
					if (p2i)
					{
						//Only could be inherited from p2
						ca.flagParent1(true);
						childGeno.phasePath = 25;
					}
					else
					{
						//Inherited from neither
						childGeno.iu_PPC_err = true;
						childGeno.iu_allAllelesUnknown = true;
						childGeno.iu_unlikelyCall = true;
						childGeno.phasePath = 26;
					}
				}
			}
		}
		
	}

	private static void phase_PC_CNV2_auto(SNPGeno parentGeno, SNPGeno childGeno)
	{
		childGeno.phasingAttempted = true;
		InheritanceMap childimap = childGeno.buildInheritanceMap(parentGeno, null);
		
		List<Integer> PC1 = canInherit_PC(childimap, true);
		
		if (PC1 == null)
		{
			//We got a PC error
			childGeno.flagAllOneParent(false);
			childGeno.iu_PC_err = true;
			childGeno.iu_unexpectedSegregation = true;
			childGeno.iu_denovoCandidate = true; //Denovo deletion of p1 allele?
			childGeno.iu_unknownGenoAssumed = true;
			//Are the remaining alleles the same?
			if(childimap.getKeySet().size() > 1)
			{
				//No
				childGeno.iu_UPHD_Candidate_p2 = true;
			}
			else
			{
				//Yes
				childGeno.iu_UPID_Candidate_p2 = true;
				childGeno.iu_inheritedDUP_candidate_p2 = true;
			}
			
			childGeno.phasePath = 27;
		}
		else
		{
			int p1i = PC1.size();
			if (p1i == 1)
			{
				//One allele kind could have been inherited from P1
				int pa = PC1.get(0);
				childGeno.flagAllNotOtherParent(pa, true);
				childGeno.flagOneAllele(pa, true);
				childGeno.iu_unknownGenoAssumed = true;
				//UPD or duplication is POSSIBLE...
				int p1c = childimap.getAllele(pa).getParent1Count();
				if (p1c > 1)
				{
					//P1 has more than one copy of heritable allele
					childGeno.iu_anyAllelesUnknown = true;
					childGeno.iu_UPID_Candidate_p1 = true;
					childGeno.phasePath = 28;
				}
				else if (p1c > 2)
				{
					//P1 has more than two copies of heritable allele
					childGeno.iu_anyAllelesUnknown = true;
					childGeno.iu_UPID_Candidate_p1 = true;
					childGeno.iu_inheritedDUP_candidate_p1 = true;
					childGeno.phasePath = 29;
				}
				else
				{
					//P1 has one copy of heritable allele
					childGeno.phasePath = 30;
				}
			}
			else if (p1i == 2)
			{
				//Two different alleles could have been inherited from P1
				//Can't tell which.
				childGeno.iu_allAllelesUnknown = true;
				if (UPHD_possible(childimap, true) != null) childGeno.iu_UPHD_Candidate_p1 = true;
				childGeno.iu_UPHD_Candidate_p2 = true;
				
				childGeno.phasePath = 31;
			}
		}
	}
	
	private static void phase_PPC_CNV2Plus_auto(SNPGeno p1Geno, SNPGeno p2Geno, SNPGeno childGeno)
	{
		childGeno.phasingAttempted = true;
		InheritanceMap childimap = childGeno.buildInheritanceMap(p1Geno, p2Geno);
		
		List<Integer> PPCE = get_uninheritableAlleles_PPC(childimap);
		
		List<Integer> PC1 = canInherit_PC(childimap, true);
		List<Integer> PC2 = canInherit_PC(childimap, false);
		List<Integer> PO1 = new ArrayList<Integer>(PC1.size());
		List<Integer> PO2 = new ArrayList<Integer>(PC2.size());
		List<Integer> PBoth = new LinkedList<Integer>();
		for (Integer i : PC1) if (!PC2.contains(i)) PO1.add(i);
		for (Integer i : PC2) if (!PC1.contains(i)) PO2.add(i);
		for (Integer i : PC1) if (PC2.contains(i)) PBoth.add(i);
		
		//Look for alleles that could have only been inherited from one or the other parent
		for(Integer i : PO1)
		{
			Allele as = childimap.getAllele(i);
			as.setFlag1Count(as.getAlleleCount());
		}
		for(Integer i : PO2)
		{
			Allele as = childimap.getAllele(i);
			as.setFlag2Count(as.getAlleleCount());
		}
		
		//See if either or both parents have been left unflagged
		boolean p1f = childimap.anyPhasedP1();
		boolean p2f = childimap.anyPhasedP2();
		//int uf = childimap.countUnphased();
		List<Allele> up = childimap.getUniqueUnphased();
		if (p1f && !p2f && (up != null))
		{
			//No allele exclusive to P2 found
			//Look for something to flag P2
			int us = up.size();
			if (us == 1){
				//Should be heritable from both. Flag one copy as p2
				Allele as = up.get(0);
				as.setFlag2Count(1);
				p2f = true;
			}
			//If if us is greater than 1, then impossible to tell which might have come from p2
			//Else if us < 1, we have no more alleles TO flag...
		}
		else if (p2f && !p1f && (up != null))
		{
			//No allele exclusive to P1 found
			int us = up.size();
			if (us == 1){
				//Should be heritable from both. Flag one copy as p2
				Allele as = up.get(0);
				as.setFlag1Count(1);
				p1f = true;
			}
		}
		
		//Recount alleles unphased and phased each way
		//Try to figure out what the deal is...
		int p1p = childimap.countPhasedP1();
		int p2p = childimap.countPhasedP2();
		int pup = childimap.countUnphased();
		
		if (PPCE == null)
		{
			//Normal, no alleles that came out of nowhere...
			
			//Count child alleles
			int acount = childGeno.getAlleleCount(); //In case mosaic
			int ccnv = childGeno.getCopyNumber();
			
			if (acount == ccnv)
			{
				//Normal, no evidence of mosaicism
				if (acount == 2)
				{
					//Again, normal.
					if (p1p == 1 && p2p == 1)
					{
						//We good.
						childGeno.phasePath = 32;
					}
					else
					{
						// 1/0, 2/0, or 0/0
						// 1/0 should not be possible due to above logic.
						if (p1p == 0 && p2p == 0)
						{
							//Neither allele could be phased, but both appear
							//to be heritable by either parent.
							childGeno.iu_allAllelesUnknown = true;
							childGeno.phasePath = 33;
						}
						else
						{
							//Check for iDEL(dnDEL)/iDUP(dnDUP) combo or UPD
							if (p1p == 2 && p2p == 0)
							{
								List<Integer> idup1 = possibleInheritedDuplications(childimap, true, ccnv, p1Geno.getCopyNumber());
								List<Integer> upid1 = UPID_possible(childimap, true);
								List<Integer> uphd1 = UPHD_possible(childimap, true);
								boolean idel2 = p2Geno.getCopyNumber() < 2;
								if (idel2) childGeno.iu_inheritedDEL_candidate_p2 = true;
								if (idup1 != null) childGeno.iu_inheritedDUP_candidate_p1 = true;
								if (upid1 != null) childGeno.iu_UPID_Candidate_p1 = true;
								if (uphd1 != null) childGeno.iu_UPHD_Candidate_p1 = true;
								if ((upid1 == null) && (uphd1 == null)) childGeno.iu_denovoCandidate = true;
								childGeno.iu_unexpectedSegregation = true;
								childGeno.phasePath = 34;
							}
							else if (p1p == 0 && p2p == 2)
							{
								List<Integer> idup2 = possibleInheritedDuplications(childimap, false, ccnv, p2Geno.getCopyNumber());
								List<Integer> upid2 = UPID_possible(childimap, false);
								List<Integer> uphd2 = UPHD_possible(childimap, false);
								boolean idel1 = p1Geno.getCopyNumber() < 2;
								if (idel1) childGeno.iu_inheritedDEL_candidate_p1 = true;
								if (idup2 != null) childGeno.iu_inheritedDUP_candidate_p2 = true;
								if (upid2 != null) childGeno.iu_UPID_Candidate_p2 = true;
								if (uphd2 != null) childGeno.iu_UPHD_Candidate_p2 = true;
								if ((upid2 == null) && (uphd2 == null)) childGeno.iu_denovoCandidate = true;
								childGeno.iu_unexpectedSegregation = true;
								childGeno.phasePath = 35;
							}
						}
					}
				}
				else
				{
					//Child cnv is >= 3
					List<Integer> idup1 = possibleInheritedDuplications(childimap, true, ccnv, p1Geno.getCopyNumber());
					List<Integer> upid1 = UPID_possible(childimap, true);
					List<Integer> uphd1 = UPHD_possible(childimap, true);
					List<Integer> idup2 = possibleInheritedDuplications(childimap, false, ccnv, p2Geno.getCopyNumber());
					List<Integer> upid2 = UPID_possible(childimap, false);
					List<Integer> uphd2 = UPHD_possible(childimap, false);
					List<Allele> p1Flagged = childimap.getUniquePhasedP1();
					List<Allele> p2Flagged = childimap.getUniquePhasedP2();
					List<Allele> unFlagged = childimap.getUniqueUnphased();
					
					if(pup > 0)
					{
						//There remain unphased alleles.
						//Do both parents have a marked allele?
						if (p1p > 0 && p2p == 0)
						{
							//P1 has allele(s), but p2 does not.
							if (p1p > 1)
							{
								//Check for each iDUP/UPD possibility, and whether each is possible with alleles flagged as they are.
								childGeno.iu_anyAllelesUnknown = true;
								//Check for idup
								if (idup1 != null) childGeno.iu_inheritedDUP_candidate_p1 = InheritanceMap.possible_idup(idup1, p1Flagged, unFlagged, true);
								if (idup2 != null) childGeno.iu_inheritedDUP_candidate_p2 = InheritanceMap.possible_idup(idup2, unFlagged);
								//Check for upid
								if (upid1 != null) childGeno.iu_UPID_Candidate_p1 = InheritanceMap.possible_upid(upid1, p1Flagged, unFlagged, true);
								if (upid2 != null) childGeno.iu_UPID_Candidate_p2 = InheritanceMap.possible_upid(upid2, unFlagged);
								//Check for uphd
								if (uphd1 != null) childGeno.iu_UPHD_Candidate_p1 = InheritanceMap.possible_uphd(uphd1, p1Flagged, unFlagged,true);
								if (uphd2 != null) childGeno.iu_UPHD_Candidate_p2 = InheritanceMap.possible_uphd(uphd2, unFlagged);
								//Assess for idel
								if (p2Geno.getCopyNumber() < 2 && childGeno.iu_inheritedDUP_candidate_p1 && childGeno.iu_UPHD_Candidate_p1) childGeno.iu_inheritedDEL_candidate_p2 = true;
								//Assess for denovo & unexpected seg
								int p1cnt = p1Geno.getAlleleCount(); //Most possible alleles to pass
								int p2cnt = p2Geno.getAlleleCount(); //Most possible alleles to pass
								if ((p1cnt + p2cnt) < ccnv)
								{
									childGeno.iu_unexpectedSegregation = true;
									childGeno.iu_denovoCandidate = true;
								}
								if (!childGeno.iu_inheritedDUP_candidate_p1 && !childGeno.iu_inheritedDUP_candidate_p2)
								{
									childGeno.iu_unexpectedSegregation = true;
									
								}
								childGeno.phasePath = 36;
							}
							else
							{
								//P1 has only one allele flagged.
								childGeno.iu_anyAllelesUnknown = true;
								//Check for idup
								if (idup1 != null) childGeno.iu_inheritedDUP_candidate_p1 = InheritanceMap.possible_idup(idup1, p1Flagged, unFlagged, true);
								if (idup2 != null) childGeno.iu_inheritedDUP_candidate_p2 = InheritanceMap.possible_idup(idup2, unFlagged);
								//Check for upid
								if (upid1 != null) childGeno.iu_UPID_Candidate_p1 = InheritanceMap.possible_upid(upid1, p1Flagged, unFlagged, true);
								if (upid2 != null) childGeno.iu_UPID_Candidate_p2 = InheritanceMap.possible_upid(upid2, unFlagged);
								//Check for uphd
								if (uphd1 != null) childGeno.iu_UPHD_Candidate_p1 = InheritanceMap.possible_uphd(uphd1, p1Flagged, unFlagged, true);
								if (uphd2 != null) childGeno.iu_UPHD_Candidate_p2 = InheritanceMap.possible_uphd(uphd2, unFlagged);
								//Assess for idel
								if (p1Geno.getCopyNumber() < 2) childGeno.iu_inheritedDEL_candidate_p1 = true;
								if (p2Geno.getCopyNumber() < 2) childGeno.iu_inheritedDEL_candidate_p2 = true;
								//Assess for denovo & unexpected seg
								int p1cnt = p1Geno.getAlleleCount(); //Most possible alleles to pass
								int p2cnt = p2Geno.getAlleleCount(); //Most possible alleles to pass
								if ((p1cnt + p2cnt) < ccnv)
								{
									childGeno.iu_unexpectedSegregation = true;
									childGeno.iu_denovoCandidate = true;
								}
								if (!childGeno.iu_inheritedDUP_candidate_p1 && !childGeno.iu_inheritedDUP_candidate_p2)
								{
									childGeno.iu_unexpectedSegregation = true;
									
								}
								childGeno.phasePath = 37;
							}
						}
						else if (p1p == 0 && p2p > 0)
						{
							//P2 has allele(s), but p1 does not.
							if (p1p > 1)
							{
								//Check for each iDUP/UPD possibility, and whether each is possible with alleles flagged as they are.
								childGeno.iu_anyAllelesUnknown = true;
								//Check for idup
								if (idup1 != null) childGeno.iu_inheritedDUP_candidate_p1 = InheritanceMap.possible_idup(idup1, unFlagged);
								if (idup2 != null) childGeno.iu_inheritedDUP_candidate_p2 = InheritanceMap.possible_idup(idup2, p2Flagged, unFlagged, false);
								//Check for upid
								if (upid1 != null) childGeno.iu_UPID_Candidate_p1 = InheritanceMap.possible_upid(upid1, unFlagged);
								if (upid2 != null) childGeno.iu_UPID_Candidate_p2 = InheritanceMap.possible_upid(upid2, p2Flagged, unFlagged, false);
								//Check for uphd
								if (uphd1 != null) childGeno.iu_UPHD_Candidate_p1 = InheritanceMap.possible_uphd(uphd1, unFlagged);
								if (uphd2 != null) childGeno.iu_UPHD_Candidate_p2 = InheritanceMap.possible_uphd(uphd2, p2Flagged, unFlagged, false);
								//Assess for idel
								if (p1Geno.getCopyNumber() < 2 && childGeno.iu_inheritedDUP_candidate_p2 && childGeno.iu_UPHD_Candidate_p2) childGeno.iu_inheritedDEL_candidate_p1 = true;
								//Assess for denovo & unexpected seg
								int p1cnt = p1Geno.getAlleleCount(); //Most possible alleles to pass
								int p2cnt = p2Geno.getAlleleCount(); //Most possible alleles to pass
								if ((p1cnt + p2cnt) < ccnv)
								{
									childGeno.iu_unexpectedSegregation = true;
									childGeno.iu_denovoCandidate = true;
								}
								if (!childGeno.iu_inheritedDUP_candidate_p1 && !childGeno.iu_inheritedDUP_candidate_p2)
								{
									childGeno.iu_unexpectedSegregation = true;
									
								}
								childGeno.phasePath = 38;
							}
							else
							{
								//P1 has only one allele flagged.
								childGeno.iu_anyAllelesUnknown = true;
								//Check for idup
								if (idup1 != null) childGeno.iu_inheritedDUP_candidate_p1 = InheritanceMap.possible_idup(idup1, unFlagged);
								if (idup2 != null) childGeno.iu_inheritedDUP_candidate_p2 = InheritanceMap.possible_idup(idup2, p2Flagged, unFlagged, false);
								//Check for upid
								if (upid1 != null) childGeno.iu_UPID_Candidate_p1 = InheritanceMap.possible_upid(upid1, unFlagged);
								if (upid2 != null) childGeno.iu_UPID_Candidate_p2 = InheritanceMap.possible_upid(upid2, p2Flagged, unFlagged, false);
								//Check for uphd
								if (uphd1 != null) childGeno.iu_UPHD_Candidate_p1 = InheritanceMap.possible_uphd(uphd1, unFlagged);
								if (uphd2 != null) childGeno.iu_UPHD_Candidate_p2 = InheritanceMap.possible_uphd(uphd2, p2Flagged, unFlagged, false);
								//Assess for idel
								if (p1Geno.getCopyNumber() < 2) childGeno.iu_inheritedDEL_candidate_p1 = true;
								if (p2Geno.getCopyNumber() < 2) childGeno.iu_inheritedDEL_candidate_p2 = true;
								//Assess for denovo & unexpected seg
								int p1cnt = p1Geno.getAlleleCount(); //Most possible alleles to pass
								int p2cnt = p2Geno.getAlleleCount(); //Most possible alleles to pass
								if ((p1cnt + p2cnt) < ccnv)
								{
									childGeno.iu_unexpectedSegregation = true;
									childGeno.iu_denovoCandidate = true;
								}
								if (!childGeno.iu_inheritedDUP_candidate_p1 && !childGeno.iu_inheritedDUP_candidate_p2)
								{
									childGeno.iu_unexpectedSegregation = true;
									
								}
								childGeno.phasePath = 39;
							}
						}
						else if (p1p == 0 && p2p == 0)
						{
							//All unphased and cnv >= 3
							//Pretty much anything is possible...
							childGeno.iu_allAllelesUnknown = true;
							//Check for idup
							if (idup1 != null) childGeno.iu_inheritedDUP_candidate_p1 = InheritanceMap.possible_idup(idup1, unFlagged);
							if (idup2 != null) childGeno.iu_inheritedDUP_candidate_p2 = InheritanceMap.possible_idup(idup2, unFlagged);
							//Check for upid
							if (upid1 != null) childGeno.iu_UPID_Candidate_p1 = InheritanceMap.possible_upid(upid1, unFlagged);
							if (upid2 != null) childGeno.iu_UPID_Candidate_p2 = InheritanceMap.possible_upid(upid2, unFlagged);
							//Check for uphd
							if (uphd1 != null) childGeno.iu_UPHD_Candidate_p1 = InheritanceMap.possible_uphd(uphd1, unFlagged);
							if (uphd2 != null) childGeno.iu_UPHD_Candidate_p2 = InheritanceMap.possible_uphd(uphd2, unFlagged);
							//Assess for idel
							if (p1Geno.getCopyNumber() < 2) childGeno.iu_inheritedDEL_candidate_p1 = true;
							if (p2Geno.getCopyNumber() < 2) childGeno.iu_inheritedDEL_candidate_p2 = true;
							//Assess for denovo & unexpected seg
							int p1cnt = p1Geno.getAlleleCount(); //Most possible alleles to pass
							int p2cnt = p2Geno.getAlleleCount(); //Most possible alleles to pass
							if ((p1cnt + p2cnt) < ccnv)
							{
								childGeno.iu_unexpectedSegregation = true;
								childGeno.iu_denovoCandidate = true;
							}
							if (!childGeno.iu_inheritedDUP_candidate_p1 && !childGeno.iu_inheritedDUP_candidate_p2)
							{
								childGeno.iu_unexpectedSegregation = true;
							}
							
							
							childGeno.phasePath = 40;
						}
						else
						{
							//At least one unphased, p1, p2
							childGeno.iu_anyAllelesUnknown = true;
							//Check for idup
							if (idup1 != null) childGeno.iu_inheritedDUP_candidate_p1 = InheritanceMap.possible_idup(idup1, p1Flagged, unFlagged, true);
							if (idup2 != null) childGeno.iu_inheritedDUP_candidate_p2 = InheritanceMap.possible_idup(idup2, p2Flagged, unFlagged, false);
							//Check for upid
							if (upid1 != null) childGeno.iu_UPID_Candidate_p1 = InheritanceMap.possible_upid(upid1, p1Flagged, unFlagged, true);
							if (upid2 != null) childGeno.iu_UPID_Candidate_p2 = InheritanceMap.possible_upid(upid2, p2Flagged, unFlagged, false);
							//Check for uphd
							if (uphd1 != null) childGeno.iu_UPHD_Candidate_p1 = InheritanceMap.possible_uphd(uphd1, p1Flagged, unFlagged, true);
							if (uphd2 != null) childGeno.iu_UPHD_Candidate_p2 = InheritanceMap.possible_uphd(uphd2, p2Flagged, unFlagged, false);
							//Skip idel
							//Assess for denovo & unexpected seg
							int p1cnt = p1Geno.getAlleleCount(); //Most possible alleles to pass
							int p2cnt = p2Geno.getAlleleCount(); //Most possible alleles to pass
							if ((p1cnt + p2cnt) < ccnv)
							{
								childGeno.iu_unexpectedSegregation = true;
								childGeno.iu_denovoCandidate = true;
							}
							if (!childGeno.iu_inheritedDUP_candidate_p1 && !childGeno.iu_inheritedDUP_candidate_p2)
							{
								childGeno.iu_unexpectedSegregation = true;
							}
							childGeno.phasePath = 41;
						}
					}
					else
					{
						//All alleles are phased. Get counts from each parent to figure out what's going on
						//Where is the duplication coming from?
						if (p1p == 0)
						{
							//None from p1
							// >1 from p2
							if (p1Geno.getCopyNumber() < 2) childGeno.iu_inheritedDEL_candidate_p1 = true;
							if (idup2 != null) childGeno.iu_inheritedDUP_candidate_p2 = true;
							if (upid2 != null) childGeno.iu_UPID_Candidate_p2 = true;
							if (uphd2 != null) childGeno.iu_UPHD_Candidate_p2 = true;
							if (ccnv > p2Geno.getAlleleCount()){
								childGeno.iu_unexpectedSegregation = true;
								childGeno.iu_denovoCandidate = true;
							}
							if (!childGeno.iu_inheritedDEL_candidate_p1)
							{
								childGeno.iu_unexpectedSegregation = true;
								childGeno.iu_denovoCandidate = true;
							}
							childGeno.phasePath = 42;
						}
						else if (p2p == 0)
						{
							//None from p2
							// >1 from p1
							if (p2Geno.getCopyNumber() < 2) childGeno.iu_inheritedDEL_candidate_p2 = true;
							if (idup1 != null) childGeno.iu_inheritedDUP_candidate_p1 = true;
							if (upid1 != null) childGeno.iu_UPID_Candidate_p1 = true;
							if (uphd1 != null) childGeno.iu_UPHD_Candidate_p1 = true;
							if (ccnv > p1Geno.getAlleleCount()){
								childGeno.iu_unexpectedSegregation = true;
								childGeno.iu_denovoCandidate = true;
							}
							if (!childGeno.iu_inheritedDEL_candidate_p2)
							{
								childGeno.iu_unexpectedSegregation = true;
								childGeno.iu_denovoCandidate = true;
							}
							childGeno.phasePath = 43;
						}
						else
						{
							//At least 1 from each.
							//Who has >1 ?
							if (p1p > 1)
							{
								if (idup1 != null) childGeno.iu_inheritedDUP_candidate_p1 = true;
								if (upid1 != null) childGeno.iu_UPID_Candidate_p1 = true;
								if (uphd1 != null) childGeno.iu_UPHD_Candidate_p1 = true;	
							}
							if (p2p > 1)
							{
								if (idup2 != null) childGeno.iu_inheritedDUP_candidate_p2 = true;
								if (upid2 != null) childGeno.iu_UPID_Candidate_p2 = true;
								if (uphd2 != null) childGeno.iu_UPHD_Candidate_p2 = true;	
							}
							if (ccnv > (p1Geno.getAlleleCount() + p2Geno.getAlleleCount())){
								childGeno.iu_unexpectedSegregation = true;
								childGeno.iu_denovoCandidate = true;
							}
							if (!childGeno.iu_inheritedDUP_candidate_p1 && !childGeno.iu_inheritedDUP_candidate_p2)
							{
								childGeno.iu_unexpectedSegregation = true;
								childGeno.iu_denovoCandidate = true;
							}
							childGeno.phasePath = 44;
						}
					}
				}
			}
			else
			{
				//Possible mosaicism? 
				//Maybe ignore for now... detect and deal with using another method?
				childGeno.iu_mosaic_candidate = true;
				childGeno.phasePath = 45;
			}
			
			childimap.flagLinkedAlleles();
			
		}
		else
		{
			//Has at least one allele that could not have come from either parent
			childGeno.iu_PPC_err = true;
			childGeno.iu_denovoCandidate = true;
			childGeno.iu_unexpectedSegregation = true;
			
			//Flag error alleles
			for (Integer i : PPCE) childimap.flagAlleleAsError(i);
			List<Allele> upp = childimap.getUniqueUnphased();
			//int pep = childimap.countErrorAlleles();
			pup = childimap.countUnphased();
			
			int ccnv = childGeno.getCopyNumber();
			
			List<Integer> idup1 = possibleInheritedDuplications(childimap, true, ccnv, p1Geno.getCopyNumber());
			List<Integer> upid1 = UPID_possible(childimap, true);
			List<Integer> uphd1 = UPHD_possible(childimap, true);
			List<Integer> idup2 = possibleInheritedDuplications(childimap, false, ccnv, p2Geno.getCopyNumber());
			List<Integer> upid2 = UPID_possible(childimap, false);
			List<Integer> uphd2 = UPHD_possible(childimap, false);
			List<Allele> p1Flagged = childimap.getUniquePhasedP1();
			List<Allele> p2Flagged = childimap.getUniquePhasedP2();
			List<Allele> unFlagged = childimap.getUniqueUnphased();
			
			//There will be unphased alleles no matter what
			if (p1p == 0 && p2p == 0)
			{
				//Unphased alleles could come from either both or neither.
				// [2] xx
				// [3] xxx x??
				// [4] x??? xx?? xxxx
				// [5] x???? xx??? xxx?? xxxxx
				
				if (upp.isEmpty())
				{
					//There are no alleles here that could have come from either parent.
					childGeno.phasePath = 46;
				}
				else
				{
					childGeno.iu_anyAllelesUnknown = true;
					//Check for idel
					if (p1Geno.getCopyNumber() < 2) childGeno.iu_inheritedDEL_candidate_p1 = true;
					if (p2Geno.getCopyNumber() < 2) childGeno.iu_inheritedDEL_candidate_p2 = true;
					if (ccnv > 2)
					{
						//If ccnv is 3+, check for idup
						if (idup1 != null) childGeno.iu_inheritedDUP_candidate_p1 = InheritanceMap.possible_idup(idup1, p1Flagged, upp, true);
						if (idup2 != null) childGeno.iu_inheritedDUP_candidate_p2 = InheritanceMap.possible_idup(idup2, p2Flagged, upp, false);
						//If ccnv is 3+, check for upd
						if (upid1 != null) childGeno.iu_UPID_Candidate_p1 = InheritanceMap.possible_upid(upid1, p1Flagged, upp, true);
						if (upid2 != null) childGeno.iu_UPID_Candidate_p2 = InheritanceMap.possible_upid(upid2, p2Flagged, upp, false);
						if (uphd1 != null) childGeno.iu_UPHD_Candidate_p1 = InheritanceMap.possible_uphd(uphd1, p1Flagged, upp, true);
						if (uphd2 != null) childGeno.iu_UPHD_Candidate_p2 = InheritanceMap.possible_uphd(uphd2, p2Flagged, upp, false);
					}
					else if (ccnv == 2)
					{
						//If child has cnv2, one parent has idel, and the other doesn't, flag upp allele as parent without idel
						if (childGeno.iu_inheritedDEL_candidate_p1 && !childGeno.iu_inheritedDEL_candidate_p2)
						{
							//P1 has idel. Flag remaining allele as p2.
							upp.get(0).setFlag2Count(1);
						}
						else if(childGeno.iu_inheritedDEL_candidate_p2 && !childGeno.iu_inheritedDEL_candidate_p1)
						{
							//P2 has idel
							upp.get(0).setFlag1Count(1);
						}
					}
					childGeno.phasePath = 47;
				}
			}
			else if (p1p > 0 && p2p == 0)
			{
				//At least 1 phased P1, none phased P2
				// [2] 1x
				// [3] 1xx 11x
				// [4] 1xxx 1x?? 11xx 111x
				// [5] 1xxxx 1x??? 1xx?? 11x?? 11xxx 111xx 1111x
				
				//Check for true unphased...
				if (pup > 0)
				{
					//Ignore error alleles, behave as it would for non-error.
					childGeno.iu_anyAllelesUnknown = true;
					//Check for idup
					if (idup1 != null) childGeno.iu_inheritedDUP_candidate_p1 = InheritanceMap.possible_idup(idup1, p1Flagged, unFlagged, true);
					if (idup2 != null) childGeno.iu_inheritedDUP_candidate_p2 = InheritanceMap.possible_idup(idup2, unFlagged);
					//Check for upid
					if (upid1 != null) childGeno.iu_UPID_Candidate_p1 = InheritanceMap.possible_upid(upid1, p1Flagged, unFlagged, true);
					if (upid2 != null) childGeno.iu_UPID_Candidate_p2 = InheritanceMap.possible_upid(upid2, unFlagged);
					//Check for uphd
					if (uphd1 != null) childGeno.iu_UPHD_Candidate_p1 = InheritanceMap.possible_uphd(uphd1, p1Flagged, unFlagged, true);
					if (uphd2 != null) childGeno.iu_UPHD_Candidate_p2 = InheritanceMap.possible_uphd(uphd2, unFlagged);
					
					childGeno.phasePath = 48;
				}
				else
				{
					if (p2Geno.getCopyNumber() < 2) childGeno.iu_inheritedDEL_candidate_p2 = true;
					if (p1p > 1)
					{
						//Check for idup
						if (idup1 != null) childGeno.iu_inheritedDUP_candidate_p1 = InheritanceMap.possible_idup(idup1, p1Flagged, unFlagged, true);
						//Check for upid
						if (upid1 != null) childGeno.iu_UPID_Candidate_p1 = InheritanceMap.possible_upid(upid1, p1Flagged, unFlagged, true);
						//Check for uphd
						if (uphd1 != null) childGeno.iu_UPHD_Candidate_p1 = InheritanceMap.possible_uphd(uphd1, p1Flagged, unFlagged, true);		
					}
					
					childGeno.phasePath = 49;
					//If phasepath 49 with no idel, deletion from other parent is assumed denovo or error
				}
				
			}
			else if (p1p == 0 && p2p > 0)
			{
				//At least 1 phased P2, none phased P1
				// [2] 2x
				// [3] 2xx 22x
				// [4] 2xxx 2x?? 22xx 222x
				// [5] 2xxxx 2x??? 2xx?? 22x?? 22xxx 222xx 2222x
				
				//Check for true unphased...
				if (pup > 0)
				{
					//Ignore error alleles, behave as it would for non-error.
					childGeno.iu_anyAllelesUnknown = true;
					//Check for idup
					if (idup2 != null) childGeno.iu_inheritedDUP_candidate_p2 = InheritanceMap.possible_idup(idup2, p2Flagged, unFlagged, false);
					if (idup1 != null) childGeno.iu_inheritedDUP_candidate_p1 = InheritanceMap.possible_idup(idup1, unFlagged);
					//Check for upid
					if (upid2 != null) childGeno.iu_UPID_Candidate_p2 = InheritanceMap.possible_upid(upid2, p2Flagged, unFlagged, false);
					if (upid1 != null) childGeno.iu_UPID_Candidate_p1 = InheritanceMap.possible_upid(upid1, unFlagged);
					//Check for uphd
					if (uphd2 != null) childGeno.iu_UPHD_Candidate_p2 = InheritanceMap.possible_uphd(uphd2, p2Flagged, unFlagged, false);
					if (uphd1 != null) childGeno.iu_UPHD_Candidate_p1 = InheritanceMap.possible_uphd(uphd1, unFlagged);
					
					childGeno.phasePath = 50;
				}
				else
				{
					if (p1Geno.getCopyNumber() < 2) childGeno.iu_inheritedDEL_candidate_p1 = true;
					if (p2p > 1)
					{
						//Check for idup
						if (idup2 != null) childGeno.iu_inheritedDUP_candidate_p2 = InheritanceMap.possible_idup(idup2, p2Flagged, unFlagged, false);
						//Check for upid
						if (upid2 != null) childGeno.iu_UPID_Candidate_p2 = InheritanceMap.possible_upid(upid2, p2Flagged, unFlagged, false);
						//Check for uphd
						if (uphd2 != null) childGeno.iu_UPHD_Candidate_p2 = InheritanceMap.possible_uphd(uphd2, p2Flagged, unFlagged, false);	
					}
										
					childGeno.phasePath = 51;
					//If phasepath 51 with no idel, deletion from other parent is assumed denovo or error
				}
				
			}
			else
			{
				//At least 1 phased P1 and 1 phased P2. Since at least one PPCE allele, CNV must be 3+
				//Make sure to remove any from unphased list that have since been flagged as P1 or P2.
				//Or generate new unphased list and remove PPCE alleles
				
				// [2] (None)
				// [3] 12x
				// [4] 12xx 112x 122x
				// [5] 12xxx 12x?? 112xx 1112x 122xx 1222x 1122x
				
				//Have to get to a pretty high CN to get true unphased... but check for anyway.
				if (pup > 0)
				{
					childGeno.iu_anyAllelesUnknown = true;
					//Check for idup
					if (idup1 != null) childGeno.iu_inheritedDUP_candidate_p1 = InheritanceMap.possible_idup(idup1, p1Flagged, unFlagged, true);
					if (idup2 != null) childGeno.iu_inheritedDUP_candidate_p2 = InheritanceMap.possible_idup(idup2, p2Flagged, unFlagged, false);
					//Check for upid
					if (upid1 != null) childGeno.iu_UPID_Candidate_p1 = InheritanceMap.possible_upid(upid1, p1Flagged, unFlagged, true);
					if (upid2 != null) childGeno.iu_UPID_Candidate_p2 = InheritanceMap.possible_upid(upid2, p2Flagged, unFlagged, false);
					//Check for uphd
					if (uphd1 != null) childGeno.iu_UPHD_Candidate_p1 = InheritanceMap.possible_uphd(uphd1, p1Flagged, unFlagged, true);
					if (uphd2 != null) childGeno.iu_UPHD_Candidate_p2 = InheritanceMap.possible_uphd(uphd2, p2Flagged, unFlagged, false);
					
					childGeno.phasePath = 52;
				}
				else
				{
					if (p1p > 1)
					{
						//Check for idup
						if (idup1 != null) childGeno.iu_inheritedDUP_candidate_p1 = InheritanceMap.possible_idup(idup1, p1Flagged, unFlagged, true);
						//Check for upid
						if (upid1 != null) childGeno.iu_UPID_Candidate_p1 = InheritanceMap.possible_upid(upid1, p1Flagged, unFlagged, true);
						//Check for uphd
						if (uphd1 != null) childGeno.iu_UPHD_Candidate_p1 = InheritanceMap.possible_uphd(uphd1, p1Flagged, unFlagged, true);	
					}
					if (p2p > 1)
					{
						//Check for idup
						if (idup2 != null) childGeno.iu_inheritedDUP_candidate_p2 = InheritanceMap.possible_idup(idup2, p2Flagged, unFlagged, false);
						//Check for upid
						if (upid2 != null) childGeno.iu_UPID_Candidate_p2 = InheritanceMap.possible_upid(upid2, p2Flagged, unFlagged, false);
						//Check for uphd
						if (uphd2 != null) childGeno.iu_UPHD_Candidate_p2 = InheritanceMap.possible_uphd(uphd2, p2Flagged, unFlagged, false);
					}
					childGeno.phasePath = 53;
				}
				
				
			}
			
			childimap.flagLinkedAlleles();
		}
		
	}
	
	private static void phase_PC_CNV3Plus_auto(SNPGeno parentGeno, SNPGeno childGeno)
	{
		childGeno.phasingAttempted = true;
		InheritanceMap childimap = childGeno.buildInheritanceMap(parentGeno, null);
		List<Integer> PC1 = canInherit_PC(childimap, true);
		List<Integer> PCE = get_uninheritableAlleles_PC(childimap);
		
		if (PC1 == null) childGeno.iu_PC_err = true;
		
		//Flag all that can't have from from P1 as P2
		for(Integer i : PCE)
		{
			Allele as = childimap.getAllele(i);
			as.setFlag2Count(as.getAlleleCount());
			childGeno.iu_unknownGenoAssumed = true;
		}
		//If only one allele type left, grab one allele and flag P1
		boolean p2f = childimap.anyPhasedP2();
		List<Allele> up = childimap.getUniqueUnphased();
		if (up != null)
		{
			int us = up.size();
			if (us == 1){
				//Should be heritable from both. Flag one copy as p2
				Allele as = up.get(0);
				as.setFlag1Count(1);
				if (!p2f) as.setFlag2Count(1);
			}
		}
		
		//Now assess results
		int p1p = childimap.countPhasedP1();
		int p2p = childimap.countPhasedP2();
		int pup = childimap.countUnphased();
		
		List<Integer> idup1 = possibleInheritedDuplications(childimap, true, childGeno.getCopyNumber(), parentGeno.getCopyNumber());
		List<Integer> upid1 = UPID_possible(childimap, true);
		List<Integer> uphd1 = UPHD_possible(childimap, true);
		List<Allele> p1Flagged = childimap.getUniquePhasedP1();
		//List<Allele> p2Flagged = childimap.getUniquePhasedP2();
		List<Allele> unFlagged = childimap.getUniqueUnphased();
		
		//Can't be more than one p1 allele, I think. And only if there is a 2 also.
		if (pup > 0)
		{
			//[3] 12? 2?? ??? 
			//[4] 12?? 122? 2??? 22?? 
			//[5] 12??? 122?? 1222? 2???? 22??? 222??
			
			//Assume uphd2 is possible be default.
			//If there are at least two identical phased 2 or two identical unphased, assume upid2 and idup2 are possible.
			if (idup1 != null) childGeno.iu_inheritedDUP_candidate_p1 = InheritanceMap.possible_idup(idup1, p1Flagged, unFlagged, true);
			if (pup == 1 && childGeno.iu_inheritedDUP_candidate_p1)
			{
				//Flag remaining up allele as p1
				Allele a = unFlagged.get(0);
				a.setFlag1Count(a.getFlag1Count() + 1);
				//Check for p2 duplications...
				if (p2p > 1)
				{
					childGeno.iu_UPHD_Candidate_p2 = true;
					for (Integer i : PCE)
					{
						if (childimap.getAllele(i).getAlleleCount() > 1)
						{
							childGeno.iu_inheritedDUP_candidate_p2 = true;
							childGeno.iu_UPID_Candidate_p2 = true;
							break;
						}
					}
				}
				p1Flagged = childimap.getUniquePhasedP1();
				if (upid1 != null) childGeno.iu_UPID_Candidate_p1 = InheritanceMap.possible_upid(upid1, p1Flagged, unFlagged, true);
				if (uphd1 != null) childGeno.iu_UPHD_Candidate_p1 = InheritanceMap.possible_uphd(uphd1, p1Flagged, unFlagged, true);
				childGeno.phasePath = 54;
			}
			else
			{
				//Leave and evaluate as is.
				childGeno.iu_UPHD_Candidate_p2 = true;
				for (Integer i : PCE)
				{
					if (childimap.getAllele(i).getAlleleCount() > 1)
					{
						childGeno.iu_inheritedDUP_candidate_p2 = true;
						childGeno.iu_UPID_Candidate_p2 = true;
						break;
					}
				}
				if (upid1 != null) childGeno.iu_UPID_Candidate_p1 = InheritanceMap.possible_upid(upid1, p1Flagged, unFlagged, true);
				if (uphd1 != null) childGeno.iu_UPHD_Candidate_p1 = InheritanceMap.possible_uphd(uphd1, p1Flagged, unFlagged, true);
				
				childGeno.phasePath = 55;
			}
			
		}
		else
		{
			//[3] 122 222
			//[4] 1222 2222
			//[5] 12222 22222
			
			//Check for p2 dup things...
			childGeno.iu_UPHD_Candidate_p2 = true;
			for (Integer i : PCE)
			{
				if (childimap.getAllele(i).getAlleleCount() > 1)
				{
					childGeno.iu_inheritedDUP_candidate_p2 = true;
					childGeno.iu_UPID_Candidate_p2 = true;
					break;
				}
			}
			
			if (p1p == 0)
			{
				childGeno.iu_unexpectedSegregation = true;
				if (parentGeno.getCopyNumber() < 2) childGeno.iu_inheritedDEL_candidate_p1 = true;
				else childGeno.iu_denovoCandidate = true;
			}
			
			childGeno.phasePath = 56;
		}
		
		childimap.flagLinkedAlleles();
	}
	
	public static void phase_PC_sxc_mammalX(SNPGeno parentGeno, SNPGeno childGeno, Sex psx, Sex csx)
	{
		childGeno.phasingAttempted = true;
		InheritanceMap childimap = childGeno.buildInheritanceMap(parentGeno, null);
		int ccnv = childGeno.getCopyNumber();
		int pcnv = parentGeno.getCopyNumber();
		boolean punsx = (psx == Sex.OTHER || psx == Sex.UNKNOWN);
		
		if (csx == Sex.MALE)
		{
			//Expected CNV 1
			if (ccnv == 1)
			{
				//Assumed to be from female parent.
				if (psx == Sex.FEMALE || (punsx && pcnv > 1))
				{
					//Check to see if could come from that parent.
					List<Integer> pc = canInherit_PC(childimap, true);
					if (pc != null)
					{
						//All is well.
						childimap.getAllele(pc.get(0)).setFlag1Count(1);
						childGeno.phasePath = 57;
					}
					else
					{
						//Weird.
						childGeno.iu_PC_err = true;
						childGeno.iu_unexpectedSegregation = true;
						childGeno.iu_allAllelesUnknown = true;
						childGeno.phasePath = 58;
					}
				}
				else if (psx == Sex.MALE || (punsx && pcnv < 2))
				{
					//Flag as other parent.
					childimap.getAllele(childGeno.getAlleles().get(0).getAllele()).setFlag2Count(1);
					childGeno.iu_unknownGenoAssumed = true;
					childGeno.phasePath = 59;
				}
			}
			else if (ccnv == 0)
			{
				//Nothing to phase, but need to flag.
				//If known parent is female, check for deletion
				//If known parent is male, assume other parent has idel
				//If known parent sex is unknown, go by CNV
				if (psx == Sex.FEMALE || (punsx && pcnv > 1))
				{
					if (pcnv < 2) childGeno.iu_inheritedDEL_candidate_p1 = true;
					else{
						childGeno.iu_denovoCandidate = true;
						childGeno.iu_unexpectedSegregation = true;
					}
					childGeno.phasePath = 60;
				}
				else if (psx == Sex.MALE || (punsx && pcnv < 2))
				{
					childGeno.iu_inheritedDEL_candidate_p2 = true;
					childGeno.iu_denovoCandidate = true;
					childGeno.phasePath = 61;
				}
			}
			else
			{
				//ccnv >= 2
				//iDUP: Local Maternal
				//UPID: Maternal
				//UPHD: Maternal or Paternal
				
				//[2] MM MF FF
				//[3] MMM MMF MFF FFF 
				//[4] MMMM MMMF MMFF MFFF FFFF
				
				//Get what could not have been inherited from known parent...
				List<Integer> PCE = get_uninheritableAlleles_PC(childimap);
				List<Integer> PC1 = canInherit_PC(childimap, true);
				
				if (psx == Sex.FEMALE || (punsx && pcnv > 1))
				{
					if (PCE != null)
					{
						for (Integer i : PCE)
						{
							Allele as = childimap.getAllele(i);
							as.setFlag2Count(as.getAlleleCount());
							if (as.getAlleleCount() > 1)
							{
								childGeno.iu_inheritedDUP_candidate_p2 = true;
								childGeno.iu_UPID_Candidate_p2 = true;
							}
						}
						childGeno.iu_unknownGenoAssumed = true;
						if (PC1 == null)
						{
							//Must be all paternal. WEIRD.
							childGeno.iu_unexpectedSegregation = true;
							childGeno.iu_PC_err = true;
							childGeno.iu_UPHD_Candidate_p2 = true;
							childGeno.phasePath = 62;
						}
						else
						{
							//PCE alleles are paternal, PC1 alleles are either
							//Assume at least one is maternal
							//[2] FM
							//[3] F?? F?M FFM
							//[4] F??? F??M FF?? FF?M FFFM
							
							//If only one left unphased, phase maternal
							if (PC1.size() == 1)
							{
								Allele a = childimap.getAllele(PC1.get(0));
								a.setFlag1Count(1);
							}
							
							int p1p = childimap.countPhasedP1();
							int p2p = childimap.countPhasedP2();
							int pup = childimap.countUnphased();
							
							List<Integer> idup1 = possibleInheritedDuplications(childimap, true, childGeno.getCopyNumber(), parentGeno.getCopyNumber());
							List<Integer> upid1 = UPID_possible(childimap, true);
							List<Integer> uphd1 = UPHD_possible(childimap, true);
							List<Allele> p1Flagged = childimap.getUniquePhasedP1();
							//List<Allele> p2Flagged = childimap.getUniquePhasedP2();
							List<Allele> unFlagged = childimap.getUniqueUnphased();
							
							if ((p1p > 0 && pup > 0) || pup > 1)
							{
								if (idup1 != null) childGeno.iu_inheritedDUP_candidate_p1 = InheritanceMap.possible_idup(idup1, p1Flagged, unFlagged, true);
								if (upid1 != null) childGeno.iu_UPID_Candidate_p1 = InheritanceMap.possible_upid(upid1, p1Flagged, unFlagged, true);
								if (uphd1 != null) childGeno.iu_UPHD_Candidate_p1 = InheritanceMap.possible_uphd(uphd1, p1Flagged, unFlagged, true);
							}
							if (p2p > 1) childGeno.iu_UPHD_Candidate_p2 = true;
							if (pup > 0) childGeno.iu_anyAllelesUnknown = true;
							childGeno.phasePath = 63;
						}
					}
					else
					{
						//All could be from either parent. Assume at least one is maternal.
						childGeno.iu_unknownGenoAssumed = true;
						if (PC1.size() == 1)
						{
							Allele a = childimap.getAllele(PC1.get(0));
							a.setFlag1Count(1);
						}
						//[2] M? ??
						//[3] M?? ???
						//[4] M??? ????
						int p1p = childimap.countPhasedP1();
						int pup = childimap.countUnphased();
						if (p1p > 0) childGeno.iu_anyAllelesUnknown = true;
						else childGeno.iu_allAllelesUnknown = true;
						if (pup > 0) childGeno.iu_UPHD_Candidate_p2 = true;
						
						//Check for maternal shenanigans
						List<Integer> idup1 = possibleInheritedDuplications(childimap, true, childGeno.getCopyNumber(), parentGeno.getCopyNumber());
						List<Integer> upid1 = UPID_possible(childimap, true);
						List<Integer> uphd1 = UPHD_possible(childimap, true);
						List<Allele> p1Flagged = childimap.getUniquePhasedP1();
						//List<Allele> p2Flagged = childimap.getUniquePhasedP2();
						List<Allele> unFlagged = childimap.getUniqueUnphased();
						
						if (idup1 != null) childGeno.iu_inheritedDUP_candidate_p1 = InheritanceMap.possible_idup(idup1, p1Flagged, unFlagged, true);
						if (upid1 != null) childGeno.iu_UPID_Candidate_p1 = InheritanceMap.possible_upid(upid1, p1Flagged, unFlagged, true);
						if (uphd1 != null) childGeno.iu_UPHD_Candidate_p1 = InheritanceMap.possible_uphd(uphd1, p1Flagged, unFlagged, true);
						
						for (Allele a : unFlagged)
						{
							if (a.getAlleleCount() > 1)
							{
								childGeno.iu_inheritedDUP_candidate_p2 = true;
								childGeno.iu_UPID_Candidate_p2 = true;
							}
						}
						
						childGeno.phasePath = 64;
					}
					
				}
				else if (psx == Sex.MALE || (punsx && pcnv < 2))
				{
					//Known parent is dad
					if (PCE != null)
					{
						for (Integer i : PCE)
						{
							Allele as = childimap.getAllele(i);
							as.setFlag2Count(as.getAlleleCount());
							if (as.getAlleleCount() > 1)
							{
								childGeno.iu_inheritedDUP_candidate_p2 = true;
								childGeno.iu_UPID_Candidate_p2 = true;
							}
						}
						childGeno.iu_unknownGenoAssumed = true;
						if (PC1 == null)
						{
							//Must be all maternal.
							childGeno.iu_UPHD_Candidate_p2 = true;
							childGeno.phasePath = 65;
						}
						else
						{
							//PCE alleles are maternal, PC1 alleles are either
							childGeno.iu_anyAllelesUnknown = true;
							childGeno.iu_UPHD_Candidate_p2 = true;
							childGeno.iu_UPHD_Candidate_p1 = true;
							if (pcnv > 1)
							{
								//Dad has more than one X
								//Not terribly likely, but definitely possible.
								List<Integer> idup1 = possibleInheritedDuplications(childimap, true, childGeno.getCopyNumber(), parentGeno.getCopyNumber());
								List<Integer> upid1 = UPID_possible(childimap, true);
								List<Allele> p1Flagged = childimap.getUniquePhasedP1();
								//List<Allele> p2Flagged = childimap.getUniquePhasedP2();
								List<Allele> unFlagged = childimap.getUniqueUnphased();
								
								if (idup1 != null) childGeno.iu_inheritedDUP_candidate_p1 = InheritanceMap.possible_idup(idup1, p1Flagged, unFlagged, true);
								if (upid1 != null) childGeno.iu_UPID_Candidate_p1 = InheritanceMap.possible_upid(upid1, p1Flagged, unFlagged, true);
							}
							childGeno.phasePath = 66;
						}
					}
					else
					{
						//All could be from either parent. Assume at least one is maternal.
						if (PC1.size() == 1)
						{
							Allele a = childimap.getAllele(PC1.get(0));
							a.setFlag2Count(1);
						}
						int p2p = childimap.countPhasedP2();
						//int pup = childimap.countUnphased();
						if (p2p > 0) childGeno.iu_anyAllelesUnknown = true;
						else childGeno.iu_allAllelesUnknown = true;
						childGeno.iu_UPHD_Candidate_p1 = true; //XY
						childGeno.iu_UPHD_Candidate_p2 = true;
						
						if (pcnv > 1)
						{
							//Dad has more than one X
							List<Integer> idup1 = possibleInheritedDuplications(childimap, true, childGeno.getCopyNumber(), parentGeno.getCopyNumber());
							List<Integer> upid1 = UPID_possible(childimap, true);
							List<Allele> p1Flagged = childimap.getUniquePhasedP1();
							List<Allele> unFlagged = childimap.getUniqueUnphased();
							
							if (idup1 != null) childGeno.iu_inheritedDUP_candidate_p1 = InheritanceMap.possible_idup(idup1, p1Flagged, unFlagged, true);
							if (upid1 != null) childGeno.iu_UPID_Candidate_p1 = InheritanceMap.possible_upid(upid1, p1Flagged, unFlagged, true);
						}
						
						childGeno.phasePath = 67;
					}
					
					
				}
			}
		}
		else
		{
			//Female or unknown sex. Treat like autosome, but remove iDEL flags from
			//paternal side.
			phase_PC_autosomal(parentGeno, childGeno);
			if (psx == Sex.FEMALE || (punsx && pcnv > 1))
			{
				childGeno.iu_inheritedDEL_candidate_p2 = false;
			}
			else if (psx == Sex.MALE || (punsx && pcnv < 2))
			{
				childGeno.iu_inheritedDEL_candidate_p1 = false;
			}
			//Don't change phase path!
		}
		
		childimap.flagLinkedAlleles();
		
	}
	
	public static void phase_PC_sxc_mammalY(SNPGeno parentGeno, SNPGeno childGeno, Sex psx, Sex csx)
	{
		childGeno.phasingAttempted = true;
		int ccnv = childGeno.getCopyNumber();
		int pcnv = parentGeno.getCopyNumber();
		boolean punsx = (psx == Sex.OTHER || psx == Sex.UNKNOWN);
		
		if (csx == Sex.FEMALE || ccnv < 1)
		{
			//Expect CNV0 !!!!
			//If not, then something wonky is afoot.
			//If CNV0, nothing happens. Just mark phase path.
			if (ccnv == 0) childGeno.phasePath = 68;
			if (ccnv > 0)
			{
				//Check parental CNV
				childGeno.iu_unexpectedSegregation = true;
				InheritanceMap childimap = childGeno.buildInheritanceMap(parentGeno, null);
				List<Integer> PCE = get_uninheritableAlleles_PC(childimap);
				List<Integer> PC1 = canInherit_PC(childimap, true);
				if (psx == Sex.FEMALE || (punsx && pcnv < 1))
				{
					childGeno.iu_UPHD_Candidate_p2 = true; //XY from dad
					if (parentGeno.getCopyNumber() > 0)
					{
						//Just... treat as autosome... Remove any maternal deletion flags
						phase_PC_autosomal(parentGeno, childGeno);
						childGeno.iu_inheritedDEL_candidate_p1 = false;
					}
					else
					{
						//Phase all as missing dad.
						for (Integer i : PCE)
						{
							Allele as = childimap.getAllele(i);
							as.setFlag2Count(as.getAlleleCount());
							if (as.getAlleleCount() > 1)
							{
								childGeno.iu_inheritedDUP_candidate_p2 = true;
								childGeno.iu_UPID_Candidate_p2 = true;
							}
						}
						//Check for dad duplications...	
						childGeno.iu_unknownGenoAssumed = true;
						childGeno.phasePath = 69;
					}
				}
				else if (psx == Sex.MALE || (punsx && pcnv > 0))
				{
					//Known parent is dad
					//Just assume all possible came from dad.
					childGeno.iu_UPHD_Candidate_p1 = true; //XY from dad
					for (Integer i : PC1)
					{
						Allele as = childimap.getAllele(i);
						as.setFlag1Count(as.getAlleleCount());
						if (as.getAlleleCount() > 1)
						{
							childGeno.iu_inheritedDUP_candidate_p1 = true;
							childGeno.iu_UPID_Candidate_p1 = true;
						}
					}
					
					if (PCE != null)
					{
						for (Integer i : PCE)
						{
							Allele as = childimap.getAllele(i);
							as.setFlag2Count(as.getAlleleCount());
							if (as.getAlleleCount() > 1)
							{
								childGeno.iu_inheritedDUP_candidate_p2 = true;
								childGeno.iu_UPID_Candidate_p2 = true;
							}
						}
						if (PCE.size() > 1) childGeno.iu_UPHD_Candidate_p2 = true;						
					}

					childGeno.phasePath = 70;
				}
				
				childimap.flagLinkedAlleles();
			}
		}
		else
		{
			//Expect CNV1
			//Any more or less and that needs to be looked into.
			InheritanceMap childimap = childGeno.buildInheritanceMap(parentGeno, null);
			if (ccnv == 1)
			{
				if (psx == Sex.FEMALE || (punsx && pcnv < 1))
				{
					//From P2. End of story.
					List<Allele> up = childimap.getUniqueUnphased();
					for (Allele a : up) a.setFlag2Count(a.getAlleleCount());
					childGeno.phasePath = 71;
				}
				else if (psx == Sex.MALE || (punsx && pcnv > 0))
				{
					//Check if from P1.
					//If not, flag as PC error. Don't bother flagging it P2.
					List<Integer> PC1 = canInherit_PC(childimap, true);
					if (PC1 != null)
					{
						for (Integer i : PC1)
						{
							Allele as = childimap.getAllele(i);
							as.setFlag1Count(as.getAlleleCount());
						}
					}
					else
					{
						childGeno.iu_PC_err = true;
						childGeno.iu_unexpectedSegregation = true;
						childGeno.iu_allAllelesUnknown = true;
					}
					childGeno.phasePath = 72;
				}
			}
			else if (ccnv == 0)
			{
				//Assume paternal deletion (either idel or denovo del)
				//Throw PC error if father is P1 and no idel.
				if (psx == Sex.FEMALE || (punsx && pcnv < 1))
				{
					childGeno.iu_inheritedDEL_candidate_p2 = true;
					childGeno.iu_denovoCandidate = true;
					childGeno.iu_unexpectedSegregation = true;
					childGeno.phasePath = 73;
				}
				else if (psx == Sex.MALE || (punsx && pcnv > 0))
				{
					childGeno.iu_inheritedDEL_candidate_p1 = (pcnv < 1);
					if (!childGeno.iu_inheritedDEL_candidate_p1)
					{
						childGeno.iu_denovoCandidate = true;
						childGeno.iu_unexpectedSegregation = true;
						childGeno.iu_PC_err = true;
					}
					childGeno.phasePath = 74;
				}
			}
			else
			{
				//Cry
				childGeno.iu_unexpectedSegregation = true;
				List<Integer> PCE = get_uninheritableAlleles_PC(childimap);
				List<Integer> PC1 = canInherit_PC(childimap, true);
				if (psx == Sex.FEMALE || (punsx && pcnv < 1))
				{
					if (pcnv > 0)
					{
						//Just go autosomal... Remove any maternal deletion flags
						phase_PC_autosomal(parentGeno, childGeno);
						childGeno.iu_inheritedDEL_candidate_p1 = false;
						childGeno.phasePath = 75;
					}
					else
					{
						if (PCE != null)
						{
							for (Integer i : PCE)
							{
								Allele as = childimap.getAllele(i);
								as.setFlag2Count(as.getAlleleCount());
								if (as.getAlleleCount() > 1)
								{
									childGeno.iu_inheritedDUP_candidate_p2 = true;
									childGeno.iu_UPID_Candidate_p2 = true;
								}
							}
						}
						childGeno.phasePath = 76;
					}
				}
				else if (psx == Sex.MALE || (punsx && pcnv > 0))
				{
					if (PCE != null)
					{
						//There are some alleles that did not come from dad...
						childGeno.iu_PC_err = true;
						//Leave these unphased.
						if (PC1 != null)
						{
							for (Integer i : PC1)
							{
								Allele as = childimap.getAllele(i);
								as.setFlag1Count(as.getAlleleCount());
								if (as.getAlleleCount() > 1)
								{
									childGeno.iu_inheritedDUP_candidate_p1 = true;
									childGeno.iu_UPID_Candidate_p1 = true;
								}
							}
						}
						int p2p = childimap.countPhasedP2();
						if (pcnv > 1 && p2p > 1) childGeno.iu_UPHD_Candidate_p2 = true;
						
						childGeno.phasePath = 77;
					}
					else
					{
						//All alleles could have come from dad. Assume they did.
						for (Integer i : PC1)
						{
							Allele as = childimap.getAllele(i);
							as.setFlag1Count(as.getAlleleCount());
							if (as.getAlleleCount() > 1)
							{
								childGeno.iu_inheritedDUP_candidate_p1 = true;
								childGeno.iu_UPID_Candidate_p1 = true;
							}
						}
						int p2p = childimap.countPhasedP2();
						if (pcnv > 1 && p2p > 1) childGeno.iu_UPHD_Candidate_p2 = true;
						childGeno.phasePath = 78;
					}
				}
			}
			
			childimap.flagLinkedAlleles();
		}
		
	}
	
	public static void phase_PPC_sxc_mammalX(SNPGeno p1Geno, SNPGeno p2Geno, SNPGeno childGeno, Sex p1sx, Sex p2sx, Sex csx)
	{
		childGeno.phasingAttempted = true;
		//InheritanceMap childimap = childGeno.buildInheritanceMap(p1Geno, p2Geno);
		int ccnv = childGeno.getCopyNumber();
		int p1cnv = p1Geno.getCopyNumber();
		//boolean p1unsx = (p1sx == Sex.OTHER || p1sx == Sex.UNKNOWN);
		int p2cnv = p2Geno.getCopyNumber();
		//boolean p2unsx = (p2sx == Sex.OTHER || p2sx == Sex.UNKNOWN);
		
		if (csx == Sex.FEMALE || ccnv > 1)
		{
			//First, phase as if autosomal...
			phase_PPC_autosomal(p1Geno, p2Geno, childGeno);
			//unmark paternal deletions if father isn't CNV0
			if (p1sx == Sex.MALE && p1cnv > 0)
			{
				childGeno.iu_inheritedDEL_candidate_p1 = false;
			}
			if (p2sx == Sex.MALE && p2cnv > 0)
			{
				childGeno.iu_inheritedDEL_candidate_p2 = false;
			}
			
		}
		else if (csx == Sex.MALE || ccnv < 2)
		{
			phase_PPC_autosomal(p1Geno, p2Geno, childGeno);
			if (p1sx == Sex.MALE)
			{
				childGeno.iu_inheritedDEL_candidate_p1 = false;
			}
			if (p2sx == Sex.MALE)
			{
				childGeno.iu_inheritedDEL_candidate_p2 = false;
			}
			//Check to see if there is one allele and it is unphased.
			if (ccnv == 1)
			{
				List<SNPAllele> call = childGeno.getAlleles();
				SNPAllele ca = call.get(0);
				if (p1sx == Sex.MALE && p2sx == Sex.FEMALE)
				{
					if(canInherit_PC(p2Geno, childGeno)) ca.flagParent2(true);
				}
				if (p2sx == Sex.MALE && p1sx == Sex.FEMALE)
				{
					if(canInherit_PC(p1Geno, childGeno)) ca.flagParent1(true);
				}
			}
		}
		else
		{
			phase_PPC_autosomal(p1Geno, p2Geno, childGeno);
		}
		
	}
	
	public static void phase_PPC_sxc_mammalY(SNPGeno p1Geno, SNPGeno p2Geno, SNPGeno childGeno, Sex p1sx, Sex p2sx, Sex csx)
	{
		childGeno.phasingAttempted = true;
		int ccnv = childGeno.getCopyNumber();
		int p1cnv = p1Geno.getCopyNumber();
		int p2cnv = p2Geno.getCopyNumber();
		
		if (csx == Sex.FEMALE || ccnv < 1)
		{
			if (ccnv > 0)
			{
				childGeno.iu_unexpectedSegregation = true;
				if (ccnv == 1)
				{
					List<SNPAllele> call = childGeno.getAlleles();
					SNPAllele ca = call.get(0);
					if (p1sx == Sex.MALE && p1cnv > 0 && p2cnv < 1)
					{
						if(canInherit_PC(p1Geno, childGeno)) ca.flagParent1(true);
						else phase_PPC_autosomal(p1Geno, p2Geno, childGeno);
					}
					if (p2sx == Sex.MALE && p2cnv > 0 && p1cnv < 1)
					{
						if(canInherit_PC(p2Geno, childGeno)) ca.flagParent2(true);
						else phase_PPC_autosomal(p1Geno, p2Geno, childGeno);
					}	
				}
				else
				{
					phase_PPC_autosomal(p1Geno, p2Geno, childGeno);
				}
				if (p1sx == Sex.FEMALE && p1cnv < 1)
				{
					childGeno.iu_inheritedDEL_candidate_p1 = false;
				}
				if (p2sx == Sex.FEMALE && p2cnv < 1)
				{
					childGeno.iu_inheritedDEL_candidate_p2 = false;
				}	
			}
			
		}
		else
		{
			//Expect ccnv 1
			if (ccnv == 1)
			{
				//See if could be inherited from father
				//If not, throw a PPC error and do not phase
				List<SNPAllele> call = childGeno.getAlleles();
				SNPAllele ca = call.get(0);
				if (p1sx == Sex.MALE && p1cnv > 0 && p2cnv < 1)
				{
					if(canInherit_PC(p1Geno, childGeno)) ca.flagParent1(true);
					else
					{
						childGeno.iu_PPC_err = true;
						childGeno.iu_unexpectedSegregation = true;
					}
				}
				if (p2sx == Sex.MALE && p2cnv > 0 && p1cnv < 1)
				{
					if(canInherit_PC(p2Geno, childGeno)) ca.flagParent2(true);
					else
					{
						childGeno.iu_PPC_err = true;
						childGeno.iu_unexpectedSegregation = true;
					}
				}
			}
			else
			{
				//Phase autosomal.
				//Remove maternal idel flags. Remove paternal idel flags if CNV>0
				phase_PPC_autosomal(p1Geno, p2Geno, childGeno);
				if (p1sx == Sex.FEMALE)
				{
					childGeno.iu_inheritedDEL_candidate_p1 = false;
				}
				if (p2sx == Sex.FEMALE)
				{
					childGeno.iu_inheritedDEL_candidate_p2 = false;
				}
				if (p1sx == Sex.MALE && p2cnv > 0)
				{
					childGeno.iu_inheritedDEL_candidate_p1 = false;
				}
				if (p2sx == Sex.MALE && p2cnv > 0)
				{
					childGeno.iu_inheritedDEL_candidate_p2 = false;
				}
			}
		}
		
		
		
	}
	
	public static void phase_PC_mtc(SNPGeno parentGeno, SNPGeno childGeno, Sex psx)
	{
		childGeno.phasingAttempted = true;
		//Just... phase all maternal. 
		//If can't inherit from mom, then leave unflagged.
		
		if (psx == Sex.FEMALE)
		{
			List<SNPAllele> call = childGeno.getAlleles();
			for (SNPAllele a : call)
			{
				if (parentGeno.countAllele(a) < 1) a.flagParent1(true);
				else
				{
					childGeno.iu_unexpectedSegregation = true;
					childGeno.iu_PC_err = true;
				}
			}
			
			childGeno.phasePath = 90;
		}
		else if (psx == Sex.MALE)
		{
			List<SNPAllele> call = childGeno.getAlleles();
			for (SNPAllele a : call)
			{
				a.flagParent2(true);
			}
			childGeno.phasePath = 91;
		}
		else
		{
			//Leave everything unphased.
			childGeno.phasePath = 92;
		}
		
	}
	
	public static void phase_PPC_mtc(SNPGeno p1Geno, SNPGeno p2Geno, SNPGeno childGeno, Sex p1sx, Sex p2sx)
	{
		childGeno.phasingAttempted = true;
		//Just... phase all maternal. 
		//If can't inherit from mom, then leave unflagged.
		if (p1sx == Sex.FEMALE)
		{
			if (p2sx == Sex.FEMALE)
			{
				//Phase autosomal...
				phase_PPC_autosomal(p1Geno, p2Geno, childGeno);
				//Remove del dup flags
				childGeno.iu_denovoCandidate = false;
				childGeno.iu_inheritedDEL_candidate_p1 = false;
				childGeno.iu_inheritedDEL_candidate_p2 = false;
				childGeno.iu_inheritedDUP_candidate_p1 = false;
				childGeno.iu_inheritedDUP_candidate_p2 = false;
				childGeno.iu_unexpectedSegregation = false;
				childGeno.iu_UPHD_Candidate_p1 = false;
				childGeno.iu_UPHD_Candidate_p2 = false;
				childGeno.iu_UPID_Candidate_p1 = false;
				childGeno.iu_UPID_Candidate_p2 = false;
				childGeno.phasePath = 93;
			}
			else
			{
				//Flag all matching as P1
				List<SNPAllele> call = childGeno.getAlleles();
				for (SNPAllele a : call)
				{
					if (p1Geno.countAllele(a) < 1) a.flagParent1(true);
					else
					{
						childGeno.iu_unexpectedSegregation = true;
						childGeno.iu_PPC_err = true;
					}
				}
				childGeno.phasePath = 94;
			}
		}
		else if (p1sx == Sex.MALE)
		{
			if (p2sx == Sex.FEMALE)
			{
				//Flag all P2
				List<SNPAllele> call = childGeno.getAlleles();
				for (SNPAllele a : call)
				{
					if (p2Geno.countAllele(a) < 1) a.flagParent2(true);
					else
					{
						childGeno.iu_unexpectedSegregation = true;
						childGeno.iu_PPC_err = true;
					}
				}
				childGeno.phasePath = 95;
			}
			else if (p2sx == Sex.MALE)
			{
				//Phase autosomal...
				phase_PPC_autosomal(p1Geno, p2Geno, childGeno);
				//Remove del dup flags
				childGeno.iu_denovoCandidate = false;
				childGeno.iu_inheritedDEL_candidate_p1 = false;
				childGeno.iu_inheritedDEL_candidate_p2 = false;
				childGeno.iu_inheritedDUP_candidate_p1 = false;
				childGeno.iu_inheritedDUP_candidate_p2 = false;
				childGeno.iu_unexpectedSegregation = false;
				childGeno.iu_UPHD_Candidate_p1 = false;
				childGeno.iu_UPHD_Candidate_p2 = false;
				childGeno.iu_UPID_Candidate_p1 = false;
				childGeno.iu_UPID_Candidate_p2 = false;
				childGeno.phasePath = 96;
			}
			else
			{
				//Flag all P2
				List<SNPAllele> call = childGeno.getAlleles();
				for (SNPAllele a : call)
				{
					if (p2Geno.countAllele(a) < 1) a.flagParent2(true);
					else
					{
						childGeno.iu_unexpectedSegregation = true;
						childGeno.iu_PPC_err = true;
					}
				}
				childGeno.phasePath = 97;
			}
		}
		else
		{
			if (p2sx == Sex.FEMALE)
			{
				//Phase autosomal...
				phase_PPC_autosomal(p1Geno, p2Geno, childGeno);
				//Remove del dup flags
				childGeno.iu_denovoCandidate = false;
				childGeno.iu_inheritedDEL_candidate_p1 = false;
				childGeno.iu_inheritedDEL_candidate_p2 = false;
				childGeno.iu_inheritedDUP_candidate_p1 = false;
				childGeno.iu_inheritedDUP_candidate_p2 = false;
				childGeno.iu_unexpectedSegregation = false;
				childGeno.iu_UPHD_Candidate_p1 = false;
				childGeno.iu_UPHD_Candidate_p2 = false;
				childGeno.iu_UPID_Candidate_p1 = false;
				childGeno.iu_UPID_Candidate_p2 = false;
				childGeno.phasePath = 98;
			}
			else if (p2sx == Sex.MALE)
			{
				//Flag all matching as P1
				List<SNPAllele> call = childGeno.getAlleles();
				for (SNPAllele a : call)
				{
					if (p1Geno.countAllele(a) < 1) a.flagParent1(true);
					else
					{
						childGeno.iu_unexpectedSegregation = true;
						childGeno.iu_PPC_err = true;
					}
				}
				childGeno.phasePath = 99;
			}
			else
			{
				//Phase autosomal...
				phase_PPC_autosomal(p1Geno, p2Geno, childGeno);
				//Remove del dup flags
				childGeno.iu_denovoCandidate = false;
				childGeno.iu_inheritedDEL_candidate_p1 = false;
				childGeno.iu_inheritedDEL_candidate_p2 = false;
				childGeno.iu_inheritedDUP_candidate_p1 = false;
				childGeno.iu_inheritedDUP_candidate_p2 = false;
				childGeno.iu_unexpectedSegregation = false;
				childGeno.iu_UPHD_Candidate_p1 = false;
				childGeno.iu_UPHD_Candidate_p2 = false;
				childGeno.iu_UPID_Candidate_p1 = false;
				childGeno.iu_UPID_Candidate_p2 = false;
				childGeno.phasePath = 100;
			}
		}
		
	}
	
		// -- Phasing Methods (General)
	
	public static void phase_PC_autosomal(SNPGeno parentGeno, SNPGeno childGeno)
	{
		int ccnv = childGeno.getCopyNumber();
		if (ccnv == 0) phase_PC_CNV0_auto(parentGeno, childGeno);
		else if (ccnv == 1) phase_PC_CNV1_auto(parentGeno, childGeno);
		else if (ccnv == 2) phase_PC_CNV2_auto(parentGeno, childGeno);
		else phase_PC_CNV3Plus_auto(parentGeno, childGeno);
	}
	
	public static void phase_PPC_autosomal(SNPGeno p1Geno, SNPGeno p2Geno, SNPGeno childGeno)
	{
		int ccnv = childGeno.getCopyNumber();
		if (ccnv == 0) phase_PPC_CNV0_auto(p1Geno, p2Geno, childGeno);
		if (ccnv == 1) phase_PPC_CNV1_auto(p1Geno, p2Geno, childGeno);
		else phase_PPC_CNV2Plus_auto(p1Geno, p2Geno, childGeno);
	}
	
	
	
}
