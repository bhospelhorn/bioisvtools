package hospelhornbg_segregation;

public class CandidateFlags {
	
	/**
	 * Manual examination  of the aligned NGS reads suggests that
	 * this variant is real.
	 */
	private boolean bBAMEv;
	
	/**
	 * There is lab evidence indicating that this variant is real.
	 */
	private boolean bLabEv;
	
	/**
	 * Variant has been verified to be real and genotypes
	 * are consistent with findings.
	 */
	private boolean bVerified;
	
	/**
	 * A field that can be used to flag the variant/candidate as
	 * interesting.
	 */
	private boolean bInteresting;
	
	/**
	 * There is potential mosaicism in one or more individuals
	 * at the locus covered by this variant.
	 */
	private boolean bMosaic;
	
	/**
	 * Candidate is a halfhet that could not be paired.
	 */
	private boolean bUnpairedHH;
	
	/**
	 * Variant and/or variant breakpoints known to fall in a poorly
	 * aligned region.
	 */
	private boolean bRPA;
	
	/**
	 * Transcript associated with this candidate is known to be
	 * difficult to align properly.
	 */
	private boolean bRPAGene;
	
	/**
	 * There is little to no evidence of this variant being real
	 * (ie. erroneous call)
	 */
	private boolean bNoEv;
	
	/**
	 * Variant is known to be common.
	 */
	private boolean bCommonVar;
	
	/**
	 * Variant is CN and covers at least 5 SNPs on a SNP chip
	 * and should be visible to PennCNV minsnp5.
	 */
	private boolean bSNPVis;
	
	/**
	 * Locus is highly variable or has a lot of similar loci.
	 * This locus is tolerant to variation and/or is subject
	 * to alignment issues.
	 * Examples: HLA, MUC, TTN
	 */
	private boolean bVariableLocus;
	
	/**
	 * Variant size falls below a specified size threshold
	 * (Default: 50bp)
	 */
	private boolean bSmall;
	
	/**
	 * Variant size falls above a specified size threshold
	 * (Default: 5Mbp)
	 */
	private boolean bLarge;
	
	/**
	 * Candidate is in a known pseudogene.
	 */
	private boolean bPseudogene;
	
	public CandidateFlags()
	{
		this.bBAMEv = false;
		this.bLabEv = false;
		this.bVerified = false;
		this.bInteresting = false;
		this.bMosaic = false;
		this.bUnpairedHH = false;
		this.bRPA = false;
		this.bRPAGene = false;
		this.bNoEv = false;
		this.bCommonVar = false;
		this.bSNPVis = false;
		this.bVariableLocus = false;
		this.bSmall = false;
		this.bLarge = false;
		this.bPseudogene = false;
	}
	
	public static CandidateFlags readVarTableField(int value)
	{
		CandidateFlags flags = new CandidateFlags();
		flags.bBAMEv = (value & 0x00000001) != 0;
		flags.bLabEv = (value & 0x00000002) != 0;
		flags.bVerified = (value & 0x00000004) != 0;
		flags.bInteresting = (value & 0x00000008) != 0;
		flags.bMosaic = (value & 0x00000010) != 0;
		flags.bUnpairedHH = (value & 0x00000020) != 0;
		flags.bRPA = (value & 0x00000040) != 0;
		flags.bRPAGene = (value & 0x00000080) != 0;
		flags.bNoEv = (value & 0x00000100) != 0;
		flags.bCommonVar = (value & 0x00000200) != 0;
		flags.bSNPVis = (value & 0x00000400) != 0;
		flags.bVariableLocus = (value & 0x00000800) != 0;
		flags.bSmall = (value & 0x00001000) != 0;
		flags.bLarge = (value & 0x00002000) != 0;
		flags.bPseudogene = (value & 0x00004000) != 0;
		return flags;
	}
	
	public static CandidateFlags readVCFBitField(String value)
	{
		//Backwards from bit field
		if (value.length() < 15) return null;
		CandidateFlags flags = new CandidateFlags();
		flags.bBAMEv = (value.charAt(0)) == '1';
		flags.bLabEv = (value.charAt(1)) == '1';
		flags.bVerified = (value.charAt(2)) == '1';
		flags.bInteresting = (value.charAt(3)) == '1';
		flags.bMosaic = (value.charAt(4)) == '1';
		flags.bUnpairedHH = (value.charAt(5)) == '1';
		flags.bRPA = (value.charAt(6)) == '1';
		flags.bRPAGene = (value.charAt(7)) == '1';
		flags.bNoEv = (value.charAt(8)) == '1';
		flags.bCommonVar = (value.charAt(9)) == '1';
		flags.bSNPVis = (value.charAt(10)) == '1';
		flags.bVariableLocus = (value.charAt(11)) == '1';
		flags.bSmall = (value.charAt(12)) == '1';
		flags.bLarge = (value.charAt(13)) == '1';
		flags.bPseudogene = (value.charAt(14)) == '1';
		return flags;
	}
	
	public int getVarTableField()
	{
		int value = 0;
		if(this.bBAMEv) value |= 0x00000001;
		if(this.bLabEv) value |= 0x00000002;
		if(this.bVerified) value |= 0x00000004;
		if(this.bInteresting) value |= 0x00000008;
		if(this.bMosaic) value |= 0x00000010;
		if(this.bUnpairedHH) value |= 0x00000020;
		if(this.bRPA) value |= 0x00000040;
		if(this.bRPAGene) value |= 0x00000080;
		if(this.bNoEv) value |= 0x00000100;
		if(this.bCommonVar) value |= 0x00000200;
		if(this.bSNPVis) value |= 0x00000400;
		if(this.bVariableLocus) value |= 0x00000800;
		if(this.bSmall) value |= 0x00001000;
		if(this.bLarge) value |= 0x00002000;
		if(this.bPseudogene) value |= 0x00004000;
		
		return value;
	}
	
	public String getVCFBitField()
	{
		StringBuilder sb = new StringBuilder(32);
		if(this.bBAMEv) sb.append("1");
		else sb.append("0");
		if(this.bLabEv) sb.append("1");
		else sb.append("0");
		if(this.bVerified) sb.append("1");
		else sb.append("0");
		if(this.bInteresting) sb.append("1");
		else sb.append("0");
		if(this.bMosaic) sb.append("1");
		else sb.append("0");
		if(this.bUnpairedHH) sb.append("1");
		else sb.append("0");
		if(this.bRPA) sb.append("1");
		else sb.append("0");
		if(this.bRPAGene) sb.append("1");
		else sb.append("0");
		if(this.bNoEv) sb.append("1");
		else sb.append("0");
		if(this.bCommonVar) sb.append("1");
		else sb.append("0");
		if(this.bSNPVis) sb.append("1");
		else sb.append("0");
		if(this.bVariableLocus) sb.append("1");
		else sb.append("0");
		if(this.bSmall) sb.append("1");
		else sb.append("0");
		if(this.bLarge) sb.append("1");
		else sb.append("0");
		if(this.bPseudogene) sb.append("1");
		else sb.append("0");
		
		return sb.toString();
	}
	
	public boolean hasBAMEvidence()
	{
		return this.bBAMEv;
	}
	
	public void setBAMEvidence(boolean b)
	{
		this.bBAMEv = b;
	}

	public boolean hasLabEvidence()
	{
		return this.bLabEv;
	}
	
	public void setLabEvidence(boolean b)
	{
		this.bLabEv = b;
	}
	
	public boolean isVerified()
	{
		return this.bVerified;
	}
	
	public void setVerified(boolean b)
	{
		this.bVerified = b;
	}
	
	public boolean flaggedInteresting()
	{
		return this.bInteresting;
	}
	
	public void setInteresting(boolean b)
	{
		this.bInteresting = b;
	}

	public boolean flaggedPotentialMosaic()
	{
		return this.bMosaic;
	}
	
	public void setMosaic(boolean b)
	{
		this.bMosaic = b;
	}
	
	public boolean isUnpairedHalfHet()
	{
		return this.bUnpairedHH;
	}
	
	public void flagUnpairedHalfHet(boolean b)
	{
		this.bUnpairedHH = b;
	}
	
	public boolean isPoorlyAligned()
	{
		return this.bRPA;
	}
	
	public void flagPoorlyAligned(boolean b)
	{
		this.bRPA = b;
	}
	
	public boolean inDifficultToAlignGene()
	{
		return this.bRPAGene;
	}
	
	public void flagInDifficultToAlignGene(boolean b)
	{
		this.bRPAGene = b;
	}
	
	public boolean hasNoEvidence()
	{
		return this.bNoEv;
	}
	
	public void flagNoEvidence(boolean b)
	{
		this.bNoEv = b;
	}
	
	public boolean isCommonVariant()
	{
		return this.bCommonVar;
	}
	
	public void flagCommonVariant(boolean b)
	{
		this.bCommonVar = b;
	}
	
	public boolean isMINSNPVisible()
	{
		return this.bSNPVis;
	}
	
	public void flagMINSNPVisible(boolean b)
	{
		this.bSNPVis = b;
	}
	
	public boolean isOnVariableLocus()
	{
		return this.bVariableLocus;
	}
	
	public void flagVariableLocus(boolean b)
	{
		this.bVariableLocus = b;
	}
	
	public boolean flaggedSmall()
	{
		return this.bSmall;
	}
	
	public void flagSmall(boolean b)
	{
		this.bSmall = b;
	}
	
	public boolean flaggedLarge()
	{
		return this.bLarge;
	}
	
	public void flagLarge(boolean b)
	{
		this.bLarge = b;
	}
	
}
