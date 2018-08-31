package hospelhornbg_genomeBuild;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import waffleoRai_Utils.FileBuffer;

/*
 * UPDATE LOG
 * 	Version 1.0.0 created April 2, 2018
 * 
 * 1.0.0 -> 1.0.1 | April 5, 2018
 * 		Created printInfo() method
 * 		Initial debug
 *
 * 1.0.1 -> 1.1.0 | April 19, 2018
 * 		Debugging
 * 
 * 1.1.0 -> 1.2.0 | August 31, 2018
 * 		Added methods to get position effect of region
 */

/**
 * A container for basic gene information - location and exons.
 * @author Blythe Hospelhorn
 * @version 1.2.0
 * @since August 31, 2018
 *
 */
public class Gene implements Comparable<Gene>{
	
	/**
	 * Distance in basepairs from exon end that is considered to be
	 * within a potential splicing region.
	 */
	public static final int SPLICE_DIST = 10;
	
	/* --- Instance Variables --- */
	
	private String name;
	private Contig chrom;
	
	private boolean strand;
	private boolean isNCRNA;
	private String ID;
	
	private int stPos; //Transcript start
	private int edPos; //Transcript end

	private int tlStart; //Start of translated region
	private int tlEnd; //End of translated region
	
	private List<Exon> exons;
	
	/* --- Inner Structures --- */
	
	/**
	 * Container for exon information.
	 * @author Blythe Hospelhorn
	 * @version 1.0.0
	 * @since April 2, 2018
	 */
	public static class Exon implements Comparable<Exon>
	{
		private int stPos;
		private int edPos;
		
		/**
		 * Construct a new exon, specifying its start and end positions
		 * relative to the contig (NOT the gene!)
		 * @param start Exon start position.
		 * @param end Exon end position.
		 */
		public Exon(int start, int end)
		{
			stPos = start;
			edPos = end;
		}
		
		/**
		 * Get the start position of the exon.
		 * @return Start position of the exon relative to the contig it lies on.
		 */
		public int getStart()
		{
			return stPos;
		}
		
		/**
		 * Get the end position of the exon.
		 * @return End position of the exon relative to the contig it lies on.
		 */
		public int getEnd()
		{
			return edPos;
		}
		
		public boolean equals(Object o)
		{
			if (o == null) return false;
			if (o == this) return true;
			if (!(o instanceof Exon)) return false;
			
			Exon e = (Exon)o;
			
			if (this.stPos != e.stPos) return false;
			if (this.edPos != e.edPos) return false;
			
			return true;
		}
		
		public int hashCode()
		{
			return Integer.hashCode(stPos) ^ Integer.hashCode(edPos);
		}
		
		public int compareTo(Exon o) 
		{
			if (o == null) return 1;
			
			int comp = this.stPos - o.stPos;
			if (comp != 0) return comp;
			
			comp = this.edPos - o.edPos;
			if (comp != 0) return comp;
			
			return 0;
		}
		
		/**
		 * Determine whether a position (relative to the same contig
		 * as this gene) falls within this exon.
		 * @param position Position relative to contig.
		 * @return True - If the position falls inside this exon, assuming it
		 * is on the same contig.
		 * <br>False - If the position falls outside this exon, regardless of contig.
		 */
		public boolean positionInExon(int position)
		{
			return (position < edPos) && (position >= stPos);
		}
		
		/**
		 * Determine whether a position (relative to the same contig as this gene)
		 * potentially falls within the splicing zone of this exon.
		 * @param position Position relative to contig.
		 * @param first Whether this exon is to be considered the first in the gene (no left flanking splicing zone)
		 * @param last Whether this exon is to be considered the last in the gene (no right flanking splicing zone)
		 * @return True - If the position falls inside this exon's potential splicing zone, assuming it
		 * is on the same contig.
		 * <br>False - If the position falls outside this exon's potential splicing zone, regardless of contig.
		 */
		public boolean positionInSpliceZone(int position, boolean first, boolean last)
		{
			if (!first && !last) return (position < (edPos + SPLICE_DIST)) && (position >= (stPos - SPLICE_DIST));
			if (first) return (position < (edPos + SPLICE_DIST)) && (position >= stPos);
			if (last) return (position < edPos) && (position >= (stPos - SPLICE_DIST));
			return false;
		}
		
		/**
		 * Determine whether a position (relative to the same contig as this gene)
		 * falls before this exon.
		 * @param position Position relative to contig.
		 * @return True - If the position falls before this exon, assuming it
		 * is on the same contig.
		 * <br>False - If the position does not fall before this exon, regardless of contig.
		 */
		public boolean isBefore(int position)
		{
			return position < stPos;
		}
		
		/**
		 * Determine whether a position (relative to the same contig as this gene)
		 * falls after this exon.
		 * @param position Position relative to contig.
		 * @return True - If the position falls after this exon, assuming it
		 * is on the same contig.
		 * <br>False - If the position does not fall after this exon, regardless of contig.
		 */
		public boolean isAfter(int position)
		{
			return position >= edPos;
		}
		
	}
	
	/* --- Construction --- */
	
	/**
	 * Construct an empty Gene container with no exons, and all instance variables set to their
	 * "unset" values.
	 */
	public Gene()
	{
		exons = new ArrayList<Exon>(); 
		setDefaults();
	}
	
	/**
	 * Construct an empty Gene with memory allocated for a known number of exons.
	 * All instance variables will be set to their "unset" values.
	 * @param numberExons Number of exons expected for this gene.
	 */
	public Gene(int numberExons)
	{
		exons = new ArrayList<Exon>(numberExons);
		setDefaults();
	}
	
	private void setDefaults()
	{
		name = null;
		chrom = null;
		strand = false;
		ID = null;
		stPos = -1;
		edPos = -1;
		tlStart = -1;
		tlEnd = -1;
		isNCRNA = false;
		exons.clear();
	}

	/* --- Getters --- */
	
	/**
	 * Get the primary name recorded for this gene.
	 * @return Gene name
	 */
	public String getName()
	{
		return name;
	}
	
	/**
	 * Get the ID string (UID) for this gene as recorded in its source database.
	 * @return Gene UID string.
	 */
	public String getID()
	{
		return ID;
	}
	
	/**
	 * Get the strand the gene is on.
	 * @return True - Plus(+) strand.
	 * <br>False - Minus(-) strand.
	 */
	public boolean getStrand()
	{
		return strand;
	}
	
	/**
	 * Get whether the gene transcript is non-coding RNA.
	 * @return True - If gene is for a non-coding RNA.
	 * <br>False - If the gene transcript is mRNA (protein coding gene).
	 */
	public boolean is_ncRNA()
	{
		return isNCRNA;
	}
	
	/**
	 * Get the chromosome/contig this gene lies on.
	 * @return Contig referenced by this gene.
	 */
	public Contig getChromosome()
	{
		return chrom;
	}
	
	/**
	 * Get the starting position of this gene. This is the lower coordinate end,
	 * regardless of which strand the gene lies on and its directionality.
	 * @return Gene start position.
	 */
	public int getTranscriptStart()
	{
		return stPos;
	}
	
	/**
	 * Get the end position of this gene. This is the higher coordinate end,
	 * regardless of which strand the gene lies on and its directionality.
	 * @return Gene end position.
	 */
	public int getTranscriptEnd()
	{
		return edPos;
	}
	
	/**
	 * Get the translation start position. This is the lower coordinate end,
	 * regardless of which strand the gene lies on and its directionality.
	 * <br>Not applicable for ncRNA, though this function will return something.
	 * @return Translation start position.
	 */
	public int getTranslationStart()
	{
		if (isNCRNA) return -1;
		return tlStart;
	}
	
	/**
	 * Get the translation end position. This is the higher coordinate end,
	 * regardless of which strand the gene lies on and its directionality.
	 * <br>Not applicable for ncRNA, though this function will return something.
	 * @return Translation end position.
	 */
	public int getTranslationStop()
	{
		if (isNCRNA) return -1;
		return tlEnd;
	}
	
	/**
	 * Get the number of exons this gene has recorded.
	 * @return Exon count for this gene.
	 */
	public int getExonCount()
	{
		return this.exons.size();
	}
	
	/**
	 * Get the Exon (as an Exon object) at the specified 0-based index.
	 * @param index Index in exon list of exon to retrieve. Keep in mind that
	 * though exon numbering generally starts with 1, the computer stores it starting
	 * at 0. Please use the 0-based coordinate.
	 * @return Exon object if index is valid. Null otherwise.
	 */
	public Exon getExon(int index)
	{
		if (index < 0) return null;
		if (index >= exons.size()) return null;
		return exons.get(index);
	}
	
	/* --- Setters --- */
	
	/**
	 * Set the gene name.
	 * @param newName String to set as gene name.
	 */
	public void setName(String newName)
	{
		name = newName;
	}
	
	/**
	 * Set the gene UID string.
	 * @param newID String to set as gene ID.
	 */
	public void setID(String newID)
	{
		ID = newID;
	}
	
	/**
	 * Set the gene strand.
	 * @param s True if plus(+) strand, false if minus(-) strand.
	 */
	public void setStrand(boolean s)
	{
		strand = s;
	}
	
	/**
	 * Set the ncRNA flag for this gene. If this flag is set, the
	 * gene is marked as ncRNA. If not, is assumed to be a coding gene
	 * with an mRNA transcript.
	 * @param b ncRNA flag.
	 */
	public void setNCRNA(boolean b)
	{
		isNCRNA = b;
	}
	
	/**
	 * Set the gene chromosome coordinate using a Contig reference.
	 * @param c Contig to link to this gene.
	 */
	public void setChromosome(Contig c)
	{
		chrom = c;
	}
	
	/**
	 * Set the start position of the gene. The start position cannot
	 * be larger than the end position, regardless of the gene directionality.
	 * @param pos New start position.
	 */
	public void setStart(int pos)
	{
		if (edPos >= 0 && pos > edPos) return;
		stPos = pos;
	}
	
	/**
	 * Set the end position of the gene. The end position cannot
	 * be smaller than the start position, regardless of the gene directionality,
	 * unless the input value is -1 (unset).
	 * @param pos New end position.
	 */
	public void setEnd(int pos)
	{
		if (edPos != -1 && edPos < stPos) return;
		edPos = pos;
	}
	
	/**
	 * Set the translation start position of the gene. This isn't 
	 * where the translation starts per se, it is the beginning of the 
	 * translated region. As a result, this value must be lower than translation
	 * end, regardless of gene directionality.
	 * @param pos New start position.
	 */
	public void setTranslationStart(int pos)
	{
		if (tlEnd >= 0 && (pos > tlEnd)) return;
		tlStart = pos;
	}
	
	/**
	 * Set the translation end position of the gene. This isn't 
	 * where the translation ends per se, it is the end of the 
	 * translated region. As a result, this value must be higher than translation
	 * start, regardless of gene directionality. The exception to this rule
	 * is if the input value is -1 (unset).
	 * @param pos New end position.
	 */
	public void setTranslationEnd(int pos)
	{
		if (pos != -1 && pos < tlStart) return;
		tlEnd = pos;
	}
	
	/**
	 * Add an exon to this gene. Its position will be determined relative
	 * to the existing exons recorded for this gene.
	 * <br>At the moment, there is no rejection of overlapping exons, so be careful.
	 * @param e Exon to add.
	 */
	public void addExon(Exon e)
	{
		exons.add(e);
		Collections.sort(exons);
	}
	
	/**
	 * Add a collection of exons to this gene. All exon numbers will be determined relative
	 * to the existing exons recorded for this gene.
	 * <br>At the moment, there is no rejection of overlapping exons, so be careful.
	 * @param eColl Collection of exons to add.
	 */
	public void addExons(Collection<Exon> eColl)
	{
		if (eColl == null) return;
		exons.addAll(eColl);
		if (!exons.isEmpty()) Collections.sort(exons);
	}
	
	/* --- Parsing/Serialization --- */
	
	/**
	 * Generate a gbdb formatted serialization of this gene for writing
	 * a gene set to gbdb.
	 * This overload does not include the chromosome ID as it is assumed
	 * to be included in a CHROM block.
	 * @return FileBuffer object (byte string wrapper) containing the serialized
	 * gene.
	 */
	public FileBuffer serializeMe()
	{
		boolean padid = false;
		boolean padn = false;
		int IDlen = ID.length() + 1;
		if (IDlen % 2 == 1){
			IDlen++;
			padid = true;
		}
		int nlen = name.length() + 1;
		if (nlen % 2 == 1){
			nlen++;
			padn = true;
		}
		int nexons = exons.size();
		
		FileBuffer myGene = new FileBuffer((4*4) + (2+2) + IDlen + nlen + (nexons * 8), true);
		myGene.addToFile(stPos);
		myGene.addToFile(edPos);
		myGene.addToFile(tlStart);
		myGene.addToFile(tlEnd);
		//Flags
		int flagraw = 0;
		if (isNCRNA) flagraw = 1;
		flagraw = flagraw << 1;
		if (strand) flagraw |= 1;
		myGene.addToFile((byte)flagraw);
		myGene.addToFile((byte)0x00); //Padding
		
		myGene.addToFile((short)nexons);
		
		for (int i = 0; i < nexons; i++)
		{
			Exon e = exons.get(i);
			myGene.addToFile(e.stPos);
			myGene.addToFile(e.edPos);
		}
		
		myGene.printASCIIToFile(ID + "\0");
		//myGene.addToFile((byte)0x00);
		if (padid) myGene.addToFile((byte)0x00);
		myGene.printASCIIToFile(name + "\0");
		//myGene.addToFile((byte)0x00);
		if (padn) myGene.addToFile((byte)0x00);
		
		return myGene;
	}
	
	/**
	 * Generate a gbdb formatted serialization of this gene for writing
	 * a gene set to gbdb.
	 * This overload includes a local chromosome UID, as provided. It may be used
	 * for a modified gbdb format.
	 * @param chromUID The file-local chromosome UID of this gene.
	 * @return FileBuffer object (byte string wrapper) containing the serialized
	 * gene.
	 */
	public FileBuffer serializeMe(int chromUID)
	{
		boolean padid = false;
		boolean padn = false;
		int IDlen = ID.length() + 1;
		if (IDlen % 2 == 1){
			IDlen++;
			padid = true;
		}
		int nlen = name.length() + 1;
		if (nlen % 2 == 1){
			nlen++;
			padn = true;
		}
		int nexons = exons.size();
		
		FileBuffer myGene = new FileBuffer(4 + (4*4) + (2+2) + IDlen + nlen + (nexons * 8), true);
		myGene.addToFile(chromUID);
		myGene.addToFile(stPos);
		myGene.addToFile(edPos);
		myGene.addToFile(tlStart);
		myGene.addToFile(tlEnd);
		//Flags
		int flagraw = 0;
		if (isNCRNA) flagraw = 1;
		flagraw = flagraw << 1;
		if (strand) flagraw |= 1;
		myGene.addToFile((byte)flagraw);
		myGene.addToFile((byte)0x00); //Padding
				
		myGene.addToFile((short)nexons);
		
		for (int i = 0; i < nexons; i++)
		{
			Exon e = exons.get(i);
			myGene.addToFile(e.stPos);
			myGene.addToFile(e.edPos);
		}
		
		myGene.printASCIIToFile(ID);
		myGene.addToFile((byte)0x00);
		if (padid) myGene.addToFile((byte)0x00);
		myGene.printASCIIToFile(name);
		myGene.addToFile((byte)0x00);
		if (padn) myGene.addToFile((byte)0x00);
		
		return myGene;
	}
	
	/* --- Utility --- */
	
	public boolean equals(Object o)
	{
		if (o == null) return false;
		if (o == this) return true;
		if (!(o instanceof Gene)) return false;
		
		Gene g = (Gene)o;
		if (!this.getID().equals(g.getID())) return false;
		if (!this.getName().equals(g.getName())) return false;
		if (this.is_ncRNA() != g.is_ncRNA()) return false;
		if (this.getStrand() != g.getStrand()) return false;
		if (!this.getChromosome().equals(g.getChromosome())) return false;
		if (this.getTranscriptStart() != g.getTranscriptStart()) return false;
		if (this.getTranscriptEnd() != g.getTranscriptEnd()) return false;
		if (this.getTranslationStart() != g.getTranslationStart()) return false;
		if (this.getTranslationStop() != g.getTranslationStop()) return false;
		if (this.getExonCount() != g.getExonCount()) return false;
		
		int eCount = this.getExonCount();
		for (int i = 0; i < eCount; i++)
		{
			if (!(this.getExon(i).equals(g.getExon(i)))) return false;
		}
		
		return true;
	}
	
	public int hashCode()
	{
		return name.hashCode() ^ ID.hashCode();
	}
	
	public int compareTo(Gene o) 
	{
		if (o == null) return 1;
		
		int chromcomp = this.getChromosome().compareTo(o.getChromosome());
		if (chromcomp != 0) return chromcomp;
		
		int comp = this.getTranscriptStart() - o.getTranscriptStart();
		if (comp != 0) return comp;
		
		comp = this.getTranscriptEnd() - o.getTranscriptEnd();
		if (comp != 0) return comp;
		
		comp = this.getTranslationStart() - o.getTranslationStart();
		if (comp != 0) return comp;
		
		comp = this.getTranslationStop() - o.getTranslationStop();
		if (comp != 0) return comp;
		
		if (this.is_ncRNA() && !o.is_ncRNA()) return 1;
		else if (!this.is_ncRNA() && o.is_ncRNA()) return -1;
		
		if (this.getStrand() && !o.getStrand()) return -1;
		else if (!this.getStrand() && o.getStrand()) return 1;
		
		if (this.getID() == null)
		{
			if (o.getID() != null) return -1;
		}
		else
		{
			if (o.getID() == null) return 1;
			else
			{
				comp = this.getID().compareTo(o.getID());
				if (comp != 0) return comp;
			}
		}
		
		if (this.getName() == null)
		{
			if (o.getName() != null) return -1;
		}
		else
		{
			if (o.getName() == null) return 1;
			else
			{
				comp = this.getName().compareTo(o.getName());
				if (comp != 0) return comp;
			}
		}
		
		return 0;
	}

	/* --- Search --- */
	
	/**
	 * Determine whether a position (assuming it is on the same Contig as
	 * this gene) falls inside this gene.
	 * @param position Position to query.
	 * @return True - If position would fall inside gene if it is on the same
	 * Contig.
	 * <br>False - If position would fall outside of this gene regardless of Contig.
	 */
	public boolean positionInGene(int position)
	{
		return (position < edPos) && (position >= stPos);
	}
	
	/**
	 * Determine whether a position (assuming it is on the same Contig as
	 * this gene) falls inside an exon in this gene, and if so, which exon
	 * it falls in.
	 * @param position Position to query.
	 * @return 0-based coordinate of exon the position falls inside of, if it
	 * is on the same contig as the gene. -1 if the position falls outside of
	 * an exon, regardless of contig.
	 */
	public synchronized int positionInExon(int position)
	{
		int ex = exons.size();
		//ZERO based exon count!!!!!
		for (int i = 0; i < ex; i++)
		{
			if (exons.get(i).positionInExon(position)) return i;
		}
		return -1;
	}
	
	/**
	 * Determine whether a position (assuming it is on the same Contig as
	 * this gene) falls inside this gene's 5' UTR region.
	 * <br>This function returns false for ncRNA.
	 * <br>This function takes strand into account.
	 * @param position Position to query.
	 * @return True - If position would fall inside 5' UTR if it is on the same
	 * Contig.
	 * <br>False - If position would fall outside of this gene's 5' UTR regardless of Contig.
	 */
	public boolean positionInUTR5(int position)
	{
		if (isNCRNA) return false;
		if (strand) return (position < tlStart) && (position >= stPos);
		return (position >= tlEnd) && (position < edPos);
	}
	
	/**
	 * Determine whether a position (assuming it is on the same Contig as
	 * this gene) falls inside this gene's 3' UTR.
	 * <br>This function returns false for ncRNA.
	 * <br>This function takes strand into account.
	 * @param position Position to query.
	 * @return True - If position would fall inside 3' UTR if it is on the same
	 * Contig.
	 * <br>False - If position would fall outside of this gene's 3' UTR regardless of Contig.
	 */
	public boolean positionInUTR3(int position)
	{
		if (isNCRNA) return false;
		if (strand) return (position < edPos) && (position >= tlEnd);
		return (position < tlStart) && (position >= stPos);
	}
	
	/**
	 * Determine whether a position (assuming it is on the same Contig as
	 * this gene) falls inside a splicing region in this gene, and if so, which exon
	 * the splicing region is associated with.
	 * @param position Position to query.
	 * @return 0-based coordinate of exon the position falls in the splicing region of, if it
	 * is on the same contig as the gene. -1 if the position falls outside of
	 * a splicing region, regardless of contig.
	 */	
	public synchronized int positionSplicing(int position)
	{
		int ex = exons.size();
		//ZERO based exon count!!!!!
		for (int i = 0; i < ex; i++)
		{
			if (i == 0){
				if (exons.get(i).positionInSpliceZone(position, true, false)) return i;
			}
			else if (i == ex - 1){
				if (exons.get(i).positionInSpliceZone(position, false, true)) return i;
			}
			else{
				if (exons.get(i).positionInSpliceZone(position, false, false)) return i;
			}
		}
		return -1;
	}

	/**
	 * Determine whether a position (assuming it is on the same Contig as
	 * this gene) falls in an intergenic area that might be considered upstream of this gene.
	 * <br>This function takes strand into account.
	 * @param position Position to query.
	 * @return True - If position would fall within a certain distance upstream of gene if it is on the same
	 * Contig.
	 * <br>False - If position would fall outside of this gene's "upstream" region regardless of Contig.
	 */
	public boolean positionIsUpstream(int position)
	{
		int dist = getTranscriptStart() - position;
		if (dist > 0 && dist <= GeneSet.NEARGENE_DIST)
		{
			if (this.getStrand()) return true;
		}
		dist = position - getTranscriptEnd();
		if (dist >= 0 && dist <= GeneSet.NEARGENE_DIST)
		{
			if (!this.getStrand()) return true;
		}
		
		return false;
	}
	
	/**
	 * Determine whether a position (assuming it is on the same Contig as
	 * this gene) falls in an intergenic area that might be considered downstream of this gene.
	 * <br>This function takes strand into account.
	 * @param position Position to query.
	 * @return True - If position would fall within a certain distance downstream of gene if it is on the same
	 * Contig.
	 * <br>False - If position would fall outside of this gene's "downstream" region regardless of Contig.
	 */
	public boolean positionIsDownstream(int position)
	{
		int dist = getTranscriptStart() - position;
		if (dist > 0 && dist <= GeneSet.NEARGENE_DIST)
		{
			if (!this.getStrand()) return true;
		}
		dist = position - getTranscriptEnd();
		if (dist >= 0 && dist <= GeneSet.NEARGENE_DIST)
		{
			if (this.getStrand()) return true;
		}
		return false;
	}
	
	/**
	 * Determine whether another gene overlaps this gene.
	 * @param o Other gene to compare.
	 * @return True - If gene o overlaps this gene.
	 * <br>False - If the other gene does not overlap this gene.
	 */
	public boolean overlaps(Gene o)
	{
		if (o.getChromosome() != this.getChromosome()) return false;
		
		int oSt = o.getTranscriptStart();
		int oEd = o.getTranscriptEnd();
		int tSt = this.getTranscriptStart();
		int tEd = this.getTranscriptEnd();
		
		if (oSt < tSt)
		{
			if (oEd >= tSt) return true;
			return false;
		}
		else if (oSt > tSt)
		{
			if (tEd >= oSt) return true;
			return false;
		}
		else return true;

		//return false;
	}

	/**
	 * Get the location effect of a position relative to this gene as a GeneFunc
	 * enum. 
	 * @param position Position to evaluate.
	 * @return GeneFunc enum describing part of the gene the position falls in.
	 */
	public synchronized GeneFunc getRelativeLocationEffect(int position)
	{
		if (positionInGene(position))
		{
			if (positionInUTR5(position)) return GeneFunc.UTR5;
			if (positionInUTR3(position)) return GeneFunc.UTR3;
			if (positionInExon(position) >= 0)
			{
				if (is_ncRNA()) return GeneFunc.NCRNA;
				else return GeneFunc.EXONIC;
			}
			else
			{
				if (positionSplicing(position) >= 0)
				{
					if (is_ncRNA()) return GeneFunc.NCRNA;
					else return GeneFunc.SPLICING;
				}
				else 
				{
					return GeneFunc.INTRONIC;
				}
			}
		}
		else
		{
			int dist = getTranscriptStart() - position;
			if (dist > 0 && dist <= GeneSet.NEARGENE_DIST)
			{
				if (this.getStrand()) return GeneFunc.UPSTREAM;
				else return GeneFunc.DOWNSTREAM;
			}
			dist = position - getTranscriptEnd();
			if (dist >= 0 && dist <= GeneSet.NEARGENE_DIST)
			{
				if (this.getStrand()) return GeneFunc.DOWNSTREAM;
				else return GeneFunc.UPSTREAM;
			}
		}
		return GeneFunc.INTERGENIC;
	}

	/**
	 * Get the 0-based index of the next exon after position in question.
	 * @param position Position (assuming same contig as gene) to query.
	 * @return Index (0-based) of next exon AFTER position (even if position falls
	 * inside of an exon).
	 * <br> -1 if position falls after final exon.
	 */
	public synchronized int nextExon(int position)
	{
		//0 based coordinates. Returns -1 if none, or after last exon (which would be 3' UTR)
		int ex = exons.size();
		//ZERO based exon count!!!!!
		//Generously assuming that the exons are ordered...
		for (int i = 0; i < ex; i++)
		{
			if (exons.get(i).isBefore(position)) return i;
		}
		return -1;
	}
	
	/**
	 * Get the location effect of a region specified by two positions (assumed to be on the
	 * same contig as the gene) as a GeneFunc enum.
	 * @param stpos Start position (lower coordinate - inclusive) of region.
	 * @param edpos End position (higher coordinate - exclusive) of region.
	 * @return GeneFunc enum describing part of the gene the position falls in.
	 */
	public synchronized GeneFunc getRelativeRegionLocationEffect(int stpos, int edpos)
	{
		//Make sure stpos is lower value
		if (stpos > edpos)
		{
			int temp = stpos;
			stpos = edpos;
			edpos = temp;
		}
		//See if any exons fall between the points
		boolean exonCaught = false;
		for (Exon e : exons)
		{
			//Is start in or before exon?
			if (stpos < e.stPos)
			{
				//Start is before exon start
				//Check end
				if (edpos <= e.stPos)
				{
					//End is before exon
					continue;
				}
				else
				{
					//End is inside or after exon
					exonCaught = true;
					break;
				}
			}
			else if (stpos >= e.stPos && stpos < e.edPos)
			{
				//Start is inside exon
				exonCaught = true;
				break;
			}
			else
			{
				//Start is after exon end
				continue;
			}
		}
		
		if (exonCaught)
		{
			//Determine if ncRNA
			if (this.is_ncRNA()) return GeneFunc.NCRNA;
			else return GeneFunc.EXONIC;
		}
		
		//Assume no exons were caught
		if (stpos >= this.edPos && edpos > this.edPos)
		{
			//Region falls after gene
			//See if start is in upstream/downstream region
			int dist = stpos - getTranscriptEnd();
			if (dist >= 0 && dist <= GeneSet.NEARGENE_DIST)
			{
				if (this.getStrand()) return GeneFunc.DOWNSTREAM;
				else return GeneFunc.UPSTREAM;
			}
			return GeneFunc.INTERGENIC;
		}
		else if (stpos< this.stPos && edpos <= this.stPos)
		{
			//Region falls before gene
			int dist = getTranscriptStart() - edpos;
			if (dist > 0 && dist <= GeneSet.NEARGENE_DIST)
			{
				if (this.getStrand()) return GeneFunc.UPSTREAM;
				else return GeneFunc.DOWNSTREAM;
			}
			return GeneFunc.INTERGENIC;
		}
		else 
		{
			if (this.is_ncRNA()) return GeneFunc.NCRNA;
			//Check for splicing/UTR. Otherwise intronic
			if (positionSplicing(stpos) >= 0) return GeneFunc.SPLICING;
			if (positionSplicing(edpos) >= 0) return GeneFunc.SPLICING;
			if (positionInUTR5(stpos) || positionInUTR5(edpos-1)) return GeneFunc.UTR5;
			if (positionInUTR3(stpos) || positionInUTR3(edpos-1)) return GeneFunc.UTR3;
			return GeneFunc.INTRONIC;
		}

		
	}
	
	/* --- View --- */
	
	/**
	 * Print information about this gene as a row of tab-delimited text
	 * @return A text line (like a table entry) representing this gene...
	 */
	public String printInfo()
	{
		String line = "";
		line += getID() + "\t";
		line += getName() + "\t";
		if (this.is_ncRNA()) line += "ncRNA\t";
		else line += "mRNA\t";
		line += getChromosome().getUDPName() + "\t";
		if (getStrand()) line += "+\t";
		else line += "-\t";
		line += getTranscriptStart() + "\t";
		line += getTranscriptEnd() + "\t";
		line += getTranslationStart() + "\t";
		line += getTranslationStop() + "\t";
		int eCount = getExonCount();
		line += eCount + "\t";
		for (int i = 0; i < eCount; i++)
		{
			Exon e = getExon(i);
			line += "[" + e.getStart() + "," + e.getEnd() + "]";
		}
		return line;
	}
	
	public String toString()
	{
		return this.getName();
	}
	
}
