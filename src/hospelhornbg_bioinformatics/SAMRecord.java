package hospelhornbg_bioinformatics;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import hospelhornbg_genomeBuild.Contig;
import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.BitStreamer;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class SAMRecord implements Comparable<SAMRecord>{
	
	/* ----- Static Variables ----- */
	
	private static SAMSortOrder sortOrder = SAMSortOrder.COORDINATE;
	
	/* ----- Instance Variables ----- */
	
	private String qname;
	private int flags;
	
	private Contig reference;
	private int position; //1-based
	private int mapq;
	private String CIGAR;
	
	private Contig refNext;
	private int pNext;
	
	private int templateLength;
	private String sequence;
	private int[] seqQuality;

	private Map<String, SAMField> alignmentFields;
	
	private WarningFlags warnings;
	
	/* ----- Inner Objects ----- */
	
	public static class FailFlags
	{
		public boolean err_syntax;
		public boolean err_nullseq_nnqual;
		public String err_bad_customFieldType;
	}
	
	public static class WarningFlags
	{
		public String err_invalid_RNAME;
		public int err_invalid_POS;
		public String err_invalid_RNEXT;
		public int err_invalid_PNEXT;
		
		public boolean err_qualstr_len_bad;
		
		public WarningFlags()
		{
			err_invalid_RNAME = null;
			err_invalid_POS = 0;
			err_invalid_RNEXT = null;
			err_invalid_PNEXT = 0;
			err_qualstr_len_bad = false;
		}
	}
	
	public static class InvalidSAMRecordException extends Exception
	{
		private static final long serialVersionUID = 6376597106929888202L;
		
		private FailFlags reason;
		
		public InvalidSAMRecordException(FailFlags flags)
		{
			reason = flags;
		}
		
		public FailFlags getFlags()
		{
			return reason;
		}
		
	}
	
	public static class ParsedRecord
	{
		private SAMRecord record;
		
		//public String raw_qname;
		//public String raw_flag;
		public String raw_rname;
		//public String raw_pos;
		//public String raw_mapq;
		//public String raw_cigar;
		public String raw_rnext;
		//public String raw_pnext;
		//public String raw_tlen;
		//public String raw_seq;
		//public String raw_qual;
		//public String[] raw_customs;
		
		public ParsedRecord(SAMRecord r)
		{
			record = r;
		}
		
		public SAMRecord getRecord()
		{
			return record;
		}
	}
	
	/* ----- Construction/Parsing ----- */
	
	public SAMRecord()
	{
		qname = "*";
		flags = 0;
		reference = null; //Serialized as *
		position = 1;
		mapq = 0;
		CIGAR = "*";
		refNext = null;
		pNext = 1;
		templateLength = 0;
		sequence = "";
		seqQuality = null;
		alignmentFields = new HashMap<String, SAMField>();
		warnings = new WarningFlags();
	}
	
	public SAMRecord(String seq)
	{
		qname = "*";
		flags = 0;
		reference = null; //Serialized as *
		position = 1;
		mapq = 0;
		CIGAR = "*";
		refNext = null;
		pNext = 1;
		templateLength = 0;
		sequence = seq;
		seqQuality = new int[seq.length()];
		alignmentFields = new HashMap<String, SAMField>();
		warnings = new WarningFlags();
	}
	
	public static ParsedRecord parseSAMRecord(String samLine, GenomeBuild gbuild, boolean verbose) throws UnsupportedFileTypeException, InvalidSAMRecordException
	{
		if (samLine == null)
		{
			if(verbose)System.err.println("SAMRecord.parseSAMRecord || ERROR: No record provided for parsing...");
			throw new FileBuffer.UnsupportedFileTypeException();
		}
		if (gbuild == null)
		{
			if(verbose)System.err.println("SAMRecord.parseSAMRecord || ERROR: Genome build required for SAM record parsing!");
			throw new FileBuffer.UnsupportedFileTypeException();
		}
		
		SAMRecord rec = new SAMRecord();
		ParsedRecord p = new ParsedRecord(rec);
		
		String[] fields = samLine.split("\t");
		
		if (fields.length < 11)
		{
			if(verbose)System.err.println("SAMRecord.parseSAMRecord || ERROR: SAM Record contains insufficient number of fields (" + fields.length + ")!");
			if(verbose)System.err.println(samLine);
			FailFlags ff = new FailFlags(); ff.err_syntax = true;
			throw new InvalidSAMRecordException(ff);
		}
		
		//QNAME
		//p.raw_qname = fields[0];
		if (fields[0].equals("*")) rec.qname = null;
		else rec.qname = fields[0];
		
		//FLAG
		//p.raw_flag = fields[1];
		try
		{
			rec.flags = Integer.parseInt(fields[1]);
		}
		catch (NumberFormatException e)
		{
			if(verbose)System.err.println("SAMRecord.parseSAMRecord || ERROR: Flag record (" + fields[1] + ") could not be parsed as a decimal integer!");
			if(verbose)System.err.println(samLine);
			FailFlags ff = new FailFlags(); ff.err_syntax = true;
			throw new InvalidSAMRecordException(ff);
		}
		
		//RNAME
		p.raw_rname = fields[2];
		if (fields[2].equals("*")) rec.reference = null;
		else
		{
			rec.reference = gbuild.getContig(fields[2]);
			if (rec.reference == null) {
				//if(verbose) System.err.println("SAMRecord.parseSAMRecord || WARNING: RNAME contig \"" + fields[2] + "\" did not match any contig in the provided build.");
				rec.warnings.err_invalid_RNAME = fields[2];
			}
		}
		
		//POS
		//p.raw_pos = fields[3];
		try
		{
			rec.position = Integer.parseInt(fields[3]);
			if (rec.reference != null && (rec.position < 0 || rec.position > rec.reference.getLength())) rec.warnings.err_invalid_POS = rec.position;
		}
		catch (NumberFormatException e)
		{
			System.err.println("SAMRecord.parseSAMRecord || ERROR: POS record (" + fields[3] + ") could not be parsed as a decimal integer!");
			System.err.println(samLine);
			FailFlags ff = new FailFlags(); ff.err_syntax = true;
			throw new InvalidSAMRecordException(ff);
		}
		
		//MAPQ
		//p.raw_mapq = fields[4];
		try
		{
			rec.mapq = Integer.parseInt(fields[4]);
		}
		catch (NumberFormatException e)
		{
			if(verbose)System.err.println("SAMRecord.parseSAMRecord || ERROR: MAPQ record (" + fields[4] + ") could not be parsed as a decimal integer!");
			if(verbose)System.err.println(samLine);
			FailFlags ff = new FailFlags(); ff.err_syntax = true;
			throw new InvalidSAMRecordException(ff);
		}
		
		//CIGAR
		//p.raw_cigar = fields[5];
		if (fields[5].equals("*")) rec.CIGAR = null;
		else rec.CIGAR = fields[5];
		
		//RNEXT
		p.raw_rnext = fields[6];
		if (fields[6].equals("*")) rec.refNext = null;
		else if (fields[6].equals("=")) rec.refNext = rec.reference;
		else
		{
			rec.refNext = gbuild.getContig(fields[6]);
			if (rec.refNext == null) {
				//if(verbose) System.err.println("SAMRecord.parseSAMRecord || WARNING: RNEXT contig \"" + fields[6] + "\" did not match any contig in the provided build.");
				rec.warnings.err_invalid_RNEXT = fields[6];
			}
		}
		
		//PNEXT
		//p.raw_pnext = fields[7];
		try
		{
			rec.pNext = Integer.parseInt(fields[7]);
			if (rec.refNext != null && (rec.pNext < 0 || rec.pNext > rec.refNext.getLength())) rec.warnings.err_invalid_PNEXT = rec.pNext;
		}
		catch (NumberFormatException e)
		{
			if(verbose)System.err.println("SAMRecord.parseSAMRecord || ERROR: PNEXT record (" + fields[7] + ") could not be parsed as a decimal integer!");
			if(verbose)System.err.println(samLine);
			FailFlags ff = new FailFlags(); ff.err_syntax = true;
			throw new InvalidSAMRecordException(ff);
		}
		
		//TLEN
		//p.raw_tlen = fields[8];
		try
		{
			rec.templateLength = Integer.parseInt(fields[8]);
		}
		catch (NumberFormatException e)
		{
			if(verbose)System.err.println("SAMRecord.parseSAMRecord || ERROR: TLEN record (" + fields[8] + ") could not be parsed as a decimal integer!");
			if(verbose)System.err.println(samLine);
			FailFlags ff = new FailFlags(); ff.err_syntax = true;
			throw new InvalidSAMRecordException(ff);
		}
		
		//SEQ
		//p.raw_seq = fields[9];
		if (fields[9].equals("*")) rec.sequence = null;
		else rec.sequence = fields[9];
		
		//QUAL
		//p.raw_qual = fields[10];
		String qualstring = fields[10];
		if (qualstring.equals("*"))
		{
			if(verbose)System.err.println("SAMRecord.parseSAMRecord || WARNING: Phred scaled base call qualities not recorded for this read!");
			rec.seqQuality = null;
		}
		else
		{
			if (rec.sequence == null)
			{
				if(verbose)System.err.println("SAMRecord.parseSAMRecord || ERROR: Non-null quality string cannot accompany null sequence string!");
				if(verbose)System.err.println(samLine);
				FailFlags ff = new FailFlags(); ff.err_nullseq_nnqual = true;
				throw new InvalidSAMRecordException(ff);
			}
			int slen = rec.sequence.length();
			int qlen = qualstring.length();
			if (slen != qlen)
			{
				if(verbose)System.err.println("SAMRecord.parseSAMRecord || WARNING: Length of sequence string does NOT match length of quality string! (slen = " + slen + " | qlen = " + qlen + ")");
			}
			if (qlen == 0)
			{
				if(verbose)System.err.println("SAMRecord.parseSAMRecord || WARNING: Quality string is empty!");
			}
			else
			{
				rec.seqQuality = new int[qlen];
				byte[] sbytes = qualstring.getBytes();
				for (int i = 0; i < qlen; i++)
				{
					int q = Byte.toUnsignedInt(sbytes[i]);
					rec.seqQuality[i] = q - 33;
				}
			}	
		}
		
		//Optional Fields
		if (fields.length > 11)
		{
			int flen = fields.length;
			//p.raw_customs = new String[flen - 11];
			for (int j = 11; j < flen; j++)
			{
				//p.raw_customs[j - 11] = fields[j];
				SAMField ofield = parseOptionalField(fields[j], verbose);
				rec.alignmentFields.put(ofield.getKey(), ofield);
			}
		}
		
		return p;
	}
	
	public static SAMField parseOptionalField(String field, boolean verbose) throws UnsupportedFileTypeException, InvalidSAMRecordException
	{
		if (field == null || field.isEmpty())
		{
			if(verbose)System.err.println("SAMRecord.parseOptionalField || ERROR: No field provided for parsing...");
			FailFlags ff = new FailFlags(); ff.err_nullseq_nnqual = true;
			throw new InvalidSAMRecordException(ff);
		}
		
		String[] fields = field.split(":");
		
		if (fields.length != 3)
		{
			if(verbose)System.err.println("SAMRecord.parseOptionalField || ERROR: Custom field is not formatted correctly: " + field);
			FailFlags ff = new FailFlags(); ff.err_nullseq_nnqual = true;
			throw new InvalidSAMRecordException(ff);
		}
		
		String tag = fields[0];
		char type = fields[1].charAt(0);
		String upval = fields[2];
		
		switch (type)
		{
		case 'A':
			return new SAMCharField(tag, upval);
		case 'i':
			return new SAMIntField(tag, upval);
		case 'f':
			return new SAMFloatField(tag, upval);
		case 'Z':
			return new SAMStringField(tag, upval);
		case 'H':
			return new SAMBytesField(tag, upval);
		case 'B':
			return new SAMNumbersField(tag, upval);
		}
		
		if(verbose)System.err.println("SAMRecord.parseOptionalField || ERROR: Custom field type not recoginized: " + field);
		FailFlags ff = new FailFlags(); ff.err_bad_customFieldType = field;
		throw new InvalidSAMRecordException(ff);
	}
	
	/* ----- Getters ----- */
	
	public String getQueryName()
	{
		return qname;
	}
	
	public boolean getFlag(int bit)
	{
		return (flags & (0x1 << bit)) != 0;
	}
	
	public boolean flaggedSegmented()
	{
		return (flags & 0x0001) != 0;
	}
	
	public boolean flaggedSegmentsAligned()
	{
		return (flags & 0x0002) != 0;
	}
	
	public boolean flaggedSegmentUnmapped()
	{
		return (flags & 0x0004) != 0;
	}
	
	public boolean flaggedNextSegmentUnmapped()
	{
		return (flags & 0x0008) != 0;
	}
	
	public boolean flaggedReverseComplemented()
	{
		return (flags & 0x0010) != 0;
	}
	
	public boolean flaggedNextReverseComplemented()
	{
		return (flags & 0x0020) != 0;
	}
	
	public boolean flaggedFirstSegment()
	{
		return (flags & 0x0040) != 0;
	}
	
	public boolean flaggedLastSegment()
	{
		return (flags & 0x0080) != 0;
	}
	
	public boolean flaggedSecondary()
	{
		return (flags & 0x0100) != 0;
	}
	
	public boolean flaggedQCFail()
	{
		return (flags & 0x0200) != 0;
	}
	
	public boolean flaggedDuplicate()
	{
		return (flags & 0x0400) != 0;
	}
	
	public boolean flaggedSupplementary()
	{
		return (flags & 0x0800) != 0;
	}

	public Contig getReferenceContig()
	{
		return reference;
	}
	
	public int getPosition()
	{
		return position;
	}
	
	public int getMapQuality()
	{
		return mapq;
	}
	
	public String getCIGAR()
	{
		return CIGAR;
	}
	
	public Contig getNextReferenceContig()
	{
		return this.refNext;
	}
	
	public int getNextPosition()
	{
		return this.pNext;
	}
	
	public int getTemplateLength()
	{
		return this.templateLength;
	}
	
	public String getSequence()
	{
		return this.sequence;
	}
	
	public int getPhredBaseQuality(int index)
	{
		if (this.seqQuality == null) return -1;
		if (index < 0 || index >= this.seqQuality.length) return -1;
		return this.seqQuality[index];
	}
	
	public SAMField getCustomField(String TAG)
	{
		return this.alignmentFields.get(TAG);
	}
	
	public List<SAMField> getAllOptionalFields()
	{
		int fcount = alignmentFields.size();
		List<SAMField> list = new ArrayList<SAMField>(fcount+1);
		
		Collection<SAMField> opvals = alignmentFields.values();
		for (SAMField f : opvals) list.add(f);
		Collections.sort(list);
		
		return list;
	}
	
	public WarningFlags getParserWarnings()
	{
		return this.warnings;
	}
	
	/* ----- Setters ----- */
	
	public void setQueryName(String s)
	{
		qname = s;
	}
	
	public void setFlag(int bit)
	{
		flags = BitStreamer.writeABit(flags, true, bit);
	}
	
	public void unsetFlag(int bit)
	{
		flags = BitStreamer.writeABit(flags, false, bit);
	}

	public void flagSegmented(boolean flag)
	{
		if (flag) flags |= 0x0001;
		else flags &= ~(0x0001);
	}
	
	public void flagSegmentsAligned(boolean flag)
	{
		if (flag) flags |= 0x0002;
		else flags &= ~(0x0002);
	}
	
	public void flagSegmentUnmapped(boolean flag)
	{
		if (flag) flags |= 0x0004;
		else flags &= ~(0x0004);
	}
	
	public void flagNextSegmentUnmapped(boolean flag)
	{
		if (flag) flags |= 0x0008;
		else flags &= ~(0x0008);
	}
	
	public void flagReverseComplemented(boolean flag)
	{
		if (flag) flags |= 0x0010;
		else flags &= ~(0x0010);
	}
	
	public void flagNextReverseComplemented(boolean flag)
	{
		if (flag) flags |= 0x0020;
		else flags &= ~(0x0020);
	}
	
	public void flagFirstSegment(boolean flag)
	{
		if (flag) flags |= 0x0040;
		else flags &= ~(0x0040);
	}
	
	public void flagLastSegment(boolean flag)
	{
		if (flag) flags |= 0x0080;
		else flags &= ~(0x0080);
	}
	
	public void flagSecondary(boolean flag)
	{
		if (flag) flags |= 0x0100;
		else flags &= ~(0x0100);
	}
	
	public void flagQCFail(boolean flag)
	{
		if (flag) flags |= 0x0200;
		else flags &= ~(0x0200);
	}
	
	public void flagDuplicate(boolean flag)
	{
		if (flag) flags |= 0x0400;
		else flags &= ~(0x0400);
	}
	
	public void flagSupplementary(boolean flag)
	{
		if (flag) flags |= 0x0800;
		else flags &= ~(0x0800);
	}

	public void setReferenceContig(Contig c)
	{
		reference = c;
	}
	
	public void setPosition(int i)
	{
		this.position = i;
	}
	
	public void setMapQuality(int q)
	{
		this.mapq = q;
	}
	
	public void setCIGAR(String cigarString)
	{
		CIGAR = cigarString;
	}
	
	public void setNextReferenceContig(Contig c)
	{
		this.refNext = c;
	}
	
	public void setNextPosition(int i)
	{
		this.pNext = i;
	}
	
	public void setTemplateLength(int len)
	{
		this.templateLength = len;
	}
	
	public void setSequence(String seq)
	{
		this.sequence = seq;
	}
	
	public void setQualities(int[] newQualities)
	{
		if(newQualities == null) seqQuality = null;
		seqQuality = new int[newQualities.length];
		for (int i = 0; i < newQualities.length; i++) seqQuality[i] = newQualities[i];
	}
	
	public void addCustomField(SAMField f)
	{
		this.alignmentFields.put(f.getKey(), f);
	}
	
	public void removeCustomField(String TAG)
	{
		this.alignmentFields.remove(TAG);
	}
	
	public void clearCustomFields()
	{
		this.alignmentFields.clear();
	}
	
	/* ----- Serialization ----- */
	
	public String writeSAMRecord(boolean use_UCSC_contig_names)
	{
		String s = "";
		
		if(qname != null) s += qname + "\t";
		else s += "*\t";
		s += flags + "\t"; //Decimal for some reason :|
		if (reference == null) s += "*\t";
		else
		{
			if (use_UCSC_contig_names) s += reference.getUCSCName() + "\t";
			else s += reference.getUDPName() + "\t";
		}
		s += position + "\t";
		s += mapq + "\t";
		if(CIGAR != null) s += CIGAR + "\t";
		else s += "*\t";
		if (refNext == null) s += "*\t";
		else if (refNext.equals(reference)) s += "=\t";
		else
		{
			if (use_UCSC_contig_names) s += refNext.getUCSCName() + "\t";
			else s += refNext.getUDPName() + "\t";
		}
		s += pNext + "\t";
		s += templateLength + "\t";
		if(sequence != null) s += sequence + "\t";
		else s += "*\t";
		
		if (seqQuality == null || seqQuality.length < 1) s += "*\t";
		else
		{
			for(int q : seqQuality)
			{
				s += (char)(q + 33);
			}
		}
		
		if (alignmentFields == null || alignmentFields.isEmpty()) return s;
		List<SAMField> opfields = this.getAllOptionalFields();
		
		for (SAMField f : opfields)
		{
			s += "\t" + f.getSAMString();
		}
		
		return s;
	}
	
	/* ----- Utility ----- */
	
	public void unmap()
	{
		reference = null;
		flagSegmentUnmapped(true);
	}
	
	public void unmapNext()
	{
		refNext = null;
		if (flaggedSegmented()) flagNextSegmentUnmapped(true);
	}
	
	/* ----- Sort ----- */
	
	public static void setSortOrder(SAMSortOrder so)
	{
		sortOrder = so;
	}
	
	public int hashCode()
	{
		return this.writeSAMRecord(false).hashCode();
	}
	
	public boolean equals(Object o)
	{
		if (o == null) return false;
		if (o == this) return true;
		if (!(o instanceof SAMRecord)) return false;
		return this.hashCode() == o.hashCode();
	}

	public boolean isNextSegment(SAMRecord rec)
	{
		if (rec == null) return false;
		
		Contig next = this.getNextReferenceContig();
		Contig oc = rec.getReferenceContig();
		if (next == null && oc != null) return false;
		if (next != null && !next.equals(oc)) return false;
		
		//Assuming same contig
		//Check pos
		int pnext = this.getNextPosition();
		int opos = rec.getPosition();
		
		if (pnext != opos) return false;
		
		return true;
	}
	
	private int compareByQueryInfo(SAMRecord other)
	{
		//Query name
		String thisqn = this.getQueryName();
		String otherqn = other.getQueryName();
		if (thisqn != null)
		{
			if(!thisqn.equals(otherqn)) return thisqn.compareTo(otherqn);
		}
		else
		{
			if (otherqn != null) return -1;
		}
		
		//Check for first/last flags
		if (this.flaggedFirstSegment()) return -1;
		if (other.flaggedFirstSegment()) return 1;
		if (this.flaggedLastSegment()) return 1;
		if (other.flaggedLastSegment()) return -1;
		
		//See if one is directly after the other
		if (this.isNextSegment(other)) return -1; //Other comes after this
		if (other.isNextSegment(this)) return 1; //This comes after other
		
		//If no other info, we have no data. Just return 0.
		
		return 0;
	}
	
	@Override
	public int compareTo(SAMRecord other) 
	{
		if (other == null) return 1;
		if (other == this) return 0;
		
		if (sortOrder == SAMSortOrder.COORDINATE)
		{
			//Coordinate mode	
			//Get contig
			Contig ct = this.getReferenceContig();
			Contig co = other.getReferenceContig();
			if (ct == null && co != null) return 1;
			else if (co == null && ct != null) return -1;
			else if (ct != null && co != null)
			{
				if (!ct.equals(co)) return ct.compareTo(co);	
			}
			//Position
			int pt = this.getPosition();
			int po = other.getPosition();
			if (pt != po) return pt - po;
			
			//Length
			int lt = this.getTemplateLength();
			int lo = other.getTemplateLength();
			if (lt != lo) return lt - lo;
			
			return compareByQueryInfo(other);
		}
		else if (sortOrder == SAMSortOrder.QUERY_NAME)
		{
			//Query group mode	
			return compareByQueryInfo(other);
		}
		
		return 0;
	}
	
}
