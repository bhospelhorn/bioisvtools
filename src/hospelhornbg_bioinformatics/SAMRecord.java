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

public class SAMRecord {
	
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
	}
	
	public static SAMRecord parseSAMRecord(String samLine, GenomeBuild gbuild) throws UnsupportedFileTypeException
	{
		if (samLine == null)
		{
			System.err.println("SAMRecord.parseSAMRecord || ERROR: No record provided for parsing...");
			throw new FileBuffer.UnsupportedFileTypeException();
		}
		if (gbuild == null)
		{
			System.err.println("SAMRecord.parseSAMRecord || ERROR: Genome build required for SAM record parsing!");
			throw new FileBuffer.UnsupportedFileTypeException();
		}
		
		SAMRecord rec = new SAMRecord();
		
		String[] fields = samLine.split("\t");
		
		if (fields.length < 11)
		{
			System.err.println("SAMRecord.parseSAMRecord || ERROR: SAM Record contains insufficient number of fields (" + fields.length + ")!");
			System.err.println(samLine);
			throw new FileBuffer.UnsupportedFileTypeException();
		}
		
		//QNAME
		if (fields[0].equals("*")) rec.qname = null;
		else rec.qname = fields[0];
		
		//FLAG
		try
		{
			rec.flags = Integer.parseInt(fields[1]);
		}
		catch (NumberFormatException e)
		{
			System.err.println("SAMRecord.parseSAMRecord || ERROR: Flag record (" + fields[1] + ") could not be parsed as a decimal integer!");
			System.err.println(samLine);
			throw new FileBuffer.UnsupportedFileTypeException();
		}
		
		//RNAME
		if (fields[2].equals("*")) rec.reference = null;
		else
		{
			rec.reference = gbuild.getContig(fields[2]);
			if (rec.reference == null) System.err.println("SAMRecord.parseSAMRecord || WARNING: RNAME contig \"" + fields[2] + "\" did not match any contig in the provided build.");
		}
		
		//POS
		try
		{
			rec.position = Integer.parseInt(fields[3]);
		}
		catch (NumberFormatException e)
		{
			System.err.println("SAMRecord.parseSAMRecord || ERROR: POS record (" + fields[3] + ") could not be parsed as a decimal integer!");
			System.err.println(samLine);
			throw new FileBuffer.UnsupportedFileTypeException();
		}
		
		//MAPQ
		try
		{
			rec.mapq = Integer.parseInt(fields[4]);
		}
		catch (NumberFormatException e)
		{
			System.err.println("SAMRecord.parseSAMRecord || ERROR: MAPQ record (" + fields[4] + ") could not be parsed as a decimal integer!");
			System.err.println(samLine);
			throw new FileBuffer.UnsupportedFileTypeException();
		}
		
		//CIGAR
		if (fields[5].equals("*")) rec.CIGAR = null;
		else rec.CIGAR = fields[5];
		
		//RNEXT
		if (fields[6].equals("*")) rec.refNext = null;
		else if (fields[6].equals("=")) rec.refNext = rec.reference;
		else
		{
			rec.refNext = gbuild.getContig(fields[6]);
			if (rec.refNext == null) System.err.println("SAMRecord.parseSAMRecord || WARNING: RNEXT contig \"" + fields[6] + "\" did not match any contig in the provided build.");
		}
		
		//PNEXT
		try
		{
			rec.pNext = Integer.parseInt(fields[7]);
		}
		catch (NumberFormatException e)
		{
			System.err.println("SAMRecord.parseSAMRecord || ERROR: PNEXT record (" + fields[7] + ") could not be parsed as a decimal integer!");
			System.err.println(samLine);
			throw new FileBuffer.UnsupportedFileTypeException();
		}
		
		//TLEN
		try
		{
			rec.templateLength = Integer.parseInt(fields[8]);
		}
		catch (NumberFormatException e)
		{
			System.err.println("SAMRecord.parseSAMRecord || ERROR: TLEN record (" + fields[8] + ") could not be parsed as a decimal integer!");
			System.err.println(samLine);
			throw new FileBuffer.UnsupportedFileTypeException();
		}
		
		//SEQ
		if (fields[9].equals("*")) rec.sequence = null;
		else rec.sequence = fields[9];
		
		//QUAL
		String qualstring = fields[10];
		if (qualstring.equals("*"))
		{
			System.err.println("SAMRecord.parseSAMRecord || WARNING: Phred scaled base call qualities not recorded for this read!");
			rec.seqQuality = null;
		}
		else
		{
			if (rec.sequence == null)
			{
				System.err.println("SAMRecord.parseSAMRecord || ERROR: Non-null quality string cannot accompany null sequence string!");
				System.err.println(samLine);
				throw new FileBuffer.UnsupportedFileTypeException();
			}
			int slen = rec.sequence.length();
			int qlen = qualstring.length();
			if (slen != qlen)
			{
				System.err.println("SAMRecord.parseSAMRecord || WARNING: Length of sequence string does NOT match length of quality string! (slen = " + slen + " | qlen = " + qlen + ")");
			}
			if (qlen == 0)
			{
				System.err.println("SAMRecord.parseSAMRecord || WARNING: Quality string is empty!");
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
			for (int j = 11; j < flen; j++)
			{
				SAMField ofield = parseOptionalField(fields[j]);
				rec.alignmentFields.put(ofield.getKey(), ofield);
			}
		}
		
		return rec;
	}
	
	public static SAMField parseOptionalField(String field) throws UnsupportedFileTypeException
	{
		if (field == null || field.isEmpty())
		{
			System.err.println("SAMRecord.parseOptionalField || ERROR: No field provided for parsing...");
			throw new FileBuffer.UnsupportedFileTypeException();
		}
		
		String[] fields = field.split(":");
		
		if (fields.length != 3)
		{
			System.err.println("SAMRecord.parseOptionalField || ERROR: Custom field is not formatted correctly: " + field);
			throw new FileBuffer.UnsupportedFileTypeException();
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
		
		System.err.println("SAMRecord.parseOptionalField || ERROR: Custom field type not recoginized: " + field);
		throw new FileBuffer.UnsupportedFileTypeException();
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
	
}
