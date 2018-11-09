package hospelhornbg_bioinformatics;

import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class SAMBytesField implements SAMField {
	
	private String key;
	private byte[] value;
	
	public SAMBytesField(String TAG, String unparsedValue) throws UnsupportedFileTypeException
	{
		key = TAG;
		if (unparsedValue == null || unparsedValue.isEmpty()) throw new FileBuffer.UnsupportedFileTypeException();
		int slen = unparsedValue.length();
		if (slen % 2 != 0) throw new FileBuffer.UnsupportedFileTypeException();
		
		value = new byte[slen/2];
		int vpos = 0;
		int spos = 0;
		while (spos < slen)
		{
			String bstr = unparsedValue.substring(spos, spos+2);
			spos += 2;
			try
			{
				value[vpos] = (byte)Integer.parseInt(bstr, 16);
			}
			catch(NumberFormatException e)
			{
				System.err.println("SAMBytesField.<init> || ERROR: Custom field value (" + unparsedValue + ") for key " + TAG + " could not be parsed as an array of hex bytes!");
				throw new FileBuffer.UnsupportedFileTypeException();
			}
			vpos++;
		}
		
	}

	@Override
	public SAMFieldType getType() {
		return SAMFieldType.HEX_BYTE_ARRAY;
	}

	@Override
	public String getKey() {
		return key;
	}

	@Override
	public String getStringValue() {
		String s = "";
		for (int i = 0; i < value.length; i++)
		{
			s += String.format("%02X ", value[i]);
		}
		return s;
	}

	@Override
	public char getCharValue() {
		throw new UnsupportedOperationException();
	}

	@Override
	public int getIntValue() {
		throw new UnsupportedOperationException();
	}

	@Override
	public double getDoubleValue() {
		throw new UnsupportedOperationException();
	}

	@Override
	public byte[] getByteArrayValue() {
		//Big endian because Java
		return value;
	}

	@Override
	public int[] getIntArrayValue() {
		int[] iarr = new int[value.length];
		for (int i = 0; i < value.length; i++)
		{
			iarr[i] = Byte.toUnsignedInt(value[i]);
		}
		return iarr;
	}
	
	public double[] getDoubleArrayValue()
	{
		double[] darr = new double[value.length];
		for (int i = 0; i < value.length; i++)
		{
			darr[i] = (double)Byte.toUnsignedInt(value[i]);
		}
		return darr;
	}
	
	public String getSAMString()
	{
		String s = "";
		for (int i = 0; i < value.length; i++)
		{
			s += String.format("%02X", value[i]);
		}
		return key + ":H:" + s;
	}

	public int compareTo(SAMField other) 
	{
		return this.getSAMString().compareTo(other.getSAMString());
	}

}
