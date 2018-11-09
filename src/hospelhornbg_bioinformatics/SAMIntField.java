package hospelhornbg_bioinformatics;

import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class SAMIntField implements SAMField {
	
	private String key;
	private int value;
	
	public SAMIntField(String TAG, String unparsedValue) throws UnsupportedFileTypeException
	{
		key = TAG;
		try
		{
			value = Integer.parseInt(unparsedValue);
		}
		catch(NumberFormatException e)
		{
			System.err.println("SAMIntField.<init> || ERROR: Custom field value (" + unparsedValue + ") for key " + TAG + " could not be parsed as a decimal integer!");
			throw new FileBuffer.UnsupportedFileTypeException();
		}
	}

	@Override
	public SAMFieldType getType() {
		return SAMFieldType.SIGNED_INT;
	}

	@Override
	public String getKey() {
		return key;
	}

	@Override
	public String getStringValue() {
		return Integer.toString(value);
	}

	@Override
	public char getCharValue() {
		return (char)value;
	}

	@Override
	public int getIntValue() {
		return value;
	}

	@Override
	public double getDoubleValue() {
		return (double)value;
	}

	@Override
	public byte[] getByteArrayValue() {
		//Big endian because Java
		byte[] ibytes = new byte[4];
		ibytes[0] = (byte)(value >>> 24);
		ibytes[1] = (byte)((value >>> 16) & 0xFF);
		ibytes[2] = (byte)((value >>> 8) & 0xFF);
		ibytes[3] = (byte)(value & 0xFF);
		return ibytes;
	}

	@Override
	public int[] getIntArrayValue() {
		int[] iarr = {value};
		return iarr;
	}

	public double[] getDoubleArrayValue()
	{
		double[] darr = {(double)value};
		return darr;
	}
	
	public String getSAMString()
	{
		return key + ":i:" + value;
	}
	
	public int compareTo(SAMField other) 
	{
		return this.getSAMString().compareTo(other.getSAMString());
	}

}
