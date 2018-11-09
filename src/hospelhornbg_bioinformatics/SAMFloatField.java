package hospelhornbg_bioinformatics;

import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class SAMFloatField implements SAMField{
	
	private String key;
	private double value;
	
	public SAMFloatField(String TAG, String unparsedValue) throws UnsupportedFileTypeException
	{
		key = TAG;
		try
		{
			value = Double.parseDouble(unparsedValue);
		}
		catch(NumberFormatException e)
		{
			System.err.println("SAMIntField.<init> || ERROR: Custom field value (" + unparsedValue + ") for key " + TAG + " could not be parsed as a decimal integer!");
			throw new FileBuffer.UnsupportedFileTypeException();
		}
	}

	@Override
	public SAMFieldType getType() {
		return SAMFieldType.SINGLE_FLOAT;
	}

	@Override
	public String getKey() {
		return key;
	}

	@Override
	public String getStringValue() {
		return Double.toString(value);
	}

	@Override
	public char getCharValue() {
		return (char)value;
	}

	@Override
	public int getIntValue() {
		return (int)Math.round(value);
	}

	@Override
	public double getDoubleValue() {
		return value;
	}

	@Override
	public byte[] getByteArrayValue() {
		long raw = Double.doubleToRawLongBits(value);
		byte[] ibytes = new byte[8];
		ibytes[0] = (byte)(raw >>> 56);
		ibytes[1] = (byte)((raw >>> 48) & 0xFF);
		ibytes[2] = (byte)((raw >>> 40) & 0xFF);
		ibytes[3] = (byte)((raw >>> 32) & 0xFF);
		ibytes[4] = (byte)((raw >>> 24) & 0xFF);
		ibytes[5] = (byte)((raw >>> 16) & 0xFF);
		ibytes[6] = (byte)((raw >>> 8) & 0xFF);
		ibytes[7] = (byte)(raw & 0xFF);
		return ibytes;
	}

	@Override
	public int[] getIntArrayValue() {
		int[] iarr = {this.getIntValue()};
		return iarr;
	}
	
	public double[] getDoubleArrayValue()
	{
		double[] darr = {value};
		return darr;
	}

	public String getSAMString()
	{
		return key + ":f:" + value;
	}
	
	public int compareTo(SAMField other) 
	{
		return this.getSAMString().compareTo(other.getSAMString());
	}

}
