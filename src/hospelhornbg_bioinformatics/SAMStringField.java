package hospelhornbg_bioinformatics;

import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class SAMStringField implements SAMField {
	
	private String key;
	private String value;
	
	public SAMStringField(String TAG, String unparsedValue) throws UnsupportedFileTypeException
	{
		key = TAG;
		value = unparsedValue;
	}

	@Override
	public SAMFieldType getType() {
		return SAMFieldType.STRING;
	}

	@Override
	public String getKey() {
		return key;
	}

	@Override
	public String getStringValue() {
		return value;
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
		return value.getBytes();
	}
	
	public double[] getDoubleArrayValue()
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public int[] getIntArrayValue() {
		byte[] bytes = value.getBytes();
		int[] iarr = new int[bytes.length];
		for (int i = 0; i < iarr.length; i++)
		{
			iarr[i] = Byte.toUnsignedInt(bytes[i]);
		}
		
		return iarr;
	}

	public String getSAMString()
	{
		return key + ":Z:" + value;
	}
	
	public int compareTo(SAMField other) 
	{
		return this.getSAMString().compareTo(other.getSAMString());
	}

}
