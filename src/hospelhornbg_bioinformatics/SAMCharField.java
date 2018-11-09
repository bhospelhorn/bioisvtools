package hospelhornbg_bioinformatics;

public class SAMCharField implements SAMField {
	
	private String key;
	private char value;
	
	public SAMCharField(String TAG, String unparsedValue)
	{
		key = TAG;
		if (unparsedValue == null || unparsedValue.isEmpty()) value = '\0';
		else value = unparsedValue.charAt(0);
	}

	@Override
	public SAMFieldType getType() {
		return SAMFieldType.CHARACTER;
	}

	@Override
	public String getKey() {
		return key;
	}

	@Override
	public String getStringValue() {
		return Character.toString(value);
	}

	@Override
	public char getCharValue() {
		return value;
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
		return Character.toString(value).getBytes();
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
		return key + ":A:" + value;
	}

	@Override
	public int compareTo(SAMField other) 
	{
		return this.getSAMString().compareTo(other.getSAMString());
	}
	
}
