package hospelhornbg_bioinformatics;

public interface SAMField extends Comparable<SAMField> {
	
	public SAMFieldType getType();
	public String getKey();
	
	public String getStringValue();
	public char getCharValue();
	public int getIntValue();
	public double getDoubleValue();
	public byte[] getByteArrayValue();
	public int[] getIntArrayValue();
	public double[] getDoubleArrayValue();
	
	public String getSAMString();

}
