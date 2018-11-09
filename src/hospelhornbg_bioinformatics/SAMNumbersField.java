package hospelhornbg_bioinformatics;

import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

public class SAMNumbersField implements SAMField {
	
	private String key;
	private int[] ivalue;
	private double[] fvalue;
	
	private char type;
	
	public SAMNumbersField(String TAG, String unparsedValue) throws UnsupportedFileTypeException
	{
		key = TAG;
		if (unparsedValue == null || unparsedValue.isEmpty()) throw new FileBuffer.UnsupportedFileTypeException();
		
		String[] fields = unparsedValue.split(",");
		char arrtype = fields[0].charAt(0);
		type = arrtype;
		if (arrtype == 'f')
		{
			//Float
			fvalue = new double[fields.length - 1];
			for (int i = 1; i < fields.length; i++)
			{
				try
				{
					fvalue[i-1] = Double.parseDouble(fields[i]);
				}
				catch(NumberFormatException e)
				{
					System.err.println("SAMBytesField.<init> || ERROR: Custom field value (" + unparsedValue + ") for key " + TAG + " could not be parsed as an array of floats!");
					throw new FileBuffer.UnsupportedFileTypeException();
				}
			}
		}
		else if (arrtype == 'c' || arrtype == 's' || arrtype == 'i')
		{
			//Signed
			ivalue = new int[fields.length - 1];
			for (int i = 1; i < fields.length; i++)
			{
				try
				{
					ivalue[i-1] = Integer.parseInt(fields[i]);
				}
				catch(NumberFormatException e)
				{
					System.err.println("SAMBytesField.<init> || ERROR: Custom field value (" + unparsedValue + ") for key " + TAG + " could not be parsed as an array of signed ints!");
					throw new FileBuffer.UnsupportedFileTypeException();
				}
			}
		}
		else if (arrtype == 'C' || arrtype == 'S' || arrtype == 'I')
		{
			//Unsigned
			ivalue = new int[fields.length - 1];
			for (int i = 1; i < fields.length; i++)
			{
				try
				{
					ivalue[i-1] = Integer.parseUnsignedInt(fields[i]);
				}
				catch(NumberFormatException e)
				{
					System.err.println("SAMBytesField.<init> || ERROR: Custom field value (" + unparsedValue + ") for key " + TAG + " could not be parsed as an array of unsigned ints!");
					throw new FileBuffer.UnsupportedFileTypeException();
				}
			}
		}
		else throw new FileBuffer.UnsupportedFileTypeException();
		
	}

	@Override
	public SAMFieldType getType() {
		if(ivalue != null) return SAMFieldType.INT_ARRAY;
		else return SAMFieldType.FLOAT_ARRAY;
	}

	@Override
	public String getKey() {
		return key;
	}

	@Override
	public String getStringValue() {
		String s = "";
		if (ivalue != null)
		{
			for (int i = 0; i < ivalue.length; i++)
			{
				s += ivalue[i];
				if (i < ivalue.length - 1) s += ",";
			}	
		}
		else
		{
			for (int i = 0; i < fvalue.length; i++)
			{
				s += fvalue[i];
				if (i < fvalue.length - 1) s += ",";
			}	
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
		throw new UnsupportedOperationException();
	}

	@Override
	public int[] getIntArrayValue() {
		return ivalue;
	}
	
	public double[] getDoubleArrayValue()
	{
		return fvalue;
	}
	
	public String getSAMString()
	{
		String s = key + ":B:" + type;
		if (ivalue != null)
		{
			if (Character.isUpperCase(type))
			{
				for (int i = 0; i < ivalue.length; i++)
				{
					s += "," + Integer.toUnsignedString(ivalue[i]);
				}	
			}
			else
			{
				for (int i = 0; i < ivalue.length; i++)
				{
					s += "," + ivalue[i];
				}	
			}
		}
		else
		{
			for (int i = 0; i < fvalue.length; i++)
			{
				s += "," + fvalue[i];
			}	
		}
		return s;
	}

	public int compareTo(SAMField other) 
	{
		return this.getSAMString().compareTo(other.getSAMString());
	}

}
