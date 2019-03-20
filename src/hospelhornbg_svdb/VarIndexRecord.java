package hospelhornbg_svdb;

import hospelhornbg_bioinformatics.SVType;

public class VarIndexRecord {
	
	private SVType svtype;
	private int line;
	
	public VarIndexRecord(SVType type, int lineNumber)
	{
		svtype = type;
		line = lineNumber;
	}
	
	public SVType getType()
	{
		return svtype;
	}
	
	public int getLine()
	{
		return line;
	}

}
