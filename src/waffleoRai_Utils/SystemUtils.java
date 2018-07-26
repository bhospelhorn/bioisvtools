package waffleoRai_Utils;

/*
 * UPDATES
 * 
 * 1.0.0 | January 18, 2018
 * 	Created
 * 	Architecture data model query written
 */

/**
 * A static method collection to simplify certain system operations.
 * @author Blythe Hospelhorn
 * @version 1.0.0
 * @since January 22, 2018
 */
public class SystemUtils {
	
	public static final int ARCH_32BIT = 32;
	public static final int ARCH_64BIT = 64;
	public static final int ARCH_UNKNOWN = 0;
	public static final int ARCH_UNSET = -1;
	
	private static int archDataModel = ARCH_UNSET;

	public static enum OperatingSystem
	{
		WINDOWS,
		OSX,
		LINUX;
	}
	
	public static int getArchModel()
	{
		if (archDataModel == ARCH_UNSET)
		{
			String dataModel = System.getProperty("sun.arch.data.model");
			if (dataModel.equals("32")) archDataModel = ARCH_32BIT;
			else if (dataModel.equals("64")) archDataModel = ARCH_64BIT;
			else archDataModel = ARCH_UNKNOWN;
		}
		return archDataModel;
	}
	
	public static int approximatePointerSize()
	{
		int arch = getArchModel();
		if (arch == ARCH_32BIT) return 4;
		if (arch == ARCH_64BIT) return 8;
		return 4;
	}
	
}
