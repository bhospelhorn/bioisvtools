package waffleoRai_Utils;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.CharBuffer;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Calendar;
import java.util.Collection;
import java.util.GregorianCalendar;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.TimeZone;



/*UPDATES
 * 2017.08.31
 * 	3.2 -> 3.2.1 | Added getMemoryBurden() function
 * 2017.09.01
 * 	3.2.1 -> 3.3.0 | Added read-only locking, documentation, java file formatting
 * 2017.09.08
 * 	3.3.0 -> 3.3.1 | Tidied up javadoc documentation. Added better referencing child management.
 * 		Switched default overload in several cases to position long overload instead of int.
 * 		Cleaned up several methods.
 * 		Deleted unsigned converters (redundant to built in Java functions) except byte->short
 * 2017.09.20
 *  3.3.1 -> 3.3.2 | Documentation minor update. Throw IndexOutOfBounds instead of ArrayIndexOutOfBounds
 * 		Minor code tidying.
 * 2017.09.24
 * 	3.3.2 -> 3.3.3 | Added framework for generating temporary paths.
 * 2017.10.14
 * 	3.3.3 -> 3.3.4 | Very minor change - added a method for getting the temp dir
 * 2017.10.14
 * 	3.3.4 -> 3.3.5 | Altered string finder functions so that if edPos is larger than file size, it is simply set to file size.
 * 2017.10.21
 * 	3.3.5 -> 3.4.0 | Added functionality for encoded strings
 * 2017.11.01
 * 	3.4.0 -> 3.4.1 | Added default date (GregorianCalendar)
 * 2017.11.04
 * 	3.4.1 -> 3.4.2 | Allow changing of size threshold.
 * 2018.01.22
 * 	3.4.2 -> 3.4.3 | Added memory approximation method.
 * 2018.01.25
 * 	3.4.3 -> 3.4.4 | Added replace method for 3 byte values
 * 2018.02.03
 * 	3.4.4 -> 3.4.5 | Fixed a completely asinine bug in the byte-order switcher for 32-bit values
 * 		Turns out the thing was setting only the LSB as all bytes (only reading byte[0]) 
 * 		Good job, you lazy copypasting dumbass.
 * 		Also fixed serious bug with Little-Endian short addition
 * 2018.04.15
 * 	3.4.5 -> 3.4.6 | Byte buffer conversions weren't being rewound. Oopsy!
 * 2018.06.21
 * 	3.4.6 -> 3.4.7 | Added a static function to write an InputStream to disk.
 * */

/**
 * Memory buffer internally implemented as a byte array to hold raw file data.
 * Includes a number of methods specific to file processing and disk I/O.
 * <br> Can handle capacity overflow.
 * <br> Standard FileBuffer class does not function via a continuous stream - it reads
 * everything it is asked to into memory at once.
 * <br> Due to byte array and byte buffer conversion procedures, maximum capacity and file size cannot exceed
 * 0x7FFFFFFF (~2GB) at a time, even with overflow.
 * @author Blythe Hospelhorn
 * @version 3.4.7*
 * @since June 21, 2018
 */
public class FileBuffer 
{
	
	public static final long DEFO_SIZE_THRESHOLD = 0x40000000; //1GB
	private static long SIZE_THRESHOLD = DEFO_SIZE_THRESHOLD;
	private static String tempDir;
	
	private byte[] contents;
	private int fSize;
	private int capacity;
	private boolean isFileFormatBE;
	private String directory;
	private String fileName;
	private String extension;
  
	private List<byte[]> overflow;
  
	private boolean readOnly;
	private List<ROSubFileBuffer> children;
	
  
	/* ----- CONSTRUCTORS ----- */
  
	/**
	 * Constructor #1:
	 * Reads in a full file (as a byte string) to the buffer.
	 * @param fPath - String indicating path of file to open and copy into memory.
	 * @throws IOException If file at specified path cannot be opened.
	 */
	public FileBuffer(String fPath) throws IOException
	{
		/*ARGS CHECK*/
		if (fPath == null || !fileExists(fPath)) throw new IOException();
	  
		/*BODY*/
		this.newReadCore(fPath, 0, FileBuffer.fileSize(fPath), true);
		this.readIn(fPath);
	}
  
	/**
	 * Constructor #2:
	 * Reads in file - everything after given offset into memory buffer
	 * @param fPath - String indicating path of file to open and copy into memory.
	 * @param stOff - Starting offset within file on disk (inclusive)
	 * @throws IOException If file at specified path cannot be opened.
	 */
	public FileBuffer(String fPath, long stOff) throws IOException
	{
		/*ARGS CHECK*/
		if (fPath == null || !fileExists(fPath) || stOff >= fileSize(fPath)) throw new IOException();
		if (stOff < 0) stOff = 0;
		  
		/*BODY*/
		this.newReadCore(fPath, stOff, FileBuffer.fileSize(fPath), true);
		this.readIn(fPath, stOff); 
  }
  
	/**
	 * Constructor #3:
	 * Reads in file - between given offset and given end offset
	 * @param fPath - String indicating path of file to open and copy into memory.
	 * @param stOff - Starting offset within file on disk (inclusive)
	 * @param edOff - Offset of the end of read chunk (exclusive)
	 * @throws IOException If file at specified path cannot be opened.
	 */
	public FileBuffer(String fPath, long stOff, long edOff) throws IOException
	{
		  /*ARGS CHECK*/
		  if (fPath == null || !fileExists(fPath) || edOff <= stOff) throw new IOException();
		  if (stOff < 0) stOff = 0;
		  if (edOff > fileSize(fPath)) edOff = (int)fileSize(fPath);
		  
		  /*BODY*/
		  this.newReadCore(fPath, stOff, edOff, true);
		  this.readIn(fPath, stOff, edOff);	  
	}
  
	/**
	 * Constructor #4:
	 * Reads in a full file (as a byte string) to the buffer
	 * @param fPath - String indicating path of file to open and copy into memory.
	 * @param isBigEndian - Whether or not multi-byte values should be read from the file BE
	 * @throws IOException If file at specified path cannot be opened.
	 */
	public FileBuffer(String fPath, boolean isBigEndian) throws IOException
	{
		/*ARGS CHECK*/
		if (fPath == null || !fileExists(fPath)) throw new IOException();
	  
		/*BODY*/
		this.newReadCore(fPath, 0, FileBuffer.fileSize(fPath), isBigEndian);
		this.readIn(fPath);
	}
  
	/**
	 * Constructor #5:
	 * Reads in a full file (as a byte string) to the buffer
	 * @param fPath - String indicating path of file to open and copy into memory.
	 * @param stOff - Starting offset within file on disk (inclusive)
	 * @param isBigEndian - Whether or not multi-byte values should be read from the file BE
	 * @throws IOException If file at specified path cannot be opened.
	 */
	public FileBuffer(String fPath, long stOff, boolean isBigEndian) throws IOException
	{
		/*ARGS CHECK*/
		if (fPath == null || !fileExists(fPath) || stOff >= fileSize(fPath)) throw new IOException();
		if (stOff < 0) stOff = 0;
		  
		/*BODY*/
		this.newReadCore(fPath, stOff, FileBuffer.fileSize(fPath), isBigEndian);
		this.readIn(fPath, stOff); 
	}

	/**
	 * Constructor #6:
	 * Reads in a full file (as a byte string) to the buffer
	 * @param fPath - String indicating path of file to open and copy into memory.
	 * @param stOff - Starting offset within file on disk (inclusive)
	 * @param edOff - Offset of the end of read chunk (exclusive)
	 * @param isBigEndian - Whether or not multi-byte values should be read from the file BE
	 * @throws IOException If file at specified path cannot be opened.
	 */
	public FileBuffer(String fPath, long stOff, long edOff, boolean isBigEndian) throws IOException
	{
		  /*ARGS CHECK*/
		  if (fPath == null || !fileExists(fPath) || edOff <= stOff) throw new IOException();
		  if (stOff < 0) stOff = 0;
		  if (edOff > fileSize(fPath)) edOff = (int)fileSize(fPath);
		  
		  /*BODY*/
		  this.newReadCore(fPath, stOff, edOff, isBigEndian);
		  this.readIn(fPath, stOff, edOff);	 
	}
  
	/**
	 * Constructor #7:
	 * Allocate memory for a new empty FileBuffer.
	 * Intended for something waiting to be written to disk.
	 * Cannot exceed 2GB. Strongly not recommended for large files.
	 * @param initialSize - Amount of memory to allocate. Must be a positive signed integer.
	 * @throws IllegalArgumentException If initialSize is invalid
	 */
	public FileBuffer(int initialSize)
	{
		this(initialSize, true);
	}
  
	/**
	 * Constructor #8:
	 * Allocate memory for a new empty FileBuffer with specified Endian-ness.
	 * Intended for something waiting to be written to disk.
	 * Cannot exceed 2GB. Strongly not recommended for large files.
	 * @param initialSize - Amount of memory to allocate. Must be a positive signed integer.
	 * @throws IllegalArgumentException If initialSize is invalid
	 */
	public FileBuffer(int initialSize, boolean isBigEndian)
	{
		 if (initialSize < 0 || initialSize > 2147483647) throw new IllegalArgumentException();
	  
		 this.newWriteCore(initialSize, isBigEndian);
	}

	/**
	 * Constructor #9:
	 * Create a new FileBuffer from an existing byte array.
	 * IMPORTANT!! It COPIES the contents of the byte array into the buffer. It does NOT reference
	 * the byte array!
	 * @param myBytes - Byte array to copy into new buffer
	 */
	public FileBuffer(byte[] myBytes)
	{
		this(myBytes, true);
	}
  
	/**
	 * Constructor #10:
	 * Create a new FileBuffer from an existing byte array, and specify the Endian-ness.
	 * IMPORTANT!! It COPIES the contents of the byte array into the buffer. It does NOT reference
	 * the byte array!
	 * @param myBytes - Byte array to copy into new buffer
	 * @param isBigEndian - Whether or not multi-byte values should be read from the file BE
	 */
	public FileBuffer(byte[] myBytes, boolean isBigEndian)
	{	  
		  this.newWriteCore(myBytes.length, true);
		  for (int i = 0; i < myBytes.length; i++) this.addToFile(myBytes[i]);
	}
  
	/**
	 * Constructor #11:
	 * FOR USE BY CHILD CLASSES ONLY!!
	 * So that the "super" parent FileBuffer has the same data as the FB it was intended to reference.
	 * Makes a copy referencing contents of the first.
	 * @param source : Source file to reference
	 */
	protected FileBuffer(FileBuffer source)
	{
		//For creating a read-only referencing FB.
		//DOES reference byte array pulled from source.
	  
		this.contents = source.contents;
		this.fSize = source.fSize;
		this.capacity = source.capacity;
		this.isFileFormatBE = source.isFileFormatBE;
		this.directory = source.directory;
		this.fileName = source.fileName;
		this.extension = source.extension;
		this.overflow = source.overflow;
	}
  
	/**
	 * Constructor #12:
	 * FOR USE BY CHILD CLASSES ONLY!! (Particularly StreamBuffer objects)
	 * Creates practically empty FileBuffer 
	 */
	protected FileBuffer()
	{
		this.setEmpty();
	}
  
	private void newReadCore(String fPath, long stOff, long edOff, boolean isBigEndian)
	{
		  this.capacity = (int)(edOff - stOff);
		  this.contents = new byte[this.capacity];
		  this.overflow = new LinkedList<byte[]>();
		  this.children = new LinkedList<ROSubFileBuffer>();
		  //File size is set in readIn
		  this.directory = chopPathToDir(fPath);
		  this.fileName = chopPathToFName(fPath);
		  this.extension = chopPathToExt(fPath);
		  this.isFileFormatBE = isBigEndian;
		  this.readOnly = false;
		  this.fSize = 0;
	}
  
	private void newWriteCore(int initialCap, boolean isBigEndian)
	{
		this.capacity = initialCap;
		this.fSize = 0;
		this.isFileFormatBE = isBigEndian;
		this.directory = null;
		this.fileName = null;
		this.extension = null;
		this.contents = new byte[this.capacity];
		this.overflow = new LinkedList<byte[]>();
		this.children = new LinkedList<ROSubFileBuffer>();
		this.readOnly = false;
	}
  
	private void setEmpty()
	{
		this.contents = null;
		this.capacity = 0;
		this.fSize = 0;
		this.isFileFormatBE = true;
		this.directory = null;
		this.fileName = null;
		this.extension = null;
		this.overflow = null;
		this.readOnly = true;
		this.children = new LinkedList<ROSubFileBuffer>();
	}
  
	/* ----- STATIC CREATORS ----- */
  
	/**
	 * Creates a FileBuffer or a StreamBuffer for data I/O.
	 * For files under 1GB, it will create and return a regular FileBuffer (load full file
	 * into memory).
	 * For files over 1GB, it will create a StreamBuffer.
	 * @param filename - Path of file to read
	 * @return FileBuffer for easy access to a file's contents either through stream or by loading into memory.
	 * @throws IOException If path is invalid or file cannot be opened
	 */
	public static FileBuffer createBuffer(String filename) throws IOException
	{
		return createBuffer(filename, 0, FileBuffer.fileSize(filename), true);
	}
  
	/**
	 * Creates a FileBuffer or a StreamBuffer for data I/O.
	 * For files under 1GB, it will create and return a regular FileBuffer (load full file
	 * into memory).
	 * For files over 1GB, it will create a StreamBuffer.
	 * @param filename - Path of file to read
	 * @param stOff - Start offset of file to read in
	 * @return FileBuffer for easy access to a file's contents either through stream or by loading into memory.
	 * @throws IOException If path is invalid or file cannot be opened
	 */
	public static FileBuffer createBuffer(String filename, long stOff) throws IOException
	{
		return createBuffer(filename, stOff, FileBuffer.fileSize(filename), true);
	}
	
	/**
	 * Creates a FileBuffer or a StreamBuffer for data I/O.
	 * For files under 1GB, it will create and return a regular FileBuffer (load full file
	 * into memory).
	 * For files over 1GB, it will create a StreamBuffer.
	 * @param filename - Path of file to read
	 * @param stOff - Start offset of file to read in
	 * @param edOff - End offset (exclusive) of file piece to read in.
	 * @return FileBuffer for easy access to a file's contents either through stream or by loading into memory.
	 * @throws IOException If path is invalid or file cannot be opened
	 */
	public static FileBuffer createBuffer(String filename, long stOff, long edOff) throws IOException
	{
		return createBuffer(filename, stOff, edOff, true);
	}
  
	/**
	 * Creates a FileBuffer or a StreamBuffer for data I/O.
	 * For files under 1GB, it will create and return a regular FileBuffer (load full file
	 * into memory).
	 * For files over 1GB, it will create a StreamBuffer.
	 * @param filename - Path of file to read
	 * @param isBE - Whether multi-byte values in the input file should be read in Big-Endian order
	 * @return FileBuffer for easy access to a file's contents either through stream or by loading into memory.
	 * @throws IOException If path is invalid or file cannot be opened
	 */
	public static FileBuffer createBuffer(String filename, boolean isBE) throws IOException
	{
		return createBuffer(filename, 0, FileBuffer.fileSize(filename), isBE);
	}
  
	/**
	 * Creates a FileBuffer or a StreamBuffer for data I/O.
	 * For files under 1GB, it will create and return a regular FileBuffer (load full file
	 * into memory).
	 * For files over 1GB, it will create a StreamBuffer.
	 * @param filename - Path of file to read
	 * @param stOff - Start offset of file to read in
	 * @param isBE - Whether multi-byte values in the input file should be read in Big-Endian order
	 * @return FileBuffer for easy access to a file's contents either through stream or by loading into memory.
	 * @throws IOException If path is invalid or file cannot be opened
	 */
	public static FileBuffer createBuffer(String filename, long stOff, boolean isBE) throws IOException
	{
		return createBuffer(filename, stOff, FileBuffer.fileSize(filename), isBE);
	}
  
	/**
	 * Creates a FileBuffer or a StreamBuffer for data I/O.
	 * For files under 1GB, it will create and return a regular FileBuffer (load full file
	 * into memory).
	 * For files over 1GB, it will create a StreamBuffer.
	 * @param filename - Path of file to read
	 * @param stOff - Start offset of file to read in
	 * @param edOff - End offset (exclusive) of file piece to read in.
	 * @param isBE - Whether multi-byte values in the input file should be read in Big-Endian order
	 * @return FileBuffer for easy access to a file's contents either through stream or by loading into memory.
	 * @throws IOException If path is invalid or file cannot be opened
	 */
	public static FileBuffer createBuffer(String filename, long stOff, long edOff, boolean isBE) throws IOException
	{
		/*Evaluates the size of the incoming file and creates an object appropriate
		 * for how much should be held in memory.
		 */
		//long fSz = FileBuffer.fileSize(filename);
		long rSz = edOff - stOff;
		if (rSz <= SIZE_THRESHOLD)
		{
			return new FileBuffer(filename, stOff, edOff, isBE);
		}
		else
		{
			return new StreamBuffer(filename, stOff, edOff, isBE);
		}	  
	}
  
	/**
	 * Create a FileBuffer of a StreamBuffer for writing.
	 * @param tempName Name stem for temporary file if StreamBuffer is required.
	 * @param size Proposed size of file
	 * @param isBE Whether the byte order of file is Big-Endian
	 * @return FileBuffer for filling.
	 * @throws IOException If there is an error creating a StreamBuffer.
	 */
	public static FileBuffer createWritableBuffer(String tempName, long size, boolean isBE) throws IOException
	{
		if (size < FileBuffer.SIZE_THRESHOLD) return new FileBuffer((int)size, isBE);
		else return new StreamBuffer(tempName, size, isBE, true);
	}
	
	/**
	 * Set the threshold for automatic buffer creation - files on disk smaller than the threshold
	 * will be read fully into memory.
	 * Files on disk larger than the threshold will be streamed from disk.
	 * This value cannot be negative and cannot exceed the highest signed positive value representable in 32 bits.
	 * (~2 GB or 0x7FFFFFFF)
	 * @param newThreshold Value to set the threshold at.
	 */
	public static void setMemoryThreshold(int newThreshold)
	{
		if (newThreshold < 0) return;
		if (newThreshold > 0x7FFFFFFF) return;
		SIZE_THRESHOLD = newThreshold;
	}
	
	/**
	 * Get the current size threshold for automatic buffer creation.
	 * Files on disk smaller than the threshold
	 * will be read fully into memory.
	 * Files on disk larger than the threshold will be streamed from disk.
	 * This value cannot be negative and cannot exceed the highest signed positive value representable in 32 bits.
	 * (~2 GB or 0x7FFFFFFF)
	 * @return A value between 0 and 0x7FFFFFFF representing the largest size a file can be to be
	 * automatically loaded into memory with the createBuffer static function.
	 */
	public static long getCurrentMemoryThreshold()
	{
		return SIZE_THRESHOLD;
	}
	
  /* ----- BASIC GETTERS ----- */
  
	/**
	 * Get the size in bytes of the file that the buffer contents would write.
	 * @return Long integer representing the number of occupied bytes in this buffer. -1 if there is an error.
	 */
	public long getFileSize()
	{
		return Integer.toUnsignedLong(this.fSize);
	}
  
	/**
	 * Get the capacity of the primary byte array (and any overflow arrays) in this buffer.
	 * @return Base capacity of FileBuffer.
	 */
	public long getBaseCapacity()
	{
		return Integer.toUnsignedLong(this.capacity);
	}
  
	/**
	 * Get whether multi-byte values in buffer are set to be read and written in Big-Endian order.
	 * @return true if Big-Endian, false if Little-Endian
	 */
	public boolean isBigEndian()
	{
		return this.isFileFormatBE;
	}
  
	/**
	 * Get a byte from this file at the provided position.
	 * @param position - Offset from FileBuffer start to get byte from.
	 * @return byte at position specified
	 * @throws IndexOutOfBoundsException If position is invalid.
	 */
	public byte getByte(int position)
	{
		if (position < 0 || position >= (int)this.getFileSize()) throw new IndexOutOfBoundsException();
		if (this.isOverflowing())
		{
			if (position >= this.capacity)
			{
				int oArr = position / this.capacity;
				oArr--;
				if (oArr >= this.overflow.size()) throw new IndexOutOfBoundsException();
				int oi = position % this.capacity;
				byte[] arr = this.overflow.get(oArr);
				if (oi < 0 || oi >= arr.length) throw new IndexOutOfBoundsException();
				byte ob = arr[oi];
				return ob;
			}
		}
		return this.contents[position];
	}
  
	/**
	 * Get a byte from this file at the provided position.
	 * @param position Offset from FileBuffer start to get byte from.
	 * @return Byte at position specified
	 * @throws IndexOutOfBoundsException If position is invalid.
	 */
	public byte getByte(long position)
	{
		int pos = (int)position;
		return this.getByte(pos);
	}
  
	/**
	 * Get the directory path stored by this FileBuffer.
	 * @return Full directory path to file opened, null if buffer not created from existing file. 
	 */
	public String getDir()
	{
		return this.directory;
	}
  
	/**
	 * Get the file name without the path specified by the FileBuffer.
	 * @return Name of file opened, null if buffer not created from existing file. 
	 */
	public String getName()
	{
		return this.fileName;
	}
  
	/**
	 * Get the extension with no file name or path as specified by the FileBuffer.
	 * @return Extension of file opened, null if buffer not created from existing file. 
	 */
	public String getExt()
	{
		return this.extension;
	}

	/**
	 * Get the full path as specified by the FileBuffer.
	 * @return Path to file opened, null if buffer not created from existing file. 
	 */
	public String getPath()
	{
		return directory + File.separator + fileName + "." + extension;
	}
  
	/**
	 * Get the minimum amount of memory it takes to hold this buffer.
	 * @return Amount of memory buffer occupies.
	 */
	public long getMemoryBurden()
	{
		return getMinimumMemoryUsage();
	}
  
	/**
   	* Overload for more accurate calculation. If this FileBuffer references another,
   	* it can scan a collection of already counted FileBuffers for its parent.
   	* If its parent has already been counted, then it will return 0.
   	* If not, it will return the amount of memory it or its parent takes up in data.
   	* @param inMem Collection (eg. List) of FileBuffer references already counted
   	* @return 0 if parent is in list, otherwise allocated size of this or parent it references.
   	*/
	public long getMemoryBurden(Collection<FileBuffer> inMem)
	{
		return getMemoryBurden();
	}
	
  	/**
  	 * Approximate the memory usage, in bytes, of this buffer. This number
  	 * is a considered a minimum possible value. It only counts the size of all
  	 * instance variables and objects supposedly unique to this instance.
  	 * It does not count the potential size of any methods or method pointers. Note
  	 * that the FileBuffer class contains a lot of methods, and it is quite possible that
  	 * this alone may increase memory usage.
  	 * <br>Also note that some estimates for the size of unique internally referenced
  	 * objects may also come in a bit low. This calculation compounds that underestimation.
  	 * <br>The object may, in reality, take up more memory than this indicates. 
  	 * Obtaining a more accurate value from the JVM is somewhat difficult, and may
  	 * not be necessary.
  	 * <br><br>IMPORTANT: This function only returns the memory usage for this instance. If
  	 * this FileBuffer is a readonly buffer that references points in other buffers, but
  	 * does not contain a byte array or structure itself, the value will NOT include
  	 * the backing data; only the approximate size the pointers used to reference the parent
  	 * might be. If the backing buffer is lost or unknown, it might be more prudent to call
  	 * the ROSubFileBuffer class's getMemoryBurden method!
  	 * @return A size representing the minimum size, in bytes, that this instance likely
  	 * uses in memory.
  	 */
  	public long getMinimumMemoryUsage()
  	{
  		//I don't know that the architecture model determines pointer size, tbh
  		//It's just an estimation
  		int estPtrSz = SystemUtils.approximatePointerSize();
  		long tot = 0;
  		//Instance variables (standard)
  		tot += 4 + 4 + 1 + 1;
  		//Instance variables (references, guess size of ptrs)
  		tot += estPtrSz * 6;
  		//Strings
  		if (directory != null) tot += directory.length();
  		if (fileName != null) tot += fileName.length();
  		if (extension != null) tot += extension.length();
  		//Primary byte array
  		if (contents != null) tot += contents.length;
  		//Overflow
  		if (overflow != null && !overflow.isEmpty())
  		{
  			for (byte[] barr : overflow) tot += barr.length;
  		}
  		//References to children
  		if (children != null) tot += estPtrSz * children.size();
  		return tot;
  	}
	
	/**
	 * Checks against another FileBuffer to see if contents are equal.
	 * @param other - FileBuffer to check against.
	 * @return true: if contain the same core byte array and overflow list reference. false: otherwise. Will return false for FB child classes as these have null-filled super parts.
	 */
	protected boolean contentsEqual(FileBuffer other)
	{
		if (other == null) return false;
		if (this.contents == other.contents && this.overflow == other.overflow) return true;
		return false;
	}
   
  /* ----- BASIC SETTERS ----- */
  
	/**
	 * Set internal directory path
	 * @param in String representing the directory path linked to this FileBuffer.
	 */
	public void setDir(String in)
	{
		this.directory = in;
	}
  
	/**
	 * Set internal file name
	 * @param in String representing the file name (without path or extension) to link to FileBuffer.
	 */
	public void setName(String in)
	{
		this.fileName = in; 
	}
  
	/**
	 * Internally set file extension string.
	 * @param in String representing the file extension string (without '.') to link to buffer.
	 */
	public void setExt(String in)
	{
		this.extension = in;
	}
  
	/**
	 * Set byte order interpretation for reading and writing multi-byte values from buffer.
	 * @param isBE Whether target byte order is Big-Endian
	 */
	public void setEndian(boolean isBE)
	{
		this.isFileFormatBE = isBE;
	}
  
	/**
	 * Set the FileBuffer to read-only mode, if applicable.
	 * While in read-only mode, FileBuffer cannot be written to in any way.
	 */
	public void setReadOnly()
  	{
		this.readOnly = true;
  	}
  
	/**
	 * If FileBuffer is in read-only mode, will release it and once again allow for writing.
	 * If FileBuffer is not in read-only mode, this will do nothing.
	 */
	public void unsetReadOnly()
	{
		this.readOnly = false;
	}
  
  /* ----- READERS ----- */
  
	private void readIn(String fPath) throws IOException
	{
		this.readIn(fPath, 0, (int)fileSize(fPath));  
	}
  
	private void readIn(String fPath, long stOff) throws IOException
	{
		this.readIn(fPath, stOff, (int)fileSize(fPath));  
	}
  
	private void readIn(String fPath, long stOff, long edOff) throws IOException
	{
		FileInputStream myStream = new FileInputStream(fPath);
	    int rBytes = (int)(edOff - stOff);
	    myStream.skip(stOff);
	    myStream.read(this.contents, 0, rBytes);
	    this.fSize = this.contents.length;
	    myStream.close();
	}

  /* ----- BYTE ORDERING ----- */
  
	/**
	 * Minimal RuntimeException child class for marking errors with Byte Ordering in file parsing
	 * and serialization.
	 * @author Blythe Hospelhorn
	 * @version 1.0
	 */
	public static class ByteOrderException extends RuntimeException
	{
		private static final long serialVersionUID = 2959602591728281548L; 
	}
  
	/**
	 * Forcibly switch the byte order of a 32-bit value (Java int)
	 * @param i32 32-bit value to reverse byte order of.
	 * @return int with bytes in opposite order as input value 
	 * <br> (ie. B3 B2 B1 B0 -> B0 B1 B2 B3)
	 */
	public static int switchByO(int i32)
	{
		byte[] byStr = numToByStr(i32);
		int b0 = Byte.toUnsignedInt(byStr[0]) << 24;
		int b1 = Byte.toUnsignedInt(byStr[1]) << 16;
		int b2 = Byte.toUnsignedInt(byStr[2]) << 8;
		int b3 = Byte.toUnsignedInt(byStr[3]);
		int out = 0x00000000;
		
		out = b0 | b1 | b2 | b3;
		  
	    return out;
	}
  
	/**
	 * Forcibly switch the byte order of a 64-bit value (Java long)
	 * @param i64 64-bit value to reverse byte order of.
	 * @return long with bytes in opposite order as input value 
	 * <br>(ie. B7 B6 B5 B4 B3 B2 B1 B0 -> 
	 * <br>B0 B1 B2 B3 B4 B5 B6 B7)
	 */
	public static long switchByO(long i64)
	{
		byte[] byStr = numToByStr(i64);
		long b0 = Byte.toUnsignedLong(byStr[0]) << 56;
		long b1 = Byte.toUnsignedLong(byStr[1]) << 48;
		long b2 = Byte.toUnsignedLong(byStr[2]) << 40;
		long b3 = Byte.toUnsignedLong(byStr[3]) << 32;
		long b4 = Byte.toUnsignedLong(byStr[4]) << 24;
		long b5 = Byte.toUnsignedLong(byStr[5]) << 16;
		long b6 = Byte.toUnsignedLong(byStr[6]) << 8;
		long b7 = Byte.toUnsignedLong(byStr[7]);
		long out = 0x0000000000000000;
		
		out = b0 | b1 | b2 | b3 | b4 | b5 | b6 | b7;
	  
		return out;
  }
  
	/**
	 * Generate a byte array representation of a multi-byte value.
	 * @param myNum Multi-byte number to break into bytes
	 * @return Byte array containing individual bytes of the input value. Lowest index corresponds to least significant byte.
	 */
	public static byte[] numToByStr(short myNum)
	{
		byte[] out = new byte[2];
	  
		int b0 = 0;
		int b1 = 0;
		b0 = (int)myNum & 0x00FF;
		b1 = (((int)myNum & 0xFF00) >> 8) & 0xFF;
	  
		out[0] = (byte)b0;
		out[1] = (byte)b1;
	  
		return out;
	}
  
	/**
	 * Generate a byte array representation of a multi-byte value.
	 * @param myNum Multi-byte number to break into bytes
	 * @return Byte array containing individual bytes of the input value. Lowest index corresponds to least significant byte.
	 */
	public static byte[] numToByStr(int myNum)
	{
		byte[] out = new byte[4];
	  
		int b0 = 0;
		int b1 = 0;
		int b2 = 0;
		int b3 = 0;
	  	b0 = myNum & 0x000000FF;
	  	b1 = ((myNum & 0x0000FF00) >> 8) & 0xFF;
	  	b2 = ((myNum & 0x00FF0000) >> 16) & 0xFF;
	  	b3 = ((myNum & 0xFF000000) >> 24) & 0xFF;
	  
	  	out[0] = (byte)b0;
	  	out[1] = (byte)b1;
	  	out[2] = (byte)b2;
	  	out[3] = (byte)b3;
	  
	  	return out;	  
	}
  
	/**
	 * Generate a byte array representation of a multi-byte value.
	 * @param myNum Multi-byte number to break into bytes
	 * @return Byte array containing individual bytes of the input value. Lowest index corresponds to least significant byte.
	 */
	public static byte[] numToByStr(long myNum)
	{
		byte[] out = new byte[8];
	  
		long b0 = 0;
		long b1 = 0;
		long b2 = 0;
		long b3 = 0;
		long b4 = 0;
		long b5 = 0;
		long b6 = 0;
		long b7 = 0;
		b0 = myNum & 0x00000000000000FFL;
		b1 = ((myNum & 0x000000000000FF00L) >> 8) & 0x00000000000000FFL;
		b2 = ((myNum & 0x0000000000FF0000L) >> 16) & 0x00000000000000FFL;
		b3 = ((myNum & 0x00000000FF000000L) >> 24) & 0x00000000000000FFL;
		b4 = ((myNum & 0x000000FF00000000L) >> 32) & 0x00000000000000FFL;
		b5 = ((myNum & 0x0000FF0000000000L) >> 40) & 0x00000000000000FFL;
		b6 = ((myNum & 0x00FF000000000000L) >> 48) & 0x00000000000000FFL;
		b7 = ((myNum & 0xFF00000000000000L) >> 56) & 0x00000000000000FFL;
	  
		out[0] = (byte)b0;
		out[1] = (byte)b1;
		out[2] = (byte)b2;
		out[3] = (byte)b3;
		out[4] = (byte)b4;
		out[5] = (byte)b5;
		out[6] = (byte)b6;
		out[7] = (byte)b7;
	  
		return out;		  
	}
  
	/**
	 * Generate a 3-byte array representation of the three LSBs of a 4-byte int.
	 * @param myNum Int to break into bytes. MSB is discarded.
	 * @return Byte array containing individual bytes of the input value. Lowest index corresponds to least significant byte.
	 */
	public static byte[] num24ToByStr(int myNum)
	{
		byte[] out = new byte[3];
	  
		int b0 = 0;
		int b1 = 0;
		int b2 = 0;
		b0 = myNum & 0x000000FF;
		b1 = ((myNum & 0x0000FF00) >> 8) & 0xFF;
		b2 = ((myNum & 0x00FF0000) >> 16) & 0xFF;
	  
		out[0] = (byte)b0;
		out[1] = (byte)b1;
		out[2] = (byte)b2;
	  
		return out;		  
	}
  
  /* ----- CAPACITY MANAGEMENT ----- */
  
	/**
	 * Change the "base capacity" to specified value.
	 * Note that the base capacity is only the capacity of the main byte array.
	 * Overflow is handled by shoving additional material into a linked list of additional byte arrays.
	 * @param newCapacity New base capacity
	 */
	public void changeBaseCapacity(int newCapacity)
	{
		if (newCapacity <= 0) return;
	  
		byte[] nArr = new byte[newCapacity];
	  
		if (this.fSize <= newCapacity)
		{
			for (int i = 0; i < this.fSize; i++){
				byte b = this.getByte(i);
				nArr[i] = b;
			}
			this.capacity = newCapacity;
			this.contents = nArr;
			this.overflow = new LinkedList<byte[]>();
		}
		else
		{
			List<byte[]> nOvr = new LinkedList<byte[]>();
			for (int i = 0; i < this.fSize; i++)
			{
				byte b = this.getByte(i);
				if (i < newCapacity) nArr[i] = b;
				else
				{
					int oai = i / newCapacity;
					int oi = i % newCapacity;
					while (nOvr.size() <= oai)
					{
						byte[] oarr = new byte[newCapacity];
						nOvr.add(oarr);
					}
					nOvr.get(oai)[oi] = b;
				}
			}
			this.capacity = newCapacity;
			this.contents = nArr;
			this.overflow = nOvr;
		}
	}
  
	/**
	 * Change the base capacity to fit the size currently occupied by the file.
	 * Keep in mind that capacity size alteration requires copying of byte arrays.
	 * Use of these functions may sacrifice speed for memory.
	 */
	public void adjustBaseCapacityToSize()
	{
		if (this.capacity == this.fSize) return;
		this.changeBaseCapacity(this.fSize);
	}
 
  /* ----- ADDITION TO FILE ----- */
  
	private void moveAllBytesUp(long position, int amount)
	{
		if (position < 0 || amount <= 0) return;
		if (position >= this.getFileSize()) return;
	  
		long topIndex = this.getFileSize() - 1;
		this.fSize += amount;
	  
		for (long i = topIndex; i >= position; i--)
		{
			this.replaceByte(this.getByte(i), i + amount);
		}
	}
  
  /**
   * Add the byte to the end of the file buffer.
   * @param i8 The byte to insert at the end of the buffer.
   * @throws UnsupportedOperationException If buffer is set to read-only
   * */
	public void addToFile(byte i8)
	{  
		if (this.readOnly()) throw new UnsupportedOperationException();
		if (this.fSize >= this.capacity)
		{
			int oai = (this.fSize / this.capacity) - 1;
			int oi = this.fSize % this.capacity;
			while (this.overflow.size() <= oai)
			{
				byte[] barr = new byte[this.capacity];
				this.overflow.add(barr);
			}
			byte[] oarr = overflow.get(oai);
			oarr[oi] = i8;
			this.fSize++;
		}
		else
		{
			this.contents[this.fSize] = i8;
			this.fSize++;
		}
 	}
  
	/**Add the byte to the file buffer at position, if there is space.
	 * Pushes everything after it down 1 byte.
	   * @param i8 - The byte to insert.
	   * @param position - The position at which to insert.
	   * @throws IndexOutOfBoundsException If position is invalid.
	   * @throws UnsupportedOperationException If buffer is set to read-only
	   * @throws BufferReferenceException If insertion would disrupt a referencing child buffer
	   * */
	public void addToFile(byte i8, int position)
	{
		if (this.readOnly()) throw new UnsupportedOperationException();
		if (position == this.fSize)
	    {
	    	addToFile(i8);
	    	return;
	    }
		if (this.children != null)
		{
			FileBuffer[] ec = this.childrenEncompassing((long)position);
			if (ec != null) throw new BufferReferenceException(ec);
		}
	    
	    if (!this.offsetValid(position)) throw new IndexOutOfBoundsException();
	    
	    this.moveAllBytesUp(position, 1);
	    this.replaceByte(i8, position);	
		this.shiftReferencesAfter((long)position, 1);
	}
  
	/**Add the byte to the file buffer at position, if there is space.
	 * Pushes everything after it down 1 byte.
	   * @param i8 - The byte to insert.
	   * @param position - The position at which to insert.
	   * @throws IndexOutOfBoundsException If position is invalid.
	   * @throws UnsupportedOperationException If buffer is set to read-only
	   * @throws BufferReferenceException If insertion would disrupt a referencing child buffer
	   * */
	public void addToFile(byte i8, long position)
	{
		int pos = (int)position;
		this.addToFile(i8, pos);
	}
  
  /**
   * Add two bytes to the end of the file buffer.
   * Order is determined by FileBuffer's set Endian-ness.
   * @param i16 The value to insert at the end of the buffer.
   * @throws UnsupportedOperationException If buffer is set to read-only
   * */
	public void addToFile(short i16)
	{
		if (this.readOnly()) throw new UnsupportedOperationException();

		byte[] myBytes = numToByStr(i16);
		if (!isFileFormatBE)
		{
			byte temp = myBytes[0];
			myBytes[0] = myBytes[1];
			myBytes[1] = temp;
		}
		for (int i = 1; i >= 0; i--) this.addToFile(myBytes[i]);
	}	
  
  /**
   * Add two bytes to the file buffer at the specified position.
   * Order is determined by FileBuffer's set Endian-ness.
	   * @param i16 The value to insert.
	   * @param position The index to insert it at.
	   * @throws IndexOutOfBoundsException If position is invalid.
	   * @throws UnsupportedOperationException If buffer is set to read-only
	   * @throws BufferReferenceException If insertion would disrupt a referencing child buffer
	   * */
	public void addToFile(short i16, int position)
	{
		if (this.readOnly()) throw new UnsupportedOperationException();
		if (position == this.fSize)
		{
			this.addToFile(i16);
			return;
		}
		if (this.children != null)
		{
			for (int i = 0; i < 2; i++)
			{
				FileBuffer[] ec = this.childrenEncompassing((long)position + i);
				if (ec != null) throw new BufferReferenceException(ec);	
			}
		}
	  
		if (position < 0 || position > this.fSize) throw new IndexOutOfBoundsException();
		this.moveAllBytesUp(position, 2);
	  
		byte[] myBytes = numToByStr(i16);
		if (!isFileFormatBE)
		{
			byte temp = myBytes[0];
			myBytes[0] = myBytes[1];
			myBytes[1] = temp;
		}
		for (int i = 1; i >= 0; i--) this.replaceByte(myBytes[i], position + (1 - i));
		this.shiftReferencesAfter((long)position, 2);
	}
  
	  /**
	   * Add two bytes to the file buffer at the specified position.
	   * Order is determined by FileBuffer's set Endian-ness.
		   * @param i16 The value to insert.
		   * @param position The index to insert it at.
		   * @throws IndexOutOfBoundsException If position is invalid.
		   * @throws UnsupportedOperationException If buffer is set to read-only
		   * @throws BufferReferenceException If insertion would disrupt a referencing child buffer
		   * */
	public void addToFile(short i16, long position)
	{
		int pos = (int)position;
		this.addToFile(i16, pos);
	}
  
	/**
	 * Add four bytes to the end of the file buffer.
	 * Order is determined by FileBuffer's set Endian-ness.
	   * @param i32 The value to insert at the end of the buffer.
	   * @throws UnsupportedOperationException If buffer is set to read-only
	   * */
	public void addToFile(int i32)
	{
		//System.err.println("FileBuffer.addToFile(int) || Called. i32 = " + String.format("%08x", i32));
		if (this.readOnly()) throw new UnsupportedOperationException();
		int myInt = 0;
		if (this.isFileFormatBE) myInt = i32;
		else myInt = switchByO(i32);
		//System.err.println("FileBuffer.addToFile(int) || myInt = " + String.format("%08x", myInt));
		byte[] myBytes = numToByStr(myInt);
		//System.err.println("FileBuffer.addToFile(int) || myBytes = " + String.format("%02x %02x %02x %02x", myBytes[3], myBytes[2], myBytes[1], myBytes[0]));
		for (int i = 3; i >= 0; i--) this.addToFile(myBytes[i]);	  
	}
  
	/**
	Add four bytes to the file buffer at the specified position.
	Order is determined by FileBuffer's set Endian-ness.
	   * @param i32 The value to insert.
	   * @param position The index to insert it at.
	   * @throws IndexOutOfBoundsException If position is invalid.
	   * @throws UnsupportedOperationException If buffer is set to read-only
	   * @throws BufferReferenceException If insertion would disrupt a referencing child buffer
	   * */
	public void addToFile(int i32, int position)
	{
		if (this.readOnly()) throw new UnsupportedOperationException();
		if (position == this.fSize)
		{
			this.addToFile(i32);
			return;
		}
		if (this.children != null)
		{
			this.checkAllBufferReferences((long)position, (long)(position + 4));
		}
	  
		if (position < 0 || position > this.fSize) throw new IndexOutOfBoundsException();
		this.moveAllBytesUp(position, 4);
	  
		int myInt = 0;
		if (this.isFileFormatBE) myInt = i32;
		else myInt = switchByO(i32);
		byte[] myBytes = numToByStr(myInt);
		for (int i = 3; i >= 0; i--) this.replaceByte(myBytes[i], position + (3 - i));		
		this.shiftReferencesAfter((long)position, 4);
	}

	/**
	Add four bytes to the file buffer at the specified position.
	Order is determined by FileBuffer's set Endian-ness.
	   * @param i32 The value to insert.
	   * @param position The index to insert it at.
	   * @throws IndexOutOfBoundsException If position is invalid.
	   * @throws UnsupportedOperationException If buffer is set to read-only
	   * @throws BufferReferenceException If insertion would disrupt a referencing child buffer
	   * */
	public void addToFile(int i32, long position)
	{
		int pos = (int)position;
		this.addToFile(i32, pos);
	}
  
  /**
   * Add eight bytes to the end of the file buffer.
   * Order is determined by FileBuffer's set Endian-ness.
	   * @param i64 The value to insert at the end of the buffer.
	   * @throws UnsupportedOperationException If buffer is set to read-only
	   * */
	public void addToFile(long i64)
	{
		if (this.readOnly()) throw new UnsupportedOperationException();
		long myLong = 0;
		if (this.isFileFormatBE) myLong = i64;
		else myLong = switchByO(i64);
		byte[] myBytes = numToByStr(myLong);
		for (int i = 7; i >= 0; i--) this.addToFile(myBytes[i]);  
	}
  
  /**
   * Add eight bytes to the file buffer at the specified position.
   * Order is determined by FileBuffer's set Endian-ness.
   * @param i64 The value to insert.
   * @param position The index to insert it at.
   * @throws IndexOutOfBoundsException If position is invalid.
   * @throws UnsupportedOperationException If buffer is set to read-only
   * @throws BufferReferenceException If insertion would disrupt a referencing child buffer
   * */
	public void addToFile(long i64, int position)
	{
		if (this.readOnly()) throw new UnsupportedOperationException();
		if (position == this.fSize)
		{
			this.addToFile(i64);
			return;
		}
		if (this.children != null)
		{
			this.checkAllBufferReferences((long)position, (long)(position + 8));
		}
	  
		if (position < 0 || position > this.fSize) throw new IndexOutOfBoundsException();
		this.moveAllBytesUp(position, 8);
	  
		long myLong = 0;
		if (this.isFileFormatBE) myLong = i64;
		else myLong = switchByO(i64);
		byte[] myBytes = numToByStr(myLong);
		for (int i = 7; i >= 0; i--) this.replaceByte(myBytes[i], position + (7 - i));
		this.shiftReferencesAfter((long)position, 8);
	}
  
	  /**
	   * Add eight bytes to the file buffer at the specified position.
	   * Order is determined by FileBuffer's set Endian-ness.
	   * @param i64 The value to insert.
	   * @param position The index to insert it at.
	   * @throws IndexOutOfBoundsException If position is invalid.
	   * @throws UnsupportedOperationException If buffer is set to read-only
	   * @throws BufferReferenceException If insertion would disrupt a referencing child buffer
	   * */
	public void addToFile(long i64, long position)
	{
		int pos = (int)position;
		this.addToFile(i64, pos);
	}
  
	/**
   	* Add three bytes to the end of the file buffer as designated by 
   	* three LSBs in provided int.
   	* Order is determined by FileBuffer's set Endian-ness.
	* @param i24 The value to insert at the end of the buffer. MSB is discarded.
	* @throws UnsupportedOperationException If buffer is set to read-only
	* */
	public void add24ToFile(int i24)
	{	 
		if (this.readOnly()) throw new UnsupportedOperationException();
		byte[] myBytes = numToByStr(i24);
		if (this.isFileFormatBE) for (int i = 2; i >= 0; i--) this.addToFile(myBytes[i]);	
		else for (int i = 0; i < 3; i++) this.addToFile(myBytes[i]);	
	}
  
	/**
   	* Add three bytes to the end of the file buffer as designated by 
   	* three LSBs in provided int.
   	* Order is determined by FileBuffer's set Endian-ness.
	* @param i24 The value to insert at the end of the buffer. MSB is discarded.
	* @param position The index to insert it at.
	* @throws IndexOutOfBoundsException If position is invalid.
	* @throws UnsupportedOperationException If buffer is set to read-only
	* @throws BufferReferenceException If insertion would disrupt a referencing child buffer
	* */
	public void add24ToFile(int i24, int position)
	{
		if (this.readOnly()) throw new UnsupportedOperationException();
		if (position == this.fSize)
		{
			this.add24ToFile(i24);
			return;
		}
		if (this.children != null)
		{
			this.checkAllBufferReferences((long)position, (long)(position + 3));
		}
	  
		if (position < 0 || position > this.fSize) throw new IndexOutOfBoundsException();
		this.moveAllBytesUp(position, 3);
	  
		byte[] myBytes = numToByStr(i24);
		if (this.isFileFormatBE) for (int i = 2; i >= 0; i--) this.replaceByte(myBytes[i], position + (2 - i));	
		else for (int i = 0; i < 3; i++) this.replaceByte(myBytes[i], position + i);
		this.shiftReferencesAfter((long)position, 3);
	}
  
	/**
   	* Add three bytes to the end of the file buffer as designated by 
   	* three LSBs in provided int.
   	* Order is determined by FileBuffer's set Endian-ness.
	* @param i24 The value to insert at the end of the buffer. MSB is discarded.
	* @param position The index to insert it at.
	* @throws IndexOutOfBoundsException If position is invalid.
	* @throws UnsupportedOperationException If buffer is set to read-only
	* @throws BufferReferenceException If insertion would disrupt a referencing child buffer
	* */
	public void add24ToFile(int i24, long position)
	{
		int pos = (int)position;
		this.add24ToFile(i24, pos);
	}
  
  /**
   * Add a piece of another FileBuffer to this FileBuffer
   * @param addition The FileBuffer to add
   * @param insertPos The index of this at which addition will be inserted.
   * @param stPos The first byte of addition to add.
   * @param edPos The final byte (exclusive) of addition to add.
   * @throws IndexOutOfBoundsException If any positions given are invalid.
   * @throws NullPointerException If addition parameter is a null reference
   * @throws UnsupportedOperationException If buffer is set to read-only
   * @throws BufferReferenceException If insertion would disrupt a referencing child buffer
   */
	public void addToFile(FileBuffer addition, int insertPos, int stPos, int edPos)
	{
		if (this.readOnly()) throw new UnsupportedOperationException();
		if (addition == null || addition.getFileSize() <= 0) throw new NullPointerException();
		if (insertPos < 0 || insertPos > (int)this.getFileSize()) throw new IndexOutOfBoundsException();
		if (stPos < 0 || stPos >= addition.getFileSize()) throw new IndexOutOfBoundsException();
		if (edPos < 0 || edPos > addition.getFileSize() || edPos <= stPos) throw new IndexOutOfBoundsException();
		if (insertPos == (int)this.getFileSize())
		{
			this.addToFile(addition, stPos, edPos);
			return;
		}
		if (this.children != null)
		{
			this.checkAllBufferReferences((long)insertPos, (long)(insertPos + (edPos - stPos)));
		}
	  
		int isize = edPos - stPos;
		this.moveAllBytesUp(insertPos, isize);
	  
		for (int i = 0; i < isize; i++)
		{
			this.replaceByte(addition.getByte(stPos + i), insertPos + i);
		}
		this.shiftReferencesAfter((long)insertPos, (edPos - stPos));
	}
  
	  /**
	   * Add a piece of another FileBuffer to this FileBuffer
	   * @param addition The FileBuffer to add
	   * @param insertPos The index of this at which addition will be inserted.
	   * @param stPos The first byte of addition to add.
	   * @param edPos The final byte (exclusive) of addition to add.
	   * @throws IndexOutOfBoundsException If any positions given are invalid.
	   * @throws NullPointerException If addition parameter is a null reference
	   * @throws UnsupportedOperationException If buffer is set to read-only
	   * @throws BufferReferenceException If insertion would disrupt a referencing child buffer
	   */
	public void addToFile(FileBuffer addition, long insertPos, long stPos, long edPos)
	{
		if (this.readOnly()) throw new UnsupportedOperationException();
		if (addition == null || addition.isEmpty()) throw new NullPointerException();
		if (!this.offsetValid(insertPos)) throw new IndexOutOfBoundsException();
		addition.checkOffsetPair(stPos, edPos);
		if (insertPos == (int)this.getFileSize())
		{
			this.addToFile(addition, stPos, edPos);
			return;
		}
		if (this.children != null)
		{
			this.checkAllBufferReferences(insertPos, insertPos + (edPos - stPos));
		}
		
		int iPos = (int)insertPos;
		long lisize = edPos - stPos;
		int isize = 0;
		if (lisize > 0x7FFFFFFFL) isize = 0x7FFFFFFF;
		else isize = (int)lisize;
		this.moveAllBytesUp(iPos, isize);
	  
		for (int i = 0; i < isize; i++)
		{
			this.replaceByte(addition.getByte(stPos + i), iPos + i);
		}
		this.shiftReferencesAfter(insertPos, (int)(edPos - stPos));
  }
 
	  /**
	   * Add a piece of another FileBuffer to this FileBuffer
	   * @param addition The FileBuffer to add
	   * @throws NullPointerException If addition parameter is a null reference
	   * @throws UnsupportedOperationException If buffer is set to read-only
	   */
	public void addToFile(FileBuffer addition)
	{
		if (this.readOnly()) throw new UnsupportedOperationException();
		if (addition == null || addition.getFileSize() <= 0) throw new NullPointerException();
		this.addToFile(addition, 0, (int)addition.getFileSize());
	}
  
	  /**
	   * Add a piece of another FileBuffer to this FileBuffer
	   * @param addition The FileBuffer to add
	   * @param stPos The first byte of addition to add.
	   * @param edPos The final byte (exclusive) of addition to add.
	   * @throws IndexOutOfBoundsException If any positions given are invalid.
	   * @throws NullPointerException If addition parameter is a null reference
	   * @throws UnsupportedOperationException If buffer is set to read-only
	   */
	public void addToFile(FileBuffer addition, int stPos, int edPos)
	{
		if (this.readOnly()) throw new UnsupportedOperationException();
		if (addition == null || addition.getFileSize() <= 0) throw new NullPointerException();
		if (addition.isEmpty()) throw new IndexOutOfBoundsException();
		if (this.capacity <= 0) throw new IndexOutOfBoundsException();
	  
		if (stPos < 0) stPos = 0;
		if (edPos > (int)addition.getFileSize()) edPos = (int)addition.getFileSize();
		if (edPos < stPos)throw new IndexOutOfBoundsException();
	  
		for (int i = stPos; i < edPos; i++)
		{
			this.addToFile(addition.getByte(i));
		}
	}
  
	  /**
	   * Add a piece of another FileBuffer to this FileBuffer
	   * @param addition The FileBuffer to add
	   * @param stPos The first byte of addition to add.
	   * @param edPos The final byte (exclusive) of addition to add.
	   * @throws IndexOutOfBoundsException If any positions given are invalid.
	   * @throws NullPointerException If addition parameter is a null reference
	   * @throws UnsupportedOperationException If buffer is set to read-only
	   */
	public void addToFile(FileBuffer addition, long stPos, long edPos)
	{
		if (this.readOnly()) throw new UnsupportedOperationException();
		if (addition == null || addition.getFileSize() <= 0) throw new NullPointerException();
		if (addition.isEmpty()) throw new IndexOutOfBoundsException();
		if (this.capacity <= 0) throw new IndexOutOfBoundsException();
	  
		if (stPos < 0) stPos = 0;
		if (edPos > addition.getFileSize()) edPos = addition.getFileSize();
		if (edPos < stPos) throw new IndexOutOfBoundsException();
		if (edPos - stPos > 0x7FFFFFFFL) edPos = stPos + 0x7FFFFFFFL;
	  
		for (long i = stPos; i < edPos; i++)
		{
			this.addToFile(addition.getByte(i));
		}
	}

  /* ----- UNSIGNED HANDLING ----- */
  
	/**
	 * To compensate for lack of a Byte.toUnsignedShort(byte) function.
	 * Converts a byte as if it was unsigned to a short (does not sign extend)
	 * @param in - Byte to convert
	 * @return short representation of byte in reading in as an unsigned value.
	 */
	public static short UByteToUShort(byte in)
	{
		short out = 0;
		int temp = 0;
		temp = (int) in;
		temp = temp & 0xFF;
		out = (short)temp;
		return out;
	}
  
  /* ----- CONTENT RETRIEVAL ----- */
  
	/**
	 * Return the next two bytes from the position (inclusive)
	 * as a short integer.
	 * If the file is Little-Endian encoded, the bytes are reversed.
	 * @param position The position in the file.
	 * @return The next two bytes from position as a BE short.
	 * @throws IndexOutOfBoundsException If position is invalid.*/
	public short shortFromFile(int position)
	{
		return this.shortFromFile((long)position);
	}
  
	/**
	 * Return the next two bytes from the position (inclusive)
	 * as a short integer.
	 * If the file is Little-Endian encoded, the bytes are reversed.
	 * @param position The position in the file.
	 * @return The next two bytes from position as a BE short.
	 * @throws IndexOutOfBoundsException If position is invalid.
	 * */
	public short shortFromFile(long position)
	{
		//return this.shortFromFile((int)position);
		short myShort = 0;
		int b0 = 0;
		int b1 = 0;
	
		if (position < 0) throw new IndexOutOfBoundsException();
	
		Byte.toUnsignedInt(this.getByte(position));
		b1 = Byte.toUnsignedInt(this.getByte(position));
	
		if (position >= this.getFileSize() - 1) b0 = 0;
		else b0 = Byte.toUnsignedInt(this.getByte(position + 1));
	
		if (this.isFileFormatBE)
		{
			b1 = (b1 << 8);
			myShort = (short)(b1 | b0);
		}
		else
		{
			b0 = (b0 << 8);
			myShort = (short) (b0 | b1);
		} 
		return myShort;
	}
  
	/**
	 * Return the next two bytes from the position (inclusive)
	 * as an integer with 0 as its MSB.
	 * If the file is Little-Endian encoded, the bytes are reversed.
	 * @param position The position in the file.
	 * @return The next two bytes from position as a BE short.
	 * @throws IndexOutOfBoundsException If position is invalid.*/
	public int shortishFromFile(int position)
	{
		return this.shortishFromFile((long)position);
	}
  
	/**
	 * Return the next two bytes from the position (inclusive)
	 * as an integer with 0 as its MSB.
	 * If the file is Little-Endian encoded, the bytes are reversed.
	 * @param position The position in the file.
	 * @return The next two bytes from position as a BE short.
	 * @throws IndexOutOfBoundsException If position is invalid.*/
	public int shortishFromFile(long position)
	{
		//return this.shortishFromFile((int)position);
		int my24 = 0;
		int b0 = 0;
		int b1 = 0;
		int b2 = 0;
	
		if (position < 0 || position >= this.getFileSize() - 2) throw new IndexOutOfBoundsException();
	
		b2 = Byte.toUnsignedInt(this.getByte(position));
		b1 = Byte.toUnsignedInt(this.getByte(position + 1));
		b0 = Byte.toUnsignedInt(this.getByte(position + 2));
	
		if (this.isFileFormatBE)
		{
			b2 = b2 << 16;
			b1 = b1 << 8;
			my24 = (b2 | b1 | b0) & (0x00FFFFFF);	
		}
		else
		{
			b0 = b0 << 16;
			b1 = b1 << 8;
			my24 = (b0 | b1 | b2) & (0x00FFFFFF);
		}
	
		return my24;
	}
  
	/**
	 * Return the next four bytes from the position (inclusive)
	 * as a standard integer.
	 * If the file is Little-Endian encoded, the bytes are reversed.
	 * @param position The position in the file.
	 * @throws IndexOutOfBoundsException If position is invalid
	 * @return The next four bytes from position as a BE short.
	 * */
	public int intFromFile(int position)
	{
		return this.intFromFile((long)position);
	}
  
	/**
	 * Return the next four bytes from the position (inclusive)
	 * as a standard integer.
	 * If the file is Little-Endian encoded, the bytes are reversed.
	 * @param position The position in the file.
	 * @throws 7IndexOutOfBoundsException If position is invalid
	 * @return The next four bytes from position as a BE short.
	 * */
	public int intFromFile(long position)
	{
		//return this.intFromFile((int)position);
		int myInt = 0;
		int b0 = 0;
		int b1 = 0;
		int b2 = 0;
		int b3 = 0;
		
		if (position < 0) throw new IndexOutOfBoundsException();
		
		Byte.toUnsignedInt(this.getByte(position));
		b3 = Byte.toUnsignedInt(this.getByte(position));
		if (!(position + 1 >= this.getFileSize())) b2 = Byte.toUnsignedInt(this.getByte(position + 1));
		if (!(position + 2 >= this.getFileSize())) b1 = Byte.toUnsignedInt(this.getByte(position + 2));
		if (!(position + 3 >= this.getFileSize())) b0 = Byte.toUnsignedInt(this.getByte(position + 3));
		
		if (this.isBigEndian())
		{
			b3 = b3 << 24; //System.out.println("b3 << 24 : " + Integer.toHexString(b3));
			b2 = b2 << 16; //System.out.println("b2 << 16 : " + Integer.toHexString(b2));
			b1 = b1 << 8; //System.out.println("b1 << 8 : " + Integer.toHexString(b1));
			myInt = b3 | b2 | b1 | b0;	
		}
		else
		{
			b0 = b0 << 24; //System.out.println("b3 << 24 : " + Integer.toHexString(b3));
			b1 = b1 << 16; //System.out.println("b2 << 16 : " + Integer.toHexString(b2));
			b2 = b2 << 8; //System.out.println("b1 << 8 : " + Integer.toHexString(b1));
			myInt = b0 | b1 | b2 | b3;	
		}	
	    return myInt;
	}
  
	/**
	 * Return the next eight bytes from the position (inclusive)
	 * as a long integer.
	 * If the file is Little-Endian encoded, the bytes are reversed.
	 * @param position The position in the file.
	 * @return The next eight bytes from position as a BE short.
	 * @throws IndexOutOfBoundsException If position is invalid
	 * */
	public long longFromFile(int position)
	{
		return longFromFile((long)position);
	}
  
	/**
	 * Return the next eight bytes from the position (inclusive)
	 * as a long integer.
	 * If the file is Little-Endian encoded, the bytes are reversed.
	 * @param position The position in the file.
	 * @return The next eight bytes from position as a BE short.
	 * @throws IndexOutOfBoundsException If position is invalid
	 * */
	public long longFromFile(long position)
	{
		long myLong = 0;
		long b0 = 0;
		long b1 = 0;
		long b2 = 0;
		long b3 = 0;
		long b4 = 0;
		long b5 = 0;
		long b6 = 0;
		long b7 = 0;
		
		if (position < 0) throw new IndexOutOfBoundsException();
		
		b7 = Byte.toUnsignedLong(this.getByte(position));
		if (position + 1 < this.getFileSize()) b6 = Byte.toUnsignedLong(this.getByte(position + 1));
		if (position + 2 < this.getFileSize()) b5 = Byte.toUnsignedLong(this.getByte(position + 2));
		if (position + 3 < this.getFileSize()) b4 = Byte.toUnsignedLong(this.getByte(position + 3));
		if (position + 4 < this.getFileSize()) b3 = Byte.toUnsignedLong(this.getByte(position + 4));
		if (position + 5 < this.getFileSize()) b2 = Byte.toUnsignedLong(this.getByte(position + 5));
		if (position + 6 < this.getFileSize()) b1 = Byte.toUnsignedLong(this.getByte(position + 6));
		if (position + 7 < this.getFileSize()) b0 = Byte.toUnsignedLong(this.getByte(position + 7));
		
		if (this.isFileFormatBE)
		{
			b7 = b7 << 56;
			b6 = b6 << 48;
			b5 = b5 << 40;
			b4 = b4 << 32;
			b3 = b3 << 24;
			b2 = b2 << 16;
			b1 = b1 << 8;
			myLong = b7 | b6 | b5 | b4 | b3 | b2 | b1 | b0;
		}
		else
		{
			b0 = b0 << 56;
			b1 = b1 << 48;
			b2 = b2 << 40;
			b3 = b3 << 32;
			b4 = b4 << 24;
			b5 = b5 << 16;
			b6 = b6 << 8;
			myLong = b0 | b1 | b2 | b3 | b4 | b5 | b6 | b7;	
		}
		
	    return myLong;
	}
  
	/**
	 * Get between 1 and 8 bits from the FileBuffer at the given bit and byte positions
	 * and return them as a byte, with all unused upper bits set to 0.
	 * @param numBits Number of bits to retrieve
	 * @param bytePos Byte offset to start bit retrieval from
	 * @param bitPos Bit offset of byte to start bit retrieval from ||
	 * [MSB] 7 6 5 4 3 2 1 0 [LSB]
	 * @return A byte containing the desired bits. If less than 8 bits are called for, upper bits of
	 * return byte will be set to 0.
	 * @throws ArrayIndexOutOfBoundsException If any parameter is invalid, including bit position or number of bits.
	 */
	public byte getBits8(int numBits, int bytePos, int bitPos)
	{
		return this.getBits8(numBits, (long)bytePos, bitPos);
	}
  
	/**
	 * Get between 1 and 8 bits from the FileBuffer at the given bit and byte positions
	 * and return them as a byte, with all unused upper bits set to 0.
	 * @param numBits Number of bits to retrieve
	 * @param bytePos Byte offset to start bit retrieval from
	 * @param bitPos Bit offset of byte to start bit retrieval from ||
	 * [MSB] 7 6 5 4 3 2 1 0 [LSB]
	 * @return A byte containing the desired bits. If less than 8 bits are called for, upper bits of
	 * return byte will be set to 0.
	 * @throws IndexOutOfBoundsException If any parameter is invalid, including bit position or number of bits.
	 */
	public byte getBits8(int numBits, long bytePos, int bitPos)
	{
		if (numBits < 1 || numBits > 8) throw new IndexOutOfBoundsException();
		if (bytePos >= this.getFileSize() || bytePos < 0) throw new IndexOutOfBoundsException();
		if (bitPos < 0 || bitPos > 7) throw new IndexOutOfBoundsException();
	  
		byte myByte = 0;
		byte sourceByte = this.getByte(bytePos);
		byte nextByte = 0;
		int readyByte = 0;
		int nByte = 0;
	  
		if (numBits == 8 && bitPos == 7) return sourceByte;
	  
		readyByte = Byte.toUnsignedInt(sourceByte);
		readyByte = (readyByte << (7 - bitPos));
	  
		readyByte = (readyByte >> (8 - numBits)) & 0xFF;
	  
		if ((bitPos + 1) < numBits && bytePos + 1 < this.getFileSize())
		{
			/*We will have to retrieve the next byte*/
			nextByte = this.getByte(bytePos + 1);
			nByte = Byte.toUnsignedInt(nextByte);
			nByte = nByte >> (8 - (numBits - (bitPos + 1)));
			nByte = nByte & 0xFF;
			readyByte = readyByte | nByte;
		}
	  
		myByte = (byte)readyByte;
		return myByte;
	}
  
	/**
	 * Get between 1 and 16 bits from the FileBuffer at the given bit and byte positions
	 * and return them as a short, with all unused upper bits set to 0.
	 * @param numBits - Number of bits to retrieve
	 * @param bytePos - Byte offset to start bit retrieval from
	 * @param bitPos - Bit offset of byte to start bit retrieval from || <br>
	 * [MSB] 15 14 13 12 11 10 9 8 | 7 6 5 4 3 2 1 0 [LSB]
	 * @return A short containing the desired bits. If less than 16 bits are called for, upper bits of
	 * return value will be set to 0.
	 * @throws IndexOutOfBoundsException If any parameter is invalid, including bit position or number of bits.
	 */
	public short getBits16(int numBits, int bytePos, int bitPos)
	{
		return getBits16(numBits, (long)bytePos, bitPos);
	}
 
	/**
	 * Get between 1 and 16 bits from the FileBuffer at the given bit and byte positions
	 * and return them as a short, with all unused upper bits set to 0.
	 * @param numBits - Number of bits to retrieve
	 * @param bytePos - Byte offset to start bit retrieval from
	 * @param bitPos - Bit offset of byte to start bit retrieval from || <br>
	 * [MSB] 15 14 13 12 11 10 9 8 | 7 6 5 4 3 2 1 0 [LSB]
	 * @return A short containing the desired bits. If less than 16 bits are called for, upper bits of
	 * return value will be set to 0.
	 * @throws IndexOutOfBoundsException If any parameter is invalid, including bit position or number of bits.
	 */
	public short getBits16(int numBits, long bytePos, int bitPos)
	{
		if (numBits < 1 || numBits > 16) throw new IndexOutOfBoundsException();
		if (bytePos >= this.getFileSize() || bytePos < 0) throw new IndexOutOfBoundsException();
		if (bitPos < 0 || bitPos > 15) throw new IndexOutOfBoundsException();
	  
		short myVal = 0;
	  	short sourceVal = 0;
	  	short nextVal = 0;
	  	int readyVal = 0;
	  	int nextReady = 0;
	  	boolean wasBE = this.isBigEndian();
	  
	  	if (!wasBE) this.isFileFormatBE = true;
	  
	  	if (bytePos < this.getFileSize()) sourceVal = this.shortFromFile(bytePos);
	  	else
	  	{
	  		if (!wasBE) this.isFileFormatBE = false;
	  		throw new IndexOutOfBoundsException();
	  	}
	  
	  	if (numBits == 16 && bitPos == 15)
	  	{
	  		if (!wasBE) this.isFileFormatBE = false;
	  		return sourceVal;
	  	}

	  	readyVal = Short.toUnsignedInt(sourceVal);
	  	readyVal = (readyVal << (15 - bitPos));
	  
	  	readyVal = (readyVal >> (16 - numBits)) & 0xFFFF;
	  
	  	if ((bitPos + 1) < numBits && bytePos + 2 < this.getFileSize())
	  	{
	  		/*We will have to retrieve the next byte*/
	  		nextVal = this.shortFromFile(bytePos + 2);
	  		nextReady = Short.toUnsignedInt(nextVal);
	  		nextReady = nextReady >> (16 - (numBits - (bitPos + 1)));
	  		nextReady = nextReady & 0xFFFF;
	  		readyVal = readyVal | nextReady;
	  	}
	  
	  	myVal = (short)readyVal;
	  	if (!wasBE) this.isFileFormatBE = false;
	  	return myVal;
	}
  
	/**
	 * Get between 1 and 32 bits from the FileBuffer at the given bit and byte positions
	 * and return them as an int, with all unused upper bits set to 0.
	 * @param numBits - Number of bits to retrieve
	 * @param bytePos - Byte offset to start bit retrieval from
	 * @param bitPos - Bit offset of byte to start bit retrieval from || <br>
	 * [MSB] 31 30 29 28 27 26 25 24 |
	 * <br> 23 22 21 20 19 18 17 16 |
	 * <br> 15 14 13 12 11 10 9 |
	 * <br> 8 7 6 5 4 3 2 1 0 [LSB]
	 * @return An int containing the desired bits. If less than 32 bits are called for, upper bits of
	 * return value will be set to 0.
	 * @throws IndexOutOfBoundsException If any parameter is invalid, including bit position or number of bits.
	 */
	public int getBits32(int numBits, int bytePos, int bitPos)
	{
		return this.getBits32(numBits, (long)bytePos, bitPos); 
	}
  
	/**
	 * Get between 1 and 32 bits from the FileBuffer at the given bit and byte positions
	 * and return them as an int, with all unused upper bits set to 0.
	 * @param numBits - Number of bits to retrieve
	 * @param bytePos - Byte offset to start bit retrieval from
	 * @param bitPos - Bit offset of byte to start bit retrieval from ||
	 * [MSB] 31 30 29 28 27 26 25 24 | 
	 * <br>23 22 21 20 19 18 17 16 | 
	 * <br>15 14 13 12 11 10 9 | 
	 * <br>8 7 6 5 4 3 2 1 0 [LSB]
	 * @return An int containing the desired bits. If less than 32 bits are called for, upper bits of
	 * return value will be set to 0.
	 * @throws IndexOutOfBoundsException If any parameter is invalid, including bit position or number of bits.
	 */
	public int getBits32(int numBits, long bytePos, int bitPos)
	{
		if (numBits < 1 || numBits > 32) throw new IndexOutOfBoundsException();
		if (bytePos >= this.getFileSize() || bytePos < 0) throw new IndexOutOfBoundsException();
		if (bitPos < 0 || bitPos > 31) throw new IndexOutOfBoundsException();
	  
		int sourceVal = 0;
		int nextVal = 0;
		int readyVal = 0;
		int nextReady = 0;
		boolean wasBE = this.isFileFormatBE;
	  
		if (!wasBE) this.isFileFormatBE = true;
	  
		if (bytePos < this.getFileSize()) sourceVal = this.intFromFile(bytePos);
		else
		{
			if (!wasBE) this.isFileFormatBE = false;
			throw new IndexOutOfBoundsException();
		}
	  
		if (numBits == 32 && bitPos == 31)
		{
			if (!wasBE) this.isFileFormatBE = false;
			return sourceVal;
		}
	  
		readyVal = sourceVal;
		readyVal = (readyVal << (31 - bitPos));
	  
		readyVal = (readyVal >> (32 - numBits));
	  
		if ((bitPos + 1) < numBits && bytePos + 4 < this.getFileSize())
		{
			/*We will have to retrieve the next byte*/
			nextVal = this.intFromFile(bytePos + 4);
			nextReady = nextVal;
			nextReady = nextReady >> (32 - (numBits - (bitPos + 1)));
			readyVal = readyVal | nextReady;
		}
	  
		if (!wasBE) this.isFileFormatBE = false;
		return readyVal; 
	}
  
	/**
	 * Get between 1 and 64 bits from the FileBuffer at the given bit and byte positions
	 * and return them as a long, with all unused upper bits set to 0.
	 * @param numBits - Number of bits to retrieve
	 * @param bytePos - Byte offset to start bit retrieval from
	 * @param bitPos - Bit offset of byte to start bit retrieval from ||
	 * <br>[MSB] 63 62 61 60 59 58 57 56 | 
	 * <br>55 54 53 52 51 50 49 48 | 
	 * <br>47 46 45 44 43 42 41 40 | 
	 * <br>39 38 37 36 35 34 33 32 | 
	 * <br>31 30 29 28 27 26 25 24 | 
	 * <br>23 22 21 20 19 18 17 16 | 
	 * <br>15 14 13 12 11 10 9 | 
	 * <br>8 7 6 5 4 3 2 1 0 [LSB]
	 * @return A long containing the desired bits. If less than 64 bits are called for, upper bits of
	 * return value will be set to 0.
	 * @throws IndexOutOfBoundsException If any parameter is invalid, including bit position or number of bits.
	 */
	public long getBits64(int numBits, int bytePos, int bitPos)
	{
		return this.getBits64(numBits, (long)bytePos, bitPos);	  
	}
  
	/**
	 * Get between 1 and 64 bits from the FileBuffer at the given bit and byte positions
	 * and return them as a long, with all unused upper bits set to 0.
	 * @param numBits - Number of bits to retrieve
	 * @param bytePos - Byte offset to start bit retrieval from
	 * @param bitPos - Bit offset of byte to start bit retrieval from ||
	 * <br>[MSB] 63 62 61 60 59 58 57 56 | 
	 * <br>55 54 53 52 51 50 49 48 | 
	 * <br>47 46 45 44 43 42 41 40 | 
	 * <br>39 38 37 36 35 34 33 32 | 
	 * <br>31 30 29 28 27 26 25 24 | 
	 * <br>23 22 21 20 19 18 17 16 | 
	 * <br>15 14 13 12 11 10 9 | 
	 * <br>8 7 6 5 4 3 2 1 0 [LSB]
	 * @return A long containing the desired bits. If less than 64 bits are called for, upper bits of
	 * return value will be set to 0.
	 * @throws IndexOutOfBoundsException If any parameter is invalid, including bit position or number of bits.
	 */
	public long getBits64(int numBits, long bytePos, int bitPos)
	{
		if (numBits < 1 || numBits > 64) throw new IndexOutOfBoundsException();
		if (bytePos >= this.getFileSize() || bytePos < 0) throw new IndexOutOfBoundsException();
		if (bitPos < 0 || bitPos > 63) throw new IndexOutOfBoundsException();
	  
		long sourceVal = 0;
		long nextVal = 0;
		long readyVal = 0;
		long nextReady = 0;
		boolean wasBE = this.isFileFormatBE;
	  
		if (!wasBE) this.isFileFormatBE = true;
	  
		if (bytePos < this.getFileSize()) sourceVal = this.intFromFile(bytePos);
		else
		{
			if (!wasBE) this.isFileFormatBE = false;
			throw new IndexOutOfBoundsException();
		}
	  
		if (numBits == 64 && bitPos == 63)
		{
			if (!wasBE) this.isFileFormatBE = false;
			return sourceVal;
		}
	  
		readyVal = sourceVal;
		readyVal = (readyVal << (63 - bitPos));
	  
		readyVal = (readyVal >> (63 - numBits));
	  
		if ((bitPos + 1) < numBits && bytePos + 8 < this.getFileSize())
		{
			/*We will have to retrieve the next byte*/
			nextVal = this.longFromFile(bytePos + 8);
			nextReady = nextVal;
			nextReady = nextReady >> (64 - (numBits - (bitPos + 1)));
			readyVal = readyVal | nextReady;
		}
	  
		if (!wasBE) this.isFileFormatBE = false;
		return readyVal; 
	}
  
	/**
	 * Get bits from the FileBuffer at the given bit and byte positions
	 * and return them as a string.
	 * @param numBits - Number of bits to retrieve
	 * @param bytePos - Byte offset to start bit retrieval from
	 * @param bitPos - Bit offset of byte to start bit retrieval from ||
	 * <br>[MSB] 63 62 61 60 59 58 57 56 | 
	 * <br>55 54 53 52 51 50 49 48 | 
	 * <br>47 46 45 44 43 42 41 40 | 
	 * <br>39 38 37 36 35 34 33 32 | 
	 * <br>31 30 29 28 27 26 25 24 | 
	 * <br>23 22 21 20 19 18 17 16 | 
	 * <br>15 14 13 12 11 10 9 | 
	 * <br>8 7 6 5 4 3 2 1 0 [LSB]
	 * @return A String length numBits representing the bits retrieved. Bits are represented by ASCII characters '0' and '1'.
	 */
	public String getBits(int numBits, int bytePos, int bitPos)
	{
		return this.getBits(numBits, (long)bytePos, bitPos);
	}
  
	/**
	 * Get bits from the FileBuffer at the given bit and byte positions
	 * and return them as a string.
	 * @param numBits - Number of bits to retrieve
	 * @param bytePos - Byte offset to start bit retrieval from
	 * @param bitPos - Bit offset of byte to start bit retrieval from ||
	 * <br>[MSB] 63 62 61 60 59 58 57 56 | 
	 * <br>55 54 53 52 51 50 49 48 | 
	 * <br>47 46 45 44 43 42 41 40 | 
	 * <br>39 38 37 36 35 34 33 32 | 
	 * <br>31 30 29 28 27 26 25 24 | 
	 * <br>23 22 21 20 19 18 17 16 | 
	 * <br>15 14 13 12 11 10 9 | 
	 * <br>8 7 6 5 4 3 2 1 0 [LSB]
	 * @return A String length numBits representing the bits retrieved. Bits are represented by ASCII characters '0' and '1'.
	 */
	public String getBits(int numBits, long bytePos, int bitPos)
	{
		String nums = "";
		if (numBits < 0 || bitPos < 0 || bytePos < 0 || bytePos >= this.getFileSize()) return nums;
	  
		if (numBits <= 8 && bitPos <= 7)
		{
			int myBits = Byte.toUnsignedInt(getBits8(numBits, bytePos, bitPos));
			for (int i = numBits - 1; i >= 0; i--)
			{
				int mask = 1 << i;
				if ((mask & myBits) != 0) nums += "1";
				else nums += "0";
			}
		}
	  
		return nums;
	}
  
  /* ----- DELETION ----- */
  
	private void moveAllBytesDown(int position, int amount)
	{
		if (position < 0 || amount <= 0) return;
		if (position >= this.fSize) return;
	  
		for (int i = position; i < this.fSize; i++)
		{
			this.replaceByte(this.getByte(i + amount), i);
		}
	  
		this.fSize -= amount;
	}
  
	/**
	 * Delete everything from the file starting (inclusive) at stOff
	 * @param stOff Position to start deletion at.
	 * @throws IndexOutOfBoundsException If offset is invalid
	 * @throws UnsupportedOperationException If buffer is read-only
	 * @throws BufferReferenceException If deletion would affect any referencing child buffers.
	 */
	public void deleteFromFile(int stOff)
	{
		if (this.readOnly()) throw new UnsupportedOperationException();
		if (stOff < 0 || stOff >= this.fSize) throw new IndexOutOfBoundsException();
		if (!this.checkAllBufferReferences((long)stOff)) this.fSize = stOff;
	}

	/**
	 * Delete everything from the file starting (inclusive) at stOff
	 * @param stOff Position to start deletion at.
	 * @throws IndexOutOfBoundsException If offset is invalid
	 * @throws UnsupportedOperationException If buffer is read-only
	 * @throws BufferReferenceException If deletion would affect any referencing child buffers.
	 */
	public void deleteFromFile(long stOff)
	{
		this.deleteFromFile((int)stOff);
	}
  
	/**
	 * Delete everything from file between two positions, and move up everything at and after edOff
	 * to fill the gap.
	 * @param stOff Starting position (inclusive) of deletion
	 * @param edOff Ending position (exclusive) of deletion
	 * @throws IndexOutOfBoundsException If offset is invalid
	 * @throws UnsupportedOperationException If buffer is read-only
	 * @throws BufferReferenceException If deletion would affect any referencing child buffers.
	 */
	public void deleteFromFile(int stOff, int edOff)
	{
		if (this.readOnly()) throw new UnsupportedOperationException();
		if (stOff < 0 || stOff >= edOff) throw new IndexOutOfBoundsException();
		if (edOff >= this.fSize) throw new IndexOutOfBoundsException();
		if (!this.checkAllBufferReferences(stOff, edOff))
		{
			int delSz = edOff - stOff;
			this.moveAllBytesDown(stOff, delSz);	
		}
	}

	/**
	 * Delete everything from file between two positions, and move up everything at and after edOff
	 * to fill the gap.
	 * @param stOff Starting position (inclusive) of deletion
	 * @param edOff Ending position (exclusive) of deletion
	 * @throws IndexOutOfBoundsException If offset is invalid
	 * @throws UnsupportedOperationException If buffer is read-only
	 * @throws BufferReferenceException If deletion would affect any referencing child buffers.
	 */
	public void deleteFromFile(long stOff, long edOff)
	{
		deleteFromFile((int)stOff, (int)edOff);
	}
  
  /* ----- REPLACEMENT ----- */
  
	/**
	 * Replace a single byte in the buffer with another byte. Cannot replace if buffer is read-only.
	 * @param b Byte to replace existing byte with
	 * @param position Position of byte to replace
	 * @return Whether the operation was successful
	 * @throws IndexOutOfBoundsException If position is invalid
	 * @throws UnsupportedOperationException If buffer is read-only
	 * @throws BufferReferenceException If deletion would affect any referencing child buffers.
	 */
	public boolean replaceByte(byte b, int position)
	{
		if (this.readOnly()) throw new UnsupportedOperationException();
		if (position < 0 || position >= (int)this.getFileSize()) throw new IndexOutOfBoundsException();
		if (children != null)
		{
			FileBuffer[] ic = this.childrenEncompassing((long)position);
			if (ic != null) throw new BufferReferenceException(ic);
		}
		if (position >= capacity)
		{
			int oai = (capacity / position) - 1;
			int oi = capacity % position;
			if (oai <= this.overflow.size()) throw new IndexOutOfBoundsException();
			if (this.overflow.get(oai).length <= oi) throw new IndexOutOfBoundsException();
			this.overflow.get(oai)[oi] = b;
		}
		else this.contents[position] = b;
		return true;
	}
  
	/**
	 * Replace a single byte in the buffer with another byte. Cannot replace if buffer is read-only.
	 * @param b Byte to replace existing byte with
	 * @param position Position of byte to replace
	 * @return Whether the operation was successful
	 * @throws IndexOutOfBoundsException If position is invalid
	 * @throws UnsupportedOperationException If buffer is read-only
	 * @throws BufferReferenceException If deletion would affect any referencing child buffers.
	 */
	public boolean replaceByte(byte b, long position)
	{
		return this.replaceByte(b, (int)position);
	}
  
	/**
	 * Replace two bytes in the buffer with a new short value. Cannot replace if buffer is read-only.
	 * Byte order according to buffer's Endian-ness.
	 * @param s Short to replace existing two bytes with
	 * @param position Position of first byte to replace
	 * @return Whether the operation was successful
	 * @throws IndexOutOfBoundsException If position is invalid
	 * @throws UnsupportedOperationException If buffer is read-only
	 * @throws BufferReferenceException If deletion would affect any referencing child buffers.
	 */
	public boolean replaceShort(short s, int position)
	{
		return replaceShort(s, (long)position);
	}
  
	/**
	 * Replace two bytes in the buffer with a new short value. Cannot replace if buffer is read-only.
	 * Byte order according to buffer's Endian-ness.
	 * @param s Short to replace existing two bytes with
	 * @param position Position of first byte to replace
	 * @return Whether the operation was successful
	 * @throws IndexOutOfBoundsException If position is invalid
	 * @throws UnsupportedOperationException If buffer is read-only
	 * @throws BufferReferenceException If deletion would affect any referencing child buffers.
	 */
	public boolean replaceShort(short s, long position)
	{
		if (this.readOnly()) throw new UnsupportedOperationException();
		if (position < 0 || position >= this.getFileSize()) throw new IndexOutOfBoundsException();
		if (children != null)
		{
			this.checkAllBufferReferences(position, position + 2);
		}
		//if (!this.isFileFormatBE) s = switchByO(s);
		byte[] sBytes = numToByStr(s);
		if (!isFileFormatBE)
		{
			byte temp = sBytes[0];
			sBytes[0] = sBytes[1];
			sBytes[1] = temp;
		}
		if (sBytes == null) return false;
		if (sBytes.length != 2) return false;
		  
		for (int p = 0; p < sBytes.length; p++)
		{
			this.replaceByte(sBytes[p], position + (1 - p));
		}
		return true;
	}
  
	/**
	 * Replace three bytes in the buffer with the least three significant bytes of an int value. 
	 * Cannot replace if buffer is read-only.
	 * Byte order according to buffer's Endian-ness.
	 * @param s Int containing bytes to replace existing three bytes with
	 * @param position Position of first byte to replace
	 * @return Whether the operation was successful
	 * @throws IndexOutOfBoundsException If position is invalid
	 * @throws UnsupportedOperationException If buffer is read-only
	 * @throws BufferReferenceException If deletion would affect any referencing child buffers.
	 */
	public boolean replaceShortish(int s, int position)
	{
		return replaceShortish(s, (long)position);
	}
	
	/**
	 * Replace three bytes in the buffer with the least three significant bytes of an int value. 
	 * Cannot replace if buffer is read-only.
	 * Byte order according to buffer's Endian-ness.
	 * @param s Int containing bytes to replace existing three bytes with
	 * @param position Position of first byte to replace
	 * @return Whether the operation was successful
	 * @throws IndexOutOfBoundsException If position is invalid
	 * @throws UnsupportedOperationException If buffer is read-only
	 * @throws BufferReferenceException If deletion would affect any referencing child buffers.
	 */
	public boolean replaceShortish(int s, long position)
	{
		if (this.readOnly()) throw new UnsupportedOperationException();
		if (position < 0 || position >= this.getFileSize()) throw new IndexOutOfBoundsException();
		if (children != null)
		{
			this.checkAllBufferReferences(position, position + 3);
		}
		
		byte[] sBytes = numToByStr(s);
		if (sBytes == null) return false;
		
		if (sBytes.length != 4) return false;
		  
		if (this.isBigEndian())
		{
			for (int p = 0; p < 3; p++)
			{
				this.replaceByte(sBytes[p], position + (2 - p));
			}	
		}
		else
		{
			for (int p = 0; p < 3; p++)
			{
				this.replaceByte(sBytes[p], position + p);
			}
		}
		return true;
	}
	
	/**
	 * Replace four bytes in the buffer with a new int value. Cannot replace if buffer is read-only.
	 * Byte order according to buffer's Endian-ness.
	 * @param i Int to replace existing four bytes with
	 * @param position Position of first byte to replace
	 * @return Whether the operation was successful
	 * @throws IndexOutOfBoundsException If position is invalid
	 * @throws UnsupportedOperationException If buffer is read-only
	 * @throws BufferReferenceException If deletion would affect any referencing child buffers.
	 */
	public boolean replaceInt(int i, int position)
	{
		return replaceInt(i, (long)position);
	}
  
	/**
	 * Replace four bytes in the buffer with a new int value. Cannot replace if buffer is read-only.
	 * Byte order according to buffer's Endian-ness.
	 * @param i Int to replace existing four bytes with
	 * @param position Position of first byte to replace
	 * @return Whether the operation was successful
	 * @throws IndexOutOfBoundsException If position is invalid
	 * @throws UnsupportedOperationException If buffer is read-only
	 * @throws BufferReferenceException If deletion would affect any referencing child buffers.
	 */
	public boolean replaceInt(int i, long position)
	{
		if (this.readOnly()) throw new UnsupportedOperationException();
		if (position < 0 || position >= this.getFileSize() - 3) throw new IndexOutOfBoundsException();
		if (children != null)
		{
			this.checkAllBufferReferences(position, position + 4);
		}
		if (!this.isFileFormatBE) i = switchByO(i);
		byte[] iBytes = numToByStr(i);
		if (iBytes == null) return false;
		if (iBytes.length != 4) return false;
	  
		for (int p = 0; p < iBytes.length; p++)
		{
			this.replaceByte(iBytes[p], position + (3 - p));
		}
		return true;
	}
  
	/**
	 * Replace eight bytes in the buffer with a new long value. Cannot replace if buffer is read-only.
	 * Byte order according to buffer's Endian-ness.
	 * @param l Long to replace existing eight bytes with
	 * @param position Position of first byte to replace
	 * @return Whether the operation was successful
	 * @throws IndexOutOfBoundsException If position is invalid
	 * @throws UnsupportedOperationException If buffer is read-only
	 * @throws BufferReferenceException If deletion would affect any referencing child buffers.
	 */
	public boolean replaceLong(long l, int position)
	{
		return replaceLong(l, (long)position);
	}
  
	/**
	 * Replace eight bytes in the buffer with a new long value. Cannot replace if buffer is read-only.
	 * Byte order according to buffer's Endian-ness.
	 * @param l Long to replace existing eight bytes with
	 * @param position Position of first byte to replace
	 * @return Whether the operation was successful
	 * @throws IndexOutOfBoundsException If position is invalid
	 * @throws UnsupportedOperationException If buffer is read-only
	 * @throws BufferReferenceException If deletion would affect any referencing child buffers.
	 */
	public boolean replaceLong(long l, long position)
	{
		if (this.readOnly()) throw new UnsupportedOperationException();
		if (position < 0 || position >= this.getFileSize() - 7) throw new IndexOutOfBoundsException();
		if (children != null)
		{
			this.checkAllBufferReferences(position, position + 8);
		}
		if (!this.isFileFormatBE) l = switchByO(l);
		byte[] lBytes = numToByStr(l);
		if (lBytes == null) return false;
		if (lBytes.length != 8) return false;
	  
		for (int p = 0; p < lBytes.length; p++)
		{
			this.replaceByte(lBytes[p], position + (7 - p));
		}
		return true;
	}
  
  /* ----- REFERENCE CHILD HANDLING ----- */
  
  	/**
  	 * Exception thrown when attempting to write to a part of the buffer
  	 * that's being referenced by another read-only child buffer in a way
  	 * that would alter the child buffer's usage.
  	 * Can obtain list of offending child buffers from exception.
  	 * @author Blythe Hospelhorn
  	 * @version 1.0
  	 */
  	public static class BufferReferenceException extends RuntimeException
  	{
  		private static final long serialVersionUID = -6374378208716127441L;
  		private FileBuffer[] referenceList;
	  
  		/**
  		 * Construct a new BufferReferenceException with list of buffers that would be
  		 * affected.
  		 * @param children - FileBuffer array containing references to child buffers whose contents
  		 * would be altered if exception was not thrown.
  		 */
  		public BufferReferenceException(FileBuffer[] children)
	  {
		  this.referenceList = children;
	  }
	  
  		/**
  		 * Get the list of child buffer references affected by the events leading up to
  		 * the throwing of this exception.
  		 * @return FileBuffer reference array containing list of child buffers affected
  		 */
  		public FileBuffer[] getRefList()
	  {
		  return this.referenceList;
	  }
  	}

  	/**
  	 * Return a list of referenced child FileBuffers that start after (but not at) a given position.
  	 * @param position - Position to find all child buffers after
  	 * @return Array of child references after position, if any found. Return null if none found.
  	 */
  	public FileBuffer[] childrenAfter(long position)
  	{
  		if (this.children == null) return null;
  		if (this.children.isEmpty()) return null;
  		List<FileBuffer> after = new LinkedList<FileBuffer>();
  		for (ROSubFileBuffer f : this.children)
  		{
  			if (f.getStartOffset(this) > position) after.add(f);
  		}
  		if (after.isEmpty()) return null;
  		FileBuffer[] aa = new FileBuffer[after.size()];
  		return after.toArray(aa);
  	}
  	
  	/**
  	 * Return a list of referenced child FileBuffers that include a byte at given position.
  	 * @param position - Position in question.
  	 * @return Array of child references containing that position, if any found. Return null if none found.
  	 */
  	public FileBuffer[] childrenEncompassing(long position)
  	{
  		if (this.children == null) return null;
  		if (this.children.isEmpty()) return null;
  		List<FileBuffer> flist = new LinkedList<FileBuffer>();
  		for (ROSubFileBuffer f : this.children)
  		{
  			if (f.getStartOffset(this) <= position && f.getEndOffset(this) > position)
  			{
  				flist.add(f);
  			}
  		}
  		if (flist.isEmpty()) return null;
  		FileBuffer[] ia = new FileBuffer[flist.size()];
  		return flist.toArray(ia);
  	}
  	
  	/**
  	 * Get a list of the child buffers referencing this buffer's contents.
  	 * @return FileBuffer array containing references to all children of this buffer.
  	 */
  	public FileBuffer[] getChildren()
  	{
  		if (this.children == null) return null;
  		if (this.children.isEmpty()) return null;
  		FileBuffer[] flist = new FileBuffer[this.children.size()];
  		return children.toArray(flist);
  	}
  	
  	/**
  	 * Check both the "after" and "encompassing" lists to see if messing with bytes at or after
  	 * stOff will affect any child referencing buffers.
  	 * @param stOff Positions at and after which bytes in this buffer may be affected.
  	 * @return If it hasn't thrown an exception, it returns false.
  	 * @throws BufferReferenceException If it finds a conflict. Exception contains list of conflicting buffers.
  	 */
  	protected boolean checkAllBufferReferences(long stOff)
  	{
		FileBuffer[] ac = this.childrenAfter(stOff);
		FileBuffer[] ic = this.childrenEncompassing(stOff);
		if (ac == null)
		{
			if (ic == null) return false;
			else throw new BufferReferenceException(ic);
		}
		else
		{
			if (ic != null)
			{
				FileBuffer[] rlist = new FileBuffer[ic.length + ac.length];
				int i = 0;
				for (FileBuffer f : ic)
				{
					rlist[i] = f;
					i++;
				}
				for (FileBuffer f : ac)
				{
					rlist[i] = f;
					i++;
				}
				throw new BufferReferenceException(rlist);
			}
			else throw new BufferReferenceException(ac);
		}
  	}
  	
  	/**
  	 * Check for referencing child buffers that include bytes between two end points.
  	 * @param stOff First position in range (inclusive)
  	 * @param edOff Last position in range (exclusive)
  	 * @return False if no child buffers will be affected.
  	 * @throws BufferReferenceException If child buffer that will be affected is found.
  	 */
  	protected boolean checkAllBufferReferences(long stOff, long edOff)
  	{
  		/*Problematic buffers...
  		 * 1. Encompassing the start position
  		 * 2. Encompassing the end position
  		 * 3. After the start, but not encompassing or after end
  		 */
		FileBuffer[] asc = this.childrenAfter(stOff);
		FileBuffer[] isc = this.childrenEncompassing(stOff);
		FileBuffer[] iec = this.childrenEncompassing(edOff - 1);
		FileBuffer[] aec = this.childrenAfter(edOff - 1);
		if (asc == null && isc == null && iec == null) return false;
		else if (asc == null && isc == null && iec != null) throw new BufferReferenceException(iec);
		else if (asc == null && isc != null  && iec == null) throw new BufferReferenceException(isc);
		else if (asc != null && isc == null  && iec == null && aec == null) throw new BufferReferenceException(asc);
		Set<FileBuffer> comp = new HashSet<FileBuffer>();
		if (isc != null)
		{
			//Automatically added
			for (FileBuffer f : isc) comp.add(f);
		}
		if (iec != null)
		{
			//Added if not already in there from isc
			for (FileBuffer f : iec) comp.add(f);
		}
		if (asc != null)
		{
			//Added if not already in set or in aec
			boolean match = false;
			for (FileBuffer f : asc)
			{
				for (FileBuffer b : aec)
				{
					if (f == b) match = true;
				}
				if (match) break;
				else comp.add(f);
			}
		}
		if (!comp.isEmpty()) 
		{
			FileBuffer[] clist = new FileBuffer[comp.size()];
			throw new BufferReferenceException(comp.toArray(clist));
		}
		return false;
  	}
  	
  	protected void shiftReferencesAfter(long pos, long amnt)
  	{
  		FileBuffer[] ac = this.childrenAfter(pos);
  		if (ac != null)
  		{
  			for (FileBuffer f : ac)
  			{
  				if (f instanceof ROSubFileBuffer || f instanceof CompositeBuffer)
  				{
  					ROSubFileBuffer s = (ROSubFileBuffer)f;
  					s.setStartOffset(this, s.getStartOffset(this) + amnt);
  					s.setEndOffset(this, s.getEndOffset(this) + amnt);
  				}
  			}
  		}
  	}
  	
  	protected boolean hasChildren()
  	{
  		if (this.children == null) return false;
  		if (this.children.size() == 0) return false;
  		return true;
  	}
  	
  /* ----- WRITING TO DISK ----- */
  
  	/**
  	 * Write a file to disk using the path in the object.
  	 * Will fail if path is invalid.
  	 * Will create any directories specified in path that do not currently exist.
  	 * @throws IOException If path is invalid or file could not be written.
  	 */
  	public void writeFile() throws IOException
  	{
  		/*If valid, uses the path in the object*/
  		String myPath = this.getPath();
  		if (this.isOverflowing()) this.adjustBaseCapacityToSize();
  		writeFile(myPath, 0, this.fSize);	   
  	}
  
  	/**
  	 * Write a file to disk using the path string provided.
  	 * Will fail if path is invalid.
  	 * Will create any directories specified in path that do not currently exist.
  	 * @param path Path of file to write object to.
  	 * @throws IOException If path is invalid or file could not be written.
  	 */
  	public void writeFile(String path) throws IOException
  	{ 
  		/*Writes the full file to the given path.*/
  		this.writeFile(path, 0, (int)this.getFileSize()); 
  	}
  
  	/**
  	 * Write a file to disk using the path in the object.
  	 * Will fail if path is invalid.
  	 * Will create any directories specified in path that do not currently exist.
  	 * @param stPos Position within buffer of first byte to write (inclusive)
  	 * @throws IOException If path is invalid or file could not be written.
  	 */
  	public void writeFile(int stPos) throws IOException
  	{
  		this.writeFile((long)stPos); 	  
  	}
  
  	/**
  	 * Write a file to disk using the path in the object.
  	 * Will fail if path is invalid.
  	 * Will create any directories specified in path that do not currently exist.
  	 * @param stPos Position within buffer of first byte to write (inclusive)
  	 * @throws IOException If path is invalid or file could not be written.
  	 */
  	public void writeFile(long stPos) throws IOException
  	{
  		/*If valid, uses the path in the object
  		 * Writes file starting at the given index*/
  		if (stPos < 0) stPos = 0;
  		if (stPos >= this.getFileSize()) throw new ArrayIndexOutOfBoundsException();
	  
  		String myPath = this.getPath();
  		if (this.isOverflowing()) this.adjustBaseCapacityToSize(); 
  		writeFile(myPath, stPos, this.getFileSize());  
  	}
  
  	/**
  	 * Write a file to disk using the path in the object.
  	 * Will fail if path is invalid.
  	 * Will create any directories specified in path that do not currently exist.
  	 * @param stPos Position within buffer of first byte to write (inclusive)
  	 * @param edPos Position within buffer of first byte to not write (exclusive)
  	 * @throws IOException If path is invalid or file could not be written.
  	 */
  	public void writeFile(int stPos, int edPos) throws IOException
  	{
  		writeFile((long)stPos, (long)edPos);
  	}
  
  	/**
  	 * Write a file to disk using the path in the object.
  	 * Will fail if path is invalid.
  	 * Will create any directories specified in path that do not currently exist.
  	 * @param stPos Position within buffer of first byte to write (inclusive)
  	 * @param edPos Position within buffer of first byte to not write (exclusive)
  	 * @throws IOException If path is invalid or file could not be written.
  	 */
  	public void writeFile(long stPos, long edPos) throws IOException
  	{
	    if (stPos < 0) stPos = 0;
	    if (stPos >= edPos) throw new ArrayIndexOutOfBoundsException();
	    if (edPos > this.getFileSize()) edPos = this.getFileSize();
	    if (this.isOverflowing()) this.adjustBaseCapacityToSize();
	  
	    String myPath = this.getPath();
	    writeFile(myPath, stPos, edPos);
  	}
  
  	/**
  	 * Write a file to disk using the path string provided.
  	 * Will fail if path is invalid.
  	 * Will create any directories specified in path that do not currently exist.
  	 * @param path Path of file to write object to.
  	 * @param stPos Position within buffer of first byte to write (inclusive)
  	 * @throws IOException If path is invalid or file could not be written.
  	 */
  	public void writeFile(String path, int stPos) throws IOException
  	{
  		this.writeFile(path, (long)stPos, this.getFileSize()); 
  	}
  
  	/**
  	 * Write a file to disk using the path string provided.
  	 * Will fail if path is invalid.
  	 * Will create any directories specified in path that do not currently exist.
  	 * @param path Path of file to write object to.
  	 * @param stPos Position within buffer of first byte to write (inclusive)
  	 * @throws IOException If path is invalid or file could not be written.
  	 */
  	public void writeFile(String path, long stPos) throws IOException
  	{
  		writeFile(path, stPos, this.getFileSize());
  	}
  
  	/**
  	 * Write a file to disk using the path string provided.
  	 * Will fail if path is invalid.
  	 * Will create any directories specified in path that do not currently exist.
  	 * @param path Path of file to write object to.
  	 * @param stPos Position within buffer of first byte to write (inclusive)
  	 * @param edPos Position within buffer of first byte to not write (exclusive)
  	 * @throws IOException If path is invalid or file could not be written.
  	 */
  	public void writeFile(String path, int stPos, int edPos) throws IOException
 	{
  		writeFile(path, (long)stPos, (long)edPos);
 	}

  	/**
  	 * Write a file to disk using the path string provided.
  	 * Will fail if path is invalid.
  	 * Will create any directories specified in path that do not currently exist.
  	 * @param path Path of file to write object to.
  	 * @param stPos Position within buffer of first byte to write (inclusive)
  	 * @param edPos Position within buffer of first byte to not write (exclusive)
  	 * @throws IOException If path is invalid or file could not be written.
  	 */
  	public void writeFile(String path, long stPos, long edPos) throws IOException
  	{
	    if (path == null) throw new IOException();
	    if (stPos < 0) stPos = 0;
	    if (stPos >= edPos) throw new IOException();
	    if (edPos > this.getFileSize()) edPos = this.getFileSize();
	    if (edPos > this.capacity) this.adjustBaseCapacityToSize();
	    
	    String dir = FileBuffer.chopPathToDir(path);
	    if(!FileBuffer.directoryExists(dir))
	    {
	    	if(!new File(dir).mkdirs()) throw new IOException();
	    }
		FileOutputStream myStream = new FileOutputStream(path);
		//These will always be in the int range if a standard FileBuffer!!
		//Therefore, for subclasses that handle larger files, this must be overridden!
		int iSt = (int)stPos;
		int iLen = (int)(edPos - stPos);
		myStream.write(this.contents, iSt, iLen);
		myStream.close();	
  	}
  	
  	/**
  	 * Append the contents of this FileBuffer to an existing file, if possible.
  	 * @param path Path of existing file to append to
  	 * @throws IOException If target file could not be opened or written to
  	 * @throws NoSuchFileException If target file does not exist
  	 */
  	public void appendToFile(String path) throws IOException, NoSuchFileException
  	{
  		this.appendToFile(path, 0, this.getFileSize());
  	}
  	
  	/**
  	 * Append the contents of this FileBuffer to an existing file, if possible.
  	 * @param path Path of existing file to append to
  	 * @param stPos Offset of first byte in buffer to write (inclusive)
  	 * @throws IOException If target file could not be opened or written to
  	 * @throws NoSuchFileException If target file does not exist
  	 */
  	public void appendToFile(String path, long stPos) throws IOException, NoSuchFileException
  	{
  		this.appendToFile(path, stPos, this.getFileSize());
  	}
  	
  	/**
  	 * Append the contents of this FileBuffer to an existing file, if possible.
  	 * @param path Path of existing file to append to
  	 * @param stPos Offset of first byte in buffer to write (inclusive)
  	 * @param edPos Offset of first byte in buffer to not write at end (exclusive)
  	 * @throws IOException If target file could not be opened or written to
  	 * @throws NoSuchFileException If target file does not exist
  	 */
  	public void appendToFile(String path, long stPos, long edPos) throws IOException, NoSuchFileException
  	{
  		if (!this.offsetValid(stPos)) throw new ArrayIndexOutOfBoundsException();
  		if (!this.offsetValid(edPos - 1)) throw new ArrayIndexOutOfBoundsException();
  		if (stPos >= edPos) throw new ArrayIndexOutOfBoundsException();
  		
  		if (stPos == 0 && edPos == this.getFileSize()) 
  		{
  			Files.write(Paths.get(path), contents, StandardOpenOption.APPEND);
  			return;
  		}
  		
  		int len = (int)(edPos - stPos);
  		byte[] b = new byte[len];
  		for (int i = 0; i < len; i++) b[i] = this.getByte(stPos + i);
  		Files.write(Paths.get(path), b, StandardOpenOption.APPEND);
  	}
  
  	/**
  	 * Copy the contents of an InputStream to a file on disk. This can be used to copy
  	 * files from inside a JAR to disk.
  	 * @param path Path of output file.
  	 * @param stream InputStream containing data to write to disk.
  	 * @throws IOException If there is an error writing to disk.
  	 */
  	public static void copyInputStreamToDisk(String path, InputStream stream) throws IOException
  	{
  		if (stream == null) return;
  		if (path == null || path.isEmpty()) return;
  		FileWriter fw = new FileWriter(path);
  		BufferedWriter bw = new BufferedWriter(fw);
  		
  		int i = -1;
		while ((i = stream.read()) != -1) bw.write(i);
  		
  		bw.close();
  		fw.close();
  	}
  	
  /* ----- STRING SEARCHING ----- */
  
  	/**
  	 * Search for a string (assuming default character set) within this buffer.
  	 * @param stPos Start offset (inclusive) of search
  	 * @param edPos End offset (exclusive) of search
  	 * @param query String to search for
  	 * @return Offset of string, if found. If not found, -1
  	 * @throws IndexOutOfBoundsException If positions are invalid
  	 */
  	public int findString(int stPos, int edPos, String query)
  	{
  		return findString(stPos, edPos, query.getBytes());
  	}
  
  	/**
  	 * Search for a string (assuming default character set) within this buffer.
  	 * @param stPos Start offset (inclusive) of search
  	 * @param edPos End offset (exclusive) of search
  	 * @param query String to search for
  	 * @return Offset of string, if found. If not found, -1
  	 * @throws IndexOutOfBoundsException If positions are invalid
  	 */
  	public long findString(long stPos, long edPos, String query)
  	{
  		return findString(stPos, edPos, query.getBytes());
  	}	
  
  	/**
  	 * Search for a byte string within this buffer.
  	 * @param stPos Start offset (inclusive) of search
  	 * @param edPos End offset (exclusive) of search
  	 * @param query String to search for
  	 * @return Offset of string, if found. If not found, -1
  	 * @throws IndexOutOfBoundsException If positions are invalid
  	 */
  	public int findString(int stPos, int edPos, byte[] query)
  	{
  		return (int)findString((long)stPos, (long)edPos, query);
  	}
  
  	/**
  	 * Search for a byte string within this buffer.
  	 * @param stPos Start offset (inclusive) of search
  	 * @param edPos End offset (exclusive) of search
  	 * @param query String to search for
  	 * @return Offset of string, if found. If not found, -1
  	 * @throws IndexOutOfBoundsException If start position is invalid
  	 */
  	public long findString(long stPos, long edPos, byte[] query)
  	{
  		long off = -1;
		int sLen = query.length;
		int c = 0;
	  
		long fsz = this.getFileSize();
		if (stPos < 0 || stPos >= fsz) throw new IndexOutOfBoundsException();
		if (edPos < 0 || edPos <= stPos) throw new IndexOutOfBoundsException();
		if (edPos > fsz) edPos = fsz;
		if (sLen > (int)fsz) return -1; //Same as can't find.
	  
		for (long i = stPos; i < edPos; i++)
		{
			if (this.getByte(i) == query[c]) c++;
			else c = 0;
			if (c >= sLen)
			{
				off = i + 1 - sLen;
				break;
			}
		}
		return off; 
  	}
  
  	/**
  	 * Count the number of times a string occurs in the file buffer (assuming default character set).
  	 * @param stPos Start offset (inclusive) of search
  	 * @param edPos End offset (exclusive) of search
  	 * @param query String to search for
  	 * @return Number of times string occurred in buffer between given offsets
  	 * @throws IndexOutOfBoundsException If positions are invalid
  	 */
  	public int countStringOcc(int stPos, int edPos, String query)
  	{
  		return this.countStringOcc((long)stPos, (long)edPos, query);
  	}
  
  	/**
  	 * Count the number of times a string occurs in the file buffer (assuming default character set).
  	 * @param stPos Start offset (inclusive) of search
  	 * @param edPos End offset (exclusive) of search
  	 * @param query String to search for
  	 * @return Number of times string occurred in buffer between given offsets
  	 * @throws IndexOutOfBoundsException If positions are invalid
  	 */
  	public int countStringOcc(long stPos, long edPos, String query)
  	{
  		return countStringOcc(stPos, edPos, query.getBytes());
  	}
  
  	/**
  	 * Count the number of times a byte string occurs in the file buffer.
  	 * @param stPos Start offset (inclusive) of search
  	 * @param edPos End offset (exclusive) of search
  	 * @param query String to search for
  	 * @return Number of times string occurred in buffer between given offsets
  	 * @throws IndexOutOfBoundsException If positions are invalid
  	 */
  	public int countStringOcc(long stPos, long edPos, byte[] query)
  	{
  		int count = 0;
  		int sLen = query.length;
  		int c = 0;
	  
  		long fsz = this.getFileSize();
  		if (stPos < 0 || stPos >= fsz) throw new IndexOutOfBoundsException();
  		if (edPos < 0 || edPos <= stPos) throw new IndexOutOfBoundsException();
  		if (edPos > fsz) edPos = fsz;
  		if (sLen > (int)fsz) return 0;
	  
  		for (long i = stPos; i < edPos; i++)
  		{
  			if (this.getByte(i) == query[c]) c++;
  			else c = 0;
  			if (c >= sLen)
  			{
  				count++;
  				c = 0;
  			}
  		}
  		return count;
  	}
  	
  	public long findEncodedString(String charset, long stPos, long edPos, String query)
  	{
  		// Convert target to byte array using charset and do usual search
  		Charset mySet = Charset.forName(charset);
  		ByteBuffer bb = mySet.encode(query);
  		byte[] queryBytes = bb.array();
  		return findString(stPos, edPos, queryBytes);
  	}
  	
  /* ----- PATH STRING HANDLING ----- */
  
  	/**
  	 * Cut a file name from a path - cut to the last slash
  	 * @param path Path string to trim
  	 * @return String representing directory path
  	 */
  	public static String chopPathToDir(String path)
  	{
  		int lastSlash = 0;
  		if (path == null) return null;
  		lastSlash = path.lastIndexOf(File.separator);
	
  		if (lastSlash < 0 || lastSlash >= path.length()) return path;
	
  		return path.substring(0, lastSlash);
  	}
  
  	/**
  	 * Cut everything before the last slash from a path string - ie. Cut out directory path
  	 * @param path Path string to trim
  	 * @return String representing only the file name and extension
  	 */
  	public static String chopDirFromPath(String path)
  	{
  		int lastSlash = 0;
		
  		if (path == null) return null;
  		lastSlash = path.lastIndexOf(File.separator);
  		if (lastSlash < 0 || lastSlash >= path.length()) return path;
  		
  		return path.substring(lastSlash + 1, path.length());
  }
  
  	/**
  	 * Cut everything but a substring after the final '.' character from a path string.
  	 * @param path Path string to trim
  	 * @return Substring of path string ideally representing the extension
  	 */
  	public static String chopPathToExt(String path)
  	{
  		int lastDot = 0;
  		if (path == null) return null;
  		lastDot = path.lastIndexOf('.');	
  		if (lastDot < 0 || lastDot >= path.length()) return path;
  		return path.substring(lastDot + 1, path.length());
  	}
  
  	/**
  	 * Cut a substring after the final '.' character from a path string - ie. cut out the extension
  	 * @param path Path string to trim
  	 * @return Substring without final '.' and characters following.
  	 */
  	public static String chopExtFromPath(String path)
  	{
  		int lastDot = 0;
		if (path == null) return null;	
		lastDot = path.lastIndexOf('.');		
		if (lastDot < 0 || lastDot >= path.length()) return path;	
		return path.substring(0, lastDot);
  }
  
  	/**
  	 * Cut path to a substring after the final slash and before the final '.'
  	 * @param path Path string to trim
  	 * @return Substring without directory path and extension.
  	 */
  	public static String chopPathToFName(String path)
  	{
  		if (path == null) return null;
  		return chopDirFromPath(chopExtFromPath(path));
  	}

  	/**
  	 * Generate a path from the stored static tempDir string for dumping temporary files.
  	 * @param name Desired name of temporary file.
  	 * @return Full file path for a temporary file.
  	 * @throws IOException If there is an error generating the path.
  	 */
  	public static String generateTemporaryPath(String name) throws IOException
  	{
  		if (tempDir == null) generateTempDir();
  		if (tempDir.isEmpty()) generateTempDir();
  		String myPath = "";
  		Random rand = new Random();
  		int rando = rand.nextInt();
  		myPath += tempDir + File.separator + ".~";
  		if (name != null) myPath += name + "_";
  		myPath += Integer.toHexString(rando) + ".tmp";
  		return myPath;
  	}
  	
  	/**
  	 * Set the tempDir String to a preferred path. This will affect the
  	 * generateTemporaryPath(String) function.
  	 * @param path The desired default temporary files path.
  	 */
  	public static void setTempDir(String path)
  	{
  		FileBuffer.tempDir = path;
  	}
  	
  	/**
  	 * Query the operating system to generate a path for temporary file
  	 * storage. This is automatically done the first time generateTemporaryPath is called, 
  	 * if the tempDir hasn't been set manually.
  	 * WARNING: There is no guarantee that this generates a valid path.
  	 * @throws IOException If there is an error finding the temp folder.
  	 */
  	public static void generateTempDir() throws IOException
  	{
  		File temp = File.createTempFile("FileBuffer_generateTempDir", ".tmp");
  		String aPath = temp.getAbsolutePath();
  		tempDir = aPath.substring(0, aPath.lastIndexOf(File.separatorChar));
  	}
  	
  	/**
  	 * Get the temp folder currently set for FileBuffer classes. If there isn't one set
  	 * by the user, this function will find one by querying the operating system.
  	 * WARNING: There is no guarantee that this returns a valid path.
  	 * @throws IOException If there is an error finding the temp folder.
  	 */
  	public static String getTempDir() throws IOException
  	{
  		if (FileBuffer.tempDir == null || FileBuffer.tempDir.isEmpty()) generateTempDir();
  		return FileBuffer.tempDir;
  	}
  	
  /* ----- STATUS CHECKERS ----- */
  
  	/**
  	 * Check whether buffer has any contents referenced.
  	 * @return True if buffer is not empty. False if buffer has no contents.
  	 */
  	public boolean isEmpty()
  	{
  		if (this.contents == null) return true;
  		else return false;
  	}
  
  	/**
  	 * Check whether input offset refers to a valid byte index within this buffer.
  	 * This EXCLUDES the file size.
  	 * @param off Offset to check
  	 * @return True if offset can be used to access a byte in this buffer. False otherwise.
  	 */
  	public boolean offsetValid(int off)
  	{
  		/*WARNING: if you intend to use the file size
  		 * for edOff (that is, read until the end of the file),
  		 * this WILL reject it as an invalid offset.
  		 * This function is to check INDEXING validity.*/
  		if (off >= this.fSize || off < 0) return false;
  		return true;
  	}
  
  	/**
  	 * Check whether input offset refers to a valid byte index within this buffer.
  	 * This EXCLUDES the file size.
  	 * @param off Offset to check
  	 * @return True if offset can be used to access a byte in this buffer. False otherwise.
  	 */
  	public boolean offsetValid(long off)
  	{
  		if (off < 0) return false;
  		if (off >= this.getFileSize()) return false;
  		return true;
  	}
  
  	/**
  	 * Check whether buffer size has exceeded base capacity, that is, the size of the 
  	 * initial internal array.
  	 * If this is the case, it may make access to overflowing values and writing to disk
  	 * slower.
  	 * @return True if buffer file size has exceeded base capacity. False otherwise.
  	 */
  	public boolean isOverflowing()
  	{
  		if (this.overflow.size() > 0 && this.fSize > this.capacity) return true;
  		return false;
  	}
  
  	/**
  	 * Check whether this buffer is locked for writing.
  	 * @return True if file buffer cannot be written to or altered. False otherwise.
  	 */
  	public boolean readOnly()
  	{
  		return this.readOnly;
  	}
  
  	/**
  	 * Check a pair of offsets (a start and an end) at the same time
  	 * to see whether they are valid for this buffer.
  	 * @param stOff Start offset to check. Must be a valid index.
  	 * @param edOff End offset to check. Must be a valid index, or the size
  	 * of the buffer.
  	 * @return True - If both are valid.
  	 * <br>False - If edOff is not a valid index, but the size of the buffer (usually valid).
  	 * @throws IndexOutOfBoundsException If either is invalid.
  	 */
  	public boolean checkOffsetPair(long stOff, long edOff)
  	{
  		long fSz = this.getFileSize();
  		if (stOff >= edOff) throw new IndexOutOfBoundsException();
  		if (stOff < 0) throw new IndexOutOfBoundsException();
  		if (edOff > fSz) throw new IndexOutOfBoundsException();
  		if (edOff == fSz) return false;
  		return true;
  	}
  	
  /* ----- STATIC FILE INFORMATION ----- */

  	/**
  	 * Retrieve the size of a file on disk at the specified path.
  	 * @param path Path to file on disk
  	 * @return Size of file if path is valid and file exists. 0 otherwise.
  	 */
  	public static long fileSize(String path)
  	{
  		if (path == null) return 0;
  		File myFile = new File(path);
  		if (myFile.isFile()) return myFile.length();
  		return 0;
  	}
  
  	/**
  	 * Check whether a file at a specified path exists on disk.
  	 * @param path Path to file on disk
  	 * @return True if file exists and is not a directory. 
  	 * False if file does not exist or is a directory.
  	 */
  	public static boolean fileExists(String path)
  	{
  		if (path == null) return false;

  		File myFile = new File(path);
  		if (myFile.isFile()) return true;
  		return false;
  	}
  
  	/**
  	 * Check whether a directory at a specified path exists on disk.
  	 * @param path Path to directory on disk
  	 * @return True if directory exists. 
  	 * False if path is invalid or directory does not exist.
  	 */
  	public static boolean directoryExists(String path)
  	{
  		if (path == null) return false;
	  
  		File myFile = new File(path);
  		if (myFile.isDirectory()) return true;
  		return false;
  	}

  /* ----- STRING FORMATTING ----- */
  
  	/**
  	 * Get a string of a byte represented by two hexadecimal digits.
  	 * @param aByte Byte to represent in hex
  	 * @return String of length two representing byte in hex
  	 */
  	public static String byteToHexString(byte aByte)
  	{
  		String hex = "";
  		int iByte = Byte.toUnsignedInt(aByte);
  		if (iByte < 0x10) hex += "0";
  		hex += Integer.toHexString(iByte);
  		return hex;
  	}

  /* ----- STRING ADDITON/ RETRIEVAL ----- */
  
  	/**
  	 * Add an ASCII string to the end of a file.
  	 * This function is intended for use by ASCII Strings. Utilizes the default character set.
  	 * @param inString The ASCII string to add.
  	 * @throws UnsupportedOperationException If buffer is read-only
  	 */
  	public void printASCIIToFile(String inString)
  	{
  		if (this.readOnly()) throw new UnsupportedOperationException();
  		byte[] s = inString.getBytes(); 
  		for (int i = 0; i < s.length; i++) this.addToFile(s[i]);
  	}
  
  	/**
  	 * Get the string of len bytes after the specified position, read as an ASCII character string,
  	 * and return a Java String.
  	 * NOTE: If a null character ('\0' or 0), or character with ASCII value less than 0x20
  	 * which is not a line break are encountered, the function will stop adding further characters
  	 * to the string and return!
  	 * @param pos Position in file to start
  	 * @param len Number of bytes to read as ASCII string
  	 * @return String which is an ASCII representation of the bytes in question
  	 * @throws IndexOutOfBoundsException If position or length are invalid
  	 */
  	public String getASCII_string(int pos, int len)
  	{
  		return getASCII_string((long)pos, len);
  	}
  
  	/**
  	 * Get the string of len bytes after the specified position, read as an ASCII character string,
  	 * and return a Java String.
  	 * NOTE: If a null character ('\0' or 0), or character with ASCII value less than 0x20
  	 * which is not a line break are encountered, the function will stop adding further characters
  	 * to the string and return!
  	 * @param pos Position in file to start
  	 * @param len Number of bytes to read as ASCII string
  	 * @return String which is an ASCII representation of the bytes in question
  	 * @throws IndexOutOfBoundsException If position or length are invalid
  	 */
  	public String getASCII_string(long pos, int len)
  	{
  		if (pos < 0 || pos >= this.getFileSize()) throw new IndexOutOfBoundsException();
  		if (len < 0) throw new IndexOutOfBoundsException();
  		if (pos + len >= this.getFileSize()) throw new IndexOutOfBoundsException();
	  
  		String s = "";
	  
  		for (int i = 0; i < len; i++)
  		{
  			byte b = this.getByte(pos + i);
  			if (b == 0) break;
  			char c = (char)(Byte.toUnsignedInt(b));
  			if (c == 0) break;
  			if (c < 32 && !(c == 0x09 || c == 0x0A || c == 0x0B || c == 0x0D)) break;
  			s += c;
  			//s += b;
  		}
	  
  		return s;
  	}
  
  	/**
  	 * Read a byte string in the file starting at position pos and ending the first
  	 * time a specific character is encountered.
  	 * NOTE: If a null character ('\0' or 0), or character with ASCII value less than 0x20
  	 * which is not a line break are encountered, the function will stop adding further characters
  	 * to the string and return, EVEN IF the character causing the termination is not the one
  	 * provided in the arguments!
  	 * @param pos Position in file to start
  	 * @param endmarker Character to terminate string read (exclusive - string will not include this character!)
  	 * @return String which is an ASCII representation of the bytes in question
  	 * @throws IndexOutOfBoundsException If position or length are invalid
  	 */
  	public String getASCII_string(int pos, char endmarker)
  	{
  		return this.getASCII_string((long)pos, endmarker);
  	}
  
  	/**
  	 * Read a byte string in the file starting at position pos and ending the first
  	 * time a specific character is encountered.
  	 * NOTE: If a null character ('\0' or 0), or character with ASCII value less than 0x20
  	 * which is not a line break are encountered, the function will stop adding further characters
  	 * to the string and return, EVEN IF the character causing the termination is not the one
  	 * provided in the arguments!
  	 * @param pos Position in file to start
  	 * @param endmarker Character to terminate string read (exclusive - string will not include this character!)
  	 * @return String which is an ASCII representation of the bytes in question
  	 * @throws IndexOutOfBoundsException If position or length are invalid
  	 */
  	public String getASCII_string(long pos, char endmarker)
  	{
  		if (pos < 0 || pos >= this.getFileSize()) throw new IndexOutOfBoundsException();
  	  
  		String s = "";
  		int i = 0;
  		char nowByte = (char)FileBuffer.UByteToUShort(this.getByte(pos));
	  
  		while (nowByte != endmarker && (pos + i) < (int)this.getFileSize())
	  	{
  			if (nowByte == 0) break;
  			if (nowByte < 32 && !(nowByte == 0x09 || nowByte == 0x0A || nowByte == 0x0B || nowByte == 0x0D)) break;
  			s += nowByte;
  			//System.out.println("Added character " + nowByte + " to string.");
  			i++;
  			nowByte = (char)(Byte.toUnsignedInt(this.getByte(pos + i)));
  			//s += b;
	  	}
	  
  		return s;
  	}
  
  	/**
  	 * Read a string of bytes from buffer between stPos and edPos and interpret
  	 * as a string encoded by a provided character set.
  	 * @param charset - Canonical String name of character set to decode string as
  	 * <p> Names for common sets: <br>
  	 * ASCII : "ASCII" <br>
  	 * UTF-8 : "UTF8" <br>
  	 * Shift_JIS : "SJIS" <br>
  	 * Shift_JISX0213 : "x-SJIS_0213" <br>
  	 * @param stPos First byte of string in buffer (inclusive)
  	 * @param edPos Byte after last byte of string in buffer (exclusive)
  	 * @return String interpreted from bytes in buffer using charset
  	 * @throws IndexOutOfBoundsException If positions are invalid.
  	 * @throws NullPointerException If charset string is null
  	 * @throws IllegalCharsetNameException If the given charset name is illegal
  	 * @throws UnsupportedCharsetException If no support for the named charset is available in this instance of the Java virtual machine
  	 * @see java.nio.charset.Charset
  	 */
  	public String readEncoded_string(String charset, int stPos, int edPos)
  	{
  		return readEncoded_string(charset, (long)stPos, (long)edPos);
  	}
  
  	/**
  	 * Read a string of bytes from buffer between stPos and edPos and interpret
  	 * as a string encoded by a provided character set.
  	 * @param charset - Canonical String name of character set to decode string as
  	 * <p> Names for common sets: <br>
  	 * ASCII : "ASCII" <br>
  	 * UTF-8 : "UTF8" <br>
  	 * Shift_JIS : "SJIS" <br>
  	 * Shift_JISX0213 : "x-SJIS_0213" <br>
  	 * @param stPos First byte of string in buffer (inclusive)
  	 * @param edPos Byte after last byte of string in buffer (exclusive)
  	 * @return String interpreted from bytes in buffer using charset
  	 * @throws IndexOutOfBoundsException If positions are invalid.
  	 * @throws NullPointerException If charset string is null
  	 * @throws IllegalCharsetNameException If the given charset name is illegal
  	 * @throws UnsupportedCharsetException If no support for the named charset is available in this instance of the Java virtual machine
  	 * @see java.nio.charset.Charset
  	 */
  	public String readEncoded_string(String charset, long stPos, long edPos)
  	{
  		if (charset == null) throw new NullPointerException();
  		if (stPos < 0) throw new IndexOutOfBoundsException();
  		if (edPos >= this.getFileSize()) throw new IndexOutOfBoundsException();
  		if (stPos >= edPos) throw new IndexOutOfBoundsException();
  		Charset mySet = Charset.forName(charset);
  		CharBuffer cb = mySet.decode(this.toByteBuffer(stPos, edPos));
  		String s = cb.toString();
  		return s;
  	}
  
  	public String readEncoded_string(String charset, long stPos, String endmarker)
  	{
  		// Convert endmarker to byte array (using charset), look for that, then read up to that.
  		Charset mySet = Charset.forName(charset);
  		ByteBuffer bb = mySet.encode(endmarker);
  		byte[] endBytes = bb.array();
  		long fSz = this.getFileSize();
  		long queryEnd = this.findString(stPos, fSz, endBytes);
  		if (queryEnd < 0) queryEnd = fSz;
  		return readEncoded_string(charset, stPos, queryEnd);
  	}
  	
  	/**
  	 * Encode a string using a given character set into a string of bytes
  	 * and add to the end of the file buffer.
  	 * @param charset Canonical String name of character set to encode string as.
  	 * <p> Names for common sets: <br>
  	 * ASCII : "ASCII" <br>
  	 * UTF-8 : "UTF8" <br>
  	 * Shift_JIS : "SJIS" <br>
  	 * Shift_JISX0213 : "x-SJIS_0213" <br>
  	 * @param myString String to encode
  	 * @throws UnsupportedOperationException If buffer is read-only.
  	 * @throws NullPointerException If charset string or provided string is null
  	 * @throws IllegalCharsetNameException If the given charset name is illegal
  	 * @throws UnsupportedCharsetException If no support for the named charset is available in this instance of the Java virtual machine
  	 * @see java.nio.charset.Charset
  	 */
  	public void addEncoded_string(String charset, String myString)
  	{
  		if (this.readOnly()) throw new UnsupportedOperationException();
  		if (charset == null || myString == null) throw new NullPointerException();
  		Charset mySet = Charset.forName(charset);
  		ByteBuffer bb = mySet.encode(myString);
  		while(bb.hasRemaining())
  		{
  			byte b = bb.get();
  			this.addToFile(b);
  		}
  	}
  
  	/**
  	 * Encode a string using a given character set into a string of bytes
  	 * and add to the end of the file buffer.
  	 * @param charset Canonical String name of character set to encode string as.
  	 * <p> Names for common sets: <br>
  	 * ASCII : "ASCII" <br>
  	 * UTF-8 : "UTF8" <br>
  	 * Shift_JIS : "SJIS" <br>
  	 * Shift_JISX0213 : "x-SJIS_0213" <br>
  	 * @param myString String to encode
  	 * @throws UnsupportedOperationException If buffer is read-only.
  	 * @throws NullPointerException If charset string or provided string is null
  	 * @throws IllegalCharsetNameException If the given charset name is illegal
  	 * @throws UnsupportedCharsetException If no support for the named charset is available in this instance of the Java virtual machine
  	 * @see java.nio.charset.Charset
  	 */
  	public void addEncoded_string(String charset, String myString, int pos)
  	{
 		if (this.readOnly()) throw new UnsupportedOperationException();
 		if (charset == null || myString == null) throw new NullPointerException();
  		this.addEncoded_string(charset, myString, (long)pos);
  	}
  
  	/**
  	 * Encode a string using a given character set into a string of bytes
  	 * and insert into the file buffer at the specified position.
  	 * @param charset Canonical String name of character set to encode string as.
  	 * <p> Names for common sets: <br>
  	 * ASCII : "ASCII" <br>
  	 * UTF-8 : "UTF8" <br>
  	 * Shift_JIS : "SJIS" <br>
  	 * Shift_JISX0213 : "x-SJIS_0213" <br>
  	 * @param myString String to encode
  	 * @throws UnsupportedOperationException If buffer is read-only.
  	 * @throws IndexOutOfBoundsException If position is invalid.
   	 * @throws NullPointerException If charset string or provided string is null
  	 * @throws IllegalCharsetNameException If the given charset name is illegal
  	 * @throws UnsupportedCharsetException If no support for the named charset is available in this instance of the Java virtual machine
  	 * @see java.nio.charset.Charset
  	 */
  	public void addEncoded_string(String charset, String myString, long pos)
  	{
  		if (this.readOnly()) throw new UnsupportedOperationException();
  		if (charset == null || myString == null) throw new NullPointerException();
  		if (pos == this.getFileSize())
  		{
  			this.addEncoded_string(charset, myString);
  			return;
  		}
  		if (!this.offsetValid(pos)) throw new IndexOutOfBoundsException();
  		Charset mySet = Charset.forName(charset);
  		ByteBuffer bb = mySet.encode(myString);
  		this.moveAllBytesUp(pos, bb.remaining());
  		while(bb.hasRemaining())
  		{
  			byte b = bb.get();
  			this.replaceByte(b, pos);
  			pos++;
  		}
  	}
  
  /* ----- CONVERSION ----- */
  
  	/**
  	 * Get a ByteBuffer version of the entire file buffer.
  	 * @return ByteBuffer with contents of file buffer
  	 * @throws NullPointerException If this buffer has no contents
  	 */
  	public ByteBuffer toByteBuffer()
  	{
  		if (this.isEmpty()) throw new NullPointerException();
  		ByteBuffer bb = ByteBuffer.allocate(this.fSize);
  		if (this.isOverflowing()) this.adjustBaseCapacityToSize();
  		bb.put(this.contents, 0, this.fSize);
  		bb.rewind();
  		return bb;
  	}
  
  	/**
  	 * Get a ByteBuffer of the bytes in this buffer between stPos and edPos
  	 * @param stPos Position of first byte in target ByteBuffer (inclusive)
  	 * @param edPos Position after last byte in target ByteBuffer (exclusive)
  	 * @return ByteBuffer containing desired bytes in file buffer.
  	 * @throws IndexOutOfBoundsException If positions are invalid.
  	 * For the standard FileBuffer, invalid positions include those greater than the maximum
  	 * positive int (0x7FFFFFFF).
  	 * @throws NullPointerException If this buffer has no contents.
  	 */
  	public ByteBuffer toByteBuffer(int stPos, int edPos)
  	{
  		if (this.isEmpty()) throw new NullPointerException();
  		if (stPos < 0) throw new IndexOutOfBoundsException();
  		if (edPos > this.fSize) throw new IndexOutOfBoundsException();
  		if (stPos >= edPos) throw new IndexOutOfBoundsException();
  		ByteBuffer bb = ByteBuffer.allocate(edPos - stPos);
  		if (edPos > this.capacity) this.adjustBaseCapacityToSize();
  		bb.put(this.contents, stPos, edPos - stPos);
  		bb.rewind();
  		return bb;
  	}
  
  	/**
  	 * Get a ByteBuffer of the bytes in this buffer between stPos and edPos
  	 * @param stPos Position of first byte in target ByteBuffer (inclusive)
  	 * @param edPos Position after last byte in target ByteBuffer (exclusive)
  	 * @return ByteBuffer containing desired bytes in file buffer.
  	 * @throws IndexOutOfBoundsException If positions are invalid.
  	 * For the standard FileBuffer, invalid positions include those greater than the maximum
  	 * positive int (0x7FFFFFFF).
  	 * @throws NullPointerException If this buffer has no contents.
  	 */
  	public ByteBuffer toByteBuffer(long stPos, long edPos)
  	{
  		long sz = edPos - stPos;
  		if (sz > 0x7FFFFFFFL) throw new IndexOutOfBoundsException();
  		return this.toByteBuffer((int)stPos, (int)edPos);
  	}
  
  	/**
  	 * Create a new copy of the file buffer. This creates a fresh copy - it does not reference
  	 * the existing file buffer. It should thus be noted that consumption of memory will
  	 * increase accordingly.
  	 * @param stPos Starting position of copy procedure (inclusive)
  	 * @param edPos Ending position of copy procedure (exclusive)
  	 * @return New FileBuffer whose contents are identical to this one.
  	 * @throws IOException If copying requires disk access and disk access failed for any reason.
  	 * @throws IndexOutOfBoundsException If positions are invalid.
  	 */
  	public FileBuffer createCopy(int stPos, int edPos) throws IOException
  	{
  		if (stPos < 0) stPos = 0;
  		if (edPos > (int)this.getFileSize()) edPos = (int)this.getFileSize();
  		if (edPos <= stPos) throw new IndexOutOfBoundsException();
	  
  		FileBuffer f = new FileBuffer(edPos - stPos, this.isBigEndian());
  		f.addToFile(this, stPos, edPos);
	  
  		return f;
  	}
  
  	/**
  	 * Create a new copy of the file buffer. This creates a fresh copy - it does not reference
  	 * the existing file buffer. It should thus be noted that consumption of memory will
  	 * increase accordingly.
  	 * @param stPos Starting position of copy procedure (inclusive)
  	 * @param edPos Ending position of copy procedure (exclusive)
  	 * @return New FileBuffer whose contents are identical to this one.
  	 * @throws IOException If copying requires disk access and disk access failed for any reason.
  	 * @throws IndexOutOfBoundsException If positions are invalid.
  	 */
  	public FileBuffer createCopy(long stPos, long edPos) throws IOException
  	{
  		if (stPos < 0) throw new IndexOutOfBoundsException();
  		return createCopy((int)stPos, (int)edPos);
  	}
  
  	/**
  	 * Create a Read-Only FileBuffer that references a piece of this buffer.
  	 * <br> WARNING: Because this derived buffer references the current buffer, it will lock 
  	 * this buffer for editing at all positions referenced by a derived buffer! This derived buffer
  	 * will need to be delinked before these regions can be edited again.
  	 * @param stPos First byte of derived buffer (inclusive)
  	 * @param edPos Byte after last byte of derived buffer (exclusive)
  	 * @return Permanently Read-Only locked FileBuffer which references calling FileBuffer
  	 * @throws IOException If copy creation requires disk access and there is an error with that disk access
  	 * @throws IndexOutOfBoundsException If positions provided are invalid.
  	 */
  	public FileBuffer createReadOnlyCopy(int stPos, int edPos) throws IOException
  	{
  		return this.createReadOnlyCopy((long)stPos, (long)edPos);
  	}
  
  	/**
  	 * Create a Read-Only FileBuffer that references a piece of this buffer.
  	 * <br>WARNING: Because this derived buffer references the current buffer, it will lock 
  	 * this buffer for editing at all positions referenced by a derived buffer! This derived buffer
  	 * will need to be delinked before these regions can be edited again.
  	 * @param stPos First byte of derived buffer (inclusive)
  	 * @param edPos Byte after last byte of derived buffer (exclusive)
  	 * @return Permanently Read-Only locked FileBuffer which references calling FileBuffer
  	 * @throws IOException If copy creation requires disk access and there is an error with that disk access
  	 * @throws IndexOutOfBoundsException If positions provided are invalid.
  	 */
  	public FileBuffer createReadOnlyCopy(long stPos, long edPos) throws IOException
  	{
  		this.checkOffsetPair(stPos, edPos);
  		ROSubFileBuffer rof = new ROSubFileBuffer(this, stPos, edPos);
  		this.children.add(rof);
  		return rof;
  	}
  	
  	/**
  	 * Delink referencing child from this buffer.
  	 * <p> WARNING: Ex-child will still reference parent contents at offsets it was created with, but
  	 * parent will now be free to edit bytes within the space the ex-child used to reference.
  	 * As a result, the "contents" of the ex-child will change along with the contents of the parent,
  	 * as the child getters are simply retrieving offset adjusted information from the parent.
  	 * Be aware that deleting large chunks from the parent can invalidate positions the ex-child
  	 * considered valid. This can lead to a lot of throwing of ArrayIndexOutOfBoundsExceptions.
  	 * <br>Delinking cannot be undone; you will have to create a new sub-buffer.
  	 * @param child Buffer referencing this to delink.
  	 * @return Length two long array. First value (index 0) is the starting offset delinked buffer was
  	 * using. Second value (index 1) was the ending offset. 
  	 * <br> Return null if child was not found linked to this buffer.
  	 */
  	public long[] delinkChild(FileBuffer child)
  	{
  		if (!(child instanceof ROSubFileBuffer)) return null;
  		if (!this.children.contains(child)) return null;
  		ROSubFileBuffer sub = (ROSubFileBuffer)child;
  		long[] offs = new long[2];
  		offs[0] = sub.getStartOffset(this);
  		offs[1] = sub.getEndOffset(this);
  		this.children.remove(child);
  		return offs;
  	}
  	
  	/**
  	 * Delink referencing child from this buffer completely so that neither buffer references the other.
  	 * <br> WARNING: Hard delinking will effectively empty the child buffer! It will have no contents
  	 * to reference. Delinking cannot be undone; you will have to create a new sub-buffer.
  	 * @param child Buffer referencing this to delink.
  	 * @return Length two long array. First value (index 0) is the starting offset delinked buffer was
  	 * using. Second value (index 1) was the ending offset. 
  	 * <br> Return null if child was not found linked to this buffer.
  	 */
  	public long[] hardDelinkChild(FileBuffer child)
  	{
  		long[] offs = delinkChild(child);
  		if (offs == null) return null;
  		ROSubFileBuffer sub = (ROSubFileBuffer)child;
  		sub.delinkParent(this);
  		return offs;
  	}
  
  	protected void addChild(ROSubFileBuffer child)
  	{
  		this.children.add(child);
  	}
  	
  	/**
  	 * Check if a buffer is a direct child of this buffer.
  	 * @param key Potential child buffer of interest.
  	 * @return Whether buffer "key" is a child of this buffer.
  	 */
  	public boolean hasChild(FileBuffer key)
  	{
  		if (this.children == null) return false;
  		if (this.children.size() <= 0) return false;
  		if (key == null) return false;
  		for (ROSubFileBuffer f : this.children)
  		{
  			if (key == f) return true;
  		}
  		return false;
  	}
  	
  	/**
  	 * Recursive function to see if the provided buffer is a child or "descendant" of
  	 * this buffer.
  	 * @param key Potential child buffer of interest.
  	 * @return Whether buffer "key" is a child or descendant of this buffer.
  	 */
  	public boolean hasDescendant(FileBuffer key)
  	{
  		if (this.children == null) return false;
  		if (this.children.size() <= 0) return false;
  		if (key == null) return false;
  		for (ROSubFileBuffer f : this.children)
  		{
  			if (key == f) return true;
  			if (f.hasDescendant(key)) return true;
  		}
  		return false;
  	}
  	
  /* ----- FILE STATISTICS/ INFO ----- */
  
  	/**
  	 * An exception to throw for parsing, conversion, and serialization errors.
  	 * @author Blythe Hospelhorn
  	 * @version 1.0.0
  	 */
  	public static class UnsupportedFileTypeException extends Exception
  	{
  		private static final long serialVersionUID = 2294917631330134570L;
  	}
  
  	/**
  	 * Get an array that counts the frequency of each possible byte value in the buffer.
  	 * Ideal for Huffman encoding.
  	 * @return Integer (int) array of length 256. Each index (00-FF) corresponds to a possible
  	 * byte value. The value is the number of times that byte occurs in the file buffer.
  	 */
  	public int[] getByteFrequencies()
  	{
  		int[] myFreqs = new int[256];
	  
  		for (int f = 0; f < 256; f++) myFreqs[f] = 0;
  		for (long i = 0; i < this.getFileSize(); i++)
  		{
  			int b = Byte.toUnsignedInt(this.getByte(i));
  			myFreqs[b]++;
  		}
  		return myFreqs;
 	}
  
  	public String toString()
  	{
  		String str = "";
	  
  		str += "Standard FileBuffer object -----\n";
  		str += "File Name: " + this.fileName + "\n";
  		str += "Directory: " + this.directory + "\n";
  		str += "Extension: " + this.extension + "\n";
  		str += "File Size: " + this.fSize + "\n";
	  	str += "Base Buffer Capacity: " + this.capacity + "\n";
	  	str += "Byte Order: ";
	  	if (this.isFileFormatBE) str += "Big-Endian \n";
	  	else str += "Little-Endian\n";
	  	str += "Overflowing: " + this.isOverflowing() + "\n";
	  	str += "First 32 Bytes: \n";
	  
	  	for (int i = 0; i < 32; i++)
	  	{
	  		if (i >= this.fSize) break;
	  		str += byteToHexString(this.contents[i]) + " ";
	  		if (i % 16 == 15) str += "\n";
	  	}

	  	return str;
  	}

  	/**
  	 * Get a String describing the type of FileBuffer this instance is.
  	 * Useful for distinguishing instances of various FileBuffer subclasses in debugging.
  	 * @return A short String describing the FileBuffer instance type.
  	 */
  	public String typeString()
  	{
  		return "Standard FileBuffer\n";
  	}
  
  	public boolean equals(Object o)
  	{
  		if (o == this) return true;
  		return false;
  	}
  
  	/* ----- TIMESTAMPING ----- */
  	
  	/**
  	 * Get a timestamp set to an arbitrarily determined "default" point in time (September 8, 1984
  	 * at 6:39:00 AM PDT). The purpose of this timestamp is to indicate an unset value.
  	 * <br>An even more "vanilla" value such as 01/01/00 might be better, but I decided to have
  	 * some fun.
  	 * @return Gregorain Calendar object stamped with an arbitrary time.
  	 */
  	public static GregorianCalendar getVanillaTimestamp()
  	{
  		int year = 1984;
		int month = Calendar.SEPTEMBER;
		int date = 8;
		int hourOfDay = 6;
		int minute = 39;
		int second = 0;
		TimeZone tz = TimeZone.getTimeZone("America/Los_Angeles");
		GregorianCalendar cal = new GregorianCalendar(year, month, date, hourOfDay, minute, second);
		cal.setTimeZone(tz);
  		return cal;
  	}
  	
  	/**
  	 * Get a string representation of a timestamp (as a GregorianCalendar object) formatted
  	 * like a full American style date.
  	 * <br>Month dd, yyyy hh:mm:ss.iii Timezone
  	 * @param timestamp Calendar object containing the time to represent.
  	 * @return Formatted string representing the provided time, if argument is null, then
  	 * this method returns "null".
  	 */
  	public static String formatTimeAmerican(GregorianCalendar timestamp)
  	{
  		if (timestamp == null) return "null";
  		String[] months = new String[]{"January", "February", "March", "April",
  										"May", "June", "July", "August",
  										"September", "October", "November", "December",
  										"MONTH"};
  		int month = timestamp.get(Calendar.MONTH);
  		if (month < 0 || month > 11) month = 12;
  		String s = months[month] + " ";
  		s += timestamp.get(Calendar.DAY_OF_MONTH) + ", ";
  		s += timestamp.get(Calendar.YEAR) + " ";
  		s += String.format("%02d", timestamp.get(Calendar.HOUR_OF_DAY)) + ":";
  		s += String.format("%02d", timestamp.get(Calendar.MINUTE)) + ":";
  		s += String.format("%02d", timestamp.get(Calendar.SECOND)) + ".";
  		s += String.format("%03d", timestamp.get(Calendar.MILLISECOND)) + " ";
  		s += timestamp.getTimeZone().getID();
  		return s;
  	}
  	
}