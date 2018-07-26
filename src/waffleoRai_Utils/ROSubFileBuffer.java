package waffleoRai_Utils;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.file.NoSuchFileException;
import java.util.Collection;

//OVERRIDE isEmpty

/*UPDATES
 * 2017.08.31
 * 	2.0 -> 2.0.1 | Added getMemoryBurden()
 * 2017.09.01
 * 	2.0.1 -> 3.0 | Altered implementation to reference parent instead of use a super class.
 * 2017.09.08
 * 	3.0 -> 3.0.1 | Updated to reflect changes in parent class 
 * 		(deletion of certain methods, internal method calling, reference implementation...)
 * 2017.09.20
 * 	3.0.1 -> 3.0.2 | Javadoc updates, change to IndexOutOfBoundsException thrown
 * 2018.01.22
 * 	3.0.2 -> 3.0.3 | Added more technical memory calculation method (getMinimumMemoryUsage())
 * */

/**
 * For creating a new file buffer from a piece of an existing one, so
 * that copying to a new array is not required.
 * Because this is basically just a wrapper for a reference to an existing buffer,
 * all write functions will throw an exception when called.
 * @author Blythe Hospelhorn
 * @version 3.0.3
 * @since January 22, 2018
 */
public class ROSubFileBuffer extends FileBuffer{

	private long posZero;
	private long posEnd;
	
	private FileBuffer parent;
	
	protected ROSubFileBuffer(FileBuffer file, int stPos, int edPos)
	{
		super();
		this.parent = file;
		if (edPos > (int)parent.getFileSize()) edPos = (int)parent.getFileSize();
		if (stPos < 0) stPos = 0;
		this.posZero = (long)stPos;
		this.posEnd = (long)edPos;
		super.setReadOnly();
	}
	
	protected ROSubFileBuffer(FileBuffer file, long stPos, long edPos)
	{
		super();
		this.parent = file;
		if (edPos > parent.getFileSize()) edPos = parent.getFileSize();
		if (stPos < 0) stPos = 0;
		this.posZero = stPos;
		this.posEnd = edPos;
		super.setReadOnly();
	}
	
	protected ROSubFileBuffer()
	{
		super();
		this.parent = null;
		this.posZero = -1;
		this.posEnd = -1;
	}
		
	private long getPosition(long locPos)
	{
		if (!this.offsetValid(locPos)) throw new ArrayIndexOutOfBoundsException();
		return locPos + posZero;
	}
	
	private long[] checkPositionPair(long stPos, long edPos)
	{
		long[] offs = new long[2];
		if (stPos >= edPos) throw new ArrayIndexOutOfBoundsException();
		if (stPos < 0) stPos = 0;
		if (edPos - stPos > this.getFileSize()) edPos = this.getFileSize() - stPos;
		stPos = getPosition(stPos);
		edPos = getPosition(edPos);
		offs[0] = stPos;
		offs[1] = edPos;
		return offs;
	}
	
	/* --- Getter/ Setter Override --- */
	
 	  public long getFileSize()
	  {
		  return posEnd - posZero;
	  }
	  
 	  /**
 	   * @throws IndexOutOfBoundsException If position is invalid.
 	   */
	  public byte getByte(int position)
	  {
		  if (!this.offsetValid(position)) throw new IndexOutOfBoundsException();
		  return parent.getByte((int)posZero + position);
	  }
	  
		/**
		 * Get a byte from this file at the provided position.
		 * @param position Offset from FileBuffer start to get byte from.
		 * @return Byte at position specified
		 * @throws IndexOutOfBoundsException If position is invalid.
		 */
	  public byte getByte(long position)
	  {
		  if (!this.offsetValid(position)) throw new IndexOutOfBoundsException();
		  return parent.getByte(posZero + position);
	  }
	  
	/* --- Reading Override --- */
	  	/*--- Anything that relies on getByte apparently doesn't need to be overridden!*/
	  
	  /**
	   * @throws NullPointerException If charset string is null.
	   * @throws IndexOutOfBoundsException If any positions are invalid.
	   * @throws IllegalCharsetNameException If the given charset name is illegal.
	   * @throws UnsupporedCharsetException If no support for the named charset is available in this instance of the Java virtual machine.
	   */
	  public String readEncoded_string(String charset, long stPos, long edPos)
	  {
		  if (charset == null) throw new NullPointerException();
		  long[] myPos = checkPositionPair(stPos, edPos);
		  if (myPos == null || myPos.length != 2) throw new IndexOutOfBoundsException();
		  return parent.readEncoded_string(charset, myPos[0], myPos[1]);
	  }
	   
	/* --- Other Override --- */
	  
	  public boolean isEmpty()
	  {
		  if (parent == null) return true;
		  return false;
	  }
	  
	  public boolean readOnly()
	  {
		  return true;
	  }
	  
	  /**
	   * This sub-type of FileBuffer can only be read-only.
	   * An exception will be thrown if this function is called.
	   * @throws UnsupportedOperationException If this function is called.
	   */
	  public void unsetReadOnly()
	  {
		  throw new UnsupportedOperationException();
	  }
	  
	  public boolean offsetValid(int off)
	  {
		  return offsetValid((long)off);
	  }
	  
	  public boolean offsetValid(long off)
	  {
		  if (parent == null) return false;
		  if (this.getFileSize() <= 0) return false;
		  if (off < 0) return false;
		  if (posZero >= parent.getFileSize()) return false;
		  if (posZero + off >= posEnd) return false;
		  if (posZero + off >= parent.getFileSize()) return false;
		  return true;
	  }
	  
	  public void writeFile() throws IOException
	  {
		  this.writeFile(this.getPath(), 0, this.getFileSize());
	  }
	  
	  public void writeFile(String path) throws IOException
	  { 
		  this.writeFile(path, 0, this.getFileSize());
	  }
	  
	  public void writeFile(int stPos) throws IOException
	  {
		  this.writeFile((long)stPos);
	  }
	  
	  public void writeFile(long stPos) throws IOException
	  {
		  this.writeFile(stPos, this.getFileSize());
	  }
	  
	  public void writeFile(int stPos, int edPos) throws IOException
	  {
		  this.writeFile((long)stPos, (long)edPos);
	  }
	  
	  public void writeFile(long stPos, long edPos) throws IOException
	  {
		  writeFile(this.getPath(), stPos, edPos);
	  }
	  
	  public void writeFile(String path, int stPos) throws IOException
	  {
		  this.writeFile(path, (long)stPos);
	  }
	  
	  public void writeFile(String path, int stPos, int edPos) throws IOException
	  {
		  this.writeFile(path, (long)stPos, (long)edPos);
	  }
	  
	  /**
	   * @throws IndexOutOfBoundsException If any position is invalid.
	   */
	  public void writeFile(String path, long stPos, long edPos) throws IOException
	  {
		  long[] myPos = checkPositionPair(stPos, edPos);
		  if (myPos == null || myPos.length != 2) throw new IndexOutOfBoundsException();
		  parent.writeFile(path, myPos[0], myPos[1]);
	  }

	  /**
	   * @throws IndexOutOfBoundsException If any position is invalid.
	   */
	  public void appendToFile(String path, long stPos, long edPos) throws IOException, NoSuchFileException
	  {
		  long[] myPos = checkPositionPair(stPos, edPos);
		  if (myPos == null || myPos.length != 2) throw new IndexOutOfBoundsException();
		  parent.appendToFile(path, myPos[0], myPos[1]);
	  }
	  
	  /**
	   * @throws IndexOutOfBoundsException If any position is invalid, or requested buffer size
	   * exceeds the maximum integer value.
	   * @throws NullPointerException If parent buffer has no contents.
	   */
	  public ByteBuffer toByteBuffer()
	  {
		  if (this.isEmpty()) throw new ArrayIndexOutOfBoundsException();
		  return parent.toByteBuffer(posZero, posEnd);
	  }
	  
	  /**
	   * @throws IndexOutOfBoundsException If any position is invalid, or requested buffer size
	   * exceeds the maximum integer value.
	   * @throws NullPointerException If parent buffer has no contents.
	   */
	  public ByteBuffer toByteBuffer(int stPos, int edPos)
	  {
		  return this.toByteBuffer((long)stPos, (long)edPos);
	  }
	  
	  /**
	   * @throws IndexOutOfBoundsException If any position is invalid, or requested buffer size
	   * exceeds the maximum integer value.
	   * @throws NullPointerException If parent buffer has no contents.
	   */
	  public ByteBuffer toByteBuffer(long stPos, long edPos)
	  {
		  long[] myPos = checkPositionPair(stPos, edPos);
		  if (myPos == null || myPos.length != 2) throw new IndexOutOfBoundsException();
		  return parent.toByteBuffer(myPos[0], myPos[1]);
	  }
	  
	  /**
	   * @throws IndexOutOfBoundsException If any position is invalid.
	   */
	  public FileBuffer createCopy(int stPos, int edPos) throws IOException
	  {
		  return this.createCopy((long)stPos, (long)edPos);
	  }
	  
	  /**
	   * @throws IndexOutOfBoundsException If any position is invalid.
	   */
	  public FileBuffer createCopy(long stPos, long edPos) throws IOException
	  {
		  long[] myPos = checkPositionPair(stPos, edPos);
		  if (myPos == null || myPos.length != 2) throw new IndexOutOfBoundsException();
		  return parent.createCopy(myPos[0], myPos[1]);
	  }
	  
	  /**
	   * @throws IndexOutOfBoundsException If any position is invalid.
	   */
	  public FileBuffer createReadOnlyCopy(int stPos, int edPos) throws IOException
	  {
		  return this.createReadOnlyCopy((long)stPos, (long)edPos);
	  }
	  
	  /**
	   * @throws IndexOutOfBoundsException If any position is invalid.
	   */
	  public FileBuffer createReadOnlyCopy(long stPos, long edPos) throws IOException
	  {
		  long[] myPos = checkPositionPair(stPos, edPos);
		  if (myPos == null || myPos.length != 2) throw new IndexOutOfBoundsException();
		  return parent.createReadOnlyCopy(myPos[0], myPos[1]);
	  }
	  
	  /**
	   * Get the minimum amount of memory taken up by buffer data.
	   * @return 0 by default as it itself takes up no memory in file data (only object info)
	   */
	  public long getMemoryBurden()
	  {
		  return getMinimumMemoryUsage();
	  }
	  
	  /**
	   * Overload for more accurate calculation.
	   * It can scan a collection of already counted FileBuffers for its parent.
	   * If its parent has already been counted, then it will return 0.
	   * If not, it will return the amount of memory its parent takes up in data.
	   * @param inMem : Collection (eg. List) of FileBuffer references already counted
	   * @return 0 if parent is in list, otherwise allocated size of parent it references.
	   */
	  public long getMemoryBurden(Collection<FileBuffer> inMem)
	  {
		  for (FileBuffer f : inMem)
		  {
			  if (f == this.parent) return getMinimumMemoryUsage();
		  }
		  inMem.add(parent);
		  return parent.getMemoryBurden(inMem) + getMinimumMemoryUsage();
	  }
	  
	  public long getMinimumMemoryUsage()
		{
			long tot = super.getMinimumMemoryUsage();
			int estPtrSz = SystemUtils.approximatePointerSize();
			tot += 8 + 8 + estPtrSz;
			return tot;
		}
	  
	  public String toString()
	  {
		  String str = "";
		  
		  str += "READ ONLY REFERENCE COPY!!\n";
		  str += "File Name: " + this.getName() + "\n";
		  str += "Directory: " + this.getDir() + "\n";
		  str += "Extension: " + this.getExt() + "\n";
		  str += "Parent File Size: " + super.getFileSize() + "\n";
		  str += "Child File Size: " + this.getFileSize() + "\n";
		  str += "Parent Base Buffer Capacity: " + super.getBaseCapacity() + "\n";
		  str += "Byte Order: ";
		  if (this.isBigEndian()) str += "Big-Endian \n";
		  else str += "Little-Endian\n";
		  str += "Parent Overflowing: " + super.isOverflowing() + "\n";
		  str += "Child Index Range: " + this.posZero + " - " + this.posEnd + "\n";
		  str += "First 32 Bytes: \n";
		  
		  for (int i = 0; i < 32; i++)
		  {
			  System.out.println("i = " + i);
			  if (i >= (int)this.getFileSize()) break;
			  str += byteToHexString(this.getByte(i)) + " ";
			  if (i % 16 == 15) str += "\n";
		  }
		  
		  return str;
	  }
	  
	  public String typeString()
	  {
		  return "Standard FileBuffer RO Reference\n";
	  }
	  
	  /* --- Parent Access --- */
	  
	  /**
	   * Returns the start offset relative to the parent.
	   * @param key Parent to retrieve offset relative to.
	   * @return Long integer representing the offset from the parent start that this buffer uses
	   * as offset 0.
	   */
	  protected long getStartOffset(FileBuffer key)
	  {
		  return this.posZero;
	  }
	  
	  /**
	   * Returns the end offset relative to the parent.
	   * @param key Parent to retrieve offset relative to.
	   * @return Long integer representing the first byte position in the parent after
	   * the end of the section this child references.
	   */
	  protected long getEndOffset(FileBuffer key)
	  {
		  return this.posEnd;
	  }
	  
	  protected void setStartOffset(FileBuffer key, long newStart)
	  {
		  if (key == this.parent)
		  {
			  if (newStart < 0) newStart = 0;
			  if (newStart > parent.getFileSize()) newStart = parent.getFileSize();
			  this.posZero = newStart;
		  }
	  }
	  
	  protected void setEndOffset(FileBuffer key, long newEnd)
	  {
		  if (key == this.parent)
		  {
			  if (newEnd < this.posZero) newEnd = posZero;
			  if (newEnd > parent.getFileSize()) newEnd = parent.getFileSize();
			  this.posEnd = newEnd;
		  }
	  }
	  
	  protected void delinkParent(FileBuffer key)
	  {
		  if (key != parent) return;
		  parent = null;
		  this.posZero = -1;
		  this.posEnd = -1;
	  }
	  
}
