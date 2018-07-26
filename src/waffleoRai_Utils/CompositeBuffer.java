package waffleoRai_Utils;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;

//Convert all ArrayIndexOutOfBoundsException to IndexOutOfBoundsException
//Multithread some of the buffered writes, including the standard FileBuffer append

/*
 * UPDATES
 * 
 * 2017.09.19 - Recorded as version 1.0.0
 * 
 * 2018.01.22
 * 		1.0.0 -> 1.0.1
 * 			Added the more instance specific memory calculation method (getMinimumMemoryUsage())
 * 			Modified getMemoryBurden to take instance variables into account.
 * 
 * 2018.02.13
 * 		1.0.1 -> 1.0.2
 * 			Fixed a bug in the writer. If the write buffer was larger than the file, then initial
 * 			buffer fill was running off the end of the file and throwing an IndexOutOfBounds
 * 
 * */

/**
 * A FileBuffer subclass that is composed of a set of sequential references to other file buffers.
 * @author Blythe Hospelhorn
 * @version 1.0.2
 * @since February 13, 2018
 */
public class CompositeBuffer extends ROSubFileBuffer
{

	private List<Entry> contents;
	
	private boolean bufferedWrite;
	private int wbuffHalfSize;
	
	private long fileSize;
	
	/* --- OBJECTS --- */

	private class Entry
	{
		private FileBuffer buffer;
		private long startPos;
		private long fileSize;
		
		public Entry(FileBuffer f, long stPos)
		{
			this.buffer = f;
			this.startPos = stPos;
			this.fileSize = buffer.getFileSize();
		}
	
		public boolean isIn(long position)
		{
			if (position < startPos) return false;
			if (position >= startPos + fileSize) return false;
			return true;
		}
		
	}
	
	private class FillThread extends Thread
	{
		private long s;
		private long e;
		private byte[] b;
		
		private int filled;
		
		public FillThread(long stPos, long edPos, byte[] buff, String name)
		{
			super(name);
			s = stPos;
			e = edPos;
			b = buff;
			filled = 0;
		}
		
		public void run()
		{
			filled = fillBuff(s, e, b);
		}
		
		public int getFillSize()
		{
			return filled;
		}
	}
	
	/* --- CONSTRUCTORS --- */
	
	/**
	 * Construct CompositeBuffer with internal ArrayList of provided size.
	 * Recommended for when the number of buffers composing it is known (ie. file serialization)
	 * Super-superclass is constructed empty with only a child list.
	 * @param initBufferNumber Initial size of internal ArrayList.
	 */
	public CompositeBuffer(int initBufferNumber)
	{
		this.contents = new ArrayList<Entry>(initBufferNumber);
		bufferedWrite = true;
		wbuffHalfSize = 0x4000000; //~64 MB
		super.setReadOnly();
	}
	
	/**
	 * Construct CompositeBuffer that utilizes a LinkedList internally.
	 * Better for composites that will add and access buffer pieces sequentially.
	 */
	public CompositeBuffer()
	{
		this.contents = new LinkedList<Entry>();
		bufferedWrite = true;
		wbuffHalfSize = 0x4000000; //~64 MB
		super.setReadOnly();
	}
	
	/* --- INTERNAL --- */
	
	private Entry getEntry(long position)
	{
		for (Entry e : this.contents)
		{
			if (e.isIn(position)) return e;
		}
		return null;
	}
	
	/* --- GETTER OVERRIDE --- */
	
	public void updateFileSize()
	{
		long tot = 0;
		for (Entry e : this.contents) tot += e.buffer.getFileSize();
		this.fileSize = tot;
	}
	
	public long getFileSize()
	{
		//updateFileSize();
		return fileSize;
	}
	
	/**
	 * Get a byte from this file at the provided position.
	 * @param position Offset from FileBuffer start to get byte from.
	 * @return Byte at position specified
	 * @throws IndexOutOfBoundsException If position is invalid.
	 */
	public byte getByte(int position)
	{
		return getByte((long)position);
	}
	
	/**
	 * Get a byte from this file at the provided position.
	 * @param position Offset from FileBuffer start to get byte from.
	 * @return Byte at position specified
	 * @throws IndexOutOfBoundsException If position is invalid.
	 */
	public byte getByte(long position)
	{
		Entry e = getEntry(position);
		if (e == null){
			/*
			System.err.println("CompositeBuffer.getByte(long) || Index out of bounds - no internal buffer found for position: 0x" + Long.toHexString(position));
			System.err.println("CompositeBuffer.getByte(long) || Number of entries: " + contents.size());
			for (int i = 0; i < contents.size(); i++)
			{
				System.err.println("\tBuffer " + i + " --- ");
				System.err.println("\t\tTotal Size: 0x" + Long.toHexString(contents.get(i).buffer.getFileSize()));
				System.err.println("\t\tRecorded Start: 0x" + Long.toHexString(contents.get(i).startPos));
				System.err.println("\t\tRecorded Size: 0x" + Long.toHexString(contents.get(i).fileSize));
			}
			*/
			throw new IndexOutOfBoundsException();
		}
		long os = position - e.startPos;
		return e.buffer.getByte(os);
	}
	
	/**
	 * Get the minimum amount of memory it takes to hold this buffer.
	 * @return Composite memory burden of all referenced buffers.
	 */
	public long getMemoryBurden()
	{
		return getMemoryBurden(null);
	}
  
	/**
	 * Get the minimum amount of memory it takes to hold this buffer.
	 * Buffers referenced in input list will not be counted.
	 * @return Composite memory burden of all referenced buffers, if not referenced in input list.
	 */
	public long getMemoryBurden(Collection<FileBuffer> inMem)
	{
		long mem = 0;
		mem += getMinimumMemoryUsage();
		for (Entry e : this.contents)
		{
			if (inMem == null) inMem = new LinkedList<FileBuffer>();
			if (!inMem.contains(e.buffer))
			{
				mem += e.buffer.getMemoryBurden(inMem);
				inMem.add(e.buffer);
			}
		}
		return mem;
	}
	
	public long getMinimumMemoryUsage()
	{
		long tot = super.getMinimumMemoryUsage();
		int estPtrSz = SystemUtils.approximatePointerSize();
		tot += 1 + 4 + 8;
		if (contents != null) tot += contents.size() * (8 + 8 + estPtrSz);
		return tot;
	}
	
	/**
	 * Determine whether the write mode is set to buffered write.
	 * @return True if write mode is set to buffered, False if set to append.
	 */
	public boolean bufferedWrite()
	{
		return this.bufferedWrite;
	}
	
	/**
	 * Determine whether the write mode is set to appended write.
	 * @return True if write mode is set to append, False if set to buffered.
	 */
	public boolean appendedWrite()
	{
		return !this.bufferedWrite;
	}

	/**
	 * Get the current size of the write buffer to be used in buffered writes, in bytes.
	 * @return The full size of the write buffer.
	 */
	public int writeBufferSize()
	{
		return 2*this.wbuffHalfSize;
	}
	
	/* --- SETTER OVERRIDE --- */
	
	private boolean insertBufferBefore(Entry e, FileBuffer a, long stPos, long edPos)
	{
		/*In retrospect, I don't know that they NEED to be in order!*/
		int lInd = contents.indexOf(e);
		if (lInd < 0) return false;
		FileBuffer add = a;
		try
		{
			if (stPos > 0 || edPos != a.getFileSize()) add = a.createReadOnlyCopy(stPos, edPos);
			Entry n = new Entry(add, e.startPos);
			contents.add(lInd, n);
			n.buffer.addChild(this);
			long sAdd = n.fileSize;
			for (int i = lInd + 1; i < contents.size(); i++) contents.get(i).startPos += sAdd;
		}
		catch (IOException ex)
		{
			ex.printStackTrace();
			return false;
		}
		
		return true;
	}
	
	private Entry splitEntry(Entry e, long localOff)
	{
		/*Returns the part after the split*/
		if (localOff <= 0 || localOff >= e.fileSize) return e;
		try
		{
			FileBuffer f1 = e.buffer.createReadOnlyCopy(0, localOff);
			FileBuffer f2 = e.buffer.createReadOnlyCopy(localOff, e.fileSize);
			Entry e1 = new Entry(f1, e.startPos);
			Entry e2 = new Entry(f2, e.startPos + e1.fileSize);
			int lInd = contents.indexOf(e);
			if (lInd < 0) return null;
			contents.remove(lInd);
			contents.add(lInd, e2);
			contents.add(lInd, e1);
			return e2;
		}
		catch (IOException ex)
		{
			ex.printStackTrace();
			return null;
		}
	}
	
	private boolean insertBufferInside(Entry e, long localOff, FileBuffer a, long stPos, long edPos)
	{
		Entry e2 = splitEntry(e, localOff);
		if (e2 == null) return false;
		if (!insertBufferBefore(e2, a, stPos, edPos)) return false;
		return true;
	}
	
	/**
	 * Insert a new buffer into the composite buffer at a position other than the end.
	 * If inserted into the middle of an existing buffer, will split buffer into two read-only
	 * pieces.
	 * @throws IndexOutOfBoundsException If any provided index is invalid.
	 * @throws NullPointerException If the given FileBuffer to add is null or empty.
	 */
	public void addToFile(FileBuffer addition, int insertPos, int stPos, int edPos)
	{
		/*In theory, you can add anywhere in the comp buffer.
		 * If in the middle of an existing reference, it will just split it into two RO's and
		 * insert this one in the middle.*/
		addToFile(addition, (long)insertPos, (long)stPos, (long)edPos);
	}
  
	/**
	 * Insert a new buffer into the composite buffer at a position other than the end.
	 * If inserted into the middle of an existing buffer, will split buffer into two read-only
	 * pieces.
	 * @throws IndexOutOfBoundsException If any provided index is invalid.
	 * @throws NullPointerException If the given FileBuffer to add is null or empty.
	 */
	public void addToFile(FileBuffer addition, long insertPos, long stPos, long edPos)
	{
		/*Check args*/
		if (addition == null) throw new NullPointerException();
		if (addition.isEmpty()) throw new NullPointerException();
		if (insertPos < 0 || stPos < 0 || edPos < 0) throw new IndexOutOfBoundsException();
		if (insertPos > this.getFileSize()) throw new IndexOutOfBoundsException();
		if (insertPos == this.getFileSize())
		{
			addToFile(addition, stPos, edPos); 
			return;
		}
		if (stPos >= edPos) throw new IndexOutOfBoundsException();
		if (edPos > addition.getFileSize()) edPos = addition.getFileSize();
		boolean superChildren = super.hasChildren();
		if (superChildren) super.checkAllBufferReferences(insertPos);
		
		/*Find breakpoint*/
		Entry t = this.getEntry(insertPos);
		long tOff = insertPos - t.startPos;
		if (tOff == 0)
		{
			/*Insertion is between*/
			if (!this.insertBufferBefore(t, addition, stPos, edPos)) throw new IndexOutOfBoundsException();
			if (superChildren) super.shiftReferencesAfter(insertPos, (edPos - stPos));
		}
		else if (tOff > 0 && tOff < t.fileSize)
		{
			/*Insertion is in the middle of buffer*/
			if (!this.insertBufferInside(t, tOff, addition, stPos, edPos)) throw new IndexOutOfBoundsException();
			if (superChildren) super.shiftReferencesAfter(insertPos, (edPos - stPos));
		}
		else
		{
			/*We have a serious error*/
			throw new IndexOutOfBoundsException();
		}
		
	}
 
	/**
	 * Insert a new buffer into the composite buffer at the end.
	 * @throws NullPointerException If the given FileBuffer to add is null or empty.
	 */
	public void addToFile(FileBuffer addition)
	{
		if (addition == null) throw new NullPointerException();
		if (addition.isEmpty()) throw new NullPointerException();
		//long fsz = addition.getFileSize();
		//System.err.println("CompositeBuffer.addToFile || Called! File size = 0x" + Long.toHexString(fsz));
		Entry e = new Entry(addition, fileSize);
		contents.add(e);
		e.buffer.addChild(this);
		fileSize += e.fileSize;
		//System.err.println("CompositeBuffer.addToFile || New total size = 0x" + Long.toHexString(fileSize));
	}
  
	/**
	 * Insert a new buffer into the composite buffer at the end.
	 * @throws NullPointerException If the given FileBuffer to add is null or empty.
	 * @throws IndexOutOfBoundsException If any provided position is invalid.
	 */
	public void addToFile(FileBuffer addition, int stPos, int edPos)
	{
		addToFile(addition, (long)stPos, (long)edPos);
	}
  
	/**
	 * Insert a new buffer into the composite buffer at the end.
	 * @throws NullPointerException If the given FileBuffer to add is null or empty.
	 * @throws IndexOutOfBoundsException If any provided position is invalid.
	 */
	public void addToFile(FileBuffer addition, long stPos, long edPos)
	{
		/*Check args*/
		if (addition == null) throw new NullPointerException();
		if (addition.isEmpty()) throw new NullPointerException();
		if (stPos >= edPos) throw new IndexOutOfBoundsException();
		long aSz = addition.getFileSize();
		if (edPos > aSz) edPos = aSz;
		
		try
		{
			FileBuffer a = addition;
			if (stPos > 0 || edPos != addition.getFileSize()) a = addition.createReadOnlyCopy(stPos, edPos);
			aSz = edPos - stPos;
			Entry e = new Entry(a, fileSize);
			contents.add(e);
			e.buffer.addChild(this);
			fileSize += aSz;
		}
		catch (IOException ex)
		{
			ex.printStackTrace();
			throw new NullPointerException();
		}
	}
	
	/**
	 * Switch write mode (to disk) to buffered write.
	 * Buffered writing utilizes two byte array buffers to write to disk - filling one
	 * in a second thread while writing the other out.
	 */
	public void setBufferedWrite()
	{
		this.bufferedWrite = true;
	}
	
	/**
	 * Switch write mode (for disk writing) to appending write.
	 * When set to appending write, each component buffer is appended to
	 * the end of the file at the target path using its append function.
	 */
	public void setAppendedWrite()
	{
		this.bufferedWrite = false;
	}
	
	/**
	 * Set total buffer size of the write buffer, if in use.
	 * @param memory Maximum amount, in bytes, of memory that this object can use
	 * as buffer space to write to disk.
	 */
	public void setWriteBufferSize(int memory)
	{
		this.wbuffHalfSize = memory / 2;
	}
	
	/* --- STATUS OVERRIDE --- */
	
  	public boolean isEmpty()
  	{
  		if (this.contents == null) return true;
  		if (this.contents.isEmpty()) return true;
  		return false;
  	}
  
  	public boolean offsetValid(int off)
  	{
  		return offsetValid((long)off);
  	}
  
  	public boolean offsetValid(long off)
  	{
  		if (off < 0) return false;
  		if (off >= this.getFileSize()) return false;
  		return true;
  	}
  
	/* --- WRITER OVERRIDE --- */
	
  	/**
  	 * @throws IndexOutOfBoundsException If any provided position is invalid.
  	 * @throws IllegalThreadStateException If there is a thread interruption when buffer filling thread
  	 * tries to join main thread. This is only possible for buffered writes.
  	 */
  	public void writeFile(String path, long stPos, long edPos) throws IOException
  	{
  		/*Arg Checks*/
  		if (path == null) throw new NullPointerException();
  		if (stPos < 0) throw new IndexOutOfBoundsException();
  		if (stPos >= edPos) throw new IndexOutOfBoundsException();
  		if (edPos > this.getFileSize()) throw new IndexOutOfBoundsException();
  		this.updateFileSize();
	    String dir = FileBuffer.chopPathToDir(path);
	    if(!FileBuffer.directoryExists(dir))
	    {
	    	if(!new File(dir).mkdirs()) throw new IOException();
	    }
  		
  		/*Delegation to other functions*/
  		if (this.bufferedWrite) bufferedWrite(path, stPos, edPos);
  		else appendingWrite(path, stPos, edPos, true);
  	}
  	
  	/**
  	 * @throws IndexOutOfBoundsException If any provided position is invalid.
  	 */
  	public void appendToFile(String path, long stPos, long edPos) throws IOException, NoSuchFileException
  	{
  		this.updateFileSize();
  		appendingWrite(path, stPos, edPos, false);
  	}
  	
  	private int fillBuff(long stPos, long edPos, byte[] buff)
  	{
  		int c = 0;
  		for (int i = 0; i < buff.length; i++)
  		{
  			long p = stPos + (long)i;
  			if (p >= edPos) return c;
  			buff[i] = this.getByte(p);
  			c++;
  		}
  		return c;
  	}
  	
  	private void bufferedWrite(String path, long stPos, long edPos) throws IOException
  	{
  		byte[] buff1 = new byte[wbuffHalfSize];
  		byte[] buff2 = new byte[wbuffHalfSize];
  		long fSize = this.getFileSize();
  		long bSt = 0;
  		long bEd = (long)wbuffHalfSize;
  		boolean done = false;
  		boolean swch = false;
  		
  		//Fill first buffer
  		if (bEd > fSize) bEd = fSize;
  		int fillSz = this.fillBuff(bSt, bEd, buff1);
  		
  		FileOutputStream myStream = new FileOutputStream(path);
  		while(!done)
  		{
  			if (bEd >= fSize) done = true;
  			//Start thread to fill other buffer
  			FillThread fillThread = null;
  			if (swch) fillThread = new FillThread(bSt, bEd, buff1, "CompositeBuffer.BufferedWrite|Buffer1Fill" + Long.toHexString(bSt) + ":" + Long.toHexString(bEd));
  			else fillThread = new FillThread(bSt, bEd, buff2, "CompositeBuffer.BufferedWrite|Buffer2Fill" + Long.toHexString(bSt) + ":" + Long.toHexString(bEd));
  			fillThread.start();
  			//Write current full buffer to disk
  			if (swch) myStream.write(buff2, 0, fillSz);
  			else myStream.write(buff1, 0, fillSz);
  			//Ensure that filler thread has completed its task
  			try 
  			{
				fillThread.join();
			} 
  			catch (InterruptedException e) 
  			{
				e.printStackTrace();
				myStream.close();
				throw new IllegalThreadStateException();
			}
  			fillSz = fillThread.getFillSize();
  			//Calculate next offsets
  			bSt += (long)fillSz;
  			bEd += (long)fillSz;
  			if (bEd > fSize) bEd = fSize;
  			swch = !swch;
  		}
  		
  		myStream.close();
  	}
  	
  	private void appendingWrite(String path, long stPos, long edPos, boolean delExisting) throws IOException
  	{
  		boolean sFound = false;
  		boolean eFound = false;
  		long sOff = 0;
  		long eOff = 0;
  		
  		if (delExisting)
  		{
  			Files.deleteIfExists(Paths.get(path));
  		}
  		
  		for (Entry e : contents)
  		{
  			//If haven't started the writing, look for start position
  			if (eFound) break;
  			if (!sFound)
  			{
  				if (e.isIn(stPos))
  				{
  					sFound = true;
  					sOff = stPos - e.startPos;
  				}
  				else continue;
  			}
  			else sOff = 0;
  			//Look for end position
  			if (e.isIn(edPos)) 
  			{
  				eOff = edPos - e.startPos;
  				eFound = true;
  			}
  			else eOff = e.fileSize;
  			//Run the append function of buffer inside e
  			e.buffer.appendToFile(path, sOff, eOff);
  		}
  	}
	
	/* --- OTHER OVERRIDE --- */
  	
  	/**
  	 * @throws NullPointerException If this or any component buffers are empty.
  	 * @throws IndexOutOfBoundsException If this buffer's size exceeds the maximum integer value.
  	 */
  	public ByteBuffer toByteBuffer()
  	{
  		long fSize = this.getFileSize();
  		if (fSize > 0x7FFFFFFFL) throw new IndexOutOfBoundsException();
  		if (this.contents == null) throw new NullPointerException();
  		if (this.contents.isEmpty()) throw new NullPointerException();
  		ByteBuffer bb = ByteBuffer.allocate((int)fSize);
  		for (Entry e : contents) bb.put(e.buffer.toByteBuffer());
  		return bb;
  	}
  	
  	/**
  	 * @throws NullPointerException If this or any component buffers are empty.
  	 * @throws IndexOutOfBoundsException If the positions given are invalid or the requested.
  	 */
  	public ByteBuffer toByteBuffer(int stPos, int edPos)
  	{
  		return toByteBuffer((long)stPos, (long)edPos);
  	}
  
  	/**
  	 * @throws NullPointerException If this or any component buffers are empty.
  	 * @throws IndexOutOfBoundsException If the positions given are invalid or the requested
  	 * region exceeds the maximum integer value in length.
  	 */
  	public ByteBuffer toByteBuffer(long stPos, long edPos)
  	{
  		if (stPos < 0) throw new IndexOutOfBoundsException();
  		if (stPos >= edPos) throw new IndexOutOfBoundsException();	
  		this.updateFileSize();
  		long fSize = this.getFileSize();
  		if (edPos > fSize) throw new IndexOutOfBoundsException();
  		
  		if (stPos == 0 && edPos == fSize) return this.toByteBuffer();
  		long sz = edPos - stPos;
  		if (sz > 0x7FFFFFFFL)throw new IndexOutOfBoundsException();
  		
  		ByteBuffer bb = ByteBuffer.allocate((int)sz);
  		
  		boolean sFound = false;
  		boolean eFound = false;
  		long sOff = 0;
  		long eOff = 0;
  		
  		for (Entry e : contents)
  		{
  			//If haven't started the writing, look for start position
  			if (eFound) break;
  			if (!sFound)
  			{
  				if (e.isIn(stPos))
  				{
  					sFound = true;
  					sOff = stPos - e.startPos;
  				}
  				else continue;
  			}
  			else sOff = 0;
  			//Look for end position
  			if (e.isIn(edPos)) 
  			{
  				eOff = edPos - e.startPos;
  				eFound = true;
  			}
  			else eOff = e.fileSize;
  			//Run the append function of buffer inside e
  			bb.put(e.buffer.toByteBuffer(sOff, eOff));
  		}
  		
  		return bb;
  	}
 
  	private void addEntry(Entry e)
  	{
  		Entry n = new Entry(e.buffer, e.startPos);
  		contents.add(n);
  	}
  	
  	/**
  	 * @throws IndexOutOfBoundsException If one or more the positions provided is invalid.
  	 * @throws NullPointerException If this buffer's content list is null.
  	 */
 	public FileBuffer createCopy(int stPos, int edPos) throws IOException
  	{
 		return createCopy((long)stPos, (long)edPos);
  	}
  
  	/**
  	 * @throws IndexOutOfBoundsException If one or more the positions provided is invalid.
  	 * @throws NullPointerException If this buffer's content list is null.
  	 */
  	public FileBuffer createCopy(long stPos, long edPos) throws IOException
  	{
  		if (this.contents == null) throw new NullPointerException();
  		if (this.contents.isEmpty()) return null;
  		this.updateFileSize();
  		CompositeBuffer myCopy = new CompositeBuffer(contents.size());
  		for (Entry e : contents) myCopy.addEntry(e);
  		return myCopy;
  	}
  	
  	public String typeString()
  	{
  		return "Composite FileBuffer";
  	}
	
	  /* --- Parent Access --- */
	  
  	/**
  	 * @throws IllegalArgumentException If key provided does not correspond to a buffer that is
  	 * a parent or ancestor of this buffer (ie. component).
  	 */
  	protected long getStartOffset(FileBuffer key)
	{
  		/*Right now, it assumes that any DIRECT references are full FileBuffers - 
  		 * any partials are already readonly buffers which should be direct children
  		 * already of the full FB's they reference and position shifts should be handled by them.*/
  		if (!key.hasDescendant(this)) throw new IllegalArgumentException();
		return 0;
	}
	  
  	/**
  	 * @throws IllegalArgumentException If key provided does not correspond to a buffer that is
  	 * a parent or ancestor of this buffer (ie. component).
  	 */
  	protected long getEndOffset(FileBuffer key)
	{
  		if (!key.hasDescendant(this)) throw new IllegalArgumentException();
		return key.getFileSize();
	}
	  
	protected void setStartOffset(FileBuffer key, long newStart)
	{
		if (key.hasChild(this))
		{
			/*Referenced buffer key has just added something outside region referenced here
			 * (Region inside this composite is locked read-only)
			 * */
			try 
			{
				for (Entry e : contents)
				{
					if (e.buffer == key)
					{
						FileBuffer sub = key.createReadOnlyCopy(newStart, key.getFileSize());
						key.delinkChild(this);	
						e.buffer = sub;
						break;
					}
				}
			} 
			catch (IOException ex) 
			{
				ex.printStackTrace();
				throw new UnsupportedOperationException();
			}
		}
	}
	  
	protected void setEndOffset(FileBuffer key, long newEnd)
	{
		if (key.hasChild(this))
		{
			try 
			{
				for (Entry e : contents)
				{
					if (e.buffer == key)
					{
						FileBuffer sub = key.createReadOnlyCopy(0, newEnd);
						key.delinkChild(this);	
						e.buffer = sub;
						break;
					}
				}
			} 
			catch (IOException ex) 
			{
				ex.printStackTrace();
				throw new UnsupportedOperationException();
			}
		}
	}
	  
	protected void delinkParent(FileBuffer key)
	{
		long upshift = 0;
		List<Entry> rlist = new LinkedList<Entry>();
		if (key.hasChild(this))
		{
			for (Entry e : contents)
			{
				e.startPos -= upshift;
				if (e.buffer == key) 
				{
					upshift += e.fileSize;
					rlist.add(e);
				}
			}
			for (Entry e : rlist) contents.remove(e);
		}
	}
	  
  	
}
