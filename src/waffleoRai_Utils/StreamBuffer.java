package waffleoRai_Utils;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Paths;
import java.nio.file.StandardCopyOption;
import java.nio.file.StandardOpenOption;
import java.util.ArrayDeque;
import java.util.Collection;
import java.util.Deque;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Random;



/*UPDATES
 * 2017.08.31
 * 	1.6 -> 1.6.1 | Added getMemoryBurden() function.
 * 2017.09.08
 * 	1.6.1 -> 1.7.0 | Updated for compatibility with parent updates, eliminate need for Stream 
 * 		specific sub buffer class.
 * 		Automatic disk dump when memory burden of insertions passes twice the memory burden of the stream buffer.
 * 2017.09.20
 * 	1.7.0 -> 1.7.1 | Throw IndexOutOfBounds instead of ArrayIndexOutOfBounds
 * 		Minor javadoc updates
 * 		Multithread the write function buffering
 * 2017.09.24
 * 	1.7.1 -> 1.7.2 | Update separator character detection.
 * 2018.01.22
 * 	1.7.2 -> 1.7.3 | Added memory calculation function more focused on instance size
 * 2018.02.03
 * 	1.7.3 -> 1.7.4 | Update for byte order switch bug
 * */

/**
 * For files larger than should be held in memory.
 * @author Blythe Hospelhorn
 * @version 1.7.4
 * @since February 3, 2018
 */
public class StreamBuffer extends FileBuffer
{
	
	//public static final int DEFO_SUBBUF_SIZE = 0x800; //2048
	//public static final int DEFO_SUBBUF_NUM = 0x8000; //32768
	//Default to 67,108,864 bytes (64 MB) loaded into memory at any given time.
	
	public static final int DEFO_SUBBUF_SIZE = 0x1000; //4096
	public static final int DEFO_SUBBUF_NUM = 0x10000; //65536
	//Default to 268,435,456 bytes (256 MB) loaded into memory at any given time.

	private int subBufferSize;
	private int subBufferNum;
	private int writeBuffHalfSize;
	
	private long streamStart;
	private long streamLength;
	
	private MiniBuffer[] contents;
	private int[] locTable; //Maps file sector (index) to buffer sector (value)
	private Deque<Integer> usageQueue; //Value is buffer sector. Most recently used at head of deque.
	
	private String filePath;
	private String tempPath;
	private boolean tempMode;
	
	//private boolean readOnly;
	private Map<Long, Change> modifications; //References RELATIVE offsets
	private int modMemThreshhold;
	
	/* ----- OBJECTS ----- */
	
	private static enum ChangeType
	{
		INSERT,
		DELETION,
		SUBSTITUTION,
		INSSUB;
	}
	
	private static interface Change
	{
		public ChangeType getType();
		public long getPosition();
		public long getLength();
		public byte getByte(long seqPos);
		public long getSize();
		public void replaceSeqByte(byte b, long seqPos);
		public void addSeqByte(byte b);
		public void addSeqFile(FileBuffer f);
		public void insertSeqByte(byte b, long seqPos);
		public void insertSeqFile(long seqPos, FileBuffer f);
		public void deleteBytes(long stPos, long len);
		public long getMemoryTax();
		public long getMinimumMemoryUsage();
	}
	
	private static class Insertion implements Change
	{
		private long fileOffset;
		//private byte[] seq;
		private Deque<FileBuffer> seq;
		
		public Insertion(long fOff, FileBuffer initSeq)
		{
			fileOffset = fOff;
			seq = new LinkedList<FileBuffer>();
			seq.add(initSeq);
		}
		
		public void replaceSeqByte(byte b, long seqPos)
		{
			if (seqPos < 0 || seqPos >= this.getLength()) return;
			long cPos = 0;
			for (FileBuffer s : this.seq)
			{
				if ((cPos + (int)s.getFileSize()) > seqPos)
				{
					//It's in this buffer
					long locOff = seqPos - cPos;
					s.replaceByte(b, locOff);
					return;
				}
				cPos += s.getFileSize();
			}
		}
		
		public void addSeqByte(byte b)
		{
			FileBuffer e = this.seq.getLast();
			if (!e.readOnly() && e.getFileSize() < e.getBaseCapacity()) e.addToFile(b);
			else
			{
				FileBuffer n = new FileBuffer(StreamBuffer.DEFO_SUBBUF_SIZE);
				n.addToFile(b);
				seq.addLast(n);
			}
		}
		
		public void insertSeqByte(byte b, long seqPos)
		{
			long l = this.getLength();
			if (seqPos < 0 || seqPos > l) return;
			if (seqPos == l)
			{
				this.addSeqByte(b);
				return;
			}
			long cPos = 0;
			boolean updated = false;
			FileBuffer[] seqArr = seq.toArray(new FileBuffer[seq.size() + 5]);
			for (int i = 0; i < seq.size(); i++)
			{
				FileBuffer s = seqArr[i];
				if ((cPos + s.getFileSize()) > seqPos)
				{
					//Position is in this buffer
					long locPos = seqPos - cPos;
					if(!s.readOnly())
					{
						//Write enabled - just try to insert it into the buffer directly
						s.addToFile(b, locPos); //Only works for full standard FB if overflow works.
						return;
					}
					else
					{
						//Read only.
						//Split the buffer into two and create a new one in between
						try 
						{
							FileBuffer FB1 = s.createReadOnlyCopy(0, locPos);
							FileBuffer FB2 = s.createReadOnlyCopy(locPos, s.getFileSize());
							FileBuffer FBn = new FileBuffer(DEFO_SUBBUF_SIZE);
							FBn.addToFile(b);
							for (int j = seq.size() - 1; j > i; j++)
							{
								seqArr[j+2] = seqArr[j];
							}
							seqArr[i] = FB1;
							seqArr[i+1] = FBn;
							seqArr[i+2] = FB2;
							//Update list
							updated = true;
						} 
						catch (IOException e) 
						{
							//Can't insert.
							e.printStackTrace();
							return;
						}
					}
				}
				else if (cPos + s.getFileSize() == seqPos)
				{
					//Position is between this buffer and the next.
					FileBuffer FBn = new FileBuffer(DEFO_SUBBUF_SIZE);
					for (int j = seq.size() - 1; j > i; j++)
					{
						seqArr[j+1] = seqArr[j];
					}
					seqArr[i+1] = FBn;
					updated = true;
				}
				if (updated)
				{
					seq.clear();
					for (int j = 0; j < seqArr.length; j++)
					{
						if (seqArr[j]!=null) seq.add(seqArr[j]);
					}
					return;
				}
				cPos += s.getFileSize();	
			}
		}
		
		public void addSeqFile(FileBuffer f)
		{
			if (f != null) seq.addLast(f);
		}
		
		public void insertSeqFile(long seqPos, FileBuffer f)
		{
			long l = this.getLength();
			if (seqPos < 0 || seqPos > l) return;
			if (seqPos == l)
			{
				this.addSeqFile(f);
				return;
			}
			long cPos = 0;
			boolean updated = false;
			int seqLen = seq.size();
			FileBuffer[] seqArr = seq.toArray(new FileBuffer[seq.size() + 5]);
			for (int i = 0; i < seqLen; i++)
			{
				FileBuffer s = seqArr[i];
				if ((cPos + s.getFileSize()) > seqPos)
				{
					//Position is in this buffer
					long locPos = seqPos - cPos;
					//Split the buffer into two and create a new one in between
					try 
					{
						FileBuffer FB1 = s.createReadOnlyCopy(0, locPos);
						FileBuffer FB2 = s.createReadOnlyCopy(locPos, s.getFileSize());
						for (int j = seq.size() - 1; j > i; j++)
						{
							seqArr[j+2] = seqArr[j];
						}
						seqArr[i] = FB1;
						seqArr[i+1] = f;
						seqArr[i+2] = FB2;
						//Update list
						updated = true;
					} 
					catch (IOException e) 
					{
						//Can't insert.
						e.printStackTrace();
						return;
					}
				}
				else if (cPos + s.getFileSize() == seqPos)
				{
					//Position is between this buffer and the next.
					for (int j = seq.size() - 1; j > i; j++)
					{
						seqArr[j+1] = seqArr[j];
					}
					seqArr[i+1] = f;
					updated = true;
				}
				if (updated)
				{
					seq.clear();
					for (int j = 0; j < seqArr.length; j++)
					{
						if (seqArr[j]!=null) seq.add(seqArr[j]);
					}
					return;
				}
				cPos += s.getFileSize();	
			}
		}
		
		public void deleteBytes(long stPos, long len)
		{
			if (stPos < 0 || len <= 0) return;
			long l = this.getLength();
			if (stPos >= l) return;
			if (stPos + len > l) return;
			long cPos = 0;
			boolean pushing = false;
			boolean update = false;
			long left = len;
			int seqSz = seq.size();
			FileBuffer[] seqArr = seq.toArray(new FileBuffer[seqSz]);
			for (int i = 0; i < seqSz; i++)
			{
				FileBuffer s = seqArr[i];
				if (!pushing)
				{
					if (cPos + s.getFileSize() > stPos)
					{
						//Position is in this buffer
						long locsp = stPos - cPos;
						if (locsp + len <= s.getFileSize())
						{
							if (locsp == 0 && locsp + len == s.getFileSize())
							{
								seq.remove(s);
								return;
							}
							else
							{
								s.deleteFromFile(locsp, locsp + len);
								return;
							}
						}
						else
						{
							left -= s.getFileSize() - locsp;
							if (locsp == 0)
							{
								seqArr[i] = null;
								update = true;
							}
							else
							{
								s.deleteFromFile(locsp);
							}
							pushing = true;
						}
					}
					cPos += s.getFileSize();
				}
				else
				{
					//Pushing
					if (left > 0)
					{
						if (left >= s.getFileSize())
						{
							left -= s.getFileSize();
							seqArr[i] = null;
							update = true;
						}
						else
						{
							s.deleteFromFile(0, left);
							left = 0;
							break;
						}
					}
					else break;
				}
			}
			if (update)
			{
				seq.clear();
				for (FileBuffer f : seqArr)
				{
					if (f != null) seq.add(f);
				}
			}
		}
		
		public ChangeType getType()
		{
			return ChangeType.INSERT;
		}
		
		public long getPosition()
		{
			return this.fileOffset;
		}
		
		public long getLength()
		{
			long t = 0;
			for (FileBuffer s : this.seq)
			{
				t += s.getFileSize();
			}
			return t;
		}
		
		public byte getByte(long seqPos)
		{
			if (seqPos < 0 || seqPos >= this.getLength()) throw new ArrayIndexOutOfBoundsException();
			long cPos = 0;
			for (FileBuffer s : this.seq)
			{
				if ((cPos + s.getFileSize()) > seqPos)
				{
					//It's in this buffer
					long locOff = seqPos - cPos;
					return s.getByte(locOff);
				}
				cPos += s.getFileSize();
			}
			return -1;
		}
		
		public long getSize()
		{
			return this.getLength();
		}
	
		public long getMemoryTax()
		{
			int estPtrSz = SystemUtils.approximatePointerSize();
			Collection<FileBuffer> counted = new LinkedList<FileBuffer>();
			long mem = 8 + estPtrSz;
			for (FileBuffer f : seq)
			{
				mem += f.getMemoryBurden(counted) + estPtrSz;
				counted.add(f);
			}
			return mem;
		}
		
		public long getMinimumMemoryUsage()
		{
			int estPtrSz = SystemUtils.approximatePointerSize();
			Collection<FileBuffer> counted = new LinkedList<FileBuffer>();
			long mem = 8 + estPtrSz;
			for (FileBuffer f : seq)
			{
				if (counted.contains(f)) mem += estPtrSz;
				else
				{
					mem += f.getMinimumMemoryUsage() + estPtrSz;
					counted.add(f);	
				}
			}
			return mem;
		}
	
	}
	
	private static class Deletion implements Change
	{
		private long fileOffset;
		private long len;
		
		public Deletion(long fOff, long len)
		{
			fileOffset = fOff;
			this.len = len;
		}
		
		public ChangeType getType()
		{
			return ChangeType.DELETION;
		}
		
		public long getPosition()
		{
			return this.fileOffset;
		}
		
		public long getLength()
		{
			return len;
		}
		
		public byte getByte(long seqPos)
		{
			return -1;
		}
	
		public long getSize()
		{
			return this.getLength() * -1;
		}

		public void replaceSeqByte(byte b, long seqPos) 
		{
			//Do nothing.
		}

		public void addSeqByte(byte b) 
		{

		}

		public void insertSeqByte(byte b, long seqPos) 
		{
			//Do nothing
		}

		public void deleteBytes(long stPos, long len) 
		{
			//Does nothing
		}

		public void addSeqFile(FileBuffer f) 
		{
			//Does nothing
		}

		public void insertSeqFile(long seqPos, FileBuffer f) 
		{
			//Does nothing
		}

		public long getMemoryTax()
		{
			return 0;
		}
		
		public long getMinimumMemoryUsage()
		{
			return 8 + 8;
		}
	
	}
	
	private static class Substitution implements Change
	{
		private long fileOffset;
		private byte subByte;
		
		public Substitution(long fOff, byte sub)
		{
			this.subByte = sub;
		}
		
		public ChangeType getType()
		{
			return ChangeType.SUBSTITUTION;
		}
		
		public long getPosition()
		{
			return this.fileOffset;
		}
		
		public long getLength()
		{
			return 0;
		}
		
		public byte getByte(long seqPos)
		{
			return subByte;
		}
	
		public long getSize()
		{
			return this.getLength();
		}

		public void replaceSeqByte(byte b, long seqPos) 
		{
			this.subByte = b;
		}

		public void addSeqByte(byte b) 
		{
			this.subByte = b;
		}

		public void insertSeqByte(byte b, long seqPos)
		{
			this.subByte = b;
		}

		public void deleteBytes(long stPos, long len) 
		{
			//Do nothing
		}

		public void addSeqFile(FileBuffer f) 
		{
			//Do nothing
		}

		public void insertSeqFile(long seqPos, FileBuffer f) 
		{
			//Do nothing
		}
	
		public long getMemoryTax()
		{
			return 1;
		}
		
		public long getMinimumMemoryUsage()
		{
			return 8 + 1;
		}
	}
	
	private static class Inssub implements Change
	{	
		private long fileOffset;
		private Substitution sub;
		private Insertion ins;
		
		public Inssub(Substitution s, Insertion i)
		{
			this.fileOffset = i.getPosition();
			this.sub = s;
			this.ins = i;
		}

		public ChangeType getType()
		{
			return ChangeType.INSSUB;
		}
		
		public long getPosition()
		{
			return this.fileOffset;
		}
		
		public long getLength()
		{
			return ins.getLength();
		}
		
		public byte getByte(long seqPos)
		{
			if (seqPos < 0 || seqPos > ins.getLength()) throw new ArrayIndexOutOfBoundsException();
			if (seqPos == ins.getLength()) return sub.getByte(0);
			return ins.getByte(seqPos);
		}
		
		public long getSize()
		{
			return ins.getSize();
		}

		public void replaceSeqByte(byte b, long seqPos) 
		{
			if (seqPos < 0 || seqPos > this.getLength()) return;
			if (seqPos == this.getLength()) sub.replaceSeqByte(b, 0);
			else ins.replaceSeqByte(b, seqPos);
		}

		public void addSeqByte(byte b) 
		{
			ins.addSeqByte(b);
		}

		public void insertSeqByte(byte b, long seqPos) 
		{
			if (seqPos < 0 || seqPos > this.getLength()) return;
			if (seqPos == this.getLength()) sub.insertSeqByte(b, 0);
			else ins.insertSeqByte(b, seqPos);
		}

		public void deleteBytes(long stPos, long len) 
		{
			ins.deleteBytes(stPos, len);
		}

		public void addSeqFile(FileBuffer f) 
		{
			ins.addSeqFile(f);
		}

		public void insertSeqFile(long seqPos, FileBuffer f) 
		{
			ins.insertSeqFile(seqPos, f);
		}
	
		public long getMemoryTax()
		{
			return ins.getMemoryTax() + 1;
		}
		
		public long getMinimumMemoryUsage()
		{
			int estPtrSz = SystemUtils.approximatePointerSize();
			long mem = 8 + (estPtrSz * 2);
			if (sub != null) mem += sub.getMinimumMemoryUsage();
			if (ins != null) mem += sub.getMinimumMemoryUsage();
			return mem;
		}
		
	}
	
	private static class TrueOffset
	{
		protected long truOff;
		protected Change cRecord;
		protected int cOff;
		
		public TrueOffset(long truOff)
		{
			this.truOff = truOff;
			this.cRecord = null;
			this.cOff = -1;
		}
		
		public TrueOffset(Change c, int cInd)
		{
			this.truOff = -1;
			this.cRecord = c;
			this.cOff = cInd;
		}
		
		public boolean isInMod()
		{
			return (this.cRecord != null);
		}
	}
	
	private class MiniBuffer
	{
		private int sectorIndex;
		private byte[] contents;
		
		private int filled; //Number of bytes from file in buffer (rest are unused)
		
		public MiniBuffer(int sector) throws IOException
		{
			sectorIndex = sector;
			contents = new byte[subBufferSize];
			//dirty = false;
			filled = 0;
			readFromDisk();
		}
		
		private void readFromDisk() throws IOException
		{
			if (!FileBuffer.fileExists(getStreamPath())) throw new IOException();
			
			long stOff = getAbsFileOffset(sectorIndex);
			long fEnd = getAbsEndOffset();
			
			FileInputStream myStream = new FileInputStream(getStreamPath());
			int rBytes = subBufferSize;
			if (stOff + Integer.toUnsignedLong(contents.length) >= fEnd)
			{
				rBytes = (int)(fEnd - stOff);
				if (rBytes <= 0) 
				{
					myStream.close();
					throw new IOException();
				}
			}
		    myStream.skip(stOff);
		    myStream.read(this.contents, 0, rBytes);
		    this.filled = rBytes;
		    myStream.close();

		}
		
		public byte readByte(int bOff)
		{
			if (bOff < 0) throw new ArrayIndexOutOfBoundsException();
			if (bOff >= this.filled) throw new ArrayIndexOutOfBoundsException();
			
			return this.contents[bOff];
		}
	
		public void printMe()
		{
			System.out.println("Sector Index: " + this.sectorIndex);
			System.out.println("Buffer Size: " + contents.length);
			System.out.println("Bytes filled: 0x" + Integer.toHexString(filled));
			for (int i = 0; i < 16; i++)
			{
				if (i >= contents.length) break;
				System.out.print(String.format("%02x ", contents[i]));
			}
			System.out.println();
		}
		
		public long estimateMemory()
		{
			int estPtrSz = SystemUtils.approximatePointerSize();
			long tot = 4 + 4 + estPtrSz;
			if (contents != null) tot += contents.length;
			return tot;
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
	
	/* ----- CONSTRUCTION ----- */
	
	/**
	 * Construct a stream buffer for the file on disk at the provided path.
	 * @param path Path of file to stream.
	 * @throws IOException If there is an error reading the file.
	 */
	public StreamBuffer(String path) throws IOException
	{
		this(path, 0, FileBuffer.fileSize(path), true);
	}
	
	/**
	 * Construct a stream buffer for the file on disk at the provided path, with offset relative
	 * to a given start position.
	 * @param path Path of file to stream.
	 * @param stOff Start offset of stream. All offsets will be relative to this.
	 * @throws IOException If there is an error reading the file.
	 */
	public StreamBuffer(String path, long stOff) throws IOException
	{
		this(path, stOff, FileBuffer.fileSize(path) - stOff, true);
	}
	
	/**
	 * Construct a stream buffer for the file on disk at the provided path, with offset relative
	 * to a given start position.
	 * @param path Path of file to stream.
	 * @param stOff Start offset of stream. All offsets will be relative to this.
	 * @param len Length, in bytes, of file from offset to include in stream.
	 * @throws IOException If there is an error reading the file.
	 */
	public StreamBuffer(String path, long stOff, long len) throws IOException
	{
		this(path, stOff, len, true);
	}
	
	/**
	 * Construct a stream buffer for the file on disk at the provided path.
	 * @param path Path of file to stream.
	 * @param isBigEndian Whether multi-byte values should be read and written in Big-Endian byte order.
	 * @throws IOException If there is an error reading the file.
	 */
	public StreamBuffer(String path, boolean isBigEndian) throws IOException
	{
		this(path, 0, FileBuffer.fileSize(path),isBigEndian);
	}
	
	/**
	 * Construct a stream buffer for the file on disk at the provided path, with offset relative
	 * to a given start position.
	 * @param path Path of file to stream.
	 * @param stOff Start offset of stream. All offsets will be relative to this.
	 * @param isBigEndian Whether multi-byte values should be read and written in Big-Endian byte order.
	 * @throws IOException If there is an error reading the file.
	 */
	public StreamBuffer(String path, long stOff, boolean isBigEndian) throws IOException
	{
		this(path, stOff, FileBuffer.fileSize(path) - stOff,isBigEndian);
	}
	
	/**
	 * Construct a stream buffer for the file on disk at the provided path, with offset relative
	 * to a given start position.
	 * @param path Path of file to stream.
	 * @param stOff Start offset of stream. All offsets will be relative to this.
	 * @param len Length, in bytes, of file from offset to include in stream.
	 * @param isBigEndian Whether multi-byte values should be read and written in Big-Endian byte order.
	 * @throws IOException If there is an error reading the file.
	 */
	public StreamBuffer(String path, long stOff, long len, boolean isBigEndian) throws IOException
	{
		super(0);
		this.subBufferSize = DEFO_SUBBUF_SIZE;
		this.subBufferNum = DEFO_SUBBUF_NUM;
		this.constructorCore(path, stOff, len, isBigEndian);
	}
	
	/**
	 * Construct a stream buffer for the file on disk at the provided path.
	 * @param path Path of file to stream.
	 * @param subBufferSize Size of each buffer page.
	 * @param subBufferNum Number of buffer pages.
	 * @throws IOException If there is an error reading the file.
	 */
	public StreamBuffer(String path, int subBufferSize, int subBufferNum) throws IOException
	{
		super(0);
		this.subBufferSize = subBufferSize;
		this.subBufferNum = subBufferNum;
		this.constructorCore(path, 0, FileBuffer.fileSize(path), true);
	}
	
	protected StreamBuffer(String tName, long potentialSize, boolean isBigEndian, boolean createNew) throws IOException
	{
		String tPath = FileBuffer.generateTemporaryPath(tName);
		commonCore(isBigEndian);
		this.subBufferSize = DEFO_SUBBUF_SIZE;
		this.subBufferNum = DEFO_SUBBUF_NUM;
		if (createNew) Files.createFile(Paths.get(tPath));
		int fileSecs = (int)(potentialSize / Integer.toUnsignedLong(this.subBufferSize)) + 1;
		locTable = new int[fileSecs];
		for (int i = 0; i < locTable.length; i++) locTable[i] = -1;
		this.filePath = tPath;
		this.streamStart = 0;
		this.streamLength = 0;
		super.unsetReadOnly();
	}
	
	private void constructorCore(String path, long stOff, long len, boolean isBigEndian) throws IOException
	{
		commonCore(isBigEndian);
		int fileSecs = (int)((len - stOff) / Integer.toUnsignedLong(this.subBufferSize)) + 1;
		locTable = new int[fileSecs];
		for (int i = 0; i < locTable.length; i++) locTable[i] = -1;
		this.filePath = path;
		this.streamStart = stOff;
		this.streamLength = len;
		this.initialRead();
		super.setReadOnly(); // DEFAULTS TO READ ONLY. HAVE TO SET IT IF WANT TO WRITE!!
	}
	
	private void commonCore(boolean isBigEndian)
	{
		contents = new MiniBuffer[this.subBufferNum];
		usageQueue = new ArrayDeque<Integer>(contents.length + 1);
		super.setEndian(isBigEndian);
		this.modifications = new HashMap<Long, Change>();
		this.modMemThreshhold = this.subBufferSize * this.subBufferNum * 3;
		this.tempPath = null;
		this.tempMode = false;
	}
	
	/* ----- READ STREAM ----- */
	
	private String getStreamPath()
	{
		if (!this.tempMode) return this.filePath;
		else return this.tempPath;
	}
	
	private void initialRead() throws IOException
	{
		if (!fileExists(getStreamPath()))
		{
			throw new IOException();
		}
		for (int i = 0; i < contents.length; i++)
		{
			if (!sectorIndexValid(i)) break;
			readSector(i);
		}
	}
	
	private void readSector(int secInd) throws IOException
	{
		//Handles everything including finding free sectors, updating lists and stacks, and all that.
		//Check to see if sector is already loaded
		if (locTable[secInd] >= 0 && locTable[secInd] < contents.length) return;
		//System.out.println("StreamBuffer.readSector || Sector is not loaded ( locTable[" + secInd + "] = " + locTable[secInd] + " )");
		MiniBuffer mb = new MiniBuffer(secInd); 
		int insSec = this.findFreeSector();
		//System.out.println("StreamBuffer.readSector || Free internal sector found: Sector " + insSec);
		if (insSec < 0) throw new IllegalStateException();
		this.contents[insSec] = mb;
		locTable[secInd] = insSec;
		//System.out.println("StreamBuffer.readSector ||  Marked in location table: External Sector " + secInd + " at Internal Sector " + insSec);
		usageQueue.push(insSec);
		//System.out.println("StreamBuffer.readSector || MiniBuffer contents: ");
		//mb.printMe();
	}
	
	/**
	 * Clear out all of the buffers, modifications, and cache records.
	 * <br>WARNING: If this buffer is write-enabled, this will delete any modifications
	 * since the last write to disk!
	 */
	public void flush()
	{
		//Deletes all modifications and clears all read buffers.
		this.usageQueue.clear();
		for (int i = 0; i < contents.length; i++) contents[i] = null;
		for (int j = 0; j < locTable.length; j++) contents[j] = null;
		modifications.clear();
	}
	
	/**
	 * Clear out all of the buffers, modifications, and cache records, and reload buffers with
	 * the first (total buffer size) bytes of the file.
	 * <br>WARNING: If this buffer is write-enabled, this will delete any modifications
	 * since the last write to disk!
	 * @throws IOException If file cannot be re-read.
	 */
	public void refresh() throws IOException
	{
		flush();
		initialRead();
	}
	
	private void cleanTemp() throws IOException
	{
		this.tempMode = false;
		Files.deleteIfExists(Paths.get(tempPath));
		this.tempPath = null;
	}
	
	/* ----- SECTOR MANAGEMENT ----- */
	
	private boolean sectorIndexValid(int secInd)
	{
		if (secInd < 0) return false;
		int lastS = getSectorIndex(this.streamLength);
		return (secInd <= lastS);
	}
	
	/**
	 * Return the index of a free buffer page to use.
	 * @return Index of the first free buffer page in array.
	 */
	public int firstFreeSector()
	{
		for (int i = 0; i < contents.length; i++)
		{
			if (contents[i] == null) return i;
		}
		if (usageQueue.size() <= 0) return -1;
		return usageQueue.getLast();
	}
	
	/**
	 * Forcible frees a buffer page.
	 * @param sectorInd Index of sector/page to free
	 */
	public void freeSector(int sectorInd)
	{
		if (sectorInd < 0) return;
		if (sectorInd >= contents.length) return;
		if (contents[sectorInd] == null) return;
		int fileSec = contents[sectorInd].sectorIndex;
		//if (!this.readOnly) contents[sectorInd].clean();
		contents[sectorInd] = null;
		locTable[fileSec] = -1;
		Integer i = sectorInd;
		usageQueue.removeLastOccurrence(i);
	}
	
	/**
	 * Find a sector/page that is free or can be freed, free it, and return its index.
	 * @return Index of freed buffer page.
	 */
	public int findFreeSector()
	{
		int myFree = -1;
		for (int i = 0; i < contents.length; i++)
		{
			if (contents[i] == null) return i;
		}
		try
		{
			myFree = usageQueue.removeLast();	
		}
		catch (NoSuchElementException e)
		{
			return -1;
		}
		if (contents[myFree] == null) return myFree;
		int fileSec = contents[myFree].sectorIndex;
		//if (!this.readOnly) contents[myFree].clean();
		contents[myFree] = null;
		locTable[fileSec] = -1;
		
		return myFree;
	}
	
	private byte ROGetByte(long relpos)
	{
		int fSec = getSectorIndex(relpos);
		int sOff = getSectorOffset(relpos);
		if (this.locTable[fSec] >= 0)
		{
			//Buffer containing this section already in memory
			int bSec = this.locTable[fSec];
			this.moveToQTop(bSec);
			return this.contents[bSec].readByte(sOff);
		}
		else
		{
			try 
			{
				this.readSector(fSec);
				int bSec = this.locTable[fSec];
				return this.contents[bSec].readByte(sOff);
			} 
			catch (IOException e) 
			{
				e.printStackTrace();
				throw new NullPointerException();
			}
		}
	}
	
	private byte WEGetByte(long pos)
	{
		TrueOffset tru = this.findTrueOffset(pos);
		if (tru.isInMod())
		{
			return tru.cRecord.getByte(tru.cOff);
		}
		return this.ROGetByte(tru.truOff);
	}
	
	private boolean moveToQTop(int sector)
	{
		if (usageQueue.getFirst() == sector) return true;
		if (usageQueue.removeLastOccurrence(sector))
		{
			usageQueue.push(sector);
			return true;
		}
		return false;
	}
	
	public void printSectorInfo(int internalSec)
	{
		if (internalSec < 0) return;
		if (internalSec >= contents.length) return;
		contents[internalSec].printMe();
	}
	
	/* ----- GETTERS ----- */
	
	  /**
	   * Get the minimum amount of memory it takes to hold this buffer.
	   * @return Amount of memory buffer occupies.
	   */
	  public long getMemoryBurden()
	  {
		  long mem = (long)this.subBufferNum * (long)this.subBufferSize;
		  if (this.readOnly()) return mem;
		  return mem + countMemoryTax();
	  }
	  
	  public long getMemoryBurden(Collection<FileBuffer> inMem)
	  {
		  return getMemoryBurden();
	  }
	  
	  public long getMinimumMemoryUsage()
	  {
	  	int estPtrSz = SystemUtils.approximatePointerSize();
	  	long tot = super.getMinimumMemoryUsage();
	  	tot += (4 * 4) + (1 * 1) + (8 * 2) + (estPtrSz * 6);
	  	//contents
	  	if (contents != null)
	  	{
	  		for (MiniBuffer buff : contents) tot += buff.estimateMemory() + estPtrSz;
	  	}
	  	//locTable
	  	if (locTable != null) tot += 4 * locTable.length;
	  	//usageQueue DOES NOT TAKE INTO ACCOUNT overhead for linked list structure!!
	  	if (usageQueue != null) tot += (8 + estPtrSz) * usageQueue.size();
	  	//Strings
	  	if (filePath != null) tot += filePath.length();
	  	if (tempPath != null) tot += tempPath.length();
	  	//Changes DOES NOT TAKE INTO ACCOUNT overhead for map structure!!
	  	if (modifications != null)
	  	{
	  		Collection<Change> changes = modifications.values();
	  		for (Change c : changes)
	  		{
	  			tot += estPtrSz;
	  			tot += 10;
	  			tot += c.getMinimumMemoryUsage();
	  		}
	  	}
	  	return tot;
	  }
		
	  
	/* ----- SETTERS ----- */
	
	  /**
	   * Enable writing for this stream.
	   * <br>WARNING: All edits will be held in memory until written to disk.
	   * Use sparingly.
	   */
	  public void setWritable()
	  {
		  super.unsetReadOnly();
	  }
	
	  /**
	   * Set the size in bytes the buffer for writing to disk will occupy in memory.
	   * <br>This is in addition to the buffer stream.
	   * @param newSize Desired size of the disk write buffer.
	   */
	  public void setWriteBufferSize(int newSize)
	  {
		  this.writeBuffHalfSize = newSize/2;
	  }
	
	  /**
	   * Set the amount of space, in bytes, that a writable buffer can hold in memory
	   * for modifications.
	   * @param newSize The desired amount of space for write modifications.
	   */
	  public void setModBufferSize(int newSize)
	  {
		  this.modMemThreshhold = newSize;
	  }
	  
	/* ----- MODIFICATIONS ----- */
	
	private long countModSize()
	{
		long c = 0;
		Collection<Change> mods = this.modifications.values();
		for (Change m : mods)
		{
			c += m.getSize();
		}
		
		return c;
	}
		
	private long countMemoryTax()
	{
		long c = 0;
		Collection<Change> mods = this.modifications.values();
		for (Change m : mods)
		{
			c += m.getMemoryTax();
		}
		
		return c;
	}
	
	private void createDeletion(long position, long length)
	{
		//Delete all modifications in the area.
		for (long i = position; i < position + length; i++)
		{
			if (this.modifications.containsKey(i)) this.modifications.remove(i);
		}
		Change c = new Deletion(position, length);
		this.modifications.put(position, c);
	}
	
	private void createSubstitution(long position, byte bSub)
	{
		if (this.modifications.containsKey(position))
		{
			Change m = this.modifications.get(position);
			if (m.getType() == ChangeType.DELETION) return;
			else if(m.getType() == ChangeType.SUBSTITUTION)
			{
				m.replaceSeqByte(bSub, 0);
				return;
			}
			else if(m.getType() == ChangeType.INSERT)
			{
				Substitution s = new Substitution(position, bSub);
				Insertion i = (Insertion)m;
				Change c = new Inssub(s, i);
				modifications.remove(position);
				modifications.put(position, c);
				return;
			}
			else if(m.getType() == ChangeType.INSSUB)
			{
				m.replaceSeqByte(bSub, (int)m.getLength());
				return;
			}
		}
		else
		{
			Change c = new Substitution(position, bSub);
			modifications.put(position, c);
		}
	}
	
	private void createInsertion(long position, FileBuffer contents)
	{
		if (this.modifications.containsKey(position))
		{
			Change m = this.modifications.get(position);
			if (m.getType() == ChangeType.DELETION) return;
			else if(m.getType() == ChangeType.SUBSTITUTION)
			{
				Substitution s = (Substitution)m;
				Insertion i = new Insertion(position, contents);
				Change c = new Inssub(s, i);
				modifications.remove(position);
				modifications.put(position, c);
			}
			else if(m.getType() == ChangeType.INSERT || m.getType() == ChangeType.INSSUB)
			{
				m.insertSeqFile(0, contents);
			}
		}
		else
		{
			Change c = new Insertion(position, contents);
			modifications.put(position, c);
		}
	}
	
	private void createInsertion(long position, byte b)
	{
		if (this.modifications.containsKey(position))
		{
			Change m = this.modifications.get(position);
			if (m.getType() == ChangeType.DELETION) return;
			else if(m.getType() == ChangeType.SUBSTITUTION)
			{
				Substitution s = (Substitution)m;
				FileBuffer f = new FileBuffer(StreamBuffer.DEFO_SUBBUF_SIZE);
				f.addToFile(b);
				Insertion i = new Insertion(position, f);
				Change c = new Inssub(s, i);
				modifications.remove(position);
				modifications.put(position, c);
			}
			else if(m.getType() == ChangeType.INSERT || m.getType() == ChangeType.INSSUB)
			{
				m.insertSeqByte(b, 0);
			}
		}
		else
		{
			FileBuffer f = new FileBuffer(StreamBuffer.DEFO_SUBBUF_SIZE);
			f.addToFile(b);
			Change c = new Insertion(position, f);
			modifications.put(position, c);
		}
	}

	/* ----- OFFSET CALCULATIONS ----- */
	
	private TrueOffset findTrueOffset(long userOffset)
	{
		if (userOffset < 0) return null;
		long uPos = 0;
		long truPos = 0;
		while (uPos < userOffset)
		{
			if (this.modifications.containsKey(truPos))
			{
				Change m = this.modifications.get(truPos);
				uPos += m.getSize();
				if ((m.getType() == ChangeType.INSERT && uPos > userOffset) || (m.getType() == ChangeType.INSSUB && uPos >= userOffset))
				{
					//It's in the insertion
					return new TrueOffset(m, (int)(userOffset - (uPos - m.getSize())));
				}
				uPos++;
				truPos++;
			}
			else
			{
				uPos++;
				truPos++;
			}
		}
		
		if (this.modifications.containsKey(truPos))
		{
			Change m = this.modifications.get(truPos);
			return new TrueOffset(m, 0);
		}
		if (truPos > this.streamLength) return null;
		return new TrueOffset(truPos);
	}
	
	private long getRelFileOffset(int sectorIndex)
	{
		long foff = Integer.toUnsignedLong(sectorIndex) * Integer.toUnsignedLong(this.subBufferSize);
		return foff;
	}
	
	private long getAbsFileOffset(int sectorIndex)
	{
		return this.getRelFileOffset(sectorIndex) + this.streamStart;
	}
	
	private int getSectorIndex(long relFileOffset)
	{
		return (int)(relFileOffset / Integer.toUnsignedLong(this.subBufferSize));
	}
	
	private int getSectorOffset(long relFileOffset)
	{
		return (int)(relFileOffset % Integer.toUnsignedLong(this.subBufferSize));
	}
 	
	/**
	 * Get the size of the source file - the "end offset" in absolute coordinates.
	 * @return Long integer representing the full size of the file being streamed.
	 */
	public long getAbsEndOffset()
	{
		return this.streamStart + this.streamLength;
	}
	
	/**
	 * Get the size in bytes of the stream. This includes only the size of the chunk of file
	 * that can be accessed by this stream buffer.
	 * @return The size of the stream.
	 */
	public long getStreamSize()
	{
		return this.streamLength;
	}
	
	/* ----- FileBuffer READ OVERRIDES ----- */
	
	/**
	 * @throws IndexOutOfBoundsException If position is invalid.
	 * @throws NullPointerException If piece of stream could not be loaded (IOException)
	 */
	public byte getByte(int position)
	{
		return this.getByte((long)position);
	}
	
	/**
	 * @throws IndexOutOfBoundsException If position is invalid.
	 * @throws NullPointerException If piece of stream could not be loaded (IOException)
	 */
	public byte getByte(long position)
	{
		if (position < 0) throw new IndexOutOfBoundsException();
		if (position > this.getFileSize()) throw new IndexOutOfBoundsException();
		if (this.readOnly()) return this.ROGetByte(position);
		else return this.WEGetByte(position);
	}
  
	/* ----- FileBuffer WRITE OVERRIDES ----- */
  	
	private void addCheck()
	{
  		if (this.readOnly()) throw new UnsupportedOperationException(); 
  		if (countMemoryTax() >= this.modMemThreshhold) dumpToDisk();
	}
	
	private int addCheck(long pos)
	{
		if (this.isEmpty()) throw new IndexOutOfBoundsException();
		int val = 0;
  		if (this.readOnly()) throw new UnsupportedOperationException(); 
  		long fSize = this.getFileSize();
  		if (pos < 0 || pos > fSize) throw new IndexOutOfBoundsException();
  		if (pos == fSize) val = 1;
  		if (countMemoryTax() >= this.modMemThreshhold) dumpToDisk();
  		return val;
	}
	
	private void dumpToDisk()
	{
		  try 
		  {
			  this.tempPath = this.generateTempPath(this.filePath);
			  this.tempMode = true;
			  this.writeMe(tempPath, 0, this.getFileSize(), false);
		  } 
		  catch (IOException e) 
		  {
			  e.printStackTrace();
			  throw new UnsupportedOperationException(); 
		  }
	}
	
    /**
     * @throws UnsupportedOperationException If buffer is read-only locked.
     * @throws IndexOutOfBoundsException If any offsets given are invalid.
     */
  	public void deleteFromFile(int stOff)
    {
  		this.deleteFromFile((long)stOff);
    }
  	
    /**
     * @throws UnsupportedOperationException If buffer is read-only locked.
     * @throws IndexOutOfBoundsException If any offsets given are invalid.
     */
    public void deleteFromFile(long stOff)
    {
    	if (this.readOnly()) throw new UnsupportedOperationException(); 
    	this.deleteFromFile(stOff, this.getFileSize());
    }
    
    /**
     * @throws UnsupportedOperationException If buffer is read-only locked.
     * @throws IndexOutOfBoundsException If any offsets given are invalid.
     */
    public void deleteFromFile(int stOff, int edOff)
    {
    	this.deleteFromFile((long)stOff, (long)edOff);
    }
    
    /**
     * @throws UnsupportedOperationException If buffer is read-only locked.
     * @throws IndexOutOfBoundsException If any offsets given are invalid.
     */
    public void deleteFromFile(long stOff, long edOff)
    {
    	if (this.readOnly()) throw new UnsupportedOperationException(); 
    	if (stOff >= edOff) throw new IndexOutOfBoundsException();
    	TrueOffset stTru = this.findTrueOffset(stOff);
    	if (stTru == null) throw new IndexOutOfBoundsException();
    	TrueOffset edTru = this.findTrueOffset(edOff);
    	if (edTru == null) throw new IndexOutOfBoundsException();
    	if (!stTru.isInMod() && !edTru.isInMod())
    	{
    		long stPos = stTru.truOff;
    		long edPos = edTru.truOff;
    		this.createDeletion(stPos, edPos - stPos);
    	}
    	else if(stTru.isInMod() && !edTru.isInMod())
    	{
    		long edPos = edTru.truOff;
    		long stend = stTru.cRecord.getPosition() + stTru.cRecord.getLength();
    		this.createDeletion(stend, edPos - stend);
    		stTru.cRecord.deleteBytes(stTru.cOff, (int)(stTru.cRecord.getLength() - stTru.cOff));
    	}
    	else if(!stTru.isInMod() && edTru.isInMod())
    	{
    		long stPos = stTru.truOff;
    		long edSt = edTru.cRecord.getPosition();
    		this.createDeletion(stPos, edSt - stPos);
    		edTru.cRecord.deleteBytes(0, edTru.cOff);
    	}
    	else if(stTru.isInMod() && edTru.isInMod())
    	{
    		long stend = stTru.cRecord.getPosition() + stTru.cRecord.getLength();
    		long edSt = edTru.cRecord.getPosition();
    		this.createDeletion(stend, edSt - stend);
    		stTru.cRecord.deleteBytes(stTru.cOff, (int)(stTru.cRecord.getLength() - stTru.cOff));
    		edTru.cRecord.deleteBytes(0, edTru.cOff);
    	}
    	
    	
    }
	
    public void addToFile(byte i8)
    {  
    	addCheck();
    	long lPos = FileBuffer.fileSize(this.getStreamPath());
    	if (modifications.containsKey(lPos))
    	{
    		Change m = modifications.get(lPos);
    		m.addSeqByte(i8);
    		return;
    	}
    	else
    	{
    		FileBuffer f = new FileBuffer(StreamBuffer.DEFO_SUBBUF_SIZE);
    		f.addToFile(i8);
    		this.createInsertion(lPos, f);
    	}
    }
    
    public void addToFile(byte i8, int position)
    {
    	if (this.readOnly()) throw new UnsupportedOperationException(); 
    	this.addToFile(i8, (long)position);
    }
    
    /**
     * @throws UnsupportedOperationException If buffer is read-only locked.
     * @throws IndexOutOfBoundsException If any offsets given are invalid.
     */
    public void addToFile(byte i8, long position)
    {
    	if (this.addCheck(position) == 1) 
    	{
    		addToFile(i8); 
    		return;
    	}
    	TrueOffset tPos = this.findTrueOffset(position);
    	if (tPos.isInMod())
    	{
    		Change c = tPos.cRecord;
    		c.insertSeqByte(i8, tPos.cOff);
    		return;
    	}
    	else
    	{
    		this.createInsertion(tPos.truOff, i8);
    	}
    } 
    
    /**
     * @throws UnsupportedOperationException If buffer is read-only locked.
     * @throws IndexOutOfBoundsException If any offsets given are invalid.
     */
    public void addToFile(short i16, int position)
    {
    	addToFile(i16, (long)position);
    }
    
    /**
     * @throws UnsupportedOperationException If buffer is read-only locked.
     * @throws IndexOutOfBoundsException If any offsets given are invalid.
     */
    public void addToFile(short i16, long position)
    {
    	if (this.addCheck(position) == 1) 
    	{
    		addToFile(i16); 
    		return;
    	}
  	  //this.moveAllBytesUp(position, 2L);
  	  
  	  byte[] myBytes = numToByStr(i16);
  	  if (!super.isBigEndian())
  	  {
			byte temp = myBytes[0];
			myBytes[0] = myBytes[1];
			myBytes[1] = temp;
  	  }
  	  for (int i = 1; i >= 0; i--) this.addToFile(myBytes[i], position + (1 - i));
  		  
  		  //this.replaceByte(myBytes[i], position + (1 - i));
  	  
    }
    
    /**
     * @throws UnsupportedOperationException If buffer is read-only locked.
     * @throws IndexOutOfBoundsException If any offsets given are invalid.
     */
    public void addToFile(int i32, int position)
    {
    	addToFile(i32, (long)position);
    }
    
    /**
     * @throws UnsupportedOperationException If buffer is read-only locked.
     * @throws IndexOutOfBoundsException If any offsets given are invalid.
     */
    public void addToFile(int i32, long position)
    {
    	if (this.addCheck(position) == 1) 
    	{
    		addToFile(i32); 
    		return;
    	}
  	  //this.moveAllBytesUp(position, 4L);
  	  
  	  int myInt = 0;
  	  if (super.isBigEndian()) myInt = i32;
  	  else myInt = switchByO(i32);
  	  byte[] myBytes = numToByStr(myInt);
  	  for (int i = 3; i >= 0; i--) this.addToFile(myBytes[i], position + (3 - i));	
    }
    
    /**
     * @throws UnsupportedOperationException If buffer is read-only locked.
     * @throws IndexOutOfBoundsException If any offsets given are invalid.
     */
    public void addToFile(long i64, int position)
    {
    	addToFile(i64, (long)position);
    }
    
    /**
     * @throws UnsupportedOperationException If buffer is read-only locked.
     * @throws IndexOutOfBoundsException If any offsets given are invalid.
     */
    public void addToFile(long i64, long position)
    {
    	if (this.addCheck(position) == 1) 
    	{
    		addToFile(i64); 
    		return;
    	}
  	  //this.moveAllBytesUp(position, 8L);
  	  
  	  long myLong = 0;
  	  if (super.isBigEndian()) myLong = i64;
  	  else myLong = switchByO(i64);
  	  byte[] myBytes = numToByStr(myLong);
  	  for (int i = 7; i >= 0; i--) this.addToFile(myBytes[i], position + (7 - i)); 
    }
    
    /**
     * @throws UnsupportedOperationException If buffer is read-only locked.
     * @throws IndexOutOfBoundsException If any offsets given are invalid.
     */
    public void add24ToFile(int i24, int position)
    {
    	add24ToFile(i24, (long)position);
    }
    
    /**
     * @throws UnsupportedOperationException If buffer is read-only locked.
     * @throws IndexOutOfBoundsException If any offsets given are invalid.
     */
    public void add24ToFile(int i24, long position)
    {    	
    	if (this.addCheck(position) == 1) 
    	{
    		add24ToFile(i24); 
    		return;
    	}
  	 
  	  //this.moveAllBytesUp(position, 3L);
  	  
  	  byte[] myBytes = numToByStr(i24);
  	  if (super.isBigEndian()) for (int i = 2; i >= 0; i--) this.replaceByte(myBytes[i], position + (2 - i));	
  	  else for (int i = 0; i < 3; i++) this.addToFile(myBytes[i], position + (2 - i));		
    }

    /**
     * @throws UnsupportedOperationException If buffer is read-only locked, or if there is an IOException.
     * @throws IndexOutOfBoundsException If any offsets given are invalid.
     */
    public void addToFile(FileBuffer addition, long insertPos, long stPos, long edPos)
    {
      if (addition == null || addition.isEmpty()) throw new NullPointerException();
      if (!(addition.offsetValid(stPos)) || !addition.offsetValid(edPos - 1)) throw new IndexOutOfBoundsException();
  		if (this.addCheck(insertPos) == 1) 
  		{
  			addToFile(addition, stPos, edPos); 
  			return;
  		}
  	  
  	  try 
  	  {
  		  FileBuffer add = addition.createReadOnlyCopy(stPos, edPos);
  		  this.createInsertion(insertPos, add);
  	  } 
  	  catch (IOException e) 
  	  {
  		  e.printStackTrace();
  		  throw new UnsupportedOperationException();
  	  }
    }
    
    /**
     * @throws UnsupportedOperationException If buffer is read-only locked, or if there is an IOException.
     * @throws IndexOutOfBoundsException If any offsets given are invalid.
     */
    public void addToFile(FileBuffer addition, long stPos, long edPos)
    {
        if (addition == null || addition.isEmpty()) throw new NullPointerException();
        if (!(addition.offsetValid(stPos)) || !addition.offsetValid(edPos - 1)) throw new IndexOutOfBoundsException();
    	this.addCheck();

  	  try 
  	  {
  		  FileBuffer add = addition.createReadOnlyCopy(stPos, edPos);
  		  if(this.modifications.containsKey(this.streamLength))
  		  {
  			  this.modifications.get(this.streamLength).addSeqFile(add);
  			  return;
  		  }
  		  else
  		  {
  			  this.createInsertion(this.streamLength, add);
  			  return;
  		  }
  	  } 
  	  catch (IOException e) 
  	  {
  		  e.printStackTrace();
  		  throw new UnsupportedOperationException();
  	  }
  	  
    }

    /**
     * @throws UnsupportedOperationException If buffer is read-only locked, or if there is an IOException.
     * @throws IndexOutOfBoundsException If any offsets given are invalid.
     */
    public void addToFile(FileBuffer addition, int insertPos, int stPos, int edPos)
    {
    	if (this.readOnly()) throw new UnsupportedOperationException(); 
    	long iPos = (long)insertPos;
    	long stOff = (long)stPos;
    	long edOff = (long)edPos;
    	this.addToFile(addition, iPos, stOff, edOff);
    }
    
    /**
     * @throws UnsupportedOperationException If buffer is read-only locked, or if there is an IOException.
     * @throws IndexOutOfBoundsException If any offsets given are invalid.
     */
    public void addToFile(FileBuffer addition)
    {
    	this.addToFile(addition, 0L, addition.getFileSize());
    }
    
    /**
     * @throws UnsupportedOperationException If buffer is read-only locked, or if there is an IOException.
     * @throws IndexOutOfBoundsException If any offsets given are invalid.
     */
    public void addToFile(FileBuffer addition, int stPos, int edPos)
    {
    	addToFile(addition, (long)stPos, (long)edPos);
    }

    /**
     * @throws UnsupportedOperationException If buffer is read-only locked.
     * @throws IndexOutOfBoundsException If position is invalid.
     */
    public boolean replaceByte(byte b, int position)
    {
    	return this.replaceByte(b, (long)position);
    }
    
    /**
     * @throws UnsupportedOperationException If buffer is read-only locked.
     * @throws IndexOutOfBoundsException If position is invalid.
     */
    public boolean replaceByte(byte b, long position)
    {
    	if (this.readOnly()) throw new UnsupportedOperationException(); 
    	//TrueOffset t = new TrueOffset(position);
    	TrueOffset t = this.findTrueOffset(position);
    	if (t == null) throw new IndexOutOfBoundsException();
    	if(!t.isInMod())
    	{
    		createSubstitution(position, b);
    		return true;
    	}
    	else
    	{
    		Change m = t.cRecord;
    		if (m.getType() == ChangeType.SUBSTITUTION)
    		{
    			m.addSeqByte(b);
    			return true;
    		}
    		else if (m.getType() == ChangeType.INSERT)
    		{
    			m.replaceSeqByte(b, t.cOff);
    		}
    		else if (m.getType() == ChangeType.INSSUB)
    		{
    			m.replaceSeqByte(b, t.cOff);
    			//Not sure if that's right.
    		}
    	}
    	return false;
    }
    
    /**
     * @throws UnsupportedOperationException If buffer is read-only locked.
     * @throws IndexOutOfBoundsException If position is invalid.
     */
    public boolean replaceShort(short s, long position)
    {
    	if (this.readOnly()) throw new UnsupportedOperationException(); 
  	  	if (position < 0 || position >= this.getFileSize() - 1) return false;
	  
  	  	byte[] sBytes = numToByStr(s);
    	if (!super.isBigEndian())
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
     * @throws UnsupportedOperationException If buffer is read-only locked.
     * @throws IndexOutOfBoundsException If position is invalid.
     */
    public boolean replaceInt(int i, long position)
    {
    	if (this.readOnly()) throw new UnsupportedOperationException(); 
    	if (position < 0 || position >= this.getFileSize() - 3) return false;
    	if (!super.isBigEndian()) i = switchByO(i);
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
     * @throws UnsupportedOperationException If buffer is read-only locked.
     * @throws IndexOutOfBoundsException If position is invalid.
     */
    public boolean replaceLong(long l, long position)
    {
    	if (this.readOnly()) throw new UnsupportedOperationException(); 
    	if (position < 0 || position >= this.getFileSize() - 7) return false;
    	if (!super.isBigEndian()) l = switchByO(l);
    	byte[] lBytes = numToByStr(l);
    	if (lBytes == null) return false;
    	if (lBytes.length != 8) return false;
  	  
    	for (int p = 0; p < lBytes.length; p++)
    	{
    		this.replaceByte(lBytes[p], position + (7 - p));
    	}
    	return true;
    }
    
    /**
     * @throws UnsupportedOperationException If buffer is read-only locked.
     * @throws IndexOutOfBoundsException If position is invalid.
     * @throws NullPointerException If charset string or provided string is null.
     * @throws IllegalCharsetNameException If the given charset name is illegal.
     * @throws UnsupportedCharsetException If no support for the named charset is available in this instance of the Java virtual machine.
     */
  	public void addEncoded_string(String charset, String myString)
  	{
  		this.addCheck();
  		super.addEncoded_string(charset, myString);
  	}
    
    /**
     * @throws UnsupportedOperationException If buffer is read-only locked.
     * @throws IndexOutOfBoundsException If position is invalid.
     * @throws NullPointerException If charset string or provided string is null.
     * @throws IllegalCharsetNameException If the given charset name is illegal.
     * @throws UnsupportedCharsetException If no support for the named charset is available in this instance of the Java virtual machine.
     */
    public void addEncoded_string(String charset, String myString, long pos)
    {
    	if (charset == null) throw new NullPointerException();
    	if (this.addCheck(pos) == 1) 
    	{
    		this.addEncoded_string(charset, myString);
    		return;
    	}

    	Charset mySet = Charset.forName(charset);
  	  	ByteBuffer bb = mySet.encode(myString);
  	  	//this.moveAllBytesUp(pos, bb.remaining());
  	  	while(bb.hasRemaining())
  	  	{
  	  		byte b = bb.get();
  	  		this.addToFile(b, pos);
  	  		pos++;
  	  	}
    }

    /* ----- FileBuffer DISK WRITE OVERRIDES ----- */
    
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
    
    public void writeFile() throws IOException
    {
    	writeFile(this.filePath, 0L, this.getFileSize());
    }
    
    public void writeFile(String path) throws IOException
    { 
    	writeFile(path, 0L, this.getFileSize());
    }
    
    public void writeFile(int stPos) throws IOException
    {
	  	  writeFile(this.filePath, stPos);
    }
    
    public void writeFile(long stPos) throws IOException
    {
    	writeFile(this.filePath, stPos, this.getFileSize());
    }
    
    public void writeFile(int stPos, int edPos) throws IOException
    {
    	writeFile(this.filePath, stPos, edPos);
    }
    
    public void writeFile(long stPos, long edPos) throws IOException
    {
    	writeFile(this.filePath, stPos, edPos);
    }
    
    public void writeFile(String path, int stPos) throws IOException
    {
    	writeFile(path, Integer.toUnsignedLong(stPos));
    }
    
    public void writeFile(String path, long stPos) throws IOException
    {
    	writeFile(path, stPos, this.getFileSize());
    }
    
    public void writeFile(String path, int stPos, int edPos) throws IOException
    {
    	long stL = Integer.toUnsignedLong(stPos);
    	long edL = Integer.toUnsignedLong(edPos);
    	writeFile(path, stL, edL);
    }

    /**
     * @throws IllegalThreadStateException If there is a threading error in writing to disk - buffer
     * fill and buffer write threads have difficulty joining.
     */
    public void writeFile(String path, long stPos, long edPos) throws IOException
    {
    	this.writeMe(path, stPos, edPos, true);
    }
    
    private void writeMe(String path, long stPos, long edPos, boolean userCall) throws IOException
    {
    	if (path == null) throw new NullPointerException();
    	if (this.readOnly() && path.equals(this.filePath)) return;
	    if (stPos < 0) stPos = 0;
	    if (stPos >= edPos) throw new IndexOutOfBoundsException();
	    if (edPos > this.getFileSize()) edPos = this.getFileSize();
	    
	    if (path.equals(this.getStreamPath()))
	    {
		    String tempPath = this.generateTempPath(path);
		    writeCore(tempPath, stPos, edPos, false);
		    Files.move(Paths.get(tempPath), Paths.get(this.getStreamPath()), StandardCopyOption.REPLACE_EXISTING);
		    this.refresh();
	    }
	    else
	    {
		    writeCore(path, stPos, edPos, false);
	    }
	    if (userCall) cleanTemp();
    }
    
    /**
     * StreamBuffer write core.
     * @param path
     * @param stPos
     * @param edPos
     * @param append
     * @throws IOException
     * @throws IllegalThreadStateException
     */
    private void writeCore(String path, long stPos, long edPos, boolean append) throws IOException
    {
	    //byte[] wb1 = new byte[this.writeBuffHalfSize];
	    //byte[] wb2 = new byte[this.writeBuffHalfSize];
    	long addWeight = this.countModSize();
  		long fSize = this.getFileSize();
  		long bSt = 0;
  		long bEd = (long)writeBuffHalfSize;
  		boolean done = false;
  		//boolean swch = false;
  		byte[] fBuff = new byte[this.writeBuffHalfSize];
  		byte[] wBuff = new byte[this.writeBuffHalfSize];
  		
  		//Fill first buffer
  		int fillSz = this.fillBuff(bSt, bEd, fBuff);
	    
	    FileOutputStream outStream = null;
	    if (!append) outStream = new FileOutputStream(path);
  		while(!done)
  		{
  			if (bEd >= fSize) done = true;
  			byte[] temp = wBuff;
  			wBuff = fBuff;
  			fBuff = temp;
  			//Start thread to fill other buffer
  			FillThread fillThread = null;
  			fillThread = new FillThread(bSt, bEd, fBuff, "StreamBuffer.writeCore|Buffer1Fill" + Long.toHexString(bSt) + ":" + Long.toHexString(bEd));
  			fillThread.start();
  			//Write current full buffer to disk
  			if (!append) outStream.write(wBuff, 0, fillSz);
  			else 
  			{
  				if (fillSz == wBuff.length) Files.write(Paths.get(path), wBuff, StandardOpenOption.APPEND);
  				else
  				{
  					byte[] tBuff = new byte[fillSz];
  					for (int i = 0; i < fillSz; i++) tBuff[i] = wBuff[i];
  					Files.write(Paths.get(path), tBuff, StandardOpenOption.APPEND);
  					
  				}
  			}
  			//Ensure that filler thread has completed its task
  			try 
  			{
				fillThread.join();
			} 
  			catch (InterruptedException e) 
  			{
				e.printStackTrace();
				outStream.close();
				throw new IllegalThreadStateException();
			}
  			fillSz = fillThread.getFillSize();
  			//Calculate next offsets
  			bSt += (long)fillSz;
  			bEd += (long)fillSz;
  			if (bEd > fSize) bEd = fSize;
  		}
	    if (!append) outStream.close();	    
	    if (!this.readOnly()) this.streamLength += addWeight;
    }
    
    private String generateTempPath(String targetPath)
    {
    	String dir = FileBuffer.chopPathToDir(targetPath);
    	String fName = FileBuffer.chopDirFromPath(targetPath);
    	char slash = File.separatorChar;
    	String t = dir + slash + "~" + fName + ".tmp";
    	while(FileBuffer.fileExists(t))
    	{
    		Random r = new Random();
    		t += Integer.toUnsignedString(r.nextInt());
    	}
    	return t;
    }
    	
    /**
     * @throws IllegalThreadStateException If there is a threading error in writing to disk - buffer
     * fill and buffer write threads have difficulty joining.
     */
  	public void appendToFile(String path, long stPos, long edPos) throws IOException, NoSuchFileException
  	{
  		long fSize = this.getFileSize();
  		if (stPos >= edPos) throw new ArrayIndexOutOfBoundsException();
  		if (stPos < 0 ) throw new ArrayIndexOutOfBoundsException();
  		if (edPos > fSize) throw new ArrayIndexOutOfBoundsException();
  		
  		writeCore(path, stPos, edPos, true);
  	}
  
	/* ----- FileBuffer OTHER OVERRIDES ----- */
	
    public boolean isEmpty()
    {
    	if (!FileBuffer.fileExists(this.filePath)) return true;
    	if (this.streamLength == 0 && this.modifications.isEmpty()) return true;
    	return false;
    }
    
	public long getFileSize()
	{
		return FileBuffer.fileSize(filePath) + this.countModSize();
	}
	
	public String getDir()
	{
		return FileBuffer.chopPathToDir(filePath);
	}
	  
	public String getName()
	{
		return FileBuffer.chopPathToFName(filePath);
	}
	  
	public String getExt()
	{
	    return FileBuffer.chopPathToExt(filePath);
	}

	public String getPath()
	{
		return this.filePath;
	}
	
	/**
	 * @throws IndexOutOfBoundsException If any position is invalid, or requested
	 * buffer size exceeds maximum integer value.
	 */
	public ByteBuffer toByteBuffer()
	{
		return this.toByteBuffer(0, this.getFileSize());
	}
	
	/**
	 * @throws IndexOutOfBoundsException If any position is invalid, or requested
	 * buffer size exceeds maximum integer value.
	 */
	public ByteBuffer toByteBuffer(int stPos, int edPos)
	{
		return this.toByteBuffer((long)stPos, (long)edPos);
	}
	
	/**
	 * @throws IndexOutOfBoundsException If any position is invalid, or requested
	 * buffer size exceeds maximum integer value.
	 */
	public ByteBuffer toByteBuffer(long stPos, long edPos)
	{
		long sz = edPos - stPos;
		if (sz > 0x7FFFFFFFL) throw new IndexOutOfBoundsException();
		  if (stPos < 0) throw new IndexOutOfBoundsException();
		  if (edPos >= this.getFileSize()) throw new IndexOutOfBoundsException();
		  if (stPos >= edPos) throw new IndexOutOfBoundsException();
		  int diff = (int)(edPos - stPos);
		  if (diff < 0 || diff > StreamBuffer.DEFO_SUBBUF_NUM * StreamBuffer.DEFO_SUBBUF_SIZE) throw new IndexOutOfBoundsException();
		  ByteBuffer bb = ByteBuffer.allocate(diff);
		  for (long i = stPos; i < edPos; i++)
		  {
			  bb.put(this.getByte(i));
		  }
		  return bb;
	}

	public FileBuffer createCopy(int stPos, int edPos) throws IOException
	{
		return createCopy((long)stPos, (long)edPos);
	}
	
	public FileBuffer createCopy(long stPos, long edPos) throws IOException
	{
		return new StreamBuffer(this.filePath, stPos, edPos - stPos);
	}

	public long getBaseCapacity()
	{
		return 0;
	}
	
	public String toString()
	{
		String s = "";
		s += "Stream Buffer\n";
		s += "-------------\n";
		s += "File Path: " + this.filePath + "\n";
		s += "Stream Start Offset: 0x" + Long.toHexString(this.streamStart) + "\n";
		s += "Stream Length: 0x" + Long.toHexString(this.streamLength) + "\n";
		s += "\n";
		s += "Sub-Buffer Size: " + this.subBufferSize + "\n";
		s += "Sub-Buffer Count: " + this.subBufferNum + "\n";
		s += "Read-Only: " + this.readOnly() + "\n";
		s += "Total Size: 0x" + Long.toHexString(this.getFileSize()) + "\n";
		
		return s;
	}
	
	public String typeString()
	{
		return "StreamBuffer";
	}
	
}
