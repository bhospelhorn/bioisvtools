package waffleoRai_Utils;

/**
 * A sleaker way to stream bits from a byte buffer (FileBuffer object).
 * Also includes static methods for individual bit reading and writing.
 * @author Blythe Hospelhorn
 * @version 1.0.1
 * @since September 8, 2017
 */
public class BitStreamer {

	private FileBuffer myFile;
	
	private long byPos;
	private int biPos;
	
	private byte tempByte;
	private boolean readMode;
	
	public BitStreamer(FileBuffer f, boolean readOnly)
	{
		this(f, 0, readOnly);
	}
	
	public BitStreamer(FileBuffer f, long stPos, boolean readOnly)
	{
		this.myFile = f;
		
		if (stPos < 0) stPos = 0;
		if (stPos >= f.getFileSize()) stPos = f.getFileSize() - 1;
		
		this.byPos = stPos;
		this.biPos = 7;
		this.tempByte = 0;
		this.readMode = readOnly;
	}
	
	/*** Getters ***/
	
	public FileBuffer getFile()
	{
		return this.myFile;
	}
	
	public long getBytePosition()
	{
		return this.byPos;
	}
	
	public int getBitPosition()
	{
		return this.biPos;
	}
	
	public boolean inReadOnlyMode()
	{
		return this.readMode;
	}
	
	/*** Setters ***/
	
	public void setFile(FileBuffer f)
	{
		this.setFile(f, 0);
	}
	
	public void setFile(FileBuffer f, long stPos)
	{
		this.myFile = f;
		
		if (stPos < 0) stPos = 0;
		if (stPos >= f.getFileSize()) stPos = f.getFileSize() - 1;
		this.byPos = stPos;
		this.biPos = 7;
	}
	
	/*** Movers ***/
	
	public void rewind()
	{
		this.byPos = 0;
		this.biPos = 7;
	}
	
	public void rewind(int bits)
	{
		int by = bits / 8;
		int bi = bits % 8;
		this.rewind(by, bi);
	}
	
	public void rewind(int bytes, int bits)
	{
		this.byPos -= bytes;
		if (this.byPos < 0) this.byPos = 0;
		
		this.biPos += bits;
		if (this.biPos >= 8)
		{
			this.biPos = this.biPos % 8;
		}
	}
	
	public void fastForward()
	{
		this.byPos = this.myFile.getFileSize() - 1;
		this.biPos = 0;
	}
	
	public void fastForward(int bits)
	{
		int by = bits/8;
		int bi = bits%8;
		this.fastForward(by,  bi);
	}
	
	public void fastForward(int bytes, int bits)
	{
		this.byPos += bytes;
		if (this.byPos >= this.myFile.getFileSize()) this.byPos = this.myFile.getFileSize();
		
		this.biPos -= bits;
		while (this.biPos < 0)
		{
			this.biPos += 8;
		}
	}
	
	public boolean canMoveForward(int bits)
	{		
		int by = bits/8;
		int bi = bits%8;
		if (this.byPos + by >= this.myFile.getFileSize()) return false;
		if (this.byPos + by == this.myFile.getFileSize() - 1 
				&& this.biPos - bi < 0) return false;
		return true;
	}
	
	public boolean canMoveBackward(int bits)
	{
		int by = bits/8;
		int bi = bits%8;
		if (this.byPos - by < 0) return false;
		if (this.byPos - by == 0 && this.biPos + bi >= 8) return false;
		return true;
	}
	
	/*** Readers ***/
	
	public static boolean readABit(byte b, int bitPos)
	{
		if (bitPos < 0 || bitPos >= 8) return false;
		return readABit((Byte.toUnsignedInt(b)), bitPos);
	}
	
	public static boolean readABit(short s, int bitPos)
	{
		if (bitPos < 0 || bitPos >= 16) return false;
		return readABit((Short.toUnsignedInt(s)), bitPos);
	}
	
	public static boolean readABit(int i, int bitPos)
	{
		if (bitPos < 0 || bitPos >= 32) return false;
		
		int mask = 1 << bitPos;
		int result = i & mask;
		if (result == 0) return false;
		else return true;
	}
	
	public static boolean readABit(long l, int bitPos)
	{
		if (bitPos < 0 || bitPos >= 64) return false;
		
		long mask = 1L << bitPos;
		long result = l & mask;
		if (result == 0) return false;
		else return true;
	}
	
	public boolean readNextBit()
	{
		byte b = this.myFile.getByte(this.byPos);
		boolean bit = readABit(b, this.biPos);
		this.fastForward(1);
		return bit;
	}
	
	public byte readToByte(int bits)
	{
		byte myByte = 0;
		myByte = this.myFile.getBits8(bits, this.byPos, this.biPos);
		this.fastForward(bits);
		return myByte;
	}
	
	public short readToShort(int bits)
	{
		short myShort = this.myFile.getBits16(bits, this.byPos, this.biPos + 8);
		this.fastForward(bits);
		return myShort;
	}
	
	public int readToInt(int bits)
	{
		int myInt = this.myFile.getBits32(bits, this.byPos, this.biPos + 24);
		this.fastForward(bits);
		return myInt;
	}
	
	public long readToLong(int bits)
	{
		long myLong = this.myFile.getBits64(bits, this.byPos, this.biPos + 56);
		this.fastForward(bits);
		return myLong;
	}
	
	/*** Writers ***/
	
	public static byte writeABit(byte target, boolean bit, int bitPos)
	{
		int ti = writeABit(Byte.toUnsignedInt(target), bit, bitPos);
		return (byte)ti;
	}
	
	public static short writeABit(short target, boolean bit, int bitPos)
	{
		int ti = writeABit(Short.toUnsignedInt(target), bit, bitPos);
		return (short)ti;
	}
	
	public static int writeABit(int target, boolean bit, int bitPos)
	{
		int i = target;
		int mask = 1 << bitPos;
		if (bit) i = target | mask;
		else i = target & (~mask);
		return i;
	}
	
	public static long writeABit(long target, boolean bit, int bitPos)
	{
		long l = target;
		long mask = 1L << bitPos;
		if (bit) l = target | mask;
		else l = target & (~mask);
		return l;
	}
	
	private void writeTempByte()
	{
		this.myFile.addToFile(this.tempByte);
		this.tempByte = 0;
		this.byPos++;
	}
	
	public void writeIncompleteTemp()
	{
		if (this.readMode) return;
		if (this.biPos == 7) return;
		else this.writeTempByte();
	}
	
	public boolean writeBits(boolean bit)
	{
		if (readMode) return false;
		
		writeABit(this.tempByte, bit, this.biPos);
		this.biPos--;
		if (this.biPos < 0)
		{
			this.writeTempByte();
			this.biPos = 7;
		}
		
		return true;
	}
	
	public boolean writeBits(byte val, int bits)
	{
		if (readMode) return false;
		if (bits < 0 || bits > 8) return false;
		
		for (int i = bits - 1; i >= 0; i--)
		{
			boolean myBit = readABit(val, i);
			this.writeBits(myBit);
		}
		
		return true;
	}
	
	public boolean writeBits(short val, int bits)
	{
		if (readMode) return false;
		if (bits < 0 || bits > 16) return false;
		
		for (int i = bits - 1; i >= 0; i--)
		{
			boolean myBit = readABit(val, i);
			this.writeBits(myBit);
		}
		
		return true;
	}
	
	public boolean writeBits(int val, int bits)
	{
		if (readMode) return false;
		if (bits < 0 || bits > 32) return false;
		
		for (int i = bits - 1; i >= 0; i--)
		{
			boolean myBit = readABit(val, i);
			this.writeBits(myBit);
		}
		
		return true;
	}
	
	public boolean writeBits(long val, int bits)
	{
		if (readMode) return false;
		if (bits < 0 || bits > 64) return false;
		
		for (int i = bits - 1; i >= 0; i--)
		{
			boolean myBit = readABit(val, i);
			this.writeBits(myBit);
		}
		
		return true;
	}
}
