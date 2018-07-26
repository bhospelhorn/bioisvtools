package waffleoRai_Utils;

import java.util.*;
import waffleoRai_Utils.HuffTable.HuffPoint;
import java.io.IOException;
import java.nio.file.Path;

/*
 * UPDATES
 * 
 * 2017.09.21
 * 	1.0.0 -> 1.1.0 | Fixed for compatibility with FileBuffer changes.
 * 		Changed all FileBuffer positions to be long instead of int
 * 
 * 2018.02.31
 * 	1.1.0 -> 1.2.0 | Finally added an overload to the static encoder that tacks an ASCII header at the top
 * 		Also fixed a dumb bug in the header adding encoder (was inserting everything instead of adding)
 * 
 */

/**
 * Huffman encoding structure. For use with the FileBuffer class.
 * @author Blythe Hospelhorn
 * @version 1.2.0
 * @since February 13, 2018
 */
public class Huffman 
{
	
	/* ~~~~~~~~~~~ Node object for tree usage ~~~~~~~~~~~ */
	private class HuffNode implements Comparable<HuffNode>
	{
		private long data;
		private long frequency;
		private long weight; 
		/*For dataless nodes, a place to store the value of the tree when
		 * weighting in the priority queue*/
		
		private boolean EOF;
		
		public HuffNode(long newFreq)
		{
			this.data = -1;
			this.frequency = newFreq;
			this.weight = -1;
			this.EOF = false;
		}
		public HuffNode(long newData, long newFreq)
		{
			this.data = newData;
			this.frequency = newFreq;
			this.weight = newData;
			this.EOF = false;
		}
		
		public long getData()
		{
			return this.data;
		}
		public long getFrequency()
		{
			return this.frequency;
		}

		public long getWeight()
		{
			return this.weight;
		}
		public void setWeight(long newWeight)
		{
			this.weight = newWeight;
		}
		
		public int hashCode()
		{
			return this.toString().hashCode();
		}
		
		public boolean equals(Object other)
		{
			if (other == null) return false;
			if (!(other instanceof HuffNode)) return false;
			if (other == this) return true;
			
			HuffNode otherHN = (HuffNode)other;
			if (this.EOF != otherHN.EOF) return false;
			if (this.data != otherHN.data) return false;
			if (this.frequency != otherHN.frequency) return false;
			
			return true;
		}
		
		public int compareTo(HuffNode other)
		{
			/*"Greater" = lower priority
			 * Higher frequency, lower weight = greater*/
			if (this.frequency > other.frequency) return 1;
			if (this.frequency < other.frequency) return -1;
			if (this.equals(other)) return 0;
			
			if (this.weight >= 0 && other.weight >= 0)
			{
				if (this.weight > other.weight) return -1;
				if (this.weight < other.weight) return 1;
				return 0;
			}
			
			if (this.weight < 0 && other.weight < 0)
			{
				if (this.weight > other.weight) return -1;
				if (this.weight < other.weight) return 1;
				return 0;	
			}
			
			if (this.weight >= 0 && other.weight < 0) return 1;
			if (this.weight < 0 && other.weight >= 0) return -1;
			
			
			return 0;
			
		}
		
		public boolean hasData()
		{
			if (this.data == -1) return false;
			return true;
		}
		public boolean isEOF()
		{
			return this.EOF;
		}
		public void makeNodeEOF()
		{
			this.data = 0x7FFFFFFFFFFFFFFFL;
			this.frequency = 1;
			this.weight = this.data;
			this.EOF = true;
		}
		
		public String toString()
		{
			String s = "";
			if (this.isEOF())
			{
				s += "EOF Node\n";
			}
			else if (!this.isEOF() && !this.hasData())
			{
				s += "Connector Node\n";
			}
			else if (!this.isEOF() && this.hasData())
			{
				s += "Data Node\n";
			}
			s += Long.toHexString(this.data) + " [" + Long.toString(this.frequency) + "]\n";
			return s;
		}
	}

	/* ~~~~~~~~~~~ Internal Variables ~~~~~~~~~~~ */
	private int bitDepth;
	private HuffTable table;
	private BinTree<HuffNode> HuffTree;
	
	/* ~~~~~~~~~~~ Private Tree Construction ~~~~~~~~~~~ */
	private void setBranchWeight(BinTree<HuffNode> target)
	{
		target.moveToRoot();
		if (target.getCurrentData().hasData())
		{
			target.getCurrentData().setWeight(target.getCurrentData().getData());
		}
		
		while (target.moveRight());
		
		long rightWeight = target.getCurrentData().getWeight();
		target.moveToRoot();
		target.getCurrentData().setWeight(rightWeight);
	}
	
	private void putInQueue(List<BinTree<HuffNode>> myQ, BinTree<HuffNode> target)
	{
		target.moveToRoot();
		
		for (int j = 0; j < myQ.size(); j++)
		{
			if (target.getRootData().compareTo(myQ.get(j).getRootData()) == -1)
			{
				/*target is less than current Q element*/
				myQ.add(j, target);
				return;
			}
			if (target.equals(myQ.get(j)))
			{
				return;
			}
		}
		myQ.add(target);
	}

	private List<BinTree<HuffNode>> constructQueue()
	{
		if (this.table == null) return null;
		
		List<BinTree<HuffNode>> myQueue = new LinkedList<BinTree<HuffNode>>();
		
		HuffNode aNode = null;
		BinTree<HuffNode> aLeaf = null;
		
		
		for (HuffTable.HuffPoint myPoint:this.table.contentsToList(false))
		{
			/*Create a single node tree and shove in queue.*/
			if (myPoint.getFreq() > 0)
			{
				aNode = new HuffNode(myPoint.getSymbol(), myPoint.getFreq());
				aLeaf = new BinTree<HuffNode>(aNode);
				
				/*Find a place to insert it at.*/
				putInQueue(myQueue, aLeaf);
			}
		}
		
		/*Add an EOF "character"*/
		HuffNode EOFNode = createEOFnode();
		aLeaf = new BinTree<HuffNode>(EOFNode);
		putInQueue(myQueue, aLeaf);
		
		//printQ(myQueue);
		/*System.out.println("Initial Queue: ");
		int cnt = 0;
		for (BinTree<HuffNode> n : myQueue)
		{
			System.out.println("----- Element " + cnt);
			System.out.println(n.toString());
			cnt++;
		}*/
		
		return myQueue;
	}
	
	/*private String printQ(List<BinTree<HuffNode>> myQ)
	{
		String s = "";
		int count = 0;
		//System.out.println("My Queue: \n");
		for (BinTree<HuffNode> n : myQ)
		{
			//System.out.println(Long.toHexString(n.getRootData().getData()) + " [" + n.getRootData().getFrequency() + "]");
			s += "----- Element " + Integer.toString(count) + "\n";
			s += n.toString();
			//s += "---------- \n";
			count++;
		}
		return s;
	}*/
	
	private void constructTree()
	{
		if (this.table == null) return;
		
		List<BinTree<HuffNode>> myQ = this.constructQueue();
		BinTree<HuffNode> aBranch = null;
		BinTree<HuffNode> newLeft = null;
		BinTree<HuffNode> newRight = null;
		
		long newFreq = 0;
		
		
		
		/*DEBUG*/
	/*	int debugCnt = 1;
		FileBuffer debugLog = new FileBuffer(50000000); 
		String logPath = "C:\\Users\\Blythe\\Documents\\Programming Projects\\Testers\\Test Case Files\\Huff\\TreeDebug.out";
		
		debugLog.printASCIIToFile("Initial Queue: \n ");
		debugLog.printASCIIToFile(printQ(myQ));
		debugLog.printASCIIToFile("=======================================\n\n");*/
		
		while (myQ.size() > 1)
		{
			/*Create new branch base from two at front of queue*/
			/*Weight new base and insert it back into the queue*/
			//debugLog.printASCIIToFile("============================== Pass # " + debugCnt + " ==============================\n");
			//debugCnt++;
			//debugLog.printASCIIToFile("Queue size before = " + myQ.size() + "\n");
			
			newLeft = myQ.get(0);
			newRight = myQ.get(1);
			myQ.remove(1);
			myQ.remove(0);
			
			newFreq = newLeft.getRootData().getFrequency();
			newFreq += newRight.getRootData().getFrequency();
			
			aBranch = new BinTree<HuffNode>(new HuffNode(newFreq));
			aBranch.insertChildNode(newLeft, false);
			aBranch.moveToRoot();
			aBranch.insertChildNode(newRight, true);
			aBranch.moveToRoot();
			this.setBranchWeight(aBranch);
			
			this.putInQueue(myQ, aBranch);
			//System.out.println("Added to Queue: ROOT = \n" + aBranch.getRootData().toString());
			aBranch = null;
			
			//debugLog.printASCIIToFile("Queue size after = " + myQ.size() + "\n");
			//debugLog.printASCIIToFile("Queue after pass: \n ");
			//debugLog.printASCIIToFile(printQ(myQ));
			//debugLog.printASCIIToFile("=======================================\n\n");
		}
		
		this.HuffTree = myQ.get(0);
		/*try {
			debugLog.writeFile(logPath);
		} catch (IOException e) {
			e.printStackTrace();
		}*/
	}
	
	private HuffNode createEOFnode()
	{
		HuffNode EOF = new HuffNode(1);
		EOF.makeNodeEOF();
		return EOF;
	}
	
	/* ~~~~~~~~~~~ Private Table Management ~~~~~~~~~~~ */
	
	private void setHuffCodesInTable()
	{
		if (this.table == null) return;
		if (this.HuffTree == null) return;
		
		List<BinTree.BinNodeInfo<HuffNode>> treeList = HuffTree.toList();
		long theData = 0;
		String aCode = null;
		
		for (int i = 0; i < treeList.size(); i++)
		{
			if (treeList.get(i).getData().hasData())
			{
				theData = treeList.get(i).getData().getData();
				aCode = treeList.get(i).getHuffCode();
				this.table.setHuffCode(theData, aCode);
			}
		}
		
	}
	
	private long highestFreq()
	{
		long max = 0;
		if (this.table == null) return 0;
		
		List<HuffTable.HuffPoint> hp = this.table.contentsToList(true);
		
		max = hp.get(hp.size() - 1).getFreq(); 
		//If sorted by frequency, should be the last entry.
		
		return max;
	}
	
	private int numberValidTblEntries()
	{
		return this.table.numberValidEntries();
	}
	
	private int sizeOfSymbolEntry()
	{
		int sz = 0;
		sz = this.bitDepth / 8;
		if (this.bitDepth % 8 != 0) sz++;
		
		return sz;
	}
	
	private int sizeOfFreqEntry()
	{
		int sz = 0;
		long topFreq = this.highestFreq();
		
		if (Math.abs(topFreq) < 256) sz = 1;
		else if (Math.abs(topFreq) < (1L << 16)) sz = 2;
		else if (Math.abs(topFreq) < (1L << 32)) sz = 4;
		else sz = 8;
		
		
		return sz;
	}
	
	private int tblEntrySize()
	{
		return this.sizeOfSymbolEntry() + this.sizeOfFreqEntry();
	}
	
	
	/* ~~~~~~~~~~~ Private Encode/ Decode Components ~~~~~~~~~~~ */
	private byte decTEMPBYTE;
	private int decTEMPBITPOS;

	private byte addBit(boolean bit, byte target, int position)
	{
		byte myByte = 0;
		int temp = Byte.toUnsignedInt(target);
		
		if (position < 0 || position > 7) return target;
		
		int mask = 1 << position;
		
		if (bit)
		{
			temp = temp | mask;
		}
		else
		{
			mask = ~mask;
			temp = temp & mask;
		}
		myByte = (byte)temp;
		
		
		return myByte;
	}
	
	private boolean getBit(byte target, int position)
	{
		if (position < 0 || position > 7) return false;
		
		int mask = 1 << position;
		int temp = Byte.toUnsignedInt(target);
		
		if ((temp & mask) != 0) return true;
		return false;
		
	}
	
	private boolean getBit(short target, int position)
	{
		if (position < 0 || position > 15) return false;
		
		int mask = 1 << position;
		int temp = Short.toUnsignedInt(target);
		
		if ((temp & mask) != 0) return true;
		return false;	
	}
	
	private void addSymbolToFileBuffer(FileBuffer decFile, long symbol)
	{
		if (this.bitDepth < 8) addSymU8(decFile, symbol);
		else if (this.bitDepth == 8) addSymE8(decFile, symbol);
		else if (this.bitDepth > 8 && this.bitDepth < 16) addSymU16(decFile, symbol);
		else if (this.bitDepth == 16) addSymE16(decFile, symbol);
		
	}
	
	private void addSymU8(FileBuffer decFile, long symbol)
	{
		byte aByte = (byte)symbol;
		for (int i = (this.bitDepth - 1); i >= 0; i--)
		{
			addBit(getBit(aByte, i), this.decTEMPBYTE, this.decTEMPBITPOS);
			this.decTEMPBITPOS--;
			if (this.decTEMPBITPOS < 0)
			{
				decFile.addToFile(this.decTEMPBYTE);
				this.decTEMPBITPOS = 7;
				this.decTEMPBYTE = 0;
			}
		}
	}
	
	private void addSymE8(FileBuffer decFile, long symbol)
	{
		byte aByte = (byte)symbol;
		decFile.addToFile(aByte);
	}
	
	private void addSymU16(FileBuffer decFile, long symbol)
	{
		short aVal = (short)symbol;
		for (int i = (this.bitDepth - 1); i >= 0; i--)
		{
			addBit(getBit(aVal, i), this.decTEMPBYTE, this.decTEMPBITPOS);
			this.decTEMPBITPOS--;
			if (this.decTEMPBITPOS < 0)
			{
				decFile.addToFile(this.decTEMPBYTE);
				this.decTEMPBITPOS = 7;
				this.decTEMPBYTE = 0;
			}
		}	
	}
	
	private void addSymE16(FileBuffer decFile, long symbol)
	{
		short aVal = (short)symbol;
		boolean dfbe = decFile.isBigEndian();
		if (!dfbe) decFile.setEndian(true);
		decFile.addToFile(aVal);
		if (!dfbe) decFile.setEndian(false);		
	}
	
	private void encodeData(FileBuffer inFile, FileBuffer encFile, long stPos, long edPos)
	{
		/*Switcher that determines which function to call.*/
		//System.err.println("Huffman.encodeData|| Called!");
		try
		{
			if (this.bitDepth <= 8) encodeDataU8(inFile, encFile, stPos, edPos);
			else if (this.bitDepth > 8 && this.bitDepth <= 16) encodeDataU16(inFile, encFile, stPos, edPos);
			else if (this.bitDepth > 16 && this.bitDepth <= 32) encodeDataU32(inFile, encFile, stPos, edPos);
			else if (this.bitDepth > 32 && this.bitDepth <= 64) encodeDataU64(inFile, encFile, stPos, edPos);	
		}
		catch (IncompatibleTableException e)
		{
			throw e;
		}
	}
	
	private void encodeDataU8(FileBuffer inFile, FileBuffer encFile, long stPos, long edPos)
	{
		//System.err.println("Huffman.encodeDataU8|| Called!");
		/*Handles encoding when bits to encode can fit in a byte*/
		
		/*1. Mark bit and byte positions*/
		/*2. Grab bits to encode*/
		
		long byPos = stPos;
		int biPos = 7;
		
		int temp = 0;
		byte tempByte = 0;
		int tempBits = 0; //How many bits in current temp byte
		
		while (byPos < edPos)
		{
			byte someBits = inFile.getBits8(this.bitDepth, byPos, biPos);
			String hCode = this.table.getHuffCode(Byte.toUnsignedLong(someBits));
			if (hCode.equals(""))
			{
				throw new IncompatibleTableException();
			}
			//if (byPos < 0x20) System.err.println("Huffman.encodeDataU8|| Bits Found at [" + byPos + "," + biPos + "]: " + String.format("%02x", someBits) + " (Huffman Code: " + hCode + ")");
			for (int i = 0; i < hCode.length(); i++)
			{
				temp = temp << 1;
				if (hCode.charAt(i) == '1')
				{
					temp = temp | 0x1;
				}
				tempBits++;
				
				if (tempBits >= 8)
				{
					tempBits = 0;
					tempByte = (byte)temp;
					temp = 0;
					encFile.addToFile(tempByte);
					//if (byPos < 0x20) System.err.println("Huffman.encodeDataU8|| Byte added to file: " + String.format("%02x", tempByte));
					tempByte = 0;
				}
			}
			
			biPos -= this.bitDepth;
			if (biPos < 0)
			{
				biPos = 8 + biPos;
				byPos++;
			}
		}
		
		/*Add EOF*/
		HuffNode EOFnode = createEOFnode();
		if (this.HuffTree.searchAndSetCurrent(EOFnode))
		{
			String EOFcode = this.HuffTree.currentLocation();
			if (EOFcode != null && EOFcode != "")
			{
				for (int i = 0; i < EOFcode.length(); i++)
				{
					temp = temp << 1;
					if (EOFcode.charAt(i) == '1')
					{
						temp = temp | 0x1;
					}
					tempBits++;
					
					if (tempBits >= 8)
					{
						tempBits = 0;
						tempByte = (byte)temp;
						temp = 0;
						encFile.addToFile(tempByte);
						tempByte = 0;
					}
				}
			}
		}
		
		/*Pad*/
		if (tempBits != 0)
		{
			while (tempBits < 8)
			{
				tempBits++;
				temp = temp << 1;
			}
			tempByte = (byte)temp;
			encFile.addToFile(tempByte);
		}
		
	}
	
	private void encodeDataU16(FileBuffer inFile, FileBuffer encFile, long stPos, long edPos)
	{
		/*Handles encoding when bits to encode cannot fit in a byte, but can fit in a short*/
		
		long byPos = stPos;
		int biPos = 15;
		
		int temp = 0;
		byte tempByte = 0;
		int tempBits = 0; //How many bits in current temp byte
		
		while (byPos < edPos)
		{
			short someBits = inFile.getBits16(this.bitDepth, byPos, biPos);
			String hCode = this.table.getHuffCode(Short.toUnsignedLong(someBits));
			if (hCode.equals(""))
			{
				throw new IncompatibleTableException();
			}
			
			for (int i = 0; i < hCode.length(); i++)
			{
				temp = temp << 1;
				if (hCode.charAt(i) == '1')
				{
					temp = temp | 0x1;
				}
				tempBits++;
				
				if (tempBits >= 8)
				{
					tempBits = 0;
					tempByte = (byte)temp;
					temp = 0;
					encFile.addToFile(tempByte);
					tempByte = 0;
				}
			}
			
			biPos -= this.bitDepth;
			if (biPos < 0)
			{
				biPos = 16 + biPos;
				byPos += 2;
			}
		}
		
		/*Add EOF*/
		HuffNode EOFnode = new HuffNode(1);
		EOFnode.makeNodeEOF();
		if (this.HuffTree.searchAndSetCurrent(EOFnode))
		{
			String EOFcode = this.HuffTree.currentLocation();
			if (EOFcode != null && EOFcode != "")
			{
				for (int i = 0; i < EOFcode.length(); i++)
				{
					temp = temp << 1;
					if (EOFcode.charAt(i) == '1')
					{
						temp = temp | 0x1;
					}
					tempBits++;
					
					if (tempBits >= 8)
					{
						tempBits = 0;
						tempByte = (byte)temp;
						temp = 0;
						encFile.addToFile(tempByte);
						tempByte = 0;
					}
				}
			}
		}
		
		/*Pad*/
		if (tempBits != 0)
		{
			while (tempBits < 8)
			{
				tempBits++;
				temp = temp << 1;
			}
			tempByte = (byte)temp;
			encFile.addToFile(tempByte);
		}	
	}
	
	private void encodeDataU32(FileBuffer inFile, FileBuffer encFile, long stPos, long edPos)
	{
		/*Handles encoding when bits to encode cannot fit in a short, but can fit in a long*/
		
		long byPos = stPos;
		int biPos = 31;
		
		int temp = 0;
		byte tempByte = 0;
		int tempBits = 0; //How many bits in current temp byte
		
		while (byPos < edPos)
		{
			int someBits = inFile.getBits32(this.bitDepth, byPos, biPos);
			String hCode = this.table.getHuffCode(Integer.toUnsignedLong(someBits));
			if (hCode.equals(""))
			{
				throw new IncompatibleTableException();
			}
			
			for (int i = 0; i < hCode.length(); i++)
			{
				temp = temp << 1;
				if (hCode.charAt(i) == '1')
				{
					temp = temp | 0x1;
				}
				tempBits++;
				
				if (tempBits >= 8)
				{
					tempBits = 0;
					tempByte = (byte)temp;
					temp = 0;
					encFile.addToFile(tempByte);
					tempByte = 0;
				}
			}
			
			biPos -= this.bitDepth;
			if (biPos < 0)
			{
				biPos = 32 + biPos;
				byPos += 4;
			}
		}
		
		/*Add EOF*/
		HuffNode EOFnode = new HuffNode(1);
		EOFnode.makeNodeEOF();
		if (this.HuffTree.searchAndSetCurrent(EOFnode))
		{
			String EOFcode = this.HuffTree.currentLocation();
			if (EOFcode != null && EOFcode != "")
			{
				for (int i = 0; i < EOFcode.length(); i++)
				{
					temp = temp << 1;
					if (EOFcode.charAt(i) == '1')
					{
						temp = temp | 0x1;
					}
					tempBits++;
					
					if (tempBits >= 8)
					{
						tempBits = 0;
						tempByte = (byte)temp;
						temp = 0;
						encFile.addToFile(tempByte);
						tempByte = 0;
					}
				}
			}
		}
		
		/*Pad*/
		if (tempBits != 0)
		{
			while (tempBits < 8)
			{
				tempBits++;
				temp = temp << 1;
			}
			tempByte = (byte)temp;
			encFile.addToFile(tempByte);
		}		
	}
	
	private void encodeDataU64(FileBuffer inFile, FileBuffer encFile, long stPos, long edPos)
	{
		/*Handles encoding when bits to encode cannot fit in a short, but can fit in a long*/
		
		long byPos = stPos;
		int biPos = 63;
		
		int temp = 0;
		byte tempByte = 0;
		int tempBits = 0; //How many bits in current temp byte
		
		while (byPos < edPos)
		{
			long someBits = inFile.getBits64(this.bitDepth, byPos, biPos);
			String hCode = this.table.getHuffCode(someBits);
			if (hCode.equals(""))
			{
				throw new IncompatibleTableException();
			}
			
			for (int i = 0; i < hCode.length(); i++)
			{
				temp = temp << 1;
				if (hCode.charAt(i) == '1')
				{
					temp = temp | 0x1;
				}
				tempBits++;
				
				if (tempBits >= 8)
				{
					tempBits = 0;
					tempByte = (byte)temp;
					temp = 0;
					encFile.addToFile(tempByte);
					tempByte = 0;
				}
			}
			
			biPos -= this.bitDepth;
			if (biPos < 0)
			{
				biPos = 64 + biPos;
				byPos += 8;
			}
		}
		
		/*Add EOF*/
		HuffNode EOFnode = new HuffNode(1);
		EOFnode.makeNodeEOF();
		if (this.HuffTree.searchAndSetCurrent(EOFnode))
		{
			String EOFcode = this.HuffTree.currentLocation();
			if (EOFcode != null && EOFcode != "")
			{
				for (int i = 0; i < EOFcode.length(); i++)
				{
					temp = temp << 1;
					if (EOFcode.charAt(i) == '1')
					{
						temp = temp | 0x1;
					}
					tempBits++;
					
					if (tempBits >= 8)
					{
						tempBits = 0;
						tempByte = (byte)temp;
						temp = 0;
						encFile.addToFile(tempByte);
						tempByte = 0;
					}
				}
			}
		}
		
		/*Pad*/
		if (tempBits != 0)
		{
			while (tempBits < 8)
			{
				tempBits++;
				temp = temp << 1;
			}
			tempByte = (byte)temp;
			encFile.addToFile(tempByte);
		}			
	}
	
	private String EOFHuffCode()
	{
		/*Gets EOF marker code*/
		if (this.HuffTree == null) return "";
		if (this.table == null) return "";
		
		HuffNode EOF = new HuffNode(1);
		EOF.makeNodeEOF();
		
		return this.HuffTree.searchAndGetLocation(EOF);
	}
	
	/* ~~~~~~~~~~~ Constructor ~~~~~~~~~~~ */
 	/**
 	 * Constructor for the Huffman object - generates a tree from the frequency
 	 * table given.
 	 * @param fTable | The frequency table used to generate the tree.
 	 */
	public Huffman(HuffTable fTable)
	{
		/*Constructor immediately constructs tree.*/
		this.bitDepth = fTable.getBitDepth();
		this.table = fTable;
		this.constructTree();
	}
	
	/* ~~~~~~~~~~~ Getters/ Setters ~~~~~~~~~~~ */
	public HuffTable getTable()
	{
		return this.table;
	}
	
 	/* ~~~~~~~~~~~ Basic Encoders/ Decoders ~~~~~~~~~~~ */
	/**Uses the tree contained in this instance to encode the file
	 * and returns the encoded version
	 * @param myFile - FileBuffer to encode
	 * @param includeTable - Whether or not to include table serialization at beginning
	 * @return new FileBuffer that contains an encoded version of myFile
	 * @throws IOException */
	public FileBuffer encodeHuff(FileBuffer myFile, boolean includeTable) throws IOException
	{	
		if (myFile == null) return null;
		if (myFile.isEmpty()) return null;
		if (this.HuffTree == null) return null;
		if (this.table == null) return null;
		
		return this.encodeHuff(myFile, includeTable, 0, myFile.getFileSize());
	}
	
	/**Uses the tree contained in this instance to encode the file
	 * and returns the encoded version
	 * @param myFile - FileBuffer to encode
	 * @param includeTable - Whether or not to include table serialization at beginning
	 * @param stOff - the position of the first byte in myFile to encode.
	 * @return new FileBuffer that contains an encoded version of myFile
	 * @throws IOException */
	public FileBuffer encodeHuff(FileBuffer myFile, boolean includeTable, long stOff) throws IOException
	{	
		if (myFile == null) return null;
		if (myFile.isEmpty()) return null;
		if (this.HuffTree == null) return null;
		if (this.table == null) return null;
		if (stOff < 0 || stOff >= myFile.getFileSize()) return null;
		
		return this.encodeHuff(myFile, includeTable, stOff, myFile.getFileSize());
	}
	
	/**Uses the tree contained in this instance to encode the file
	 * and returns the encoded version
	 * @param myFile - FileBuffer to encode
	 * @param includeTable - Whether or not to include table serialization at beginning\
	 * @param stOff - the position of the first byte in myFile to encode
	 * @param edOff - the position in myFile to stop encoding
	 * @return new FileBuffer that contains an encoded version of myFile
	 * @throws IOException */
	public FileBuffer encodeHuff(FileBuffer myFile, boolean includeTable, long stOff, long edOff) throws IOException
	{
		return encodeHuff(myFile, includeTable, stOff, edOff, null);
	}
	
	/**Uses the tree contained in this instance to encode the file
	 * and returns the encoded version. This one can also add a string at the beginning.
	 * @param myFile - FileBuffer to encode
	 * @param includeTable - Whether or not to include table serialization at beginning\
	 * @param stOff - the position of the first byte in myFile to encode
	 * @param edOff - the position in myFile to stop encoding
	 * @param head - A string to add to the start of the file (if wanted for tagging)
	 * @return new FileBuffer that contains an encoded version of myFile
	 * @throws IOException */
	public FileBuffer encodeHuff(FileBuffer myFile, boolean includeTable, long stOff, long edOff, String head) throws IOException
	{
		//System.err.println("Huffman.encodeHuff || Entering header add overload!");
		//System.err.println("Huffman.encodeHuff(FileBuffer,bool,long,long,String) || Called!");
		/*Args check*/
		if (myFile == null) return null;
		if (myFile.isEmpty()) return null;
		if (this.HuffTree == null) return null;
		if (this.table == null) return null;
		if (edOff > myFile.getFileSize() || edOff < 0) return null;
		if (stOff < 0 || stOff >= edOff || stOff >= myFile.getFileSize()) return null;
		//System.err.println("Huffman.encodeHuff(FileBuffer,bool,long,long,String) || Args check passed!");
		
		/*Set Huff Codes*/
		this.setHuffCodesInTable();
		
		/*Set up buffer to dump info in*/
		long rSize = (edOff - stOff) + (edOff - stOff)/2;
		FileBuffer encFile = FileBuffer.createWritableBuffer(myFile.getName() + ".huff", rSize, true);
		
			/*DEBUG CHECK*/
			//System.out.println("Serialized Table Size: " + this.serializedTableSize(true));
			//System.out.println("Initial Encoded File Buffer Capacity: " + initBufferCap);
		
		/*Add tag if specified*/
		if (head != null)
		{
			//System.err.println("Huffman.encodeHuff || head = " + head);
			encFile.printASCIIToFile(head);
		}
		
		
		/*Add serialized table*/
		if (includeTable)
		{
			FileBuffer serialTable = serializeTable(edOff - stOff);
			encFile.addToFile(serialTable, 0, serialTable.getFileSize());
			/*DEBUG CHECK*/
			//System.out.println("Serialized Table: " + "\n" + serialTable.toString());
			//System.out.println("");
			//System.out.println("Checking File Addition Success... encFile: \n" + encFile.toString());
			//String DEBUGTABLEPATH = DEBUGDIR + "hufftableencode.bin";
			//System.err.println("Huffman.encodeHuff(FileBuffer,bool,long,long,String) || Writing table alone to " + DEBUGTABLEPATH);
			//serialTable.writeFile(DEBUGTABLEPATH);
			
		}
		
		try
		{
			this.encodeData(myFile, encFile, stOff, edOff);
		}
		catch (IncompatibleTableException e)
		{
			System.out.println("Table is incompatible with file! (File contains entries not in table)");
			return null;
		}
		
		//System.err.println("Huffman.encodeHuff || Leaving overload!");
		return encFile;		
	}

	public FileBuffer encodeHuff(FileBuffer myFile, String headerAdd) throws IOException
	{
		//System.err.println("Huffman.encodeHuff(FileBuffer,String) || Called!");
		return this.encodeHuff(myFile, true, 0, myFile.getFileSize(), headerAdd);
	}
	
	/**Uses the tree contained in this instance to decode the file
	 * and returns the decoded version
	 * @param myFile - FileBuffer to decode
	 * @param stOff - the position of the first byte in myFile to encode
	 * @return new FileBuffer that contains an encoded version of myFile
	 * @throws IOException */
	public FileBuffer decodeHuff(FileBuffer myFile, long stOff) throws IOException
	{	
		//System.err.println("Huffman.decodeHuff(FileBuffer,long) || Called!");
		if (myFile == null) return null;
		if (myFile.isEmpty()) return null;
		if (this.HuffTree == null) return null;
		if (this.table == null) return null;
		if (stOff < 0 || stOff >= myFile.getFileSize()) return null;	
		
		//System.err.println("Huffman.decodeHuff(FileBuffer,long) || args check passed!");
		return this.decodeHuff(myFile, stOff, myFile.getFileSize());
	}
	
	/**Uses the tree contained in this instance to decode the file
	 * and returns the decoded version
	 * @param myFile - FileBuffer to decode
	 * @param stOff - the position of the first byte in myFile to encode
	 * @param edOff - the position in myFile to stop encoding
	 * @return new FileBuffer that contains an encoded version of myFile
	 * @throws IOException */
	public FileBuffer decodeHuff(FileBuffer myFile, long stOff, long edOff) throws IOException
	{	
		//System.err.println("Huffman.decodeHuff(FileBuffer,long,long) || Called!");
		if (myFile == null) return null;
		if (myFile.isEmpty()) return null;
		if (this.HuffTree == null) return null;
		if (this.table == null) return null;
		if (edOff > myFile.getFileSize() || edOff < 0) return null;
		if (stOff < 0 || stOff >= edOff || stOff >= myFile.getFileSize()) return null;
		//System.err.println("Huffman.decodeHuff(FileBuffer,long,long) || Args check passed!");
		//System.err.println("Huffman.decodeHuff(FileBuffer,long,long) || Start position: 0x" + Long.toHexString(stOff));
		
		long rSize = (edOff - stOff) + (edOff - stOff)/2;
		FileBuffer decFile = FileBuffer.createWritableBuffer(myFile.getName() + ".huffdec", rSize, true);
		
		/*Get a byte from the file*/
		/*Move around tree until we get a symbol*/
		/*Determine size of symbol and encode bits into the appropriate amount of bytes*/
		
		//System.err.println("Huffman.decodeHuff(FileBuffer,long,long) || TABLE ----");
		//this.printTableToSTDERR();
		
		this.HuffTree.moveToRoot();
		byte eByte = 0;
		long sym = 0;
		this.decTEMPBYTE = 0;
		this.decTEMPBITPOS = 7;
		
		//String codeDEBUG = "";
		
		boolean EOFhit = false;
		
		for (long i = stOff; i < edOff; i++)
		{
			eByte = myFile.getByte(i);
			//if (i - stOff < 0x20) System.err.println("Huffman.decodeHuff(FileBuffer,long,long) || Byte found: " + String.format("%02x", eByte));
			for (int b = 7; b >= 0; b--)
			{
				if (this.getBit(eByte, b))
				{
					//codeDEBUG += "1";
					if (!this.HuffTree.moveRight())
					{
						//Get value
						if (this.HuffTree.getCurrentData().isEOF())
						{
							EOFhit = true;
							break;
						}
						sym = this.HuffTree.getCurrentData().getData();
						if (this.HuffTree.getCurrentData().hasData()) {
							//if (i - stOff < 0x20) System.err.println("Huffman.decodeHuff(FileBuffer,long,long) || Adding Symbol: " + String.format("%02x", sym) + " (" + codeDEBUG + ")");
							this.addSymbolToFileBuffer(decFile, sym);
							//codeDEBUG = "";
						}
						//Return to root and move right
						this.HuffTree.moveToRoot();
						this.HuffTree.moveRight();
					}
				}
				else
				{
					if (!this.HuffTree.moveLeft())
					{
						//codeDEBUG += "0";
						//Get value
						if (this.HuffTree.getCurrentData().isEOF())
						{
							EOFhit = true;
							break;
						}
						sym = this.HuffTree.getCurrentData().getData();
						if (this.HuffTree.getCurrentData().hasData()) {
							//if (i - stOff < 0x20) System.err.println("Huffman.decodeHuff(FileBuffer,long,long) || Adding Symbol: " + String.format("%02x", sym) + " (" + codeDEBUG + ")");
							this.addSymbolToFileBuffer(decFile, sym);
							//codeDEBUG = "";
						}
						//Return to root and move left
						this.HuffTree.moveToRoot();
						this.HuffTree.moveLeft();
					}		
				}
			}
			if (EOFhit) break;
		}
		
		/*Terminate if at risk of exceeding original file size so that padding is not
		 * read as data.*/
		/*Check the current on Hufftree - if it has symbol, then add.*/
		
		sym = this.HuffTree.getCurrentData().getData();
		if (this.HuffTree.getCurrentData().hasData() && !this.HuffTree.getCurrentData().isEOF()) this.addSymbolToFileBuffer(decFile, sym);
		this.HuffTree.moveToRoot();
		
		return decFile;	
	}

	
	/* ~~~~~~~~~~~ Table Serialization ~~~~~~~~~~~ */
	/**
	 * Returns the number of bytes it would take to serialize the table contained
	 * in the object.
	 * @param includeFSize | Whether or not to include file size in the byte string.
	 * @return Size of serialized table
	 */
	public int serializedTableSize(boolean includeFSize)
	{
		int tot = 0;
		tot += 12; /*Standard*/
		
		if (includeFSize) tot += 4;
		
		int nEntries = this.table.numberValidEntries();
		
		int entrySize = this.tblEntrySize();
		
		tot += (nEntries * entrySize);
		
		return tot;
	}
	
	/**
	 * Converts the table object contained in the Huffman to a serialized byte string.
	 * @param includeFSize | Whether or not to include table's stored file size in serialization.
	 * @return FileBuffer style string of bytes representing the table.
	 */
	public FileBuffer serializeTable(boolean includeFSize)
	{
		return serializeTable(includeFSize, this.table.getFileSize());
	}
	
	/**
	 * Converts the table object contained in the Huffman to a serialized byte string.
	 * @param fileSize | The file size to include in the serialized data (overriding stored size in table object)
	 * @return FileBuffer style string of bytes representing the table.
	 */
	public FileBuffer serializeTable(long fileSize)
	{
		//System.err.println("Huffman.serializedTable(long) || Called!");
		return serializeTable(true, fileSize);
	}
	
	private FileBuffer serializeTable(boolean includeFSize, long fileSize)
	{
		//System.err.println("Huffman.serializeTable(bool,long) || Called!");
		byte hasSize = 0;
		byte bDepth = (byte)this.bitDepth;
		byte fRecBits = 0;
		byte reserved = 0x7F;
		//System.err.println("Huffman.serializeTable(bool,long) || Bit Depth: " + bDepth);
		
		if (includeFSize && fileSize <= 0xFFFFFFFFL) hasSize = 4;
		else if (includeFSize && fileSize > 0xFFFFFFFFL) hasSize = 8;
		//System.err.println("Huffman.serializeTable(bool,long) || Size of file size record: " + fileSize);
		
		int nEntries = this.numberValidTblEntries();
		
		int sSize = this.serializedTableSize(true);
		int fEntrySize = this.sizeOfFreqEntry();
		int sEntrySize = this.sizeOfSymbolEntry();
		
		fRecBits = (byte)(fEntrySize * 8);
		int datarecSize = sEntrySize * 8;
		
		//System.err.println("Huffman.serializeTable(bool,long) || Number of entries: " + nEntries);
		//System.err.println("Huffman.serializeTable(bool,long) || Serialized size: 0x" + Integer.toHexString(sSize));
		//System.err.println("Huffman.serializeTable(bool,long) || Freq Record Size: " + fEntrySize);
		//System.err.println("Huffman.serializeTable(bool,long) || Symbol Record Size: " + sEntrySize);
		
		FileBuffer myTable = new FileBuffer(sSize, false);

		myTable.addToFile((byte)0x14);
		myTable.addToFile((byte)0x02);
		myTable.addToFile((byte)0x0C);
		myTable.addToFile((byte)0x05);
		myTable.addToFile(hasSize);
		myTable.addToFile(bDepth);
		myTable.addToFile(fRecBits);
		myTable.addToFile(reserved);
		myTable.addToFile(nEntries);
		
		if (includeFSize)
		{
			if (hasSize == 4) myTable.addToFile((int)fileSize);
			else if (hasSize == 8) myTable.addToFile(fileSize);	
		}
		
		for (HuffTable.HuffPoint e:this.table.contentsToList(false))
		{
			if (e.getFreq() > 0)
			{
				switch (datarecSize)
				{
				case 8:
					myTable.addToFile((byte)e.getSymbol());
					break;
				case 16:
					myTable.addToFile((short)e.getSymbol());
					break;
				case 32:
					myTable.addToFile((int)e.getSymbol());
					break;
				case 64:
					myTable.addToFile(e.getSymbol());
					break;
				default:
					break;
				}
				
				switch (fRecBits)
				{
				case 8:
					myTable.addToFile((byte)e.getFreq());
					break;
				case 16:
					myTable.addToFile((short)e.getFreq());
					break;
				case 32:
					myTable.addToFile((int)e.getFreq());
					break;
				case 64:
					myTable.addToFile(e.getFreq());
					break;
				default:
					break;
				}
			}
		}
		
		return myTable;
	}

	/**
	 * Parses a serialized HuffTable byte string stored in a file that has been
	 * written by this library.
	 * @param myFile | File to read table from
	 * @param stOff | Offset to start looking for table
	 * @return HuffTable object created from serialized data
	 */
	public static HuffTable readHuffTable(FileBuffer myFile, long stOff)
	{
		//System.err.println("Huffman.readHuffTable || Called!");
		/*Parses a hufftable serialized in a file.*/
		if (stOff < 0 || stOff >= myFile.getFileSize()) return null;
		//System.out.println("DEBUG | Table read args check passed.");
		
		/*1. Look for the first tble magic past the start offset*/
		byte[] magic = {0x14, 0x02, 0x0C, 0x05};
		long tStart = myFile.findString(stOff, myFile.getFileSize(), magic);
		//System.err.println("Huffman.readHuffTable || Table magic found at: 0x" + Long.toHexString(tStart));
		//System.out.println("DEBUG | Magic found at offset: " + Integer.toHexString(tStart));
		
		if (tStart < 0) return null;
		
		/*2. Read 4 info bytes*/
		boolean hasSize = false;
		byte sizeBytes = myFile.getByte(tStart + 4);
		if (sizeBytes > 0) hasSize = true;
		//System.err.println("Huffman.readHuffTable || Size of size record: " + sizeBytes);
		//System.out.println("DEBUG | Size of size record: " + sizeBytes);
		
		byte bitD = myFile.getByte(tStart + 5);
		if (bitD < 2 || bitD > 16) return null;
		//System.err.println("Huffman.readHuffTable || Bit Depth: " + bitD);
		//System.out.println("DEBUG | Bit depth recorded: " + bitD);
		
		byte freqDepth = myFile.getByte(tStart + 6);
		//System.err.println("Huffman.readHuffTable || Size of frequency records: " + freqDepth + " bits");
		//System.out.println("DEBUG | Size of frequency records: " + freqDepth);
		
		/*3. Read number entries*/
		
		boolean fileEnd = myFile.isBigEndian();
		if (fileEnd) myFile.setEndian(false);
		
		int nEntries = myFile.intFromFile(tStart + 8);
		//System.err.println("Huffman.readHuffTable || Number of entries: " + nEntries);
		//System.out.println("DEBUG | Number of entries: " + nEntries);
		
		/*4. Read size, if applicable.*/
		long fSize = 0;
		if (hasSize)
		{
			if (sizeBytes == 4) fSize = myFile.intFromFile(tStart + 12);
			else if (sizeBytes == 8) fSize = myFile.longFromFile(tStart + 12);
			else return null;
			//System.err.println("Huffman.readHuffTable || Recorded file size: " + fSize);
			//System.out.println("DEBUG | Recorded file size: " + fSize);
		}
		
		/*5. Read entries*/
		long eStart = tStart + 12 + sizeBytes;
		//System.err.println("Huffman.readHuffTable || Entry offset: 0x" + Long.toHexString(eStart));
		//System.out.println("DEBUG | Entry offset: " + Integer.toHexString(eStart));
		HuffTable myTable = new HuffTable(bitD);
		if (hasSize) myTable.setFileSize((int)fSize);
		
		int entrySize = (bitD + freqDepth) / 8;
		//System.out.println("DEBUG | Entry Size: " + entrySize);
		//System.err.println("Huffman.readHuffTable || Entry Size: " + entrySize + " bytes");
		int symbol = 0;
		long frequency = 0;
		int bytes = bitD / 8;
		if (bitD % 8 != 0)
		{
			bytes++;
			entrySize++;
		}
		
		for (int j = 0; j < nEntries; j++)
		{
			long k = eStart + (j * entrySize);
			if (bitD <= 8)
			{
				symbol = Byte.toUnsignedInt(myFile.getByte(k));
			}
			else if (bitD > 8 && bitD <= 16)
			{
				symbol = Short.toUnsignedInt(myFile.shortFromFile(k));
			}
			
			if (freqDepth <= 8)
			{
				frequency = Byte.toUnsignedLong(myFile.getByte(k + bytes));
			}
			else if (freqDepth > 8 && freqDepth <= 16)
			{
				frequency = Short.toUnsignedLong(myFile.shortFromFile(k + bytes));
			}
			else if (freqDepth > 16 && freqDepth <= 32)
			{
				frequency = Integer.toUnsignedLong(myFile.intFromFile(k + bytes));
			}
			else if (freqDepth > 32 && freqDepth <= 64)
			{
				frequency = myFile.longFromFile(k + bytes);
			}
			
			myTable.setFrequency(symbol, frequency);
			//System.err.println("Huffman.readHuffTable || Symbol #" + j + ": " + String.format("%02x", symbol));
			//System.err.println("Huffman.readHuffTable || Frequency #" + j + ": " + frequency);
			//System.out.println("DEBUG | Symbol #" + j + ": " + Integer.toHexString(symbol));
			//System.out.println("DEBUG | Frequency #" + j + ": " + frequency);
		}
		
		myFile.setEndian(fileEnd);
		
		//System.err.println("Huffman.readHuffTable || Function returning...");
		return myTable;
	}
	
	public void printTableToSTDERR()
	{
		setHuffCodesInTable();
		List<HuffPoint> tbl = this.table.contentsToList(false);
		for (HuffPoint p : tbl)
		{
			System.err.println(Long.toHexString(p.getSymbol()) + "\t" + p.getFreq() + "\t" + p.getHuffCode());
		}
	}
	
	/* ~~~~~~~~~~~ De Novo Table Generation ~~~~~~~~~~~ */
	/**
	 * Generates a new HuffTable from the frequencies of bit chunks (of width
	 * specified) in given file.
	 * @param myFile | File to generate table from.
	 * @param bitDepth | Bits to analyze at once
	 * @return HuffTable object generated
	 */
	public static HuffTable freqTableFromFile(FileBuffer myFile, int bitDepth)
	{
		/*Generates a table of frequencies of each bit set (usually byte)*/
		HuffTable myTable = new HuffTable(bitDepth);
		myTable.setFileSize(myFile.getFileSize());
		
		int byPos = 0;
		int biPos = 0;
		
		if (bitDepth < 2 || bitDepth > 16) return null;
		
		if (bitDepth >= 2 && bitDepth < 8)
		{
			//Unusual bit number > less than a byte
			byte aByte = 0;
			
			while (byPos < myFile.getFileSize())
			{
				aByte = myFile.getBits8(bitDepth, byPos, biPos);
				myTable.incrementFrequency(Byte.toUnsignedLong(aByte));
				biPos += bitDepth;
				if (biPos > 7)
				{
					biPos = biPos % 8;
					byPos++;
				}
			}
			
		}
		else if (bitDepth == 8)
		{
			//Even byte
			byte aByte = 0;
			for (int i = 0; i < myFile.getFileSize(); i++)
			{
				aByte = myFile.getByte(i);
				myTable.incrementFrequency(Byte.toUnsignedLong(aByte));
			}
		}
		else if (bitDepth > 8 && bitDepth < 16)
		{
			//Unusual bit number, between byte and short
			short aShort = 0;
			while (byPos < myFile.getFileSize())
			{
				aShort = myFile.getBits16(bitDepth, byPos, biPos);
				myTable.incrementFrequency(Short.toUnsignedLong(aShort));
				biPos += bitDepth;
				if (biPos > 15)
				{
					biPos = biPos % 16;
					byPos += 2;
				}
			}
		}
		else if (bitDepth == 16)
		{
			//Even short
			short aShort = 0;
			for (int i = 0; i < myFile.getFileSize(); i+=2)
			{
				aShort = myFile.shortFromFile(i);
				myTable.incrementFrequency(Short.toUnsignedLong(aShort));
			}
		}
		
		return myTable;
	}
	
	/* ~~~~~~~~~~~ Error Handling ~~~~~~~~~~~ */
	public static class IncompatibleTableException extends RuntimeException
	{

		/**
		 * 
		 */
		private static final long serialVersionUID = -2554826591558733363L;
		
	}
	
	/* ~~~~~~~~~~~ Printers ~~~~~~~~~~~ */
	/**
	 * Prints frequency table contained in object to console.
	 */
	public void printTable()
	{
		if (this.table == null)
		{
			System.out.println("This Huffman object does not contain a table!");
			return;
		}
		
		List<HuffPoint> sortedList = this.table.contentsToList(false);
		System.out.println("Symbol Code\t\tFrequency\t\tHuff Code");
		System.out.println("-----------\t\t---------\t\t---------");
		
		for (HuffPoint p:sortedList)
		{
			System.out.print(Long.toHexString(p.getSymbol()) + "\t\t");
			System.out.print(Long.toString(p.getFreq()) + "\t\t");
			System.out.print(p.getHuffCode() + "\n");
		}
		
	}
	
	/**
	 * Prints frequency table contained in object to specified path.
	 * @param printPath | Path to write text to
	 */
	public void printTableToFile(String printPath)
	{
		if (this.table == null)
		{
			System.out.println("This Huffman object does not contain a table!");
			return;
		}
		
		int estLen = 0;
		/*Estimated line length is 90.*/
		int estLines = this.numberValidTblEntries() + 4;
		estLen = estLines * 90;
		
		FileBuffer printTbl = new FileBuffer(estLen);
		printTbl.printASCIIToFile("Symbol Code\t\tFrequency\t\tHuff Code\n");
		printTbl.printASCIIToFile("-----------\t\t---------\t\t---------\n");
		
		List<HuffPoint> eList = this.table.contentsToList(false);
		for (HuffPoint p : eList)
		{
			printTbl.printASCIIToFile(Long.toHexString(p.getSymbol()) + "\t\t");
			printTbl.printASCIIToFile(Long.toString(p.getFreq()) + "\t\t");
			printTbl.printASCIIToFile(p.getHuffCode() + "\n");
		}
		printTbl.printASCIIToFile("EOF Marker: ");
		printTbl.printASCIIToFile(this.EOFHuffCode() + "\n");
		
		printTbl.printASCIIToFile("\nHuffman Tree: \n");
		printTbl.printASCIIToFile(this.HuffTree.toString() + "\n");
		
		try
		{
			printTbl.writeFile(printPath);
		}
		catch (IOException e)
		{
			System.out.println("Error writing to file \"" + printPath + "\" !");
			printTable();
		}
		
		
	}
	
	/* ~~~~~~~~~~~ Quick (Static) Encode/ Decode ~~~~~~~~~~~ */
	
	/**
	 * Encodes the given file using the default table settings.
	 * Generates the frequency table, serializes it, and sticks it at the beginning
	 * of the file.
	 * @param myFile: File to encode
	 * @param bDepth: Bit depth of Huffman encoding
	 * @param stPos: First byte to encode
	 * @param edPos: First byte not encoded (position after final byte)
	 * @return Encoded file.
	 * @throws IOException 
	 */
	public static FileBuffer HuffEncodeFile(FileBuffer myFile, int bDepth, long stPos, long edPos) throws IOException
	{
		return HuffEncodeFile(myFile, bDepth, stPos, edPos, null);
	}
	
	/**
	 * Encodes the given file using the default table settings.
	 * Generates the frequency table, serializes it, and sticks it at the beginning
	 * of the file.
	 * This overload also prints the frequency table to an output text file.
	 * @param myFile | File to encode
	 * @param bDepth | Bit depth
	 * @param stPos | First byte in file to encode
	 * @param edPos | Offset to encode to
	 * @param fTablePath | Path to put frequency table text file
	 * @return FileBuffer containing encoded file with table at the beginning
	 * @throws IOException 
	 */
	public static FileBuffer HuffEncodeFile(FileBuffer myFile, int bDepth, long stPos, long edPos, String fTablePath) throws IOException
	{
		/*Does the same as its overloaded partner, except it prints the frequency table to a file.*/
		HuffTable fTable = freqTableFromFile(myFile, bDepth);
		//System.out.println("DEBUG || General Encoder - Table Generated");
		Huffman huff = new Huffman(fTable);
		//System.out.println("DEBUG || General Encoder - Huffman Generated");
		FileBuffer enc = huff.encodeHuff(myFile, true, stPos, edPos);
		//System.out.println("DEBUG || General Encoder - Encode Attempted");
		if (fTablePath != null) huff.printTableToFile(fTablePath);
		//System.out.println("DEBUG || General Encoder - Table Printed");
		
		return enc;
	}
	
	public static FileBuffer HuffEncodeFile(FileBuffer myFile, int bDepth) throws IOException
	{
		return HuffEncodeFile(myFile, bDepth, 0, myFile.getFileSize(), null);
	}
	
	public static FileBuffer HuffEncodeFile(FileBuffer myFile, int bDepth, Path fTablePath) throws IOException
	{
		return HuffEncodeFile(myFile, bDepth, 0, myFile.getFileSize(), fTablePath.toString());
	}
	
	public static FileBuffer HuffEncodeFile(FileBuffer myFile, int bDepth, String header) throws IOException
	{
		//System.err.println("Huffman.HuffEncodeFile(FileBuffer,int,String) || Called");
		HuffTable fTable = freqTableFromFile(myFile, bDepth);
		//System.err.println("Huffman.HuffEncodeFile(FileBuffer,int,String) || Frequency table generated.");
		Huffman huff = new Huffman(fTable);
		//System.err.println("Huffman.HuffEncodeFile(FileBuffer,int,String) || Huffman tree generated.");
		FileBuffer enc = huff.encodeHuff(myFile, header);
		//System.err.println("Huffman.HuffEncodeFile(FileBuffer,int,String) || Encoding complete.");
		
		return enc;
	}
	
	/**
	 * Decodes the given file assuming it has been encoded by the above quick
	 * encoder.
	 * THIS METHOD WILL RETURN NULL IF IT DOES NOT FIND A PROPERLY FORMATTED FREQ
	 * TABLE AT THE BEGINNING OF THE FILE BUFFER OR IT CANNOT READ IT.
	 * @param myHuff: File to decode
	 * @return Decoded file
	 * @throws IOException 
	 */
	public static FileBuffer HuffDecodeFile(FileBuffer myHuff) throws IOException
	{
		return HuffDecodeFile(myHuff, 0);
	}
	
	/**
	 * Decodes the given file assuming it has been encoded by the above quick
	 * encoder.
	 * THIS METHOD WILL RETURN NULL IF IT DOES NOT FIND A PROPERLY FORMATTED FREQ
	 * TABLE AT THE BEGINNING OF THE FILE BUFFER OR IT CANNOT READ IT.
	 * @param myHuff: File to decode
	 * @return Decoded file
	 * @throws IOException 
	 */
	public static FileBuffer HuffDecodeFile(FileBuffer myHuff, long stPos) throws IOException
	{
		//System.err.println("Huffman.HuffEDecodeFile(FileBuffer,long) || Called!");
		HuffTable fTable = readHuffTable(myHuff, stPos);
		//System.err.println("Huffman.HuffEDecodeFile(FileBuffer,long) || Table read.");
		if (fTable == null) throw new IncompatibleTableException();
		//System.err.println("Huffman.HuffEDecodeFile(FileBuffer,long) || stPos = 0x" + Long.toHexString(stPos));
		
		long fEnd = 0;
		
		Huffman huff = new Huffman(fTable);
		if (huff.getTable().getFileSize() > 0) fEnd = huff.serializedTableSize(true);
		else fEnd = huff.serializedTableSize(false);
		//System.err.println("Huffman.HuffEDecodeFile(FileBuffer,long) || Frequency Table End: 0x" + Long.toHexString(fEnd));
		
		FileBuffer dec  = huff.decodeHuff(myHuff, fEnd + stPos);
		//System.err.println("Huffman.HuffEDecodeFile(FileBuffer,long) || Decode complete!");
		
		return dec;
	}

	
	/* ~~~~~~~~~~~ Utility Overrides ~~~~~~~~~~~ */
	public String toString()
	{
		String s = "";
		
		s += "Huffman Object: \n";
		s += "TABLE\n";
		s += this.table.toString();
		s += "TREE\n";
		s += this.HuffTree.toString();
		
		return s;
	}
	
	public int hashCode()
	{
		return this.toString().hashCode();
	}
	
	public boolean equals(Object other)
	{
		if (other == null) return false;
		if (!(other instanceof Huffman)) return false;
		if (this == other) return true;
		
		Huffman h = (Huffman)other;
		
		return this.toString().equals(h.toString());
	}
}
