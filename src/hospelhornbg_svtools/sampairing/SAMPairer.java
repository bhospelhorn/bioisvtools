package hospelhornbg_svtools.sampairing;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import hospelhornbg_bioinformatics.SAMRecord;
import hospelhornbg_bioinformatics.SAMRecord.InvalidSAMRecordException;
import hospelhornbg_genomeBuild.GenomeBuild;
import waffleoRai_Utils.FileBuffer;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;
import waffleoRai_Utils.Huffman;

public class SAMPairer {
	
	public static final int RECORDS_PER_FILE = 10000;

	private BufferedWriter[] openStreams;
	private int[] fileCounts;
	private int[] recordCounts;
	
	private String tempdir;
	
	public SAMPairer(String tempDirPath)
	{
		tempdir = tempDirPath;
		openStreams = new BufferedWriter[256];
		fileCounts = new int[256];
		recordCounts = new int[256];
	}
	
	private String getFilePath(int group, int fileNumber) throws IOException
	{
		String gnum = "g_" + String.format("%02X", group);
		String dir = tempdir + File.separator + gnum;
		if (!FileBuffer.directoryExists(dir))
		{
			Files.createDirectories(Paths.get(dir));
		}
		return dir + File.separator + gnum + "_" + fileNumber + ".tmp";
	}
	
	public synchronized void saveRecord(SAMRecord r) throws IOException
	{
		//Hash the qname
		int qhash = r.getQueryName().hashCode();
		
		//LSB of the hash determines the group
		int group = qhash & 0xFF;
		
		//See if there is already an open stream
		BufferedWriter bw = openStreams[group];
		
		//If not, open one
		if (bw == null)
		{
			String fpath = getFilePath(group, fileCounts[group]);
			FileWriter fw = new FileWriter(fpath);
			bw = new BufferedWriter(fw);
			openStreams[group] = bw;
		}
		
		//Write record
		bw.write(r.writeSAMRecord(false) + "\n");
		recordCounts[group]++;
		
		//See if record count is at max for this file
		if (recordCounts[group] >= RECORDS_PER_FILE)
		{
			//Close stream
			bw.close();
			openStreams[group] = null;
			//Compress result
			String fpath = getFilePath(group, fileCounts[group]);
			compressFile(fpath);
			//Increment file count
			fileCounts[group]++;
			recordCounts[group] = 0;
		}
		
	}
	
	private void compressFile(String path) throws IOException
	{
		String huffpath = path + ".huff";
		FileBuffer in = FileBuffer.createBuffer(path);
		FileBuffer out = Huffman.HuffEncodeFile(in, 8);
		out.writeFile(huffpath);
		Files.deleteIfExists(Paths.get(path));
	}
	
	private void decompressFile(String inpath, String outpath) throws IOException
	{
		FileBuffer in = FileBuffer.createBuffer(inpath);
		FileBuffer out = Huffman.HuffDecodeFile(in);
		out.writeFile(outpath);
		Files.deleteIfExists(Paths.get(inpath));
	}
	
	public synchronized void closeOpenStreams() throws IOException
	{
		for (int i = 0; i < 256; i++)
		{
			if (openStreams[i] != null)
			{
				openStreams[i].close();
				openStreams[i] = null;
				//Compress
				String path = getFilePath(i, fileCounts[i]);
				compressFile(path);
				fileCounts[i]++;
				recordCounts[i] = 0;
			}
		}
	}
	
	public String sortAllReads(int group, GenomeBuild gb) throws IOException
	{
		//Will close all open streams
		//Returns the path of the new dir with the sorted reads
		closeOpenStreams();
		
		//Read in records, sort into new files based on hash MSB
		//Don't need to parse whole record. Just read the first field to hash
		//Delete old files as they are processed
		String gnum = "g_" + String.format("%02X", group);
		String dir = tempdir + File.separator + gnum;
		
		int fcount = fileCounts[group];
		for (int i = 0; i < fcount; i++)
		{
			String inpath = dir + File.separator + gnum + "_" + i + ".tmp";
			//Decompress
			decompressFile(inpath + ".huff", inpath);
			//Read in
			FileReader fr = new FileReader(inpath);
			BufferedReader br = new BufferedReader(fr);
			
			String line = null;
			while((line = br.readLine()) != null)
			{
				String[] fields = line.split("\t");
				if (fields.length >= 1)
				{
					int sg = fields[0].hashCode() >>> 24;
					String outpath = dir + File.separator + "sg_" + sg + ".tmp";
					line += "\n";
					Files.write(Paths.get(outpath), line.getBytes(), StandardOpenOption.CREATE, StandardOpenOption.APPEND);
				}
			}
			
			br.close();
			Files.deleteIfExists(Paths.get(inpath));
		}
		
		//Go to the file for each subgroup and read the file into a record list
			//Sort the record list by query name
			//Output to a new file
			//Delete old file for the subgroup
		for (int i = 0; i < 256; i++)
		{
			String inpath = dir + File.separator + "sg_" + i + ".tmp";
			FileReader fr = new FileReader(inpath);
			BufferedReader br = new BufferedReader(fr);
			
			List<SAMRecord> list = new LinkedList<SAMRecord>();
			String line = null;
			while((line = br.readLine()) != null)
			{
				try 
				{
					list.add(SAMRecord.parseSAMRecord(line, gb, false).getRecord());
				} 
				catch (UnsupportedFileTypeException e) 
				{
					e.printStackTrace();
				} 
				catch (InvalidSAMRecordException e) 
				{
					e.printStackTrace();
				}
			}
			
			br.close();
			
			Collections.sort(list);
			
			//Write back out
			String outpath = dir + File.separator + "sg_" + i + "_sorted.tmp";
			FileWriter fw = new FileWriter(outpath);
			BufferedWriter bw = new BufferedWriter(fw);
			for (SAMRecord r : list) bw.write(r.writeSAMRecord(false) + "\n");
			bw.close();
			
			Files.deleteIfExists(Paths.get(inpath));
		}
		
		
		return dir;
	}
	
}
