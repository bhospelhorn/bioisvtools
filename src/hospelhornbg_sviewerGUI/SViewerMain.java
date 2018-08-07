package hospelhornbg_sviewerGUI;

import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.IOException;

import javax.swing.SwingUtilities;

import waffleoRai_Utils.FileBuffer;

public class SViewerMain {

	private static String lastVCF;
	private static String lastBED;
	private static String lastExport;
	private static String lastTruth;
	private static String lastTable;
	
	public static void writeLastPaths()
	{
		try
		{
			String tempdir = FileBuffer.getTempDir();
			String tPath = tempdir + File.separator + ".SViewerpaths.tmp";
			FileBuffer temp = new FileBuffer(2048);
			temp.printASCIIToFile(lastVCF + "\n");
			temp.printASCIIToFile(lastBED + "\n");
			temp.printASCIIToFile(lastExport + "\n");
			temp.printASCIIToFile(lastTruth + "\n");
			temp.printASCIIToFile(lastTable + "\n");
			temp.writeFile(tPath);
		}
		catch (IOException e)
		{
			e.printStackTrace();
			System.out.println("Path saving failed...");
		}
	}
	
	public static void readLastPaths()
	{
		try
		{
			String tempdir = FileBuffer.getTempDir();
			String tPath = tempdir + File.separator + ".SViewerpaths.tmp";
			if (!FileBuffer.fileExists(tPath)) return;
			FileBuffer temp = new FileBuffer(tPath);
			long cPos = 0;
			lastVCF = temp.getASCII_string(cPos, '\n'); cPos += lastVCF.length() + 1;
			if (lastVCF.isEmpty() || lastVCF.equals("null")) lastVCF = null;
			lastBED = temp.getASCII_string(cPos, '\n'); cPos += lastBED.length() + 1;
			if (lastBED.isEmpty() || lastBED.equals("null")) lastBED = null;
			lastExport = temp.getASCII_string(cPos, '\n'); cPos += lastExport.length() + 1;
			if (lastExport.isEmpty() || lastExport.equals("null")) lastExport = null;
			lastTruth = temp.getASCII_string(cPos, '\n'); cPos += lastTruth.length() + 1;
			if (lastTruth.isEmpty() || lastTruth.equals("null")) lastTruth = null;
			lastTable = temp.getASCII_string(cPos, '\n'); cPos += lastTable.length() + 1;
			if (lastTable.isEmpty() || lastTable.equals("null")) lastTable = null;
			
		}
		catch (IOException e)
		{
			e.printStackTrace();
			System.out.println("Path reading failed...");
		}
	}

	public static void main(String[] args) {
		
		readLastPaths();
		SwingUtilities.invokeLater(new Runnable()
		{
			public void run()
			{
				VarViewer GUI = new VarViewer();
				GUI.setLastVCFPath(lastVCF);
				GUI.setLastBEDPath(lastBED);
				GUI.setLastExportPath(lastExport);
				GUI.setLastTruthPath(lastTruth);
				GUI.setLastTablePath(lastTable);
				System.out.println("SViewerMain.main || lastVCF = " + lastVCF);
				GUI.render();
				GUI.addWindowListener(new WindowAdapter() {
					public void windowClosing(WindowEvent e)
					{
						lastVCF = GUI.getLastVCFPath();
						System.out.println("SViewerMain.main || lastVCF set to = " + lastVCF);
						lastBED = GUI.getLastBEDPath();
						lastExport = GUI.getLastExportPath();
						lastTruth = GUI.getLastTruthPath();
						lastTable = GUI.getLastTablePath();
						writeLastPaths();
						System.exit(0);
					}
				});
			}
		});

	}

}
