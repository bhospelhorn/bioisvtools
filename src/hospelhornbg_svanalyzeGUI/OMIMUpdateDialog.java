package hospelhornbg_svanalyzeGUI;

import java.awt.Frame;
import java.io.IOException;

import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JProgressBar;
import javax.swing.SwingWorker;

import hospelhornbg_genomeBuild.GeneSet;
import hospelhornbg_genomeBuild.GenomeBuild;
import hospelhornbg_genomeBuild.GenomeBuildUID;
import hospelhornbg_svproject.CommonLoader;
import hospelhornbg_svproject.OMIMUpdateListener;

import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Font;

//Right now, only looks for GRCh37 & GRCh38

public class OMIMUpdateDialog extends JDialog implements OMIMUpdateListener{

	/* --- Constants --- */
	
	private static final long serialVersionUID = -7962828078633935537L;
	
	public static final int WIDTH = 280;
	public static final int HEIGHT = 160;
	
	public static final GenomeBuildUID[] GENOME_LIST = {GenomeBuildUID.GRCh37, GenomeBuildUID.GRCh38};
	
	/* --- Instance Variables --- */
	
	private JLabel lblGenome;
	
	private JLabel lblStage;
	private JLabel lblProgress;
	private JProgressBar progressBar;
	
	private int geneCount;
	
	/* --- Construction --- */
	
	public OMIMUpdateDialog(Frame parent)
	{
		super(parent, true);
		initGUI();
		geneCount = 0;
	}
	
	private void initGUI()
	{
		setTitle("Update OMIM Gene List");
		setResizable(false);
		setMinimumSize(new Dimension(WIDTH, HEIGHT));
		setPreferredSize(new Dimension(WIDTH, HEIGHT));
		getContentPane().setLayout(null);
		
		lblStage = new JLabel("Downloading...");
		lblStage.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblStage.setBounds(22, 44, 220, 14);
		getContentPane().add(lblStage);
		
		lblProgress = new JLabel("0/0 bytes");
		lblProgress.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblProgress.setBounds(22, 93, 220, 14);
		getContentPane().add(lblProgress);
		
		progressBar = new JProgressBar();
		progressBar.setBounds(22, 62, 220, 20);
		getContentPane().add(progressBar);
		
		lblGenome = new JLabel("Working on: GENOME (0/0)");
		lblGenome.setFont(new Font("Tahoma", Font.BOLD, 11));
		lblGenome.setBounds(22, 19, 247, 14);
		getContentPane().add(lblGenome);
	}
	
	public void renderAndRun()
	{
		this.pack();
		this.setVisible(true);
		start();
	}
	
	/* --- Listening --- */
	
	public void onDownloadStart()
	{
		lblStage.setText("Downloading...");
		lblStage.repaint();
		lblProgress.setText("0 of 0 bytes downloaded");
		lblProgress.repaint();
		progressBar.setMinimum(0);
		progressBar.setMaximum(0);
		progressBar.setValue(0);
		progressBar.repaint();
	}
	
	public void onDownloadProgressUpdate(int readBytes, int appTotalBytes)
	{
		lblProgress.setText(readBytes + " of " + appTotalBytes + " bytes downloaded");
		lblProgress.repaint();
		progressBar.setMaximum(appTotalBytes);
		progressBar.setValue(readBytes);
		progressBar.repaint();
	}
	
	public void onDownloadFail()
	{
		showError("OMIM update for current genome build failed!");
	}
	
	public void onDownloadComplete()
	{
		lblStage.setText("Download complete!");
		lblStage.repaint();
	}
	
	public void onInvalidGeneSetFound()
	{
		showError("GeneSet for this build is invalid! Transcript list could not be generated...");
	}
	
	public void onTableReadStart()
	{
		lblStage.setText("Reading Table...");
		lblStage.repaint();
		lblProgress.setText("0 lines read");
		lblProgress.repaint();
		progressBar.setMinimum(0);
		progressBar.setMaximum(0);
		progressBar.setValue(0);
		progressBar.repaint();
	}
	
	public void onReadTableLine(int lineNumber)
	{
		lblProgress.setText(lineNumber + " lines read");
		lblProgress.repaint();
	}
	
	public void onTableReadComplete()
	{
		lblStage.setText("Table read complete!");
		lblStage.repaint();
	}
	
	public void onWritePrepareStart()
	{
		lblStage.setText("Preparing to write transcript list...");
		lblStage.repaint();
		lblProgress.setText("");
		lblProgress.repaint();
	}
	
	public void onWritePrepareComplete(int geneCount)
	{
		lblStage.setText("Write prepare complete! " + geneCount + " genes found!");
		lblStage.repaint();
		lblProgress.setText("");
		lblProgress.repaint();
		this.geneCount = geneCount;
	}
	
	public void onWriteStart()
	{
		lblStage.setText("Writing transcript list...");
		lblStage.repaint();
		lblProgress.setText("0 of " + geneCount + " records written...");
		lblProgress.repaint();
		progressBar.setMinimum(0);
		progressBar.setMaximum(geneCount);
		progressBar.setValue(0);
		progressBar.repaint();
	}
	
	public void onWriteGeneIDs(int geneNumber)
	{
		lblProgress.setText(geneNumber + " of " + geneCount + " records written...");
		lblProgress.repaint();
		progressBar.setValue(geneNumber);
		progressBar.repaint();
	}
	
	public void onWriteComplete()
	{
		lblStage.setText("Writing complete!");
		lblStage.repaint();
	}
	
	/* --- Function --- */
	
	public void start()
	{
		//Set wait
		setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
		OMIMUpdateDialog me = this;
		
		SwingWorker<Void, Void> task = new SwingWorker<Void, Void>(){

			protected Void doInBackground() throws Exception 
			{
				try
				{
					for(GenomeBuildUID gbuid : GENOME_LIST)
					{
						GenomeBuild gb = CommonLoader.loadGenomeBuild(gbuid.getUID());
						GeneSet gs = CommonLoader.loadGeneSet(gb, true);
						CommonLoader.updateOMIMList(gs, me);
					}
				}
				catch (IOException e)
				{
					e.printStackTrace();
					showError("I/O Error! There was an error reading or writing a file!\n"
							+ "See stderr for details.");
				}
				catch (Exception e)
				{
					e.printStackTrace();
					showError("Unknown Error! Update for current genome failed!");
				}
				return null;
			}
			
			public void done()
			{
				setCursor(null);
			}
			
		};
		task.execute();
		
	}
	
	/* ---- Error ---- */
	
	public void showError(String message)
	{
		JOptionPane.showMessageDialog(this, message, "Error", JOptionPane.ERROR_MESSAGE);
	}
	
	
}
