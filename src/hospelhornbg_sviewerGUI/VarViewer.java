package hospelhornbg_sviewerGUI;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenuBar;
import javax.swing.JMenu;
import java.awt.GridBagLayout;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;

import java.awt.Cursor;
import java.awt.GridBagConstraints;
import javax.swing.border.BevelBorder;
import javax.swing.filechooser.FileFilter;
import javax.swing.table.DefaultTableModel;

import hospelhornbg_bioinformatics.UCSCGVBED;
import hospelhornbg_bioinformatics.VCF;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_segregation.Candidate;
import hospelhornbg_segregation.Pedigree;
import waffleoRai_GUITools.ComponentGroup;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.IOException;

import javax.swing.JTable;
import javax.swing.SwingWorker;
import javax.swing.JPanel;
import javax.swing.JLabel;
import java.awt.Dimension;
import javax.swing.JRadioButtonMenuItem;

public class VarViewer extends JFrame{

	/* --- Constants --- */
	
	private static final long serialVersionUID = -4843504164505648686L;
	
	/* --- Instance Variables --- */
	
	private ViewManager manager;
	
	private ComponentGroup group_always;
	private ComponentGroup main_loaded;
	private ComponentGroup truth_loaded;
	private ComponentGroup pedigree_loaded;
	
	private JScrollPane scrollPane;
	private JTable table;
	
	private String lastVCF;
	private String lastBED;
	private String lastExport;
	private String lastTruthset;
	private String lastTable;
	private String lastPedigree;
	
	private JLabel lblVarCount;
	
	private JRadioButtonMenuItem rbmHG19;
	private JRadioButtonMenuItem rbmHG38;
	
	/* --- Construction/Parsing --- */
	
	public VarViewer()
	{
		manager = new ViewManager("hg19");
		lastVCF = "";
		lastBED = "";
		lastExport = "";
		lastTruthset = "";
		group_always = new ComponentGroup();
		main_loaded = new ComponentGroup();
		truth_loaded = new ComponentGroup();
		pedigree_loaded = new ComponentGroup();
		initGUI();
		restoreEnabled();
	}
	
	private void initGUI()
	{
		setTitle("Variant Viewer");
		
		JMenuBar menuBar = new JMenuBar();
		setJMenuBar(menuBar);
		
		JMenu mnFile = new JMenu("File");
		menuBar.add(mnFile);
		group_always.addComponent("mnFile", mnFile);
		
		JMenuItem mntmOpenVcf = new JMenuItem("Open VCF...");
		mnFile.add(mntmOpenVcf);
		group_always.addComponent("mntmOpenVcf", mntmOpenVcf);
		mntmOpenVcf.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) 
			{
				loadVCF();
			}
			
		});
		
		JMenuItem mntmOpenBed = new JMenuItem("Open BED...");
		mnFile.add(mntmOpenBed);
		group_always.addComponent("mntmOpenBed", mntmOpenBed);
		mntmOpenBed.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) 
			{
				loadBED();
			}
			
		});
		
		JMenuItem mntmExportVcf = new JMenuItem("Export VCF...");
		mnFile.add(mntmExportVcf);
		main_loaded.addComponent("mntmExportVcf", mntmExportVcf);
		mntmExportVcf.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) 
			{
				writeVCF();
			}
			
		});
		
		JMenuItem mntmExportTrack = new JMenuItem("Export UCSC Genome Viewer BED...");
		mnFile.add(mntmExportTrack);
		main_loaded.addComponent("mntmExportTrack", mntmExportTrack);
		mntmExportTrack.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) 
			{
				exportViewerTrackBED();
			}
			
		});

		JMenu mnFilters = new JMenu("Filters");
		menuBar.add(mnFilters);
		group_always.addComponent("mnFilters", mnFilters);
		
		JMenuItem mntmEditFilters = new JMenuItem("Edit Filters...");
		mnFilters.add(mntmEditFilters);
		main_loaded.addComponent("mntmEditFilters", mntmEditFilters);
		mntmEditFilters.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) 
			{
				editFilters();
			}
			
		});
		
		JMenuItem mntmLoadTruthSet = new JMenuItem("Load Truth Set...");
		mnFilters.add(mntmLoadTruthSet);
		main_loaded.addComponent("mntmLoadTruthSet", mntmLoadTruthSet);
		mntmLoadTruthSet.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) 
			{
				loadTruthSet();
			}
			
		});
		
		JMenuItem mntmResetTruthSet = new JMenuItem("Reset Truth Set");
		mnFilters.add(mntmResetTruthSet);
		truth_loaded.addComponent("mntmResetTruthSet", mntmResetTruthSet);
		mntmResetTruthSet.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) 
			{
				resetTruthSet();
			}
			
		});
		
		JMenuItem mntmClearTruthSet = new JMenuItem("Clear Truth Set");
		mnFilters.add(mntmClearTruthSet);
		truth_loaded.addComponent("mntmClearTruthSet", mntmClearTruthSet);
		mntmClearTruthSet.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) 
			{
				clearTruthSet();
			}
			
		});
		
		JMenuItem mntmEditColumns = new JMenuItem("Edit Columns...");
		mnFilters.add(mntmEditColumns);
		main_loaded.addComponent("mntmEditColumns", mntmEditColumns);
		mntmEditColumns.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) 
			{
				editColumns();
			}
			
		});
		
		
		JMenu mnStatistics = new JMenu("Statistics");
		menuBar.add(mnStatistics);
		group_always.addComponent("mnStatistics", mnStatistics);
		
		JMenuItem mntmGenerateCustomSv = new JMenuItem("Generate Custom SV Distribution Table...");
		mnStatistics.add(mntmGenerateCustomSv);
		main_loaded.addComponent("mntmGenerateCustomSv", mntmGenerateCustomSv);
		mntmGenerateCustomSv.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) 
			{
				exportViewerTrackBED();
			}
			
		});
		
		JMenu mnPedigree = new JMenu("Pedigree");
		main_loaded.addComponent("mnPedigree", mnPedigree);
		menuBar.add(mnPedigree);
		
		JMenuItem mntmLoadPedigree = new JMenuItem("Load Pedigree...");
		mnPedigree.add(mntmLoadPedigree);
		main_loaded.addComponent("mntmLoadPedigree", mntmLoadPedigree);
		mntmLoadPedigree.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) 
			{
				openPedigreeForm();
			}
			
		});
		
		JMenuItem mntmAdvancedFilters = new JMenuItem("Advanced Filters...");
		mnPedigree.add(mntmAdvancedFilters);
		pedigree_loaded.addComponent("mntmAdvancedFilters", mntmAdvancedFilters);
		mntmAdvancedFilters.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) 
			{
				editAdvancedFilters();
			}
			
		});
		
		JMenu mnGenome = new JMenu("Genome");
		menuBar.add(mnGenome);
		group_always.addComponent("mnGenome", mnGenome);
		
		rbmHG19 = new JRadioButtonMenuItem("GRCh37");
		mnGenome.add(rbmHG19);
		group_always.addComponent("rbmHG19", rbmHG19);
		rbmHG19.setSelected(true);
		rbmHG19.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) 
			{
				if (rbmHG19.isSelected()) return;
				rbmHG19.setSelected(true);
				rbmHG38.setSelected(false);
				setGenome("hg19");
			}
			
		});
		
		rbmHG38 = new JRadioButtonMenuItem("GRCh38");
		mnGenome.add(rbmHG38);
		group_always.addComponent("rbmHG38", rbmHG38);
		rbmHG38.setSelected(false);
		rbmHG38.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) 
			{
				if (rbmHG38.isSelected()) return;
				rbmHG19.setSelected(false);
				rbmHG38.setSelected(true);
				setGenome("hg38");
			}
			
		});
		
		GridBagLayout gridBagLayout = new GridBagLayout();
		gridBagLayout.columnWidths = new int[]{0, 0};
		gridBagLayout.rowHeights = new int[]{167, 18, 0};
		gridBagLayout.columnWeights = new double[]{1.0, Double.MIN_VALUE};
		gridBagLayout.rowWeights = new double[]{1.0, 1.0, Double.MIN_VALUE};
		getContentPane().setLayout(gridBagLayout);
		
		scrollPane = new JScrollPane();
		scrollPane.setViewportBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
		GridBagConstraints gbc_scrollPane = new GridBagConstraints();
		gbc_scrollPane.weighty = 1.0;
		gbc_scrollPane.weightx = 1.0;
		gbc_scrollPane.insets = new Insets(10, 10, 10, 10);
		gbc_scrollPane.fill = GridBagConstraints.BOTH;
		gbc_scrollPane.gridx = 0;
		gbc_scrollPane.gridy = 0;
		getContentPane().add(scrollPane, gbc_scrollPane);
		
		table = new JTable();
		table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		table.setFillsViewportHeight(true);
		scrollPane.setViewportView(table);
		table.addMouseListener(new MouseAdapter(){
			
			private int lastr = -1;
			private int lastl = -1;
			
			@Override
			public void mouseClicked(MouseEvent e) {
				
				if (!manager.pedigreeLoaded()) return;
				
				int r = table.getSelectedRow();
				int l = table.getSelectedColumn();
				if (lastr == r && lastl == l)
				{
					lastr = -1;
					lastl = -1;
					openCandidateInfo(r);
				}
				else
				{
					lastr = r;
					lastl = l;
				}
			}
		});
		
		JPanel panel = new JPanel();
		panel.setMaximumSize(new Dimension(32767, 50));
		panel.setMinimumSize(new Dimension(10, 20));
		panel.setLayout(null);
		GridBagConstraints gbc_panel = new GridBagConstraints();
		gbc_panel.anchor = GridBagConstraints.NORTHWEST;
		gbc_panel.fill = GridBagConstraints.BOTH;
		gbc_panel.gridx = 0;
		gbc_panel.gridy = 1;
		getContentPane().add(panel, gbc_panel);
		
		JLabel lblVariantCount = new JLabel("Variant Count:");
		lblVariantCount.setBounds(10, 11, 91, 14);
		panel.add(lblVariantCount);
		
		lblVarCount = new JLabel("0");
		lblVarCount.setBounds(108, 11, 46, 14);
		panel.add(lblVarCount);
	}

	public void render()
	{
		this.pack();
		this.setVisible(true);
	}
	
	/* --- Enabling --- */
	
	public void setWait()
	{
		setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
		disableAll();
	}
	
	public void unsetWait()
	{
		setCursor(null);
		restoreEnabled();
	}
	
	public void disableAll()
	{
		group_always.setEnabling(false);
		main_loaded.setEnabling(false);
		truth_loaded.setEnabling(false);
		pedigree_loaded.setEnabling(true);
		table.setEnabled(false);
		scrollPane.setEnabled(false);
	}
	
	public void enableAll()
	{
		group_always.setEnabling(true);
		main_loaded.setEnabling(true);
		truth_loaded.setEnabling(true);
		pedigree_loaded.setEnabling(true);
		table.setEnabled(true);
		scrollPane.setEnabled(true);
	}
	
	public void restoreEnabled()
	{
		group_always.setEnabling(true);
		table.setEnabled(true);
		scrollPane.setEnabled(true);
		main_loaded.setEnabling(manager.poolLoaded());
		truth_loaded.setEnabling(manager.truthSetLoaded());
		pedigree_loaded.setEnabling(manager.pedigreeLoaded());
	}
	
	public void repaintAll()
	{
		group_always.repaint();
		main_loaded.repaint();
		truth_loaded.repaint();
		table.repaint();
		scrollPane.repaint();
		lblVarCount.repaint();
		pedigree_loaded.repaint();
	}
	
	/* --- Syncing --- */
	
	public void updateForm()
	{
		updateTable();
		lblVarCount.setText(Integer.toString(manager.getFilteredSet().size()));
		repaintAll();
	}
	
	public void updateTable()
	{
		if (!manager.poolLoaded()) table.setModel(new DefaultTableModel());
		else table.setModel(new DefaultTableModel(manager.getTable(), manager.getColumnHeader()));
	}
	
	/* --- Action --- */
	
	public void setGenome(String gname)
	{
		SwingWorker<Void, Void> task = new SwingWorker<Void, Void>(){

			protected Void doInBackground() throws Exception 
			{
				setWait();
				try
				{
					manager.setGenome(gname);
					updateForm();
				}
				catch (IllegalArgumentException e)
				{
					e.printStackTrace();
					showError("Unknown Error - Build (" + gname + ") could not be loaded.");
				}
				return null;
			}
			
			public void done()
			{
				unsetWait();
			}
			
		};
		task.execute();

	}
	
	public void loadVCF()
	{
		//TODO: Update for candidate compatibility
		JFileChooser fc = new JFileChooser(this.lastVCF);
		fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
		fc.addChoosableFileFilter(new FileFilter()
				{
					public boolean accept(File f)
					{
						if (f.isDirectory()) return true;
						String n = f.getName();
						int dot = n.lastIndexOf('.');
						if (dot < 0) return false;
						n = n.substring(dot + 1);
						if (n.equalsIgnoreCase("vcf")) return true;
						return false;
					}
					
					public String getDescription()
					{
						return "Variant Call Format (.vcf)";
					}
				});
		int retVal = fc.showOpenDialog(this);
		if (retVal == JFileChooser.APPROVE_OPTION)
		{
			File f = fc.getSelectedFile();
			String p = f.getAbsolutePath();
			boolean sv = false;
			String question = "Would you like to attempt to parse structural variants from this VCF?";
			String title = "Open as SV File";
			int n = JOptionPane.showConfirmDialog(this, question, title, JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE);
			if (n == JOptionPane.YES_OPTION) sv = true;
			boolean readsv = sv;
			SwingWorker<Void, Void> task = new SwingWorker<Void, Void>(){

				protected Void doInBackground() throws Exception 
				{
					setWait();
					try
					{
						VariantPool pool = VCF.readVCF(p, readsv);
						manager.loadVariantPool(pool, readsv);
						updateForm();
					}
					catch (IOException e)
					{
						e.printStackTrace();
						showError("I/O Error - File (" + p + ") could not be read.");
					}
					catch (UnsupportedFileTypeException e)
					{
						e.printStackTrace();
						showError("Parsing Error - File (" + p + ") could not be read.");
					}
					catch (Exception e)
					{
						e.printStackTrace();
						showError("Unknown Error - File (" + p + ") could not be read.");
					}
					return null;
				}
				
				public void done()
				{
					unsetWait();
				}
				
			};
			task.execute();
			this.lastVCF = p;
		}
	}
	
	public void loadBED()
	{
		
	}
	
	public void writeVCF()
	{
		
	}
	
	public void editFilters()
	{
		if (manager == null)
		{
			showError("Fatal Error: No view manager present!");
			System.exit(1);
		}
		if (!manager.poolLoaded())
		{
			showError("Please load a VCF or BED file!");
			return;
		}
		FilterForm fDialog = new FilterForm(this, manager);
		fDialog.addApplyListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				fDialog.disableAll();
				fDialog.updateManager();
				fDialog.setVisible(false);
				fDialog.dispose();
				updateForm();
			}
			
		});
		fDialog.addCancelListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				fDialog.setVisible(false);
				fDialog.dispose();
			}
		});
		fDialog.setLocationRelativeTo(this);
		fDialog.pack();
		fDialog.setVisible(true);
	}
	
	public void loadTruthSet()
	{
		
	}
	
	public void resetTruthSet()
	{
		
	}
	
	public void clearTruthSet()
	{
		
	}
	
	public void editColumns()
	{
		if (manager == null)
		{
			showError("Fatal Error: No view manager present!");
			System.exit(1);
		}
		if (!manager.poolLoaded())
		{
			showError("Please load a VCF or BED file!");
			return;
		}
		ColumnForm cDialog = new ColumnForm(this, manager);
		cDialog.addApplyListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) {
				cDialog.setWait();
				cDialog.updateManager();
				cDialog.setVisible(false);
				cDialog.dispose();
				updateForm();
			}
			
		});
		cDialog.addCancelListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				cDialog.setVisible(false);
				cDialog.dispose();
			}
		});
		cDialog.setLocationRelativeTo(this);
		cDialog.pack();
		cDialog.setVisible(true);
				
	}
	
	public void editSVStats()
	{
		
	}
	
	public void exportViewerTrackBED()
	{
		JFileChooser fc = new JFileChooser(this.lastExport);
		fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
		fc.addChoosableFileFilter(new FileFilter()
				{
					public boolean accept(File f)
					{
						if (f.isDirectory()) return true;
						String n = f.getName();
						int dot = n.lastIndexOf('.');
						if (dot < 0) return false;
						n = n.substring(dot + 1);
						if (n.equalsIgnoreCase("bed")) return true;
						return false;
					}
					
					public String getDescription()
					{
						return "Browser Extensible Data Format (.bed)";
					}
				});
		int retVal = fc.showSaveDialog(this);
		if (retVal == JFileChooser.APPROVE_OPTION)
		{
			File f = fc.getSelectedFile();
			String p = f.getAbsolutePath();
			if (fc.getFileFilter().getDescription().equals("Browser Extensible Data Format (.bed)"))
			{
				int dot = p.lastIndexOf('.');
				if (dot >= 0) p = p.substring(0, dot);
				p += ".bed";
			}
			
			ExportPrompt ep = new ExportPrompt(this, manager.getAllSamples());
			ep.setLocationRelativeTo(this);
			ep.setVisible(true);
			
			String name = ep.getNameField();
			String desc = ep.getDescriptionField();
			String sample = ep.getSample();
			String path = p;
			
			SwingWorker<Void, Void> task = new SwingWorker<Void, Void>(){

				protected Void doInBackground() throws Exception 
				{
					setWait();
					try
					{
						UCSCGVBED outfile = new UCSCGVBED(name, manager.getFilteredSet());
						outfile.setDescription(desc);
						outfile.write(path, sample);
					}
					catch (IOException e)
					{
						e.printStackTrace();
						showError("I/O Error - Path (" + path + ") could not be written to.");
					}
					catch (Exception e)
					{
						e.printStackTrace();
						showError("Unknown Error - Current set could not be exported!");
					}
					return null;
				}
				
				public void done()
				{
					unsetWait();
				}
				
			};
			task.execute();
			this.lastExport = p;
		}
	}
	
	public void openCandidateInfo(int i)
	{
		if (!manager.pedigreeLoaded()) return;
		
		//Get candidate and pedigree
		Candidate c = manager.getLinkedCandidate(i);
		Pedigree fam = manager.getPedigree();
		if (c == null) return;
		if (fam == null) return;
		
		CandidateInfoForm dialog = new CandidateInfoForm(c, fam);
		dialog.setLocationRelativeTo(this);
		SwingWorker<Void, Void> task = new SwingWorker<Void, Void>(){

			protected Void doInBackground() throws Exception 
			{
				dialog.render();
				return null;
			}
		
			
		};
		task.execute();

	}
	
	public void openPedigreeForm()
	{
		if (manager == null)
		{
			showError("Fatal Error: No view manager present!");
			System.exit(1);
		}
		if (!manager.poolLoaded())
		{
			showError("Please load a VCF or BED file!");
			return;
		}
		PedigreeForm dialog = new PedigreeForm(this, manager);
		dialog.setLocationRelativeTo(this);
		dialog.setLastPath(lastPedigree);
		dialog.addWindowListener(new WindowAdapter(){
			
			public void windowClosing(WindowEvent e)
			{
				lastPedigree = dialog.getLastPath();
				dialog.dispose();
			}
		});
		dialog.render();
	}
	
	public void editAdvancedFilters()
	{
		if (manager == null)
		{
			showError("Fatal Error: No view manager present!");
			System.exit(1);
		}
		if (!manager.poolLoaded())
		{
			showError("Please load a VCF or BED file!");
			return;
		}
		if (!manager.pedigreeLoaded())
		{
			showError("Please load pedigree to edit advanced filters!");
			return;
		}
		
		MoreFilterForm dialog = new MoreFilterForm(this, manager);
		dialog.setLocationRelativeTo(this);
		dialog.addWindowListener(new WindowAdapter(){
			
			public void windowClosing(WindowEvent e)
			{
				dialog.dispose();
			}
		});
		dialog.pack();
		dialog.setVisible(true);
	}
	
	/* --- Path Bookmarking --- */
	
	public String getLastVCFPath()
	{
		return lastVCF;
	}
	
	public String getLastBEDPath()
	{
		return lastBED;
	}
	
	public String getLastTruthPath()
	{
		return lastTruthset;
	}
	
	public String getLastExportPath()
	{
		return lastExport;
	}
	
	public String getLastTablePath()
	{
		return lastTable;
	}
	
	public void setLastVCFPath(String path)
	{
		lastVCF = path;
	}
	
	public void setLastBEDPath(String path)
	{
		lastBED = path;
	}
	
	public void setLastTruthPath(String path)
	{
		lastTruthset = path;
	}
	
	public void setLastExportPath(String path)
	{
		lastExport = path;
	}
	
	public void setLastTablePath(String path)
	{
		lastTable = path;
	}
	
	/* --- Error --- */
	
	public void showError(String message)
	{
		JOptionPane.showMessageDialog(this, message, "Error", JOptionPane.ERROR_MESSAGE);
	}
}
