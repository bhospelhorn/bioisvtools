package hospelhornbg_sviewerGUI;

import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JTextField;
import javax.swing.SwingWorker;
import javax.swing.JButton;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.border.BevelBorder;
import javax.swing.filechooser.FileFilter;
import javax.swing.table.DefaultTableModel;

import hospelhornbg_bioinformatics.VCF;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_segregation.Individual;
import hospelhornbg_segregation.Pedigree;
import waffleoRai_GUITools.ComponentGroup;
import waffleoRai_Utils.FileBuffer.UnsupportedFileTypeException;

import javax.swing.JLabel;
import javax.swing.JOptionPane;

import java.awt.Font;
import java.awt.Frame;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.awt.Cursor;
import java.awt.Dimension;

public class PedigreeForm extends JDialog {
	
	private static final long serialVersionUID = -7873302436199578305L;
	
	private JTextField txtPedPath;
	private JTable table;
	private JLabel lblFamilyName;
	
	private ComponentGroup disable_group;
	
	private ViewManager manager;
	private Pedigree loadedFam;
	
	private String lastped;
	
	public PedigreeForm(Frame parent, ViewManager m) 
	{
		super(parent, true);
		disable_group = new ComponentGroup();
		loadedFam = null;
		manager = m;
		initGUI();
		if (m != null) loadedFam = manager.getPedigree();
		updateForm();
	}
	
	private void initGUI()
	{
		setPreferredSize(new Dimension(450, 350));
		setMinimumSize(new Dimension(450, 350));
		setTitle("Load Pedigree");
		getContentPane().setLayout(null);
		
		txtPedPath = new JTextField();
		txtPedPath.setBounds(10, 11, 414, 20);
		getContentPane().add(txtPedPath);
		txtPedPath.setColumns(10);
		disable_group.addComponent("txtPedPath", txtPedPath);
		
		JButton btnLoad = new JButton("Load");
		btnLoad.setBounds(335, 42, 89, 23);
		getContentPane().add(btnLoad);
		disable_group.addComponent("btnLoad", btnLoad);
		
		JButton btnBrowse = new JButton("Browse...");
		btnBrowse.setBounds(241, 42, 89, 23);
		getContentPane().add(btnBrowse);
		disable_group.addComponent("btnBrowse", btnBrowse);
		
		JButton btnDone = new JButton("Done");
		btnDone.setBounds(335, 277, 89, 23);
		getContentPane().add(btnDone);
		disable_group.addComponent("btnDone", btnDone);
		
		JScrollPane scrollPane = new JScrollPane();
		scrollPane.setViewportBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
		scrollPane.setBounds(10, 101, 414, 165);
		getContentPane().add(scrollPane);
		disable_group.addComponent("scrollPane", scrollPane);
		
		table = new JTable();
		table.setFillsViewportHeight(true);
		scrollPane.setViewportView(table);
		disable_group.addComponent("table", table);
		
		JLabel lblFamily = new JLabel("Family:");
		lblFamily.setFont(new Font("Tahoma", Font.BOLD, 11));
		lblFamily.setBounds(10, 74, 46, 14);
		getContentPane().add(lblFamily);
		
		lblFamilyName = new JLabel("[]");
		lblFamilyName.setBounds(61, 74, 163, 14);
		getContentPane().add(lblFamilyName);
		
		
	}

	public void render()
	{
		this.pack();
		this.setVisible(true);
	}
	
	/* ---- Enable/Disable ---- */
	
	public void enableAll()
	{
		disable_group.setEnabling(true);
	}
	
	public void disableAll()
	{
		disable_group.setEnabling(false);
	}
	
	public void repaintAll()
	{
		disable_group.repaint();
	}
	
	public void setWait()
	{
		setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
		disableAll();
	}
	
	public void unsetWait()
	{
		setCursor(null);
		enableAll();
	}
	
	/* ---- Table ---- */
	
	public String[] getColumnHeaders()
	{
		String[] cheaders = {"Name", "Affected", "Sex", "Father", "Mother"};
		return cheaders;
	}
	
	public String[][] makePEDTable()
	{
		if (loadedFam == null) return null;
		List<Individual> ilist = loadedFam.getAllMembers();
		String[][] table = new String[ilist.size()][5]; 
		int r = 0;
		for (Individual i : ilist)
		{
			table[r][0] = i.getName();
			table[r][1] = i.getENGString_affected();
			table[r][2] = i.getENGString_sex();
			Individual dad = i.getFather();
			if (dad != null) table[r][3] = dad.getName();
			else table[r][3] = "[Unknown]";
			Individual mom = i.getMother();
			if (mom != null) table[r][4] = mom.getName();
			else table[r][4] = "[Unknown]";
			r++;
		}
		
		return table;
	}
	
	
	/* ---- Form Updates ---- */
	
	public void updateForm()
	{
		txtPedPath.setText("");
		if(loadedFam !=  null)
		{
			lblFamilyName.setText(loadedFam.getFamilyName());
			String[] colheaders = getColumnHeaders();
			String[][] tbl = makePEDTable();
			table.setModel(new DefaultTableModel(tbl, colheaders));
		}
		else
		{
			lblFamilyName.setText("[No Family Loaded]");
			table.setModel(new DefaultTableModel());
		}
		repaintAll();
	}
	
	/* ---- Actions ---- */
	
	public void browseForPed()
	{
		JFileChooser fc = new JFileChooser(lastped);
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
						if (n.equalsIgnoreCase("ped")) return true;
						return false;
					}
					
					public String getDescription()
					{
						return "Pedigree File (.ped)";
					}
				});
		int retVal = fc.showOpenDialog(this);
		if (retVal == JFileChooser.APPROVE_OPTION)
		{
			File f = fc.getSelectedFile();
			String p = f.getAbsolutePath();
			
			this.txtPedPath.setText(p);
			
		}
	}
	
	public void loadPed()
	{
		String p = this.txtPedPath.getText();
		
		SwingWorker<Void, Void> task = new SwingWorker<Void, Void>(){

			protected Void doInBackground() throws Exception 
			{
				setWait();
				try
				{
					loadedFam = new Pedigree(p);
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
				updateForm();
				unsetWait();
			}
			
		};
		task.execute();
		this.lastped = p;
	}
	
	public void loadPedigreeIntoManager()
	{
		manager.loadPedigree(loadedFam);
	}
	
	/* ---- Save Paths ---- */
	
	public String getLastPath()
	{
		return lastped;
	}
	
	public void setLastPath(String path)
	{
		lastped = path;
	}
	
	/* ---- Error ---- */
	
	public void showError(String message)
	{
		JOptionPane.showMessageDialog(this, message, "Error", JOptionPane.ERROR_MESSAGE);
	}

}
