package hospelhornbg_sviewerGUI;

import javax.swing.JDialog;
import javax.swing.JTextField;
import javax.swing.JButton;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.border.BevelBorder;

import hospelhornbg_segregation.Pedigree;
import waffleoRai_GUITools.ComponentGroup;

import javax.swing.JLabel;
import javax.swing.JOptionPane;

import java.awt.Font;
import java.awt.Frame;
import java.awt.Dimension;

public class PedigreeForm extends JDialog {
	
	private static final long serialVersionUID = -7873302436199578305L;
	
	private JTextField txtPedPath;
	private JTable table;
	private JLabel lblFamilyName;
	
	private ComponentGroup disable_group;
	
	private ViewManager manager;
	private Pedigree loadedFam;
	
	public PedigreeForm(Frame parent, ViewManager m) 
	{
		super(parent, true);
		disable_group = new ComponentGroup();
		loadedFam = null;
		manager = m;
		initGUI();
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

	/* ---- Enable/Disable ---- */
	
	
	
	/* ---- Form Updates ---- */
	
	/* ---- Actions ---- */
	
	public void browseForPed()
	{
		
	}
	
	public void loadPed()
	{
		
	}
	
	public void finishAndClose()
	{
		
	}
	
	/* ---- Error ---- */
	
	public void showError(String message)
	{
		JOptionPane.showMessageDialog(this, message, "Error", JOptionPane.ERROR_MESSAGE);
	}

}
