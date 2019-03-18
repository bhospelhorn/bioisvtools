package hospelhornbg_svanalyzeGUI;

import javax.swing.JDialog;
import javax.swing.JOptionPane;
import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.JList;
import javax.swing.border.BevelBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import hospelhornbg_genomeBuild.GenomeBuildUID;
import hospelhornbg_svproject.CommonLoader;
import waffleoRai_GUITools.InExListPair;

import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.List;

import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JSeparator;

public class InstallGenomeDialog extends JDialog{

	/* --- Constants --- */
	
	private static final long serialVersionUID = 1777481056208879121L;
	
	public static final int WIDTH = 410;
	public static final int HEIGHT = 410;
	
	/* --- Instance Variables --- */
	
	private InExListPair<GenomeBuildUID> listPair;
	
	/* --- Construction --- */
	
	public InstallGenomeDialog()
	{
		initGUI();
	}
	
	private void initGUI()
	{
		setTitle("Install Genome Build");
		setResizable(false);
		setMinimumSize(new Dimension(WIDTH, HEIGHT));
		setPreferredSize(new Dimension(WIDTH, HEIGHT));
		getContentPane().setLayout(null);
		
		JLabel lblBuildsInJar = new JLabel("Builds In JAR:");
		lblBuildsInJar.setBounds(10, 140, 87, 14);
		getContentPane().add(lblBuildsInJar);
		
		JLabel lblInstalledBuilds = new JLabel("Installed Builds:");
		lblInstalledBuilds.setBounds(10, 11, 87, 14);
		getContentPane().add(lblInstalledBuilds);
		
		JLabel lblBuildsToInstall = new JLabel("Builds to Install:");
		lblBuildsToInstall.setBounds(209, 140, 87, 14);
		getContentPane().add(lblBuildsToInstall);
		
		JScrollPane spInstalled = new JScrollPane();
		spInstalled.setViewportBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
		spInstalled.setBounds(10, 36, 181, 89);
		getContentPane().add(spInstalled);
		
		JList<GenomeBuildUID> lstInstall = new JList<GenomeBuildUID>();
		spInstalled.setViewportView(lstInstall);
		
		JButton btnInstall = new JButton("Install");
		btnInstall.setBounds(207, 347, 89, 23);
		getContentPane().add(btnInstall);
		btnInstall.setEnabled(false);
		btnInstall.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) 
			{
				installGenomes();
			}
			
		});
		
		JButton btnCancel = new JButton("Cancel");
		btnCancel.setBounds(306, 347, 89, 23);
		getContentPane().add(btnCancel);
		btnCancel.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) 
			{
				setVisible(false);
				dispose();
			}
			
		});
		
		JScrollPane spJarBuilds = new JScrollPane();
		spJarBuilds.setViewportBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
		spJarBuilds.setBounds(20, 164, 158, 125);
		getContentPane().add(spJarBuilds);
		
		JList<GenomeBuildUID> lstJarBuilds = new JList<GenomeBuildUID>();
		spJarBuilds.setViewportView(lstJarBuilds);
		
		JScrollPane spToInstall = new JScrollPane();
		spToInstall.setViewportBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
		spToInstall.setBounds(221, 165, 158, 125);
		getContentPane().add(spToInstall);
		
		JList<GenomeBuildUID> lstToInstall = new JList<GenomeBuildUID>();
		spToInstall.setViewportView(lstToInstall);
		
		JButton btnAdd = new JButton("Add");
		btnAdd.setBounds(30, 297, 89, 23);
		getContentPane().add(btnAdd);
		btnAdd.setEnabled(false);
		btnAdd.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) 
			{
				listPair.includeSelected();
				listPair.updateSourceLists();
				listPair.updateGraphicLists();
				btnInstall.setEnabled(!listPair.getSourceIncludeList().isEmpty());
			}
			
		});
		lstJarBuilds.addListSelectionListener(new ListSelectionListener() {

			@Override
			public void valueChanged(ListSelectionEvent e) 
			{
				btnAdd.setEnabled(!lstJarBuilds.isSelectionEmpty());
			}
			
		});
		
		JButton btnRemove = new JButton("Remove");
		btnRemove.setBounds(231, 297, 89, 23);
		getContentPane().add(btnRemove);
		btnRemove.setEnabled(false);
		btnRemove.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) 
			{
				listPair.excludeSelected();
				listPair.updateSourceLists();
				listPair.updateGraphicLists();
			}
			
		});
		lstToInstall.addListSelectionListener(new ListSelectionListener() {

			@Override
			public void valueChanged(ListSelectionEvent e) 
			{
				btnRemove.setEnabled(!lstToInstall.isSelectionEmpty());
			}
			
		});
		
		JSeparator separator = new JSeparator();
		separator.setBounds(10, 330, 384, 2);
		getContentPane().add(separator);
		
		//Prepare lists
		GenomeBuildUID[] all = GenomeBuildUID.values();
		List<GenomeBuildUID> gblist = new ArrayList<GenomeBuildUID>(all.length + 1);
		for (GenomeBuildUID gb : all) gblist.add(gb);
		listPair = new InExListPair<GenomeBuildUID>(lstJarBuilds, lstToInstall, null, gblist);
		
		//Installed list
		List<GenomeBuildUID> installed = CommonLoader.getInstalledGenomes();
		populateList(lstInstall, installed);
		
		listPair.updateGraphicLists();
		
	}
	
	private void populateList(JList<GenomeBuildUID> guiList, List<GenomeBuildUID> source)
	{
		DefaultListModel<GenomeBuildUID> model = new DefaultListModel<GenomeBuildUID>();
		for (GenomeBuildUID id : source) model.addElement(id);
		guiList.setModel(model);
		guiList.repaint();
	}
	
	public void render()
	{
		pack();
		setVisible(true);
	}
	
	/* --- Function --- */
	
	public void installGenomes()
	{
		//TODO: Write
	}
	
	/* ---- Error ---- */
	
	public void showError(String message)
	{
		JOptionPane.showMessageDialog(this, message, "Error", JOptionPane.ERROR_MESSAGE);
	}
	
}
