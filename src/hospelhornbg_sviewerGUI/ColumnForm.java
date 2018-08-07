package hospelhornbg_sviewerGUI;

import java.awt.Cursor;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JDialog;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.ListModel;
import javax.swing.JList;
import javax.swing.border.BevelBorder;

import hospelhornbg_sviewerGUI.ViewManager.ViewColumn;
import java.awt.Dimension;

public class ColumnForm extends JDialog {
	
	/* --- Constants --- */
	
	private static final long serialVersionUID = 1680332470401824721L;
	
	/* --- Instance Variables --- */
	
	private ViewManager manager;
	
	private JButton btnApply;
	private JButton btnCancel;
	
	private JScrollPane spStandardOut;
	private JScrollPane spStandardIn;
	private JScrollPane spInfoOut;
	private JScrollPane spInfoIn;
	private JScrollPane spSmplOut;
	private JScrollPane spSmplIn;
	
	private JList<ViewColumn> lstStOut;
	private JList<ViewColumn> lstStIn;
	private JList<String> lstInfoOut;
	private JList<String> lstInfoIn;
	private JList<String> lstSmplOut;
	private JList<String> lstSmplIn;
	
	private JButton btnInStd;
	private JButton btnExStd;
	private JButton btnInInfo;
	private JButton btnExInfo;
	private JButton btnInSmpl;
	private JButton btnExSmpl;
	
	/* --- Construction/Parsing --- */
	
	public ColumnForm(Frame parent, ViewManager manager)
	{
		super(parent, true);
		if (manager == null) throw new IllegalArgumentException();
		this.manager = manager;
		this.setLocationRelativeTo(parent);
		initGUI();
		setWait();
		updateForm();
		unsetWait();
	}
	
	private void initGUI()
	{
		setTitle("Edit Columns");
		setResizable(false);
		getContentPane().setLayout(null);
		
		setPreferredSize(new Dimension(450, 622));
		setMinimumSize(new Dimension(450, 622));
		
		
		btnApply = new JButton("Apply");
		btnApply.setBounds(246, 552, 89, 35);
		getContentPane().add(btnApply);
		
		btnCancel = new JButton("Cancel");
		btnCancel.setBounds(345, 552, 89, 35);
		getContentPane().add(btnCancel);
		
		JLabel lblStandardColumns = new JLabel("Standard Columns");
		lblStandardColumns.setBounds(10, 12, 116, 14);
		getContentPane().add(lblStandardColumns);
		
		spStandardOut = new JScrollPane();
		spStandardOut.setViewportBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
		spStandardOut.setBounds(10, 37, 149, 150);
		getContentPane().add(spStandardOut);
		
		lstStOut = new JList<ViewColumn>();
		spStandardOut.setViewportView(lstStOut);
		
		JLabel lblInfoFields = new JLabel("INFO Fields");
		lblInfoFields.setBounds(10, 198, 63, 14);
		getContentPane().add(lblInfoFields);
		
		JLabel lblSampleGenotypes = new JLabel("Sample Genotypes");
		lblSampleGenotypes.setBounds(10, 407, 133, 14);
		getContentPane().add(lblSampleGenotypes);
		
		btnInStd = new JButton("Include");
		btnInStd.setBounds(169, 74, 89, 23);
		getContentPane().add(btnInStd);
		btnInStd.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) {
				includeStandard();
			}
			
		});
		
		btnExStd = new JButton("Exclude");
		btnExStd.setBounds(169, 113, 89, 23);
		getContentPane().add(btnExStd);
		btnExStd.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) {
				excludeStandard();
			}
			
		});
		
		spStandardIn = new JScrollPane();
		spStandardIn.setBounds(268, 37, 149, 150);
		getContentPane().add(spStandardIn);
		
		lstStIn = new JList<ViewColumn>();
		lstStIn.setBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
		spStandardIn.setViewportView(lstStIn);
		
		spInfoOut = new JScrollPane();
		spInfoOut.setViewportBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
		spInfoOut.setBounds(10, 220, 149, 173);
		getContentPane().add(spInfoOut);
		
		lstInfoOut = new JList<String>();
		spInfoOut.setViewportView(lstInfoOut);
		
		spInfoIn = new JScrollPane();
		spInfoIn.setViewportBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
		spInfoIn.setBounds(268, 220, 149, 173);
		getContentPane().add(spInfoIn);
		
		lstInfoIn = new JList<String>();
		spInfoIn.setViewportView(lstInfoIn);
		
		spSmplOut = new JScrollPane();
		spSmplOut.setViewportBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
		spSmplOut.setBounds(10, 426, 149, 115);
		getContentPane().add(spSmplOut);
		
		lstSmplOut = new JList<String>();
		spSmplOut.setViewportView(lstSmplOut);
		
		spSmplIn = new JScrollPane();
		spSmplIn.setViewportBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
		spSmplIn.setBounds(268, 426, 149, 115);
		getContentPane().add(spSmplIn);
		
		lstSmplIn = new JList<String>();
		spSmplIn.setViewportView(lstSmplIn);
		
		btnInInfo = new JButton("Include");
		btnInInfo.setBounds(169, 281, 89, 23);
		getContentPane().add(btnInInfo);
		btnInInfo.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) {
				includeInfo();
			}
			
		});
		
		btnExInfo = new JButton("Exclude");
		btnExInfo.setBounds(169, 315, 89, 23);
		getContentPane().add(btnExInfo);
		btnExInfo.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) {
				excludeInfo();
			}
			
		});
		
		btnInSmpl = new JButton("Include");
		btnInSmpl.setBounds(169, 453, 89, 23);
		getContentPane().add(btnInSmpl);
		btnInSmpl.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) {
				includeSample();
			}
			
		});
		
		btnExSmpl = new JButton("Exclude");
		btnExSmpl.setBounds(169, 487, 89, 23);
		getContentPane().add(btnExSmpl);
		btnExSmpl.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) {
				excludeSample();
			}
			
		});
		
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
		spStandardOut.setEnabled(false);
		spStandardIn.setEnabled(false);
		spInfoOut.setEnabled(false);
		spInfoIn.setEnabled(false);
		spSmplOut.setEnabled(false);
		spSmplIn.setEnabled(false);
		
		lstStOut.setEnabled(false);
		lstStIn.setEnabled(false);
		lstInfoOut.setEnabled(false);
		lstInfoIn.setEnabled(false);
		lstSmplOut.setEnabled(false);
		lstSmplIn.setEnabled(false);
		
		btnInStd.setEnabled(false);
		btnExStd.setEnabled(false);
		btnInInfo.setEnabled(false);
		btnExInfo.setEnabled(false);
		btnInSmpl.setEnabled(false);
		btnExSmpl.setEnabled(false);
	}
	
	public void restoreEnabled()
	{
		spStandardOut.setEnabled(true);
		spStandardIn.setEnabled(true);
		spInfoOut.setEnabled(true);
		spInfoIn.setEnabled(true);
		spSmplOut.setEnabled(true);
		spSmplIn.setEnabled(true);
		
		lstStOut.setEnabled(true);
		lstStIn.setEnabled(true);
		lstInfoOut.setEnabled(true);
		lstInfoIn.setEnabled(true);
		lstSmplOut.setEnabled(true);
		lstSmplIn.setEnabled(true);
		
		btnInStd.setEnabled(lstStOut.getModel().getSize() != 0);
		btnExStd.setEnabled(lstStIn.getModel().getSize() != 0);
		btnInInfo.setEnabled(lstInfoOut.getModel().getSize() != 0);
		btnExInfo.setEnabled(lstInfoIn.getModel().getSize() != 0);
		btnInSmpl.setEnabled(lstSmplOut.getModel().getSize() != 0);
		btnExSmpl.setEnabled(lstSmplIn.getModel().getSize() != 0);
	}
	
	public void repaintAll()
	{
		lstStOut.repaint();
		lstStIn.repaint();
		lstInfoOut.repaint();
		lstInfoIn.repaint();
		lstSmplOut.repaint();
		lstSmplIn.repaint();
		
		spStandardOut.repaint();
		spStandardIn.repaint();
		spInfoOut.repaint();
		spInfoIn.repaint();
		spSmplOut.repaint();
		spSmplIn.repaint();
		
		btnInStd.repaint();
		btnExStd.repaint();
		btnInInfo.repaint();
		btnExInfo.repaint();
		btnInSmpl.repaint();
		btnExSmpl.repaint();
	}
	
	/* --- Syncing --- */
	
	public static ListModel<String> generateListModel(Collection<String> source)
	{
		DefaultListModel<String> model = new DefaultListModel<String>();
		for (String s : source) model.addElement(s);
		return model;
	}
	
	public static ListModel<ViewColumn> generateListModelColumn(Collection<ViewColumn> source)
	{
		DefaultListModel<ViewColumn> model = new DefaultListModel<ViewColumn>();
		for (ViewColumn c : source) model.addElement(c);
		return model;
	}
	
	public void updateForm()
	{
		List<ViewColumn> exColList = manager.getAllPossibleStandardFields();
		List<ViewColumn> inColList = manager.getIncludedStandardFields();
		inColList.remove(ViewColumn.OTHERINFO);
		inColList.remove(ViewColumn.SAMPLEGENO);
		exColList.removeAll(inColList);
		
		lstStOut.setModel(generateListModelColumn(exColList));
		lstStIn.setModel(generateListModelColumn(inColList));
		
		List<String> exStrList = manager.getAllPossibleInfoFields();
		List<String> inStrList = manager.getIncludedInfoFields();
		exStrList.removeAll(inStrList);
		
		lstInfoOut.setModel(generateListModel(exStrList));
		lstInfoIn.setModel(generateListModel(inStrList));
		
		exStrList = manager.getAllSamples();
		inStrList = manager.getIncludedSamples();
		exStrList.removeAll(inStrList);
		
		lstSmplOut.setModel(generateListModel(exStrList));
		lstSmplIn.setModel(generateListModel(inStrList));
		
		restoreEnabled();
		repaintAll();
	}
	
	public void updateManager()
	{
		manager.clearColumns();
		//Info fields
		boolean info = (lstInfoIn.getModel().getSize() > 0);
		if (info)
		{
			int msz = lstInfoIn.getModel().getSize();
			for (int i = 0; i < msz; i++)
			{
				//System.out.println("ColumnForm.updateManager || info field included: " + lstInfoIn.getModel().getElementAt(i));
				manager.includeInfoColumn(lstInfoIn.getModel().getElementAt(i));
			}
		}
		
		//Sample geno
		boolean geno = (lstSmplIn.getModel().getSize() > 0);
		if (geno)
		{
			int msz = lstSmplIn.getModel().getSize();
			for (int i = 0; i < msz; i++) manager.includeSampleColumn(lstSmplIn.getModel().getElementAt(i));
		}
		
		//Standard columns
		int msz = lstStIn.getModel().getSize();
		if (msz > 0)
		{
			for (int i = 0; i < msz; i++) manager.includeStandardColumn(lstStIn.getModel().getElementAt(i));
		}
		if (info) manager.includeStandardColumn(ViewColumn.OTHERINFO);
		if (geno) manager.includeStandardColumn(ViewColumn.SAMPLEGENO);
		
		
		updateForm();
	}
	
	/* --- Action --- */
	
	public static void moveBetweenListsString(JList<String> source, JList<String> target)
	{
		List<String> tempt = new LinkedList<String>();
		List<String> temps = new LinkedList<String>();
		int tSz = target.getModel().getSize();
		int sSz = source.getModel().getSize();
		for (int i = 0; i < tSz; i++) tempt.add(target.getModel().getElementAt(i));
		for (int i = 0; i < sSz; i++) temps.add(source.getModel().getElementAt(i));

		tempt.addAll(source.getSelectedValuesList());
		temps.removeAll(source.getSelectedValuesList());
		ListModel<String> modelt = generateListModel(tempt);
		ListModel<String> models = generateListModel(temps);
		target.setModel(modelt);
		source.setModel(models);
		
	}
	
	public static void moveBetweenListsColumn(JList<ViewColumn> source, JList<ViewColumn> target)
	{
		List<ViewColumn> tempt = new LinkedList<ViewColumn>();
		List<ViewColumn> temps = new LinkedList<ViewColumn>();
		int tSz = target.getModel().getSize();
		int sSz = source.getModel().getSize();
		for (int i = 0; i < tSz; i++) tempt.add(target.getModel().getElementAt(i));
		for (int i = 0; i < sSz; i++) temps.add(source.getModel().getElementAt(i));

		tempt.addAll(source.getSelectedValuesList());
		temps.removeAll(source.getSelectedValuesList());
		ListModel<ViewColumn> modelt = generateListModelColumn(tempt);
		ListModel<ViewColumn> models = generateListModelColumn(temps);
		target.setModel(modelt);
		source.setModel(models);
		
	}
	
	public void addApplyListener(ActionListener l)
	{
		btnApply.addActionListener(l);
	}
	
	public void addCancelListener(ActionListener l)
	{
		btnCancel.addActionListener(l);
	}
	
	public void includeStandard()
	{
		moveBetweenListsColumn(lstStOut, lstStIn);
		restoreEnabled();
		repaintAll();
	}
	
	public void excludeStandard()
	{
		moveBetweenListsColumn(lstStIn, lstStOut);
		restoreEnabled();
		repaintAll();
	}
	
	public void includeInfo()
	{
		moveBetweenListsString(lstInfoOut, lstInfoIn);
		restoreEnabled();
		repaintAll();
	}
	
	public void excludeInfo()
	{
		moveBetweenListsString(lstInfoIn, lstInfoOut);
		restoreEnabled();
		repaintAll();
	}
	
	public void includeSample()
	{
		moveBetweenListsString(lstSmplOut, lstSmplIn);
		restoreEnabled();
		repaintAll();
	}
	
	public void excludeSample()
	{
		moveBetweenListsString(lstSmplIn, lstSmplOut);
		restoreEnabled();
		repaintAll();
	}
	

}
