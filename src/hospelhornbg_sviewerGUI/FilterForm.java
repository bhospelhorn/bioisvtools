package hospelhornbg_sviewerGUI;

import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JCheckBox;
import java.awt.Font;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import javax.swing.JScrollPane;
import javax.swing.JButton;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.border.BevelBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import hospelhornbg_bioinformatics.FlagInfoFilter;
import hospelhornbg_bioinformatics.FloatInfoFilter;
import hospelhornbg_bioinformatics.IntegerInfoFilter;
import hospelhornbg_bioinformatics.SVType;
import hospelhornbg_bioinformatics.StringInfoFilter;
import hospelhornbg_bioinformatics.VariantFilter;
import hospelhornbg_bioinformatics.VariantPool;
import hospelhornbg_genomeBuild.Contig;
import waffleoRai_GUITools.CheckBoxGroup;
import waffleoRai_GUITools.ComponentGroup;
import waffleoRai_GUITools.InExListPair;
import waffleoRai_GUITools.ModelManager;
import waffleoRai_GUITools.RadioButtonGroup;

import javax.swing.JSeparator;
import javax.swing.SwingConstants;
import javax.swing.JTextField;
import javax.swing.JComboBox;
import javax.swing.JRadioButton;
import java.awt.Dimension;

public class FilterForm extends JDialog{
	
	/* --- Constants --- */
	
	private static final long serialVersionUID = 4077851671466248942L;
	
	public static final int SVTYPE_INDEX_DEL = 0;
	public static final int SVTYPE_INDEX_DUP = 1;
	public static final int SVTYPE_INDEX_INS = 2;
	public static final int SVTYPE_INDEX_INV = 3;
	public static final int SVTYPE_INDEX_CNV = 4;
	public static final int SVTYPE_INDEX_BND = 5;
	public static final int SVTYPE_INDEX_OTHER = 6;
	
	public static final int MULTIMODE_INDEX_ALL = 0;
	public static final int MULTIMODE_INDEX_ANY = 1;
	public static final int MULTIMODE_INDEX_NONE = 2;
	
	public static final int COMPMODE_INDEX_EQ = 0;
	public static final int COMPMODE_INDEX_NE = 1;
	public static final int COMPMODE_INDEX_GT = 2;
	public static final int COMPMODE_INDEX_GE = 3;
	public static final int COMPMODE_INDEX_LT = 4;
	public static final int COMPMODE_INDEX_LE = 5;
	
	/* --- Instance Variables --- */
	
	private ViewManager manager;
	
	private ComponentGroup painter;
	private JButton btnApply;
	private JButton btnCancel;
	
	private JCheckBox cbConfirmed;
	private JCheckBox cbSVIn;
	private JCheckBox cbSNVIn;
	
	private JTextField txtMinSize;
	private JTextField txtMaxSize;
	private JTextField txtQual;
	private JCheckBox cbQual;
	private JCheckBox cbEXqual;
	
	private JList<VariantFilter> lstCustom;
	private JComboBox<String> cmbxField;
	private JLabel lblField;
	private JTextField txtThreshold;
	private JButton btnSave;
	
	private InExListPair<Contig> chromosomeLists;
	
	private RadioButtonGroup multiMode;
	private RadioButtonGroup compareMode;
	private RadioButtonGroup BNDMode;
	private CheckBoxGroup SVTYPES;
	
	//private VariantFilter myFilter;
	private int fieldType; //Bookkeeping for quick enabling, checking
	
	/* --- Construction/Parsing --- */
	
	public FilterForm(Frame parent, ViewManager viewmanager)
	{
		super(parent, true);
		manager = viewmanager;
		painter = new ComponentGroup();
		multiMode = new RadioButtonGroup(3);
		BNDMode = new RadioButtonGroup(2);
		compareMode = new RadioButtonGroup(6);
		SVTYPES = new CheckBoxGroup(7);
		initGUI();
		//myFilter = null;
		fieldType = VariantPool.INFODEF_UNK;
		syncToManager();
	}
	
	private void initGUI()
	{
		setPreferredSize(new Dimension(840, 525));
		setResizable(false);
		setMinimumSize(new Dimension(840, 525));
		setTitle("Variant Filters");
		getContentPane().setLayout(null);
		
		JLabel lblStandard = new JLabel("Standard");
		lblStandard.setFont(new Font("Tahoma", Font.BOLD, 13));
		lblStandard.setBounds(10, 11, 72, 14);
		getContentPane().add(lblStandard);
		
		cbConfirmed = new JCheckBox("Confirmed Variants Only");
		cbConfirmed.setFont(new Font("Tahoma", Font.PLAIN, 11));
		cbConfirmed.setBounds(20, 32, 147, 23);
		getContentPane().add(cbConfirmed);
		
		JScrollPane spChromIn = new JScrollPane();
		spChromIn.setViewportBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
		spChromIn.setBounds(40, 104, 172, 142);
		getContentPane().add(spChromIn);
		painter.addComponent("spChromIn", spChromIn);
		
		JList<Contig> lstChromIn = new JList<Contig>();
		spChromIn.setViewportView(lstChromIn);
		
		JLabel lblChromosomes = new JLabel("Chromosomes");
		lblChromosomes.setFont(new Font("Tahoma", Font.BOLD, 11));
		lblChromosomes.setBounds(30, 64, 81, 14);
		getContentPane().add(lblChromosomes);
		
		JLabel lblIncluded = new JLabel("Included");
		lblIncluded.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblIncluded.setBounds(40, 89, 46, 14);
		getContentPane().add(lblIncluded);
		
		JScrollPane spChromEx = new JScrollPane();
		spChromEx.setViewportBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
		spChromEx.setBounds(272, 104, 172, 142);
		getContentPane().add(spChromEx);
		painter.addComponent("spChromEx", spChromEx);
		
		JList<Contig> lstChromEx = new JList<Contig>();
		spChromEx.setViewportView(lstChromEx);
		
		JLabel lblExcluded = new JLabel("Excluded");
		lblExcluded.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblExcluded.setBounds(271, 89, 46, 14);
		getContentPane().add(lblExcluded);
		
		chromosomeLists = new InExListPair<Contig>(lstChromEx, lstChromIn, manager.getChromosomeList(), manager.getAllChromosomes());
		
		JButton btnChromEx = new JButton("->");
		btnChromEx.setFont(new Font("Tahoma", Font.PLAIN, 11));
		btnChromEx.setBounds(216, 138, 46, 23);
		getContentPane().add(btnChromEx);
		painter.addComponent("btnChromEx", btnChromEx);
		btnChromEx.setEnabled(!lstChromIn.isSelectionEmpty());
		lstChromIn.addListSelectionListener(new ListSelectionListener(){

			public void valueChanged(ListSelectionEvent e) {
				btnChromEx.setEnabled(!lstChromIn.isSelectionEmpty());
			}
			
		});
		btnChromEx.addActionListener(new ActionListener(){

			public void actionPerformed(ActionEvent e) {
				excludeChrom();
			}
			
		});
		
		JButton btnChromIn = new JButton("<-");
		btnChromIn.setFont(new Font("Tahoma", Font.PLAIN, 11));
		btnChromIn.setBounds(216, 172, 46, 23);
		getContentPane().add(btnChromIn);
		painter.addComponent("btnChromIn", btnChromIn);
		btnChromIn.setEnabled(!lstChromEx.isSelectionEmpty());
		lstChromEx.addListSelectionListener(new ListSelectionListener(){

			public void valueChanged(ListSelectionEvent e) {
				btnChromIn.setEnabled(!lstChromEx.isSelectionEmpty());
			}
			
		});
		btnChromIn.addActionListener(new ActionListener(){

			public void actionPerformed(ActionEvent e) {
				includeChrom();
			}
			
		});
		
		JLabel lblSvTypes = new JLabel("SV Types");
		lblSvTypes.setFont(new Font("Tahoma", Font.BOLD, 11));
		lblSvTypes.setBounds(30, 257, 59, 14);
		getContentPane().add(lblSvTypes);
		
		JCheckBox cbDEL = new JCheckBox("DEL");
		cbDEL.setFont(new Font("Tahoma", Font.PLAIN, 11));
		cbDEL.setBounds(40, 274, 46, 23);
		getContentPane().add(cbDEL);
		SVTYPES.addCheckBox(cbDEL, SVTYPE_INDEX_DEL);
		
		JCheckBox cbDUP = new JCheckBox("DUP");
		cbDUP.setFont(new Font("Tahoma", Font.PLAIN, 11));
		cbDUP.setBounds(88, 274, 46, 23);
		getContentPane().add(cbDUP);
		SVTYPES.addCheckBox(cbDUP, SVTYPE_INDEX_DUP);
		
		JCheckBox cbINV = new JCheckBox("INV");
		cbINV.setFont(new Font("Tahoma", Font.PLAIN, 11));
		cbINV.setBounds(181, 274, 46, 23);
		getContentPane().add(cbINV);
		SVTYPES.addCheckBox(cbINV, SVTYPE_INDEX_INV);
		
		JCheckBox cbINS = new JCheckBox("INS");
		cbINS.setFont(new Font("Tahoma", Font.PLAIN, 11));
		cbINS.setBounds(136, 274, 46, 23);
		getContentPane().add(cbINS);
		SVTYPES.addCheckBox(cbINS, SVTYPE_INDEX_INS);
		
		JCheckBox cbCNV = new JCheckBox("CNV");
		cbCNV.setFont(new Font("Tahoma", Font.PLAIN, 11));
		cbCNV.setBounds(229, 274, 46, 23);
		getContentPane().add(cbCNV);
		SVTYPES.addCheckBox(cbCNV, SVTYPE_INDEX_CNV);
		
		JCheckBox cbBND = new JCheckBox("BND");
		cbBND.setFont(new Font("Tahoma", Font.PLAIN, 11));
		cbBND.setBounds(277, 274, 46, 23);
		getContentPane().add(cbBND);
		SVTYPES.addCheckBox(cbBND, SVTYPE_INDEX_BND);
		
		JCheckBox cbSVOther = new JCheckBox("Other");
		cbSVOther.setFont(new Font("Tahoma", Font.PLAIN, 11));
		cbSVOther.setBounds(325, 274, 59, 23);
		getContentPane().add(cbSVOther);
		SVTYPES.addCheckBox(cbSVOther, SVTYPE_INDEX_OTHER);
		
		JLabel lblVariantSize = new JLabel("Variant Size");
		lblVariantSize.setFont(new Font("Tahoma", Font.BOLD, 11));
		lblVariantSize.setBounds(30, 304, 72, 14);
		getContentPane().add(lblVariantSize);
		
		JSeparator separator = new JSeparator();
		separator.setOrientation(SwingConstants.VERTICAL);
		separator.setBounds(479, 11, 2, 451);
		getContentPane().add(separator);
		
		btnApply = new JButton("Apply");
		btnApply.setBounds(272, 8, 89, 29);
		getContentPane().add(btnApply);
		painter.addComponent("btnApply", btnApply);
		
		btnCancel = new JButton("Cancel");
		btnCancel.setBounds(368, 8, 89, 29);
		getContentPane().add(btnCancel);
		painter.addComponent("btnCancel", btnCancel);
		
		JLabel lblMinimum = new JLabel("Minimum:");
		lblMinimum.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblMinimum.setBounds(40, 325, 49, 14);
		getContentPane().add(lblMinimum);
		
		txtMinSize = new JTextField();
		txtMinSize.setBounds(88, 322, 86, 20);
		getContentPane().add(txtMinSize);
		txtMinSize.setColumns(10);
		
		JLabel lblBp = new JLabel("bp");
		lblBp.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblBp.setBounds(177, 325, 18, 14);
		getContentPane().add(lblBp);
		
		JLabel lblMaximuml = new JLabel("Maximum:");
		lblMaximuml.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblMaximuml.setBounds(40, 350, 49, 14);
		getContentPane().add(lblMaximuml);
		
		txtMaxSize = new JTextField();
		txtMaxSize.setColumns(10);
		txtMaxSize.setBounds(88, 347, 86, 20);
		getContentPane().add(txtMaxSize);
		
		JLabel label_1 = new JLabel("bp");
		label_1.setFont(new Font("Tahoma", Font.PLAIN, 11));
		label_1.setBounds(177, 350, 18, 14);
		getContentPane().add(label_1);
		
		JLabel lblTypes = new JLabel("Types");
		lblTypes.setFont(new Font("Tahoma", Font.BOLD, 11));
		lblTypes.setBounds(30, 379, 46, 14);
		getContentPane().add(lblTypes);
		
		cbSVIn = new JCheckBox("Include Structural Variants");
		cbSVIn.setFont(new Font("Tahoma", Font.PLAIN, 11));
		cbSVIn.setBounds(40, 397, 155, 23);
		getContentPane().add(cbSVIn);
		
		cbSNVIn = new JCheckBox("Include Small Variants and SNVs");
		cbSNVIn.setFont(new Font("Tahoma", Font.PLAIN, 11));
		cbSNVIn.setBounds(40, 423, 187, 23);
		getContentPane().add(cbSNVIn);
		
		JLabel lblCustom = new JLabel("Custom");
		lblCustom.setFont(new Font("Tahoma", Font.BOLD, 13));
		lblCustom.setBounds(491, 12, 72, 14);
		getContentPane().add(lblCustom);
		
		JScrollPane spCustom = new JScrollPane();
		spCustom.setViewportBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
		spCustom.setBounds(513, 32, 281, 142);
		getContentPane().add(spCustom);
		painter.addComponent("spCustom", spCustom);
		
		lstCustom = new JList<VariantFilter>();
		spCustom.setViewportView(lstCustom);
		
		JButton btnRemove = new JButton("Remove");
		btnRemove.setBounds(705, 182, 89, 23);
		getContentPane().add(btnRemove);
		painter.addComponent("btnRemove", btnRemove);
		lstCustom.clearSelection();
		btnRemove.setEnabled(!lstCustom.isSelectionEmpty());
		lstCustom.addListSelectionListener(new ListSelectionListener(){

			public void valueChanged(ListSelectionEvent e) {
				btnRemove.setEnabled(!lstCustom.isSelectionEmpty());
			}
			
		});
		btnRemove.addActionListener(new ActionListener(){

			public void actionPerformed(ActionEvent e) {
				removeFilters();
			}
			
		});

		JLabel lblNewFilter = new JLabel("New Filter");
		lblNewFilter.setFont(new Font("Tahoma", Font.BOLD, 11));
		lblNewFilter.setBounds(502, 211, 59, 14);
		getContentPane().add(lblNewFilter);
		
		JLabel lblInfoFieldKey = new JLabel("INFO Field Key:");
		lblInfoFieldKey.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblInfoFieldKey.setBounds(517, 232, 81, 14);
		getContentPane().add(lblInfoFieldKey);
		
		cmbxField = new JComboBox<String>();
		cmbxField.setBounds(605, 229, 172, 20);
		getContentPane().add(cmbxField);
		cmbxField.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				selectFieldKey();
			}
		});
		
		JRadioButton rbALL = new JRadioButton("ALL");
		rbALL.setFont(new Font("Tahoma", Font.PLAIN, 11));
		rbALL.setBounds(513, 284, 50, 23);
		getContentPane().add(rbALL);
		multiMode.addButton(rbALL, MULTIMODE_INDEX_ALL);
		rbALL.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				selectMultiMode(MULTIMODE_INDEX_ALL);
			}
		});
		
		JRadioButton rbANY = new JRadioButton("ANY");
		rbANY.setFont(new Font("Tahoma", Font.PLAIN, 11));
		rbANY.setBounds(563, 284, 49, 23);
		getContentPane().add(rbANY);
		multiMode.addButton(rbANY, MULTIMODE_INDEX_ANY);
		rbANY.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				selectMultiMode(MULTIMODE_INDEX_ANY);
			}
		});
		
		JRadioButton rbNONE = new JRadioButton("NONE");
		rbNONE.setBounds(614, 284, 59, 23);
		getContentPane().add(rbNONE);
		multiMode.addButton(rbNONE, MULTIMODE_INDEX_NONE);
		rbNONE.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				selectMultiMode(MULTIMODE_INDEX_NONE);
			}
		});
		
		JLabel lblFieldType = new JLabel("Field Type:");
		lblFieldType.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblFieldType.setBounds(517, 257, 59, 14);
		getContentPane().add(lblFieldType);
		
		lblField = new JLabel("null");
		lblField.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblField.setBounds(605, 257, 172, 14);
		getContentPane().add(lblField);
		
		JRadioButton rbEQ = new JRadioButton("==");
		rbEQ.setFont(new Font("Tahoma", Font.PLAIN, 11));
		rbEQ.setBounds(513, 321, 46, 23);
		getContentPane().add(rbEQ);
		compareMode.addButton(rbEQ, COMPMODE_INDEX_EQ);
		rbEQ.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				selectCompareMode(COMPMODE_INDEX_EQ);
			}
		});
		
		JRadioButton rbNE = new JRadioButton("!=");
		rbNE.setBounds(563, 321, 46, 23);
		getContentPane().add(rbNE);
		compareMode.addButton(rbNE, COMPMODE_INDEX_NE);
		rbNE.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				selectCompareMode(COMPMODE_INDEX_NE);
			}
		});
		
		JRadioButton rbGT = new JRadioButton(">");
		rbGT.setBounds(613, 321, 46, 23);
		getContentPane().add(rbGT);
		compareMode.addButton(rbGT, COMPMODE_INDEX_GT);
		rbGT.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				selectCompareMode(COMPMODE_INDEX_GT);
			}
		});
		
		JRadioButton rbGE = new JRadioButton(">=");
		rbGE.setBounds(663, 321, 46, 23);
		getContentPane().add(rbGE);
		compareMode.addButton(rbGE, COMPMODE_INDEX_GE);
		rbGE.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				selectCompareMode(COMPMODE_INDEX_GE);
			}
		});
		
		JRadioButton rbLT = new JRadioButton("<");
		rbLT.setBounds(713, 321, 46, 23);
		getContentPane().add(rbLT);
		compareMode.addButton(rbLT, COMPMODE_INDEX_LT);
		rbLT.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				selectCompareMode(COMPMODE_INDEX_LT);
			}
		});
		
		JRadioButton rbLE = new JRadioButton("<=");
		rbLE.setBounds(763, 321, 46, 23);
		getContentPane().add(rbLE);
		compareMode.addButton(rbLE, COMPMODE_INDEX_LE);
		rbLE.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				selectCompareMode(COMPMODE_INDEX_LE);
			}
		});
		
		JLabel lblComparisonValue = new JLabel("Comparison Value:");
		lblComparisonValue.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblComparisonValue.setBounds(513, 361, 99, 14);
		getContentPane().add(lblComparisonValue);
		
		txtThreshold = new JTextField();
		txtThreshold.setBounds(526, 376, 268, 20);
		getContentPane().add(txtThreshold);
		txtThreshold.setColumns(10);
		
		btnSave = new JButton("Save");
		btnSave.setBounds(705, 407, 89, 23);
		getContentPane().add(btnSave);
		painter.addComponent("btnSave", btnSave);
		
		cbQual = new JCheckBox("Minimum Quality:");
		cbQual.setFont(new Font("Tahoma", Font.BOLD, 11));
		cbQual.setBounds(229, 321, 132, 23);
		getContentPane().add(cbQual);
		cbQual.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				txtQual.setEnabled(cbQual.isSelected());
				cbEXqual.setEnabled(cbQual.isSelected());
			}
		});
		
		txtQual = new JTextField();
		txtQual.setBounds(368, 322, 86, 20);
		getContentPane().add(txtQual);
		txtQual.setColumns(10);
		
		cbEXqual = new JCheckBox("Exclude Variants with No Quality Score");
		cbEXqual.setFont(new Font("Tahoma", Font.PLAIN, 11));
		cbEXqual.setBounds(241, 346, 216, 23);
		getContentPane().add(cbEXqual);
		
		JLabel lblBndChromosomeExclusion = new JLabel("BND Chromosome Exclusion Mode");
		lblBndChromosomeExclusion.setFont(new Font("Tahoma", Font.BOLD, 11));
		lblBndChromosomeExclusion.setBounds(229, 379, 196, 14);
		getContentPane().add(lblBndChromosomeExclusion);
		
		JRadioButton rbBNDChrModeIn = new JRadioButton("Inclusive");
		rbBNDChrModeIn.setFont(new Font("Tahoma", Font.PLAIN, 11));
		rbBNDChrModeIn.setBounds(252, 397, 72, 23);
		getContentPane().add(rbBNDChrModeIn);
		BNDMode.addButton(rbBNDChrModeIn, 0);
		
		JRadioButton rbBNDChrModeEx = new JRadioButton("Exclusive");
		rbBNDChrModeEx.setFont(new Font("Tahoma", Font.PLAIN, 11));
		rbBNDChrModeEx.setBounds(338, 397, 81, 23);
		getContentPane().add(rbBNDChrModeEx);
		BNDMode.addButton(rbBNDChrModeEx, 1);
		
		rbBNDChrModeIn.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				selectBNDMode(0);
			}
		});
		rbBNDChrModeEx.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				selectBNDMode(1);
			}
		});
		
		BNDMode.select(0);
	}

	/* --- Enabling --- */
	
	public void disableAll()
	{
		painter.setEnabling(false);
		cbConfirmed.setEnabled(false);
		cbSVIn.setEnabled(false);
		cbSNVIn.setEnabled(false);
		txtMinSize.setEnabled(false);
		txtMaxSize.setEnabled(false);
		lstCustom.setEnabled(false);
		cmbxField.setEnabled(false);
		txtThreshold.setEnabled(false);
		chromosomeLists.disable();
		multiMode.disableAll();
		compareMode.disableAll();
		BNDMode.disableAll();
		SVTYPES.disableAll();
		btnSave.setEnabled(false);
		cbQual.setEnabled(false);
		txtQual.setEnabled(false);
		cbEXqual.setEnabled(false);
	}
	
	public void syncEnabled()
	{
		painter.setEnabling(true);
		cbConfirmed.setEnabled(true);
		cbSVIn.setEnabled(true);
		cbSNVIn.setEnabled(true);
		txtMinSize.setEnabled(true);
		txtMaxSize.setEnabled(true);
		lstCustom.setEnabled(true);
		lstCustom.clearSelection();
		chromosomeLists.enable();
		SVTYPES.enableAll();
		cbQual.setEnabled(true);
		txtQual.setEnabled(cbQual.isSelected());
		cbEXqual.setEnabled(cbQual.isSelected());
		BNDMode.enableAll();
		
		cmbxField.setEnabled(true);
		setFilterEditorEnabled(cmbxField.getSelectedIndex() >= 0);
		
	}

	public void repaintAll()
	{
		cbConfirmed.repaint();
		cbSVIn.repaint();
		cbSNVIn.repaint();
		txtMinSize.repaint();
		txtMaxSize.repaint();
		lstCustom.repaint();
		cmbxField.repaint();
		txtThreshold.repaint();
		chromosomeLists.repaint();
		multiMode.repaintAll();
		compareMode.repaintAll();
		SVTYPES.repaintAll();
		painter.repaint();
		lblField.repaint();
		cbQual.repaint();
		txtQual.repaint();
		cbEXqual.repaint();
		BNDMode.repaintAll();
	}
	
	private void setFilterEditorEnabled(boolean b)
	{
		txtThreshold.setEnabled(b && (fieldType != VariantPool.INFODEF_FLAG));
		multiMode.setEnabledAll(b);
		compareMode.setEnabledAll(b);
		btnSave.setEnabled(b);
	}
	
	/* --- Syncing --- */
	
	public void syncToManager()
	{
		//Chromosomes
		chromosomeLists.updateGraphicLists();
		
		if (manager.compositeChromosomeExclusive()) BNDMode.select(1);
		else BNDMode.select(0);
		
		//Initial checkboxes
		cbConfirmed.setSelected(manager.confirmedOnly());
		cbSNVIn.setSelected(manager.includeNonSV());
		cbSVIn.setSelected(manager.includeSV());
		
		//SVTypes
		SVTYPES.select(SVTYPE_INDEX_BND, manager.SVTypeIncluded(SVType.BND));
		SVTYPES.select(SVTYPE_INDEX_DEL, manager.SVTypeIncluded(SVType.DEL));
		SVTYPES.select(SVTYPE_INDEX_DUP, manager.SVTypeIncluded(SVType.DUP));
		SVTYPES.select(SVTYPE_INDEX_INS, manager.SVTypeIncluded(SVType.INS));
		SVTYPES.select(SVTYPE_INDEX_INV, manager.SVTypeIncluded(SVType.INV));
		SVTYPES.select(SVTYPE_INDEX_CNV, manager.SVTypeIncluded(SVType.CNV));
		SVTYPES.select(SVTYPE_INDEX_OTHER, manager.SVTypeIncluded(SVType.OTHER));
		
		//Lengths
		txtMinSize.setText(Integer.toString(manager.getMinimumSize()));
		txtMaxSize.setText(Integer.toString(manager.getMaximumSize()));
		
		//Quality
		cbQual.setSelected(manager.filterByQuality());
		txtQual.setText(Double.toString(manager.getMinimumQuality()));
		cbEXqual.setSelected(manager.excludeNoQuality());
		
		//Custom filters
		lstCustom.setModel(manager.getCustomFilters_Swing());
		List<String> infolist = manager.getAllPossibleInfoFields();
		for (String i : infolist) this.cmbxField.addItem(i);
		clearFilterEditor();
	
		
		syncEnabled();
		repaintAll();
	}
	
	private void updateFieldType(String key)
	{
		if (key == null || key.isEmpty()){
			fieldType = -2;
			return;
		}
		fieldType = manager.getInfoFieldType(key);
	}
	
	private void updateFieldTypeLabel()
	{
		switch (fieldType)
		{
		case VariantPool.INFODEF_UNK: lblField.setText("unknown"); break;
		case VariantPool.INFODEF_FLAG: lblField.setText("Boolean"); break;
		case VariantPool.INFODEF_FLOAT: lblField.setText("Float"); break;
		case VariantPool.INFODEF_INT: lblField.setText("Integer"); break;
		case VariantPool.INFODEF_STRING: lblField.setText("String"); break;
		default: lblField.setText("null"); break;
		}
		lblField.repaint();
	}

	public void clearFilterEditor()
	{
		cmbxField.setSelectedIndex(-1);
		fieldType = -2;
		updateFieldTypeLabel();
		txtThreshold.setText("");
		multiMode.select(0);
		compareMode.select(0);
		setFilterEditorEnabled(false);
	}
	
	public int getFilterMultiMode()
	{
		int sel = multiMode.getSelectedIndex();
		if (sel < 0)
		{
			multiMode.select(0);
			sel = 0;
		}
		return sel;
	}
	
	public int getFilterCompareMode()
	{
		int sel = compareMode.getSelectedIndex();
		if (sel < 0)
		{
			compareMode.select(0);
			sel = 0;
		}
		return sel;
	}
	
	public IntegerInfoFilter createIntFilter(String field)
	{
		int t = 0;
		try
		{
			t = Integer.parseInt(txtThreshold.getText());
		}
		catch (NumberFormatException e)
		{
			showError("Please enter a valid integer for threshold value!");
			txtThreshold.setText("0");
			return null;
		}
		IntegerInfoFilter filter = new IntegerInfoFilter(field, getFilterCompareMode(), t);
		int multi = getFilterMultiMode();
		switch (multi)
		{
		case VariantFilter.MULTIFIELD_ALL: filter.set_multi_AND(); break;
		case VariantFilter.MULTIFIELD_ANY: filter.set_multi_OR(); break;
		case VariantFilter.MULTIFIELD_NONE: filter.set_multi_NOT(); break;
		}
		return filter;
	}
	
	public FlagInfoFilter createBoolFilter(String field)
	{
		txtThreshold.setText("");
		FlagInfoFilter filter = new FlagInfoFilter(field);
		return filter;
	}
	
	public FloatInfoFilter createFloatFilter(String field)
	{
		double t = 0;
		try
		{
			t = Double.parseDouble(txtThreshold.getText());
		}
		catch (NumberFormatException e)
		{
			showError("Please enter a valid float for threshold value!");
			txtThreshold.setText("0.00");
			return null;
		}
		FloatInfoFilter filter = new FloatInfoFilter(field, getFilterCompareMode(), t);
		int multi = getFilterMultiMode();
		switch (multi)
		{
		case VariantFilter.MULTIFIELD_ALL: filter.set_multi_AND(); break;
		case VariantFilter.MULTIFIELD_ANY: filter.set_multi_OR(); break;
		case VariantFilter.MULTIFIELD_NONE: filter.set_multi_NOT(); break;
		}
		return filter;
	}
	
	public StringInfoFilter createStringFilter(String field)
	{
		StringInfoFilter filter = new StringInfoFilter(field, getFilterCompareMode(), txtThreshold.getText());
		int multi = getFilterMultiMode();
		switch (multi)
		{
		case VariantFilter.MULTIFIELD_ALL: filter.set_multi_AND(); break;
		case VariantFilter.MULTIFIELD_ANY: filter.set_multi_OR(); break;
		case VariantFilter.MULTIFIELD_NONE: filter.set_multi_NOT(); break;
		}
		return filter;
	}
	
	public VariantFilter createNewFilter()
	{
		//Nab field and field type
		if (cmbxField.getSelectedIndex() < 0)
		{
			showError("Please enter a field to filter!");
			return null;
		}
		String field = (String)cmbxField.getSelectedItem();
		if (fieldType == VariantPool.INFODEF_UNK || fieldType == -2) updateFieldType(field);
		if (field == null || field.isEmpty())
		{
			showError("INFO field string is empty!");
			return null;
		}
		switch (fieldType)
		{
		case VariantPool.INFODEF_FLAG: return createBoolFilter(field);
		case VariantPool.INFODEF_FLOAT: return createFloatFilter(field);
		case VariantPool.INFODEF_INT: return createIntFilter(field);
		case VariantPool.INFODEF_STRING: return createStringFilter(field);
		case VariantPool.INFODEF_UNK: return createStringFilter(field);
		default:
			showError("Field type could not be recognized!\nPlease select a different INFO field!");
			break;
		}
		return null;
	}
	
	public void addFilterToList(VariantFilter myFilter)
	{
		ModelManager<VariantFilter> m = new ModelManager<VariantFilter>();
		List<VariantFilter> fList = m.listModelToList(lstCustom.getModel());
		fList.add(myFilter);
		Collections.sort(fList);
		lstCustom.setModel(m.collectionToListModel(fList));
	}
	
	public boolean SVTYPE_selected(SVType t)
	{
		switch (t)
		{
		case BED_REGION: return SVTYPES.isSelected(SVTYPE_INDEX_OTHER);
		case BND: return SVTYPES.isSelected(SVTYPE_INDEX_BND);
		case CNV: return SVTYPES.isSelected(SVTYPE_INDEX_CNV);
		case DEL: return SVTYPES.isSelected(SVTYPE_INDEX_DEL);
		case DELME: return SVTYPES.isSelected(SVTYPE_INDEX_DEL);
		case DUP: return SVTYPES.isSelected(SVTYPE_INDEX_DUP);
		case INS: return SVTYPES.isSelected(SVTYPE_INDEX_INS);
		case INSME: return SVTYPES.isSelected(SVTYPE_INDEX_INS);
		case INV: return SVTYPES.isSelected(SVTYPE_INDEX_INV);
		case OTHER: return SVTYPES.isSelected(SVTYPE_INDEX_OTHER);
		case TANDEM: return SVTYPES.isSelected(SVTYPE_INDEX_DUP);
		default: return false;
		}
	}
	
	public void updateManager()
	{
		//Custom filters
		manager.clearCustomFilters();
		ModelManager<VariantFilter> m = new ModelManager<VariantFilter>();
		List<VariantFilter> fList = m.listModelToList(lstCustom.getModel());
		for (VariantFilter f : fList) manager.addCustomFilter(f);
		
		//Chromosomes
		chromosomeLists.updateSourceLists();
		manager.clearIncludedChromosomes();
		Collection<Contig> cList = chromosomeLists.getSourceIncludeList();
		for (Contig c : cList) manager.includeChromosome(c);
		
		manager.setCompositeChrosomsomeExclusionMode(BNDMode.getSelectedIndex() == 1);
		
		//Confirmed
		manager.setConfirmedOnly(cbConfirmed.isSelected());
		
		//SV Types
		manager.clearIncludedSVTypes();
		for (SVType t : SVType.allTypes())
		{
			if (SVTYPE_selected(t)) manager.includeSVType(t);
		}
		
		//Var Types
		manager.setIncludeNonSV(cbSNVIn.isSelected());
		manager.setIncludeSV(cbSVIn.isSelected());
		
		//Quality
		manager.setFilterByQuality(cbQual.isSelected());
		manager.setExcludeNoQuality(cbEXqual.isSelected());
		double q = 0.0;
		try
		{
			q = Double.parseDouble(txtQual.getText());
		}
		catch (NumberFormatException e)
		{
			showError("Quality is invalid. Setting values to default...");
			q = 0.0;
		}
		manager.setMinimumQuality(q);
		
		//Lengths
		int min = 0;
		int max = Integer.MAX_VALUE;
		
		try
		{
			min = Integer.parseInt(this.txtMinSize.getText());
			max = Integer.parseInt(this.txtMaxSize.getText());
		}
		catch (NumberFormatException e)
		{
			showError("Minimum or maximum size is invalid. Setting values to default...");
			min = 0;
			max = Integer.MAX_VALUE;
		}
		if (min < 0) min = 0;
		if (max < 0) max = 0;
		if (min > max) max = min;
		manager.setMinimumSize(min);
		manager.setMaximumSize(max);
		
	}
	
	/* --- Action --- */
	
	public void addApplyListener(ActionListener l)
	{
		btnApply.addActionListener(l);
	}
	
	public void addCancelListener(ActionListener l)
	{
		btnCancel.addActionListener(l);
	}
	
	public void includeChrom()
	{
		chromosomeLists.includeSelected();
		repaintAll();
	}
	
	public void excludeChrom()
	{
		chromosomeLists.excludeSelected();
		repaintAll();
	}
	
	public void removeFilters()
	{
		ModelManager<VariantFilter> m = new ModelManager<VariantFilter>();
		List<VariantFilter> contents = m.listModelToList(lstCustom.getModel());
		List<VariantFilter> selected = lstCustom.getSelectedValuesList();
		contents.removeAll(selected);
		lstCustom.setModel(m.collectionToListModel(contents));
		repaintAll();
	}
	
	public void saveFilter()
	{
		VariantFilter f = createNewFilter();
		if (f == null)
		{
			showError("Filter could not be created!");
			return;
		}
		addFilterToList(f);
		clearFilterEditor();
		repaintAll();
	}
	
	public void selectFieldKey()
	{
		String key = (String)cmbxField.getSelectedItem();
		updateFieldType(key);
		updateFieldTypeLabel();
	}
	
	public void selectMultiMode(int index)
	{
		multiMode.select(index);
	}
	
	public void selectCompareMode(int index)
	{
		compareMode.select(index);
	}
	
	public void selectBNDMode(int index)
	{
		BNDMode.select(index);
		BNDMode.repaintAll();
	}
	
	/* --- Error --- */
	
	public void showError(String message)
	{
		JOptionPane.showMessageDialog(this, message, "Error", JOptionPane.ERROR_MESSAGE);
	}
}
