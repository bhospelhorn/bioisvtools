package hospelhornbg_sviewerGUI;

import javax.swing.JDialog;
import javax.swing.JSeparator;
import javax.swing.JButton;
import javax.swing.JLabel;
import java.awt.Font;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Set;

import javax.swing.JCheckBox;
import javax.swing.JSlider;

import hospelhornbg_genomeBuild.GeneFunc;
import hospelhornbg_segregation.Inheritance;
import waffleoRai_GUITools.CheckBoxGroup;

import java.awt.Cursor;
import java.awt.Dimension;

public class MoreFilterForm extends JDialog{
	
	/* --- Constants --- */
	
	private static final long serialVersionUID = -477841198420014140L;
	
	public static final int POS_EXONIC = 0;
	public static final int POS_SPLICE = 1;
	public static final int POS_UTR5 = 2;
	public static final int POS_UTR3 = 3;
	public static final int POS_NCRNA = 4;
	public static final int POS_INTRONIC = 5;
	public static final int POS_UPSTREAM = 6;
	public static final int POS_DOWNSTREAM = 7;
	public static final int POS_INTERGENIC = 8;
	
	public static final int INH_HOMREC = 0;
	public static final int INH_DOM = 1;
	public static final int INH_HALFHET = 2;
	public static final int INH_MVIOL = 3;
	public static final int INH_DNS = 4;
	public static final int INH_HHUP = 5;
	
	/* --- Instance Variables --- */
	
	private ViewManager manager;
	
	private CheckBoxGroup cbgPosition;
	private CheckBoxGroup cbgSegregation;
	
	private JSlider sldMinMap;
	private JButton btnApply;
	
	/* --- Construction/Parsing --- */

	public MoreFilterForm(Frame parent, ViewManager vm)
	{
		super(parent, true);
		setPreferredSize(new Dimension(450, 420));
		setMinimumSize(new Dimension(450, 420));
		cbgPosition = new CheckBoxGroup(9);
		cbgSegregation = new CheckBoxGroup(6);
		manager = vm;
		initGUI();
		updateForm();
	}

	private void initGUI()
	{
		setResizable(false);
		setTitle("More Filters");
		getContentPane().setLayout(null);
		
		JSeparator separator = new JSeparator();
		separator.setBounds(10, 115, 424, 2);
		getContentPane().add(separator);
		
		JSeparator separator_1 = new JSeparator();
		separator_1.setBounds(10, 230, 424, 2);
		getContentPane().add(separator_1);
		
		JSeparator separator_2 = new JSeparator();
		separator_2.setBounds(10, 341, 424, 2);
		getContentPane().add(separator_2);
		
		btnApply = new JButton("Apply");
		btnApply.setBounds(345, 354, 89, 23);
		getContentPane().add(btnApply);
		btnApply.addActionListener(new ActionListener(){

			@Override
			public void actionPerformed(ActionEvent e) {
				apply();
				
			}
			
		});
		
		JLabel lblCodingEffect = new JLabel("Coding Effect");
		lblCodingEffect.setFont(new Font("Tahoma", Font.PLAIN, 13));
		lblCodingEffect.setBounds(10, 11, 82, 16);
		getContentPane().add(lblCodingEffect);
		
		JLabel lblSegregation = new JLabel("Segregation");
		lblSegregation.setFont(new Font("Tahoma", Font.PLAIN, 13));
		lblSegregation.setBounds(10, 125, 82, 16);
		getContentPane().add(lblSegregation);
		
		JLabel lblMappability = new JLabel("Mappability");
		lblMappability.setFont(new Font("Tahoma", Font.PLAIN, 13));
		lblMappability.setBounds(10, 242, 82, 16);
		getContentPane().add(lblMappability);
		
		JCheckBox chckbxExonic = new JCheckBox("Exonic");
		chckbxExonic.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxExonic.setBounds(20, 32, 97, 23);
		getContentPane().add(chckbxExonic);
		cbgPosition.addCheckBox(chckbxExonic, POS_EXONIC);
		
		JCheckBox chckbxSplicing = new JCheckBox("Splicing");
		chckbxSplicing.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxSplicing.setBounds(20, 57, 97, 23);
		getContentPane().add(chckbxSplicing);
		cbgPosition.addCheckBox(chckbxSplicing, POS_SPLICE);
		
		JCheckBox chckbxUtr = new JCheckBox("UTR5");
		chckbxUtr.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxUtr.setBounds(20, 83, 97, 23);
		getContentPane().add(chckbxUtr);
		cbgPosition.addCheckBox(chckbxUtr, POS_UTR5);
		
		JCheckBox chckbxUtr_1 = new JCheckBox("UTR3");
		chckbxUtr_1.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxUtr_1.setBounds(127, 32, 97, 23);
		getContentPane().add(chckbxUtr_1);
		cbgPosition.addCheckBox(chckbxUtr_1, POS_UTR3);
		
		JCheckBox chckbxNcrna = new JCheckBox("ncRNA");
		chckbxNcrna.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxNcrna.setBounds(127, 57, 97, 23);
		getContentPane().add(chckbxNcrna);
		cbgPosition.addCheckBox(chckbxNcrna, POS_NCRNA);
		
		JCheckBox chckbxIntronic = new JCheckBox("Intronic");
		chckbxIntronic.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxIntronic.setBounds(127, 83, 97, 23);
		getContentPane().add(chckbxIntronic);
		cbgPosition.addCheckBox(chckbxIntronic, POS_INTRONIC);
		
		JCheckBox chckbxIntergenic = new JCheckBox("Intergenic");
		chckbxIntergenic.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxIntergenic.setBounds(226, 83, 97, 23);
		getContentPane().add(chckbxIntergenic);
		cbgPosition.addCheckBox(chckbxIntergenic, POS_INTERGENIC);
		
		JCheckBox chckbxHomozygousRecessive = new JCheckBox("Homozygous Recessive");
		chckbxHomozygousRecessive.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxHomozygousRecessive.setBounds(20, 148, 144, 23);
		getContentPane().add(chckbxHomozygousRecessive);
		cbgSegregation.addCheckBox(chckbxHomozygousRecessive, INH_HOMREC);
		
		JCheckBox chckbxDominant = new JCheckBox("Dominant");
		chckbxDominant.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxDominant.setBounds(20, 174, 71, 23);
		getContentPane().add(chckbxDominant);
		cbgSegregation.addCheckBox(chckbxDominant, INH_DOM);
		
		JCheckBox chckbxHalfhet = new JCheckBox("Half-Het (Paired)");
		chckbxHalfhet.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxHalfhet.setBounds(20, 200, 117, 23);
		getContentPane().add(chckbxHalfhet);
		cbgSegregation.addCheckBox(chckbxHalfhet, INH_HALFHET);
		
		JCheckBox chckbxMendelianViolation = new JCheckBox("Mendelian Violation");
		chckbxMendelianViolation.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxMendelianViolation.setBounds(185, 148, 117, 23);
		getContentPane().add(chckbxMendelianViolation);
		cbgSegregation.addCheckBox(chckbxMendelianViolation, INH_MVIOL);
		
		JCheckBox chckbxInheritanceUnresolved = new JCheckBox("Inheritance Unresolved");
		chckbxInheritanceUnresolved.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxInheritanceUnresolved.setBounds(185, 174, 144, 23);
		getContentPane().add(chckbxInheritanceUnresolved);
		cbgSegregation.addCheckBox(chckbxInheritanceUnresolved, INH_DNS);
		
		sldMinMap = new JSlider();
		sldMinMap.setBounds(10, 269, 200, 26);
		getContentPane().add(sldMinMap);
		
		JLabel lblNewLabel = new JLabel("Minimum Mappability:");
		lblNewLabel.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblNewLabel.setBounds(10, 306, 107, 14);
		getContentPane().add(lblNewLabel);
		
		JLabel lblMinMap = new JLabel("0");
		lblMinMap.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblMinMap.setBounds(127, 306, 46, 14);
		getContentPane().add(lblMinMap);
		
		JCheckBox chckbxHalfhetunpaired = new JCheckBox("Half-Het (Unpaired)");
		chckbxHalfhetunpaired.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxHalfhetunpaired.setBounds(185, 200, 135, 23);
		getContentPane().add(chckbxHalfhetunpaired);
		cbgSegregation.addCheckBox(chckbxHalfhetunpaired, INH_HHUP);
		
		JCheckBox chckbxUpstream = new JCheckBox("Upstream");
		chckbxUpstream.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxUpstream.setBounds(226, 32, 97, 23);
		getContentPane().add(chckbxUpstream);
		cbgPosition.addCheckBox(chckbxUpstream, POS_UPSTREAM);
		
		JCheckBox chckbxDownstream = new JCheckBox("Downstream");
		chckbxDownstream.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxDownstream.setBounds(226, 57, 97, 23);
		getContentPane().add(chckbxDownstream);
		cbgPosition.addCheckBox(chckbxDownstream, POS_DOWNSTREAM);
		
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
		enableAll();
	}
	
	public void disableAll()
	{
		cbgPosition.disableAll();
		cbgSegregation.disableAll();
		btnApply.setEnabled(false);
		sldMinMap.setEnabled(false);
	}
	
	public void enableAll()
	{
		cbgPosition.enableAll();
		cbgSegregation.enableAll();
		btnApply.setEnabled(true);
		sldMinMap.setEnabled(true);
	}
	
	/* --- Syncing --- */
	
	public void updateManager()
	{
		//Update position effect filters
		manager.clearCodingEffects();
		if (cbgPosition.isSelected(POS_EXONIC)) manager.addCodingEffect(GeneFunc.EXONIC);
		if (cbgPosition.isSelected(POS_SPLICE)) manager.addCodingEffect(GeneFunc.SPLICING);
		if (cbgPosition.isSelected(POS_UTR5)) manager.addCodingEffect(GeneFunc.UTR5);
		if (cbgPosition.isSelected(POS_UTR3)) manager.addCodingEffect(GeneFunc.UTR3);
		if (cbgPosition.isSelected(POS_NCRNA)) manager.addCodingEffect(GeneFunc.NCRNA);
		if (cbgPosition.isSelected(POS_INTRONIC)) manager.addCodingEffect(GeneFunc.INTRONIC);
		if (cbgPosition.isSelected(POS_UPSTREAM)) manager.addCodingEffect(GeneFunc.UPSTREAM);
		if (cbgPosition.isSelected(POS_DOWNSTREAM)) manager.addCodingEffect(GeneFunc.DOWNSTREAM);
		if (cbgPosition.isSelected(POS_INTERGENIC)) manager.addCodingEffect(GeneFunc.INTERGENIC);
		
		//Update inheritance pattern filters
		manager.clearInheritancePatterns();
		if (cbgSegregation.isSelected(INH_HOMREC)) manager.addInheritancePattern(Inheritance.HOM_REC);
		if (cbgSegregation.isSelected(INH_DOM)) manager.addInheritancePattern(Inheritance.DOMINANT);
		if (cbgSegregation.isSelected(INH_HALFHET)) manager.addInheritancePattern(Inheritance.HALF_HET);
		if (cbgSegregation.isSelected(INH_HALFHET)) manager.addInheritancePattern(Inheritance.HALF_HET_SV);
		if (cbgSegregation.isSelected(INH_HALFHET)) manager.addInheritancePattern(Inheritance.COMP_HET);
		if (cbgSegregation.isSelected(INH_HALFHET)) manager.addInheritancePattern(Inheritance.COMP_HET_SV);
		if (cbgSegregation.isSelected(INH_MVIOL)) manager.addInheritancePattern(Inheritance.DENOVO_DOM);
		if (cbgSegregation.isSelected(INH_MVIOL)) manager.addInheritancePattern(Inheritance.DENOVO_HET);
		if (cbgSegregation.isSelected(INH_MVIOL)) manager.addInheritancePattern(Inheritance.DENOVO_HET_SV);
		if (cbgSegregation.isSelected(INH_MVIOL)) manager.addInheritancePattern(Inheritance.DENOVO_REC);
		if (cbgSegregation.isSelected(INH_DNS)) manager.addInheritancePattern(Inheritance.UNRESOLVED);
		manager.setIncludeUnpairedHalfHets(cbgSegregation.isSelected(INH_HHUP));
		
		//Update mappability filter 
		//TODO: Implement when mappability files are ready
	}
	
	public void updateForm()
	{
		//Updates form to match manager. Needed when form is first loaded.
		if (manager == null) throw new NullPointerException();
		
		//Update positions
		Set<GeneFunc> ce = manager.getPassedCodingEffects();
		cbgPosition.select(POS_EXONIC, ce.contains(GeneFunc.EXONIC));
		cbgPosition.select(POS_SPLICE, ce.contains(GeneFunc.SPLICING));
		cbgPosition.select(POS_UTR5, ce.contains(GeneFunc.UTR5));
		cbgPosition.select(POS_UTR3, ce.contains(GeneFunc.UTR3));
		cbgPosition.select(POS_NCRNA, ce.contains(GeneFunc.NCRNA));
		cbgPosition.select(POS_INTRONIC, ce.contains(GeneFunc.INTRONIC));
		cbgPosition.select(POS_UPSTREAM, ce.contains(GeneFunc.UPSTREAM));
		cbgPosition.select(POS_DOWNSTREAM, ce.contains(GeneFunc.DOWNSTREAM));
		cbgPosition.select(POS_INTERGENIC, ce.contains(GeneFunc.INTERGENIC));
		
		//Update inheritance patterns
		Set<Inheritance> ipset = manager.getPassedInheritancePatterns();
		cbgSegregation.select(INH_HOMREC, ipset.contains(Inheritance.HOM_REC));
		cbgSegregation.select(INH_DOM, ipset.contains(Inheritance.DOMINANT));
		cbgSegregation.select(INH_HALFHET, ipset.contains(Inheritance.HALF_HET) || ipset.contains(Inheritance.HALF_HET_SV) || ipset.contains(Inheritance.COMP_HET) || ipset.contains(Inheritance.COMP_HET_SV));
		cbgSegregation.select(INH_MVIOL, ipset.contains(Inheritance.DENOVO_DOM) || ipset.contains(Inheritance.DENOVO_HET) || ipset.contains(Inheritance.DENOVO_HET_SV) || ipset.contains(Inheritance.DENOVO_REC));
		cbgSegregation.select(INH_DNS, ipset.contains(Inheritance.UNRESOLVED));
		cbgSegregation.select(INH_HHUP, manager.includeUnpairedHalfHets());
		
		//Update mappability slider
		//TODO: Implement when mappability files are ready
		
		cbgPosition.repaintAll();
		cbgSegregation.repaintAll();
	}
	
	/* --- Action --- */
	
	public void apply()
	{
		setWait();
		updateManager();
		unsetWait();
		setVisible(false); //Close?
	}
	
}
