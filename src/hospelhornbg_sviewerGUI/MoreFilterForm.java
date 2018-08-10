package hospelhornbg_sviewerGUI;

import javax.swing.JDialog;
import javax.swing.JSeparator;
import javax.swing.JButton;
import javax.swing.JLabel;
import java.awt.Font;
import javax.swing.JCheckBox;
import javax.swing.JSlider;
import java.awt.Dimension;

public class MoreFilterForm extends JDialog{
	
	private static final long serialVersionUID = -477841198420014140L;

	public MoreFilterForm()
	{
		setPreferredSize(new Dimension(450, 420));
		setMinimumSize(new Dimension(450, 420));
		initGUI();
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
		
		JButton btnApply = new JButton("Apply");
		btnApply.setBounds(345, 354, 89, 23);
		getContentPane().add(btnApply);
		
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
		
		JCheckBox chckbxSplicing = new JCheckBox("Splicing");
		chckbxSplicing.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxSplicing.setBounds(20, 57, 97, 23);
		getContentPane().add(chckbxSplicing);
		
		JCheckBox chckbxUtr = new JCheckBox("UTR5");
		chckbxUtr.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxUtr.setBounds(20, 83, 97, 23);
		getContentPane().add(chckbxUtr);
		
		JCheckBox chckbxUtr_1 = new JCheckBox("UTR3");
		chckbxUtr_1.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxUtr_1.setBounds(127, 32, 97, 23);
		getContentPane().add(chckbxUtr_1);
		
		JCheckBox chckbxNcrna = new JCheckBox("ncRNA");
		chckbxNcrna.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxNcrna.setBounds(127, 57, 97, 23);
		getContentPane().add(chckbxNcrna);
		
		JCheckBox chckbxIntronic = new JCheckBox("Intronic");
		chckbxIntronic.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxIntronic.setBounds(127, 83, 97, 23);
		getContentPane().add(chckbxIntronic);
		
		JCheckBox chckbxIntergenic = new JCheckBox("Intergenic");
		chckbxIntergenic.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxIntergenic.setBounds(223, 32, 97, 23);
		getContentPane().add(chckbxIntergenic);
		
		JCheckBox chckbxHomozygousRecessive = new JCheckBox("Homozygous Recessive");
		chckbxHomozygousRecessive.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxHomozygousRecessive.setBounds(20, 148, 144, 23);
		getContentPane().add(chckbxHomozygousRecessive);
		
		JCheckBox chckbxDominant = new JCheckBox("Dominant");
		chckbxDominant.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxDominant.setBounds(20, 174, 71, 23);
		getContentPane().add(chckbxDominant);
		
		JCheckBox chckbxHalfhet = new JCheckBox("Half-Het (Paired)");
		chckbxHalfhet.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxHalfhet.setBounds(20, 200, 117, 23);
		getContentPane().add(chckbxHalfhet);
		
		JCheckBox chckbxMendelianViolation = new JCheckBox("Mendelian Violation");
		chckbxMendelianViolation.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxMendelianViolation.setBounds(185, 148, 117, 23);
		getContentPane().add(chckbxMendelianViolation);
		
		JCheckBox chckbxInheritanceUnresolved = new JCheckBox("Inheritance Unresolved");
		chckbxInheritanceUnresolved.setFont(new Font("Tahoma", Font.PLAIN, 11));
		chckbxInheritanceUnresolved.setBounds(185, 174, 144, 23);
		getContentPane().add(chckbxInheritanceUnresolved);
		
		JSlider sldMinMap = new JSlider();
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
		
	}
}
