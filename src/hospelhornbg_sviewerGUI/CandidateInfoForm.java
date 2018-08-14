package hospelhornbg_sviewerGUI;

import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.border.BevelBorder;

import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_genomeBuild.Gene;
import hospelhornbg_segregation.Candidate;
import hospelhornbg_segregation.Individual;
import hospelhornbg_segregation.Pedigree;

import javax.swing.JSeparator;
import java.awt.Font;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class CandidateInfoForm extends JDialog{
	
	private static final long serialVersionUID = -2491579039576303676L;
	
	private Candidate iCandidate;
	private Pedigree iFamily;
	private Map<Integer, Candidate> iPartnerLink;
	
	private JTable tblPartners;
	private JTable tblGeno;
	
	public CandidateInfoForm(Candidate c, Pedigree f) 
	{
		iPartnerLink = new HashMap<Integer, Candidate>();
		if (c == null) throw new IllegalArgumentException();
		if (f == null) throw new IllegalArgumentException();
		iCandidate = c;
		iFamily = f;
		initGUI();
	}
	
	private void initGUI()
	{
		setResizable(false);
		setTitle("Candidate Details");
		getContentPane().setLayout(null);
		
		JLabel lblId = new JLabel("ID:");
		lblId.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblId.setBounds(10, 11, 15, 14);
		getContentPane().add(lblId);
		
		JLabel lblVarID = new JLabel("[ID]");
		lblVarID.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblVarID.setBounds(35, 11, 164, 14);
		getContentPane().add(lblVarID);
		
		JLabel lblSvType = new JLabel("SV Type:");
		lblSvType.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblSvType.setBounds(10, 31, 46, 14);
		getContentPane().add(lblSvType);
		
		JLabel lblStartPosition = new JLabel("Start Position:");
		lblStartPosition.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblStartPosition.setBounds(10, 55, 70, 14);
		getContentPane().add(lblStartPosition);
		
		JLabel lblEndPosition = new JLabel("End Position:");
		lblEndPosition.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblEndPosition.setBounds(10, 120, 70, 14);
		getContentPane().add(lblEndPosition);
		
		JLabel lblGene = new JLabel("Gene:");
		lblGene.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblGene.setBounds(230, 10, 35, 14);
		getContentPane().add(lblGene);
		
		JLabel lblVarType = new JLabel("[Type]");
		lblVarType.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblVarType.setBounds(59, 31, 140, 14);
		getContentPane().add(lblVarType);
		
		JLabel lblStart90 = new JLabel("[CI90]");
		lblStart90.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblStart90.setBounds(90, 95, 122, 14);
		getContentPane().add(lblStart90);
		
		JLabel lblStart95 = new JLabel("[CI95]");
		lblStart95.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblStart95.setBounds(90, 75, 122, 14);
		getContentPane().add(lblStart95);
		
		JLabel lblStartCall = new JLabel("[Call]");
		lblStartCall.setBounds(90, 55, 122, 14);
		getContentPane().add(lblStartCall);
		
		JLabel lblGeneName = new JLabel("[Gene]");
		lblGeneName.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblGeneName.setBounds(267, 10, 136, 14);
		getContentPane().add(lblGeneName);
		
		JLabel lblCodingEffect = new JLabel("Coding Effect:");
		lblCodingEffect.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblCodingEffect.setBounds(230, 50, 77, 14);
		getContentPane().add(lblCodingEffect);
		
		JLabel lblInheritance = new JLabel("Inheritance:");
		lblInheritance.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblInheritance.setBounds(230, 70, 59, 14);
		getContentPane().add(lblInheritance);
		
		JLabel lblComphetPartners = new JLabel("Comp-Het Partners");
		lblComphetPartners.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblComphetPartners.setBounds(10, 187, 95, 14);
		getContentPane().add(lblComphetPartners);
		
		JScrollPane spPartners = new JScrollPane();
		spPartners.setViewportBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
		spPartners.setBounds(10, 208, 393, 77);
		getContentPane().add(spPartners);
		
		tblPartners = new JTable();
		tblPartners.setFillsViewportHeight(true);
		spPartners.setViewportView(tblPartners);
		
		JLabel lblGenotypeCalls = new JLabel("Genotype Calls");
		lblGenotypeCalls.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblGenotypeCalls.setBounds(10, 296, 77, 14);
		getContentPane().add(lblGenotypeCalls);
		
		JScrollPane spGeno = new JScrollPane();
		spGeno.setViewportBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
		spGeno.setBounds(10, 316, 393, 98);
		getContentPane().add(spGeno);
		
		tblGeno = new JTable();
		tblGeno.setFillsViewportHeight(true);
		spGeno.setViewportView(tblGeno);
		
		JLabel lblSize = new JLabel("Size:");
		lblSize.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblSize.setBounds(230, 103, 23, 14);
		getContentPane().add(lblSize);
		
		JLabel lblEndCall = new JLabel("[Call]");
		lblEndCall.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblEndCall.setBounds(90, 120, 122, 14);
		getContentPane().add(lblEndCall);
		
		JLabel lblEnd95 = new JLabel("[CI95]");
		lblEnd95.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblEnd95.setBounds(90, 140, 122, 14);
		getContentPane().add(lblEnd95);
		
		JLabel lblEnd90 = new JLabel("[CI90]");
		lblEnd90.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblEnd90.setBounds(90, 160, 122, 14);
		getContentPane().add(lblEnd90);
		
		JLabel lblEffect = new JLabel("[Effect]");
		lblEffect.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblEffect.setBounds(305, 50, 98, 14);
		getContentPane().add(lblEffect);
		
		JLabel lblIP = new JLabel("[ip]");
		lblIP.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblIP.setBounds(294, 70, 98, 14);
		getContentPane().add(lblIP);
		
		JSeparator separator = new JSeparator();
		separator.setBounds(5, 50, 200, 2);
		getContentPane().add(separator);
		
		JSeparator separator_1 = new JSeparator();
		separator_1.setBounds(5, 115, 200, 2);
		getContentPane().add(separator_1);
		
		JSeparator separator_2 = new JSeparator();
		separator_2.setBounds(5, 177, 200, 2);
		getContentPane().add(separator_2);
		
		JLabel lblSVSize = new JLabel("[sv size]");
		lblSVSize.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblSVSize.setBounds(261, 103, 131, 14);
		getContentPane().add(lblSVSize);
		
		JLabel lblTranscript = new JLabel("Transcript:");
		lblTranscript.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblTranscript.setBounds(230, 30, 52, 14);
		getContentPane().add(lblTranscript);
		
		JLabel lblTranscriptName = new JLabel("[transcript]");
		lblTranscriptName.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblTranscriptName.setBounds(287, 30, 116, 14);
		getContentPane().add(lblTranscriptName);
	}
	
	public void render()
	{
		this.pack();
		this.setVisible(true);
	}
	
	/* --- Partner Table --- */
	
	private String[] generatePartnerTableColumnHeaders()
	{
		String[] sarr = {"Transcript", "Gene", "Type", "Chromosome", "Position", "End", "Inheritance Pattern"};
		return sarr;
	}
	
	private String[][] generateAndLinkPartnerTable()
	{
		iPartnerLink.clear();
		List<Candidate> plist = iCandidate.getAllPartners();
		if (plist == null) return null;
		if (plist.isEmpty()) return null;
		int sz = plist.size();
		String[][] tbl = new String[sz][7];
		int i = 0;
		for (Candidate p : plist)
		{
			iPartnerLink.put(i, p);
			Gene g = p.getGene();
			tbl[i][0] = g.getID();
			tbl[i][1] = g.getName();
			Variant v = p.getVariant();
			if (v instanceof StructuralVariant)
			{
				StructuralVariant sv = (StructuralVariant)v;
				tbl[i][2] = sv.getType().getString();
				tbl[i][3] = sv.getChromosome().getUDPName();
				tbl[i][4] = Integer.toString(sv.getPosition());
				tbl[i][5] = Integer.toString(sv.getEndPosition());
				tbl[i][6] = p.getInheritancePattern().toString();
			}
			else
			{
				tbl[i][2] = "SNV or Small InDel";
				tbl[i][3] = v.getChromosome().getUDPName();
				tbl[i][4] = Integer.toString(v.getPosition());
				tbl[i][5] = "";
				tbl[i][6] = p.getInheritancePattern().toString();
			}
			i++;
		}
		return tbl;
	}
	
	private void selectTableRecord(int recIndex)
	{
		
	}
	
	/* --- Genotype Table --- */
	
	private String[] generateGenoTableColumnHeaders()
	{
		String[] sarr = {"Sample Name", "Relation", "Affected", "Sex", "Genotype"};
		return sarr;
	}
	
	private String[][] generateGenoTable()
	{
		List<Individual> ilist = iFamily.getAllMembers();
		if (ilist == null) return null;
		if (ilist.isEmpty()) return null;
		int inum = ilist.size();
		return null;
	}
	
	

}
