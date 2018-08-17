package hospelhornbg_sviewerGUI;

import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.border.BevelBorder;
import javax.swing.table.DefaultTableModel;

import hospelhornbg_bioinformatics.AffectedStatus;
import hospelhornbg_bioinformatics.Genotype;
import hospelhornbg_bioinformatics.Sex;
import hospelhornbg_bioinformatics.StructuralVariant;
import hospelhornbg_bioinformatics.Variant;
import hospelhornbg_genomeBuild.Gene;
import hospelhornbg_genomeBuild.GeneFunc;
import hospelhornbg_segregation.Candidate;
import hospelhornbg_segregation.Individual;
import hospelhornbg_segregation.Inheritance;
import hospelhornbg_segregation.Pedigree;
import hospelhornbg_segregation.Relationship;

import javax.swing.JSeparator;
import java.awt.Font;
import java.awt.Frame;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.awt.Dimension;
import java.awt.Color;

public class CandidateInfoForm extends JDialog{
	
	private static final long serialVersionUID = -2491579039576303676L;
	
	private Candidate iCandidate;
	private Pedigree iFamily;
	private Map<Integer, Candidate> iPartnerLink;
	
	private JScrollPane spPartners;
	private JScrollPane spGeno;
	private JTable tblPartners;
	private JTable tblGeno;
		
	/**
	 * @wbp.parser.constructor
	 */
	public CandidateInfoForm(Candidate c, Pedigree f) 
	{
		setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
		//iDescendants = new LinkedList<CandidateInfoForm>();
		setPreferredSize(new Dimension(420, 455));
		setMinimumSize(new Dimension(420, 455));
		iPartnerLink = new HashMap<Integer, Candidate>();
		if (c == null) throw new IllegalArgumentException();
		if (f == null) throw new IllegalArgumentException();
		iCandidate = c;
		iFamily = f;
		initGUI();
		updatePartnerTable();
		updateGenotypeTable();
	}
	
	public CandidateInfoForm(Candidate c, Pedigree f, Frame parent) 
	{
		super(parent, true);
		//iDescendants = new LinkedList<CandidateInfoForm>();
		setPreferredSize(new Dimension(420, 455));
		setMinimumSize(new Dimension(420, 455));
		iPartnerLink = new HashMap<Integer, Candidate>();
		if (c == null) throw new IllegalArgumentException();
		if (f == null) throw new IllegalArgumentException();
		iCandidate = c;
		iFamily = f;
		initGUI();
		updatePartnerTable();
		updateGenotypeTable();
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
		
		Variant v = iCandidate.getVariant();
		
		JLabel lblVarID = new JLabel("[ID]");
		if (v != null) lblVarID.setText(v.getVarID());
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
		
		StructuralVariant sv = null;
		if (v instanceof StructuralVariant) sv = (StructuralVariant)v;
		
		JLabel lblVarType = new JLabel("[Type]");
		if (sv != null) lblVarType.setText(sv.getType().getString());
		lblVarType.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblVarType.setBounds(59, 31, 140, 14);
		getContentPane().add(lblVarType);
		
		JLabel lblStart90 = new JLabel("[CI90]");
		lblStart90.setForeground(Color.BLUE);
		if (sv != null){
			int st = sv.getCIPosition(false, false, false);
			int ed = sv.getCIPosition(false, false, true);
			lblStart90.setText(sv.getChromosome().getUDPName() + ":" + st + "-" + ed);
		}
		lblStart90.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblStart90.setBounds(90, 95, 122, 14);
		getContentPane().add(lblStart90);
		
		JLabel lblStart95 = new JLabel("[CI95]");
		lblStart95.setForeground(Color.GREEN);
		if (sv != null){
			int st = sv.getCIPosition(false, true, false);
			int ed = sv.getCIPosition(false, true, true);
			lblStart95.setText(sv.getChromosome().getUDPName() + ":" + st + "-" + ed);
		}
		lblStart95.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblStart95.setBounds(90, 75, 122, 14);
		getContentPane().add(lblStart95);
		
		JLabel lblStartCall = new JLabel("[Call]");
		if (sv != null){
			int st = sv.getPosition();
			lblStartCall.setText(sv.getChromosome().getUDPName() + ":" + st);
		}
		lblStartCall.setBounds(90, 55, 122, 14);
		getContentPane().add(lblStartCall);
		
		Gene gene = iCandidate.getGene();
		JLabel lblGeneName = new JLabel("[Gene]");
		if (gene != null) lblGeneName.setText(gene.getName());
		else lblGeneName.setText("[None]");
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
		
		spPartners = new JScrollPane();
		spPartners.setViewportBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
		spPartners.setBounds(10, 208, 393, 77);
		getContentPane().add(spPartners);
		
		tblPartners = new JTable();
		tblPartners.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		tblPartners.setFillsViewportHeight(true);
		spPartners.setViewportView(tblPartners);
		tblPartners.addMouseListener(new MouseAdapter(){

			private int lastr = -1;
			private int lastl = -1;
			
			@Override
			public void mouseClicked(MouseEvent e) {
				int r = tblPartners.getSelectedRow();
				int l = tblPartners.getSelectedColumn();
				if (lastr == r && lastl == l)
				{
					lastr = -1;
					lastl = -1;
					selectTableRecord();
				}
				else
				{
					lastr = r;
					lastl = l;
				}
			}
			
		});
		
		JLabel lblGenotypeCalls = new JLabel("Genotype Calls");
		lblGenotypeCalls.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblGenotypeCalls.setBounds(10, 296, 77, 14);
		getContentPane().add(lblGenotypeCalls);
		
		spGeno = new JScrollPane();
		spGeno.setViewportBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
		spGeno.setBounds(10, 316, 393, 98);
		getContentPane().add(spGeno);
		
		tblGeno = new JTable();
		tblGeno.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		tblGeno.setFillsViewportHeight(true);
		spGeno.setViewportView(tblGeno);
		
		JLabel lblSize = new JLabel("Size:");
		lblSize.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblSize.setBounds(230, 103, 23, 14);
		getContentPane().add(lblSize);
		
		JLabel lblEndCall = new JLabel("[Call]");
		if (sv != null){
			int ed = sv.getEndPosition();
			lblEndCall.setText(sv.getEndChromosome().getUDPName() + ":" + ed);
		}
		lblEndCall.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblEndCall.setBounds(90, 120, 122, 14);
		getContentPane().add(lblEndCall);
		
		JLabel lblEnd95 = new JLabel("[CI95]");
		if (sv != null){
			int st = sv.getCIPosition(true, true, false);
			int ed = sv.getCIPosition(true, true, true);
			lblEnd95.setText(sv.getEndChromosome().getUDPName() + ":" + st + "-" + ed);
		}
		lblEnd95.setForeground(Color.GREEN);
		lblEnd95.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblEnd95.setBounds(90, 140, 122, 14);
		getContentPane().add(lblEnd95);
		
		JLabel lblEnd90 = new JLabel("[CI90]");
		if (sv != null){
			int st = sv.getCIPosition(true, false, false);
			int ed = sv.getCIPosition(true, false, true);
			lblEnd90.setText(sv.getEndChromosome().getUDPName() + ":" + st + "-" + ed);
		}
		lblEnd90.setForeground(Color.BLUE);
		lblEnd90.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblEnd90.setBounds(90, 160, 122, 14);
		getContentPane().add(lblEnd90);
		
		GeneFunc func = null;
		if (v != null) func = v.getGeneFunction();
		JLabel lblEffect = new JLabel("[Effect]");
		if (func != null) lblEffect.setText(func.toString());
		lblEffect.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblEffect.setBounds(305, 50, 98, 14);
		getContentPane().add(lblEffect);
		
		Inheritance ip = iCandidate.getInheritancePattern();
		JLabel lblIP = new JLabel("[ip]");
		if (ip != null)
		{
			lblIP.setText(ip.toString());
		}
		else lblIP.setText("[Unknown]");
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
		if (sv != null)
		{
			int svlen = sv.getAbsoluteSVLength();
			lblSVSize.setText(Integer.toString(svlen) + "bp");
		}
		lblSVSize.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblSVSize.setBounds(261, 103, 131, 14);
		getContentPane().add(lblSVSize);
		
		JLabel lblTranscript = new JLabel("Transcript:");
		lblTranscript.setFont(new Font("Tahoma", Font.PLAIN, 11));
		lblTranscript.setBounds(230, 30, 52, 14);
		getContentPane().add(lblTranscript);
		
		JLabel lblTranscriptName = new JLabel("[transcript]");
		if (gene != null) lblTranscriptName.setText(gene.getID());
		else lblTranscriptName.setText("[None]");
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
	
	private void selectTableRecord()
	{
		if (!tblPartners.isEnabled()) return;
		int i = tblPartners.getSelectedRow();
		Candidate p = iPartnerLink.get(i);
		if (p == null) return;
		
		CandidateInfoForm pdialog = new CandidateInfoForm(p, iFamily);
		pdialog.setLocationRelativeTo(this);
		pdialog.render();
	}
	
	public void updatePartnerTable()
	{
		if (iCandidate == null) return;
		if (!iCandidate.hasPartners())
		{
			tblPartners.setModel(new DefaultTableModel());
			tblPartners.setEnabled(false);
			spPartners.setEnabled(false);
			tblPartners.repaint();
			spPartners.repaint();
			return;
		}
		
		String[] headers = generatePartnerTableColumnHeaders();
		String[][] data = generateAndLinkPartnerTable();
		tblPartners.setModel(new DefaultTableModel(data, headers));
		tblPartners.setEnabled(true);
		spPartners.setEnabled(true);
		tblPartners.repaint();
		spPartners.repaint();

	}
	
	/* --- Genotype Table --- */
	
	private String[] generateGenoTableColumnHeaders()
	{
		String[] sarr = {"Sample Name", "Relation", "Affected", "Sex", "Genotype", "Copy Number"};
		return sarr;
	}
	
	private String[][] generateGenoTable()
	{
		List<Individual> ilist = iFamily.getAllMembers();
		if (ilist == null) return null;
		if (ilist.isEmpty()) return null;
		Collections.sort(ilist);
		Individual proband = iFamily.getProband();
		int inum = ilist.size();
		String[][] tbl = new String[inum][6];
		int i = 0;
		for (Individual indiv : ilist)
		{
			tbl[i][0] = indiv.getName();
			if (proband != null)
			{
				Relationship r = proband.getRelationship(indiv);
				if (r != null) tbl[i][1] = r.toString_English();
				else tbl[i][1] = "Unrelated";
			}
			else tbl[i][1] = "Unrelated";
			AffectedStatus aff = indiv.getAffectedStatus();
			if (aff == null) tbl[i][2] = "UNK";
			else
			{
				switch(aff)
				{
				case AFFECTED:
					tbl[i][2] = "YES";
					break;
				case PARTIALLY_AFFECTED:
					tbl[i][2] = "PAR";
					break;
				case POSSIBLY_AFFECTED:
					tbl[i][2] = "PBL";
					break;
				case UNAFFECTED:
					tbl[i][2] = "NO";
					break;
				case UNKNOWN:
					tbl[i][2] = "UNK";
					break;
				default:
					tbl[i][2] = "UNK";
					break;
				}	
			}
			Sex sx = indiv.getSex();
			if (aff == null) tbl[i][3] = "U";
			else
			{
				switch(sx)
				{
				case FEMALE:
					tbl[i][3] = "F";
					break;
				case MALE:
					tbl[i][3] = "M";
					break;
				case OTHER:
					tbl[i][3] = "O";
					break;
				case UNKNOWN:
					tbl[i][3] = "U";
					break;
				default:
					tbl[i][3] = "U";
					break;
				}
			}
			//Genotype
			Variant v = iCandidate.getVariant();
			if (v == null){
				tbl[i][4] = "./.";
				continue;
			}
			Genotype g = v.getSampleGenotype(indiv.getName());
			if (g == null){
				tbl[i][4] = "./.";
				continue;
			}
			//int[] alleles = g.getAlleles();
			int cn = g.getCopyNumber();
			String astr = g.getField(Genotype.INFODEF_GT.getKey());
			tbl[i][4] = astr;
			tbl[i][5] = Integer.toString(cn);
		}
		
		return tbl;
	}
	
	public void updateGenotypeTable()
	{
		if (iFamily == null){
			tblGeno.setModel(new DefaultTableModel());
			tblGeno.setEnabled(false);
			spGeno.setEnabled(false);
			tblGeno.repaint();
			spGeno.repaint();
			return;
		}
		
		String[] headers = generateGenoTableColumnHeaders();
		String[][] data = generateGenoTable();
		tblGeno.setModel(new DefaultTableModel(data, headers));
		tblGeno.setEnabled(true);
		spGeno.setEnabled(true);
		tblGeno.repaint();
		spGeno.repaint();
	}

	
}
