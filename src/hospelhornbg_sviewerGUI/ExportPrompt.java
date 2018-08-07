package hospelhornbg_sviewerGUI;

import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Collection;

import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JTextField;
import java.awt.Dimension;

import javax.swing.JButton;
import javax.swing.JTextPane;
import javax.swing.JScrollPane;
import javax.swing.border.BevelBorder;
import javax.swing.JComboBox;

public class ExportPrompt extends JDialog{

	private static final long serialVersionUID = -5944832251955645029L;
	private JTextField textField;
	private JTextPane textPane;
	private JComboBox<String> cmbxSample;
	
	public ExportPrompt(Frame parent, Collection<String> samples)
	{
		super(parent, true);
		setPreferredSize(new Dimension(450, 300));
		setMinimumSize(new Dimension(450, 300));
		setTitle("Export UCSC Genome Viewer Track");
		setResizable(false);
		getContentPane().setLayout(null);
		
		JLabel lblTrackName = new JLabel("Track Name:");
		lblTrackName.setBounds(10, 11, 75, 14);
		getContentPane().add(lblTrackName);
		
		textField = new JTextField();
		textField.setBounds(20, 31, 414, 20);
		getContentPane().add(textField);
		textField.setColumns(10);
		
		JButton btnOkay = new JButton("Okay");
		btnOkay.setBounds(333, 224, 101, 36);
		getContentPane().add(btnOkay);
		btnOkay.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				setVisible(false);
				dispose();
			}
			
		});
		
		JLabel lblDescription = new JLabel("Description:");
		lblDescription.setBounds(10, 62, 95, 14);
		getContentPane().add(lblDescription);
		
		JScrollPane scrollPane = new JScrollPane();
		scrollPane.setViewportBorder(new BevelBorder(BevelBorder.LOWERED, null, null, null, null));
		scrollPane.setBounds(20, 87, 414, 124);
		getContentPane().add(scrollPane);
		
		textPane = new JTextPane();
		scrollPane.setViewportView(textPane);
		
		JLabel lblSample = new JLabel("Sample:");
		lblSample.setBounds(10, 224, 46, 14);
		getContentPane().add(lblSample);
		
		cmbxSample = new JComboBox<String>();
		cmbxSample.setBounds(20, 240, 292, 20);
		getContentPane().add(cmbxSample);
		for (String sample : samples) cmbxSample.addItem(sample);
		cmbxSample.setSelectedIndex(0);
		
	}
	
	public String getNameField()
	{
		return this.textField.getText();
	}
	
	public String getDescriptionField()
	{
		return this.textPane.getText();
	}
	
	public String getSample()
	{
		return (String)cmbxSample.getSelectedItem();
	}
	
}
