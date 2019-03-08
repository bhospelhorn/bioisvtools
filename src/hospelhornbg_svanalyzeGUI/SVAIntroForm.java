package hospelhornbg_svanalyzeGUI;

import java.awt.Point;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JFrame;
import javax.swing.JSeparator;
import javax.swing.JButton;
import javax.swing.SwingConstants;

import hospelhornbg_svproject.CommonLoader;

import java.awt.Dimension;
import java.awt.Font;

public class SVAIntroForm extends JFrame{
	
	/* --- Constants --- */
	
	private static final long serialVersionUID = 8809439691696199613L;
	
	public static final int RETURN_VALUE_CANCEL = -1;
	public static final int RETURN_VALUE_LAUNCH = 0;
	public static final int RETURN_VALUE_INSTALL = 1;
	public static final int RETURN_VALUE_ADD_GENOME = 2;
	public static final int RETURN_VALUE_UPDATE_OMIM = 3;
	
	public static final int WIDTH = 250;
	public static final int HEIGHT = 250;
	
	/* --- Instance Variables --- */
	
	private int returnValue;
	
	/* --- Construction --- */
	
	public SVAIntroForm()
	{
		returnValue = RETURN_VALUE_CANCEL;
		//Determine if installed - ASSUMES that loading has already been attempted!
		initGUI(CommonLoader.getCommonDirPath() != null);
	}
	
	private void initGUI(boolean installed)
	{
		setTitle("SVAnalyze - Welcome!");
		setResizable(false);
		this.setMinimumSize(new Dimension(WIDTH, HEIGHT));
		this.setPreferredSize(new Dimension(WIDTH, HEIGHT));
		getContentPane().setLayout(null);
		
		JSeparator separator = new JSeparator();
		separator.setBounds(15, 123, 220, 2);
		getContentPane().add(separator);
		
		JButton btnLaunch = new JButton("Launch");
		btnLaunch.setFont(new Font("Tahoma", Font.BOLD, 12));
		btnLaunch.setToolTipText("Launch the SV viewer.");
		btnLaunch.setBounds(55, 26, 140, 50);
		getContentPane().add(btnLaunch);
		btnLaunch.setEnabled(installed);
		btnLaunch.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) 
			{
				returnValue = RETURN_VALUE_LAUNCH;
				setVisible(false);
				dispose();
			}
			
		});
		
		JButton btnInstall = new JButton("Install");
		btnInstall.setFont(new Font("Tahoma", Font.BOLD, 11));
		btnInstall.setToolTipText("Install directories for SVAnalyze");
		btnInstall.setBounds(55, 87, 140, 25);
		getContentPane().add(btnInstall);
		btnInstall.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) 
			{
				returnValue = RETURN_VALUE_INSTALL;
				setVisible(false);
				dispose();
			}
			
		});
		
		JButton btnGenomes = new JButton("Manage Genomes");
		btnGenomes.setFont(new Font("Tahoma", Font.PLAIN, 11));
		btnGenomes.setToolTipText("Manage the installed genome builds that can be used for projects.");
		btnGenomes.setBounds(55, 136, 140, 25);
		getContentPane().add(btnGenomes);
		btnGenomes.setEnabled(installed);
		btnGenomes.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) 
			{
				returnValue = RETURN_VALUE_ADD_GENOME;
				setVisible(false);
				dispose();
			}
			
		});
		
		JButton btnOmim = new JButton("Update OMIM");
		btnOmim.setFont(new Font("Tahoma", Font.PLAIN, 11));
		btnOmim.setToolTipText("Update the OMIM whitelist so variants falling in OMIM genes can be flagged.");
		btnOmim.setHorizontalTextPosition(SwingConstants.CENTER);
		btnOmim.setBounds(55, 172, 140, 25);
		getContentPane().add(btnOmim);
		btnOmim.setEnabled(installed);
		btnOmim.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) 
			{
				returnValue = RETURN_VALUE_UPDATE_OMIM;
				setVisible(false);
				dispose();
			}
			
		});
		
		
		Point centering = getCenteringCoordinates();
		setLocation(centering.x, centering.y);
	}
	
	public void render()
	{
		this.pack();
		this.setVisible(true);
	}
	
	private static Point getCenteringCoordinates()
	{
		Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
		
		int xc = screenSize.width / 2;
		int yc = screenSize.height / 2;
		
		int x = xc - (WIDTH/2);
		int y = yc - (HEIGHT/2);
		
		return new Point(x,y);
	}
	
	/* --- Getters --- */
	
	public int getReturnValue()
	{
		return returnValue;
	}
	
}
