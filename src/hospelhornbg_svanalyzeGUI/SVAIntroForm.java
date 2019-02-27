package hospelhornbg_svanalyzeGUI;

import java.awt.Point;

import javax.swing.JFrame;

public class SVAIntroForm extends JFrame{
	
	/* --- Constants --- */
	
	private static final long serialVersionUID = 8809439691696199613L;
	
	public static final int RETURN_VALUE_CANCEL = -1;
	public static final int RETURN_VALUE_LAUNCH = 0;
	public static final int RETURN_VALUE_INSTALL = 1;
	public static final int RETURN_VALUE_ADD_GENOME = 2;
	public static final int RETURN_VALUE_UPDATE_OMIM = 3;
	
	public static final int WIDTH = 500;
	public static final int HEIGHT = 500;
	
	/* --- Instance Variables --- */
	
	private int returnValue;
	
	/* --- Construction --- */
	
	public SVAIntroForm()
	{
		
	}
	
	private void initGUI()
	{
		
	}
	
	public void render()
	{
		this.pack();
		this.setVisible(true);
	}
	
	private static Point getCenteringCoordinates()
	{
		return null;
	}
	
	/* --- Getters --- */

}
