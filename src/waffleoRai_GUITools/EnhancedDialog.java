package waffleoRai_GUITools;

import java.awt.Cursor;
import java.awt.Frame;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import javax.swing.JDialog;
import javax.swing.JTextField;
import javax.swing.JToggleButton;

public class EnhancedDialog extends JDialog{

	private static final long serialVersionUID = 2338202784019770645L;
	
	private boolean waiting;
	
	protected Map<String, ComponentGroup> groups;
	protected Map<String, String> paths; //For last used paths in JFileChooser
	
	private Map<Integer, JToggleButton> boolMap;
	private Map<Integer, JTextField> intMap;
	private Map<Integer, JTextField> dblMap;
	private Map<Integer, JTextField> strMap;
	
	
	protected EnhancedDialog()
	{
		super();
		constructorCore();
	}
	
	protected EnhancedDialog(Frame owner, boolean modal)
	{
		super(owner, modal);
		constructorCore();
	}
	
	private void constructorCore()
	{
		waiting = false;
		this.groups = new HashMap<String, ComponentGroup>();
		this.paths = new HashMap<String, String>();
		boolMap = new HashMap<Integer, JToggleButton>();
		intMap = new HashMap<Integer, JTextField>();
		dblMap = new HashMap<Integer, JTextField>();
		strMap = new HashMap<Integer, JTextField>();
	}
	
	public void setGroupEnabled(String groupName, boolean enabled)
	{
		ComponentGroup g = groups.get(groupName);
		if (g == null) return;
		g.setEnabling(enabled);
	}
	
	public void repaintGroup(String groupName)
	{
		ComponentGroup g = groups.get(groupName);
		if (g == null) return;
		g.repaint();
	}
	
	public void disableAll()
	{
		Collection<ComponentGroup> allGroups = groups.values();
		for (ComponentGroup g : allGroups) g.setEnabling(false);
	}
	
	public void enableAll()
	{
		Collection<ComponentGroup> allGroups = groups.values();
		for (ComponentGroup g : allGroups) g.setEnabling(true);
	}
	
	public void repaintAll()
	{
		Collection<ComponentGroup> allGroups = groups.values();
		for (ComponentGroup g : allGroups) g.repaint();
	}
	
	public void restoreEnabled()
	{
		enableAll();
		repaintAll();
	}

	public void setWait()
	{
		if (waiting) return;
		waiting = true;
		setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
		disableAll();
	}
	
	public void unsetWait()
	{
		waiting = false;
		setCursor(null);
		restoreEnabled();
	}
	
	public void pairFieldAndCheckbox(int field, JToggleButton myCheckbox)
	{
		boolMap.put(field, myCheckbox);
	}
	
	public void pairFieldAndTextfield(int field, JTextField myTextbox)
	{
		strMap.put(field, myTextbox);
	}
	
	public void pairFieldAndNumericalTextfield(int field, JTextField myTextbox)
	{
		intMap.put(field, myTextbox);
	}
	
	public void pairFieldAndFloatTextfield(int field, JTextField myTextbox)
	{
		dblMap.put(field, myTextbox);
	}
	
	public JToggleButton getTiedCheckbox(int field)
	{
		return boolMap.get(field);
	}
	
	public JTextField getTiedTextField(int field)
	{
		return strMap.get(field);
	}
	
	public JTextField getTiedNumericalTextField(int field)
	{
		return intMap.get(field);
	}
	
	public JTextField getTiedFloatTextField(int field)
	{
		return dblMap.get(field);
	}

	public boolean getToggle_FieldValue(int field)
	{
		JToggleButton cb = getTiedCheckbox(field);
		if (cb == null) throw new IllegalArgumentException();
		return cb.isSelected();
	}
	
	public void setToggle_FieldValue(int field, boolean value)
	{
		JToggleButton cb = getTiedCheckbox(field);
		if (cb == null) throw new IllegalArgumentException();
		cb.setSelected(value);
	}
	
	public String getText_FieldValue(int field)
	{
		JTextField txt = getTiedTextField(field);
		if (txt == null) throw new IllegalArgumentException();
		return txt.getText();
	}
	
	public void setText_FieldValue(int field, String value)
	{
		JTextField txt = getTiedTextField(field);
		if (txt == null) throw new IllegalArgumentException();
		txt.setText(new String(value));
	}
	
	public int getNumerical_TextFieldValue(int field)
	{
		JTextField txt = getTiedNumericalTextField(field);
		if (txt == null) throw new IllegalArgumentException();
		try
		{
			int i = Integer.parseInt(txt.getText());
			return i;
		}
		catch (NumberFormatException ex)
		{
			txt.setText("0");
			return 0;
		}
	}
	
	public void setNumerical_TextFieldValue(int field, int value)
	{
		JTextField txt = getTiedNumericalTextField(field);
		if (txt == null) throw new IllegalArgumentException();
		txt.setText(Integer.toString(value));
	}
	
	public double getFloat_TextFieldValue(int field)
	{
		JTextField txt = getTiedFloatTextField(field);
		if (txt == null) throw new IllegalArgumentException();
		try
		{
			double f = Double.parseDouble(txt.getText());
			return f;
		}
		catch (NumberFormatException| NullPointerException ex)
		{
			txt.setText("0.0");
			return 0.0;
		}
	}
	
	public void setFloat_TextFieldValue(int field, double value)
	{
		JTextField txt = getTiedFloatTextField(field);
		if (txt == null) throw new IllegalArgumentException();
		txt.setText(Double.toString(value));
	}

	
}
