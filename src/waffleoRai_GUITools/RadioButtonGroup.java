package waffleoRai_GUITools;

import javax.swing.JRadioButton;

public class RadioButtonGroup {
	
	private JRadioButton[] buttons;
	private int selected;
	
	public RadioButtonGroup(int buttonNumber)
	{
		if (buttonNumber <= 0) buttonNumber = 16;
		buttons = new JRadioButton[buttonNumber];
	}
	
	public void addButton(JRadioButton button, int index)
	{
		if (index < 0 || index >= buttons.length) throw new IndexOutOfBoundsException();
		buttons[index] = button;
		selected = 0;
	}
	
	public void select(int index)
	{
		if (index < 0 || index >= buttons.length) return;
		selected = index;
		for (int i = 0; i < buttons.length; i++)
		{
			if (buttons[i] != null)
			{
				if (!buttons[i].isEnabled())
				{
					buttons[i].setSelected(false);
					continue;
				}
				buttons[i].setSelected(i == index);	
			}
		}
	}
	
	public void select(JRadioButton button)
	{
		if (button == null) return;
		if (!button.isEnabled()) return;
		boolean found = false;
		for (int i = 0; i < buttons.length; i++)
		{
			if (buttons[i] == null) continue;
			if (button == buttons[i] && !found) //Yes, they must be identical REFERENCES
			{
				selected = i;
				found = true;
				buttons[i].setSelected(true);
			}
			else buttons[i].setSelected(false);
		}
	}
	
	public int getSelectedIndex()
	{
		return selected;
	}
	
	public JRadioButton getSelectedButton()
	{
		return buttons[selected];
	}

	private int firstEnabledButton()
	{
		for (int i = 0; i < buttons.length; i++)
		{
			if (buttons[i] != null)
			{
				if (buttons[i].isSelected()) return i;
			}
		}
		return -1;
	}
	
	public void disable(int index)
	{
		if (index < 0 || index >= buttons.length) return;
		if (buttons[index] != null)
		{
			buttons[index].setEnabled(false);
			if (selected == index) selected = firstEnabledButton();
		}
	}
	
	public void enable(int index)
	{
		if (index < 0 || index >= buttons.length) return;
		if (buttons[index] != null)
		{
			buttons[index].setEnabled(true);
		}
	}
	
	public void disableAll()
	{
		setEnabledAll(false);
	}
	
	public void enableAll()
	{
		setEnabledAll(true);
	}
	
	public void setEnabledAll(boolean b)
	{
		for (JRadioButton rb : buttons) rb.setEnabled(b);
	}
	
	public void repaintAll()
	{
		for (JRadioButton rb : buttons) rb.repaint();
	}
	
}
