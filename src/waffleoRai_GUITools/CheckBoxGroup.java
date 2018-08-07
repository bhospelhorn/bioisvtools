package waffleoRai_GUITools;

import java.util.ArrayList;
import java.util.Collection;

import javax.swing.JCheckBox;

public class CheckBoxGroup {
	
	private JCheckBox[] boxes;
	
	public CheckBoxGroup(int numberBoxes)
	{
		if (numberBoxes <= 0) numberBoxes = 16;
		boxes = new JCheckBox[numberBoxes];
	}
	
	public void addCheckBox(JCheckBox box, int index)
	{
		if (index < 0 || index >= boxes.length) throw new IndexOutOfBoundsException();
		boxes[index] = box;
	}
	
	public void select(int index)
	{
		if (index < 0 || index >= boxes.length) return;
		JCheckBox cb = boxes[index];
		if (cb == null) return;
		if (!cb.isEnabled()){
			cb.setSelected(false); 
			return;
		}
		cb.setSelected(!cb.isSelected());
	}
	
	public void select(JCheckBox button)
	{
		if (button == null) return;
		if (!button.isEnabled()) return;
		for (int i = 0; i < boxes.length; i++)
		{
			if (boxes[i] == button)
			{
				button.setSelected(!button.isSelected());
			}
		}
	}
	
	public void select(int index, boolean isSelected)
	{
		if (index < 0 || index >= boxes.length) return;
		JCheckBox cb = boxes[index];
		if (cb == null) return;
		if (!cb.isEnabled()){
			cb.setSelected(false); 
			return;
		}
		cb.setSelected(isSelected);
	}
	
	public boolean isSelected(int index)
	{
		if (index < 0 || index >= boxes.length) return false;
		if (boxes[index] == null) return false;
		return boxes[index].isSelected();
	}
	
	public Collection<Integer> getSelectedIndices()
	{
		Collection<Integer> selected = new ArrayList<Integer>(boxes.length);
		for (int i = 0; i < boxes.length; i++)
		{
			if (boxes[i] != null)
			{
				if (boxes[i].isEnabled() && boxes[i].isSelected()) selected.add(i);
			}
		}
		return selected;
	}
	
	public Collection<JCheckBox> getSelectedBoxes()
	{
		Collection<JCheckBox> selected = new ArrayList<JCheckBox>(boxes.length);
		for (int i = 0; i < boxes.length; i++)
		{
			if (boxes[i] != null)
			{
				if (boxes[i].isEnabled() && boxes[i].isSelected()) selected.add(boxes[i]);
			}
		}
		return selected;
	}
	
	public void disable(int index)
	{
		if (index < 0 || index >= boxes.length) return;
		if (boxes[index] != null) boxes[index].setEnabled(false);
	}
	
	public void enable(int index)
	{
		if (index < 0 || index >= boxes.length) return;
		if (boxes[index] != null) boxes[index].setEnabled(true);
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
		for (JCheckBox cb : boxes) cb.setEnabled(b);
	}
	
	public void repaintAll()
	{
		for (JCheckBox cb : boxes) cb.repaint();
	}
	

}
