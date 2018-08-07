package waffleoRai_GUITools;

import java.awt.Component;
import java.awt.event.ActionListener;
import java.awt.event.WindowListener;

public interface ObjectEditorDialog<T> {

	public void updateObject();
	public void updateForm();
	public void addApproveButtonListener(ActionListener l);
	public void addCancelButtonListener(ActionListener l);
	public T retrieveContents();
	public void pack();
	public void setLocationRelativeTo(Component parent);
	public void addWindowListener(WindowListener l);
	public void setVisible(boolean b);
	public void dispose();
	
}
