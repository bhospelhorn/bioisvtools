package waffleoRai_GUITools;

import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

public abstract class ObjectEditor<T> {

	public static final int OBJECT_EDITOR_APPROVE_OPTION = 1;
	public static final int OBJECT_EDITOR_CANCEL_OPTION = 0;
	public static final int OBJECT_EDITOR_ERROR_OPTION = -1;
	
	protected ObjectEditorDialog<T> myDialog = null;
	protected T myObject = null;
	
	protected int returnOption;
	
	protected abstract ObjectEditorDialog<T> instantiateDialog(Frame parent, T target);
	
	public int showMe(Frame parent, T target)
	{
		if (myDialog != null) throw new IllegalStateException();
		returnOption = OBJECT_EDITOR_ERROR_OPTION;
		myDialog = createDialog(parent, target);
		myDialog.addWindowListener(new WindowAdapter() 
		{
			public void windowClosing(WindowEvent e)
			{
				returnOption = OBJECT_EDITOR_CANCEL_OPTION;
			}
		});
		myDialog.setVisible(true);
		
		myDialog.dispose();
		myDialog = null;
		return 0;
	}
	
	public T getTarget()
	{
		return myObject;
	}
	
	public ObjectEditorDialog<T> createDialog(Frame parent, T target)
	{
		ObjectEditorDialog<T> d = instantiateDialog(parent, target);
		d.addApproveButtonListener(new ActionListener() 
		{
			public void actionPerformed(ActionEvent e)
			{
				approveChanges();
			}
		});
		d.addCancelButtonListener(new ActionListener() 
		{
			public void actionPerformed(ActionEvent e)
			{
				rejectChanges();
			}
		});
		d.pack();
		d.setLocationRelativeTo(parent);
		return d;
	}
	
	public void approveChanges()
	{
		if (myDialog == null) return;
		myDialog.updateObject();
		myObject = myDialog.retrieveContents();
		returnOption = OBJECT_EDITOR_APPROVE_OPTION;
		myDialog.setVisible(false);
	}
	
	public void rejectChanges()
	{
		myObject = null;
		returnOption = OBJECT_EDITOR_CANCEL_OPTION;
		myDialog.setVisible(false);
	}
	
}
