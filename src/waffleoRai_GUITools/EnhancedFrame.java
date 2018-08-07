package waffleoRai_GUITools;

import java.awt.Cursor;
import java.util.Collection;
import java.util.Map;

import javax.swing.JFrame;

public abstract class EnhancedFrame extends JFrame{

	private static final long serialVersionUID = 4947483465018400352L;
	
	protected Map<String, ComponentGroup> groups;
	protected Map<String, String> paths; //For last used paths in JFileChooser
	
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
	}

	public void setWait()
	{
		setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
		disableAll();
	}
	
	public void unsetWait()
	{
		setCursor(null);
		restoreEnabled();
	}
	

}
