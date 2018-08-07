package waffleoRai_GUITools;

import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;

import javax.swing.JList;

public class InExListPair<T> {
	
	private JList<T> exList;
	private JList<T> inList;
	
	private Collection<T> sourceIn;
	private Collection<T> sourceAll;
	
	public InExListPair(JList<T> excludedList, JList<T> includedList, Collection<T> included, Collection<T> all)
	{
		if (excludedList != null) exList = excludedList;
		else exList = new JList<T>();
		if (includedList != null) inList = includedList;
		else inList = new JList<T>();
		if (included != null) sourceIn = included;
		else sourceIn = new LinkedList<T>();
		if (all != null) sourceAll = all;
		else sourceAll = new LinkedList<T>();
	}
	
	public JList<T> getExcludeList()
	{
		return exList;
	}
	
	public JList<T> getIncludeList()
	{
		return inList;
	}
	
	public Collection<T> getSourceList()
	{
		return sourceAll;
	}
	
	public Collection<T> getSourceIncludeList()
	{
		return sourceIn;
	}
	
	protected void setSourceList(Collection<T> source)
	{
		sourceAll = source;
	}
	
	protected void setSourceIncludeList(Collection<T> source)
	{
		sourceIn = source;
	}
	
	public void repaint()
	{
		exList.repaint();
		inList.repaint();
	}
	
	public void disable()
	{
		exList.setEnabled(false);
		inList.setEnabled(false);
	}
	
	public void enable()
	{
		exList.setEnabled(true);
		inList.setEnabled(true);
	}
	
	private Collection<T> generateExcludeList()
	{
		Collection<T> sourceEx = new ArrayList<T>(sourceAll.size() - sourceIn.size());
		for (T o : sourceAll)
		{
			if (!(sourceIn.contains(o))) sourceEx.add(o);
		}
		return sourceEx;
	}
	
	public void updateGraphicLists()
	{
		ModelManager<T> modeler = new ModelManager<T>();
		Collection<T> sourceEx = generateExcludeList();
		exList.setModel(modeler.collectionToListModel(sourceEx));
		inList.setModel(modeler.collectionToListModel(sourceIn));
		repaint();
	}
	
	public void updateSourceLists()
	{
		ModelManager<T> modeler = new ModelManager<T>();
		sourceIn.clear();
		sourceIn.addAll(modeler.listModelToList(inList.getModel()));
	}
	
	public void includeSelected()
	{
		ModelManager<T> modeler = new ModelManager<T>();
		modeler.moveSelectedToAnotherList(exList, inList);
		repaint();
	}
	
	public void excludeSelected()
	{
		ModelManager<T> modeler = new ModelManager<T>();
		modeler.moveSelectedToAnotherList(inList, exList);
		repaint();
	}
	
}
