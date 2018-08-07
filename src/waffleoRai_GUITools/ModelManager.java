package waffleoRai_GUITools;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import javax.swing.DefaultListModel;
import javax.swing.JList;
import javax.swing.ListModel;

public class ModelManager<T> {

	public ListModel<T> arrayToListModel(T[] array)
	{
		DefaultListModel<T> model = new DefaultListModel<T>();
		for (T e : array) model.addElement(e);
		return model;
	}
	
	public ListModel<T> collectionToListModel(Collection<T> collection)
	{
		DefaultListModel<T> model = new DefaultListModel<T>();
		for (T e : collection) model.addElement(e);
		return model;
	}
	
	public List<T> listModelToList(ListModel<T> model)
	{
		int size = model.getSize();
		List<T> list = new ArrayList<T>(size);
		for (int i = 0; i < size; i++) list.add(model.getElementAt(i));
		return list;
	}
	
	public void moveSelectedToAnotherList(JList<T> source, JList<T> target)
	{
		if (source.isSelectionEmpty()) return;
		List<T> s = listModelToList(source.getModel());
		List<T> t = listModelToList(target.getModel());
		List<T> v = source.getSelectedValuesList();
		t.addAll(v);
		s.removeAll(v);
		source.setModel(collectionToListModel(s));
		target.setModel(collectionToListModel(t));
	}
	
}
