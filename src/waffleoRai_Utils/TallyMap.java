package waffleoRai_Utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class TallyMap {

	private Map<Integer, Integer> coremap;
	
	public TallyMap()
	{
		coremap = new HashMap<Integer, Integer>();
	}
	
	public synchronized void increment(int value)
	{
		Integer t = coremap.get(value);
		if (t == null) coremap.put(value, 1);
		else coremap.put(value, t + 1);
	}
	
	public synchronized void decrement(int value)
	{
		Integer t = coremap.get(value);
		if (t == null) coremap.put(value, -1);
		else coremap.put(value, t - 1);
	}
	
	public synchronized List<Integer> getAllValues()
	{
		List<Integer> list = new ArrayList<Integer>(coremap.size());
		list.addAll(coremap.keySet());
		Collections.sort(list);
		
		return list;
	}
	
	public synchronized int getCount(int value)
	{
		Integer t = coremap.get(value);
		if (t == null) return 0;
		else return t;
	}
	
}
