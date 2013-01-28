package edu.psu.compbio.seqcode.gse.tools.genenames;

import java.util.*;

public class PropertiesLoader {
	
	private ResourceBundle bundle;
	private HashSet<String> bundleKeys;
	
	public PropertiesLoader(String pn) { 
		bundle = ResourceBundle.getBundle("edu.psu.compbio.seqcode.gse.tools.genenames." + pn);
		bundleKeys = new HashSet<String>();
		Enumeration<String> keys = bundle.getKeys();
		while(keys.hasMoreElements()) { 
			bundleKeys.add(keys.nextElement());
		}
	}
	
	public Set<String> getKeys() { 
		return bundleKeys;
	}
	
	public boolean hasKey(String k) { 
		return bundleKeys.contains(k); 
	}
	
	public String getValue(String k) { 
		return bundle.getString(k);
	}
}
