package org.seqcode.projects.shaun.dataexplorer;

import java.util.ArrayList;
import java.util.Collection;

import org.seqcode.genome.location.Region;


public class DataCollectionFilter {

	private String filterRule;
	
	public DataCollectionFilter(String rule){filterRule=rule;}
	
	private boolean rowPasses(DataCollection collection, Region r){
		boolean accept=false;
		
		int passCount=0;
		for(DataSource d : collection.getSources()){
			boolean pass=false;
			double[] v = collection.getDataVal(r, d);
			for(int i=0; i<v.length; i++)
				if(v[i]>=d.getThreshold())
					pass=true;
			if(pass)
				passCount++;
		}
		if(filterRule.equals("AND") && passCount==collection.getSources().size())
			accept=true;
		else if(filterRule.equals("OR") && passCount>0)
			accept=true;
		
		return(accept);
	}
	
	public void filter(DataCollection dc){
		ArrayList<Region> removeRegs = new ArrayList<Region>();
		for(Region r : dc.getRegions()){
			if(!rowPasses(dc,r))
				removeRegs.add(r);
		}
		for(Region r : removeRegs)
			dc.removeRow(r);
	}
}
