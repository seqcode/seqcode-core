package org.seqcode.data.seqdata;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;


/**
 * Simple pairing of a SeqLocator with the SeqExpt common to all replicates
 * @author mahony
 *
 */
public class SeqLocatorMatchedAligns  implements Comparable<SeqLocatorMatchedAligns>{

	public List<SeqExpt> expts=null;
	public Collection<SeqAlignment> aligns=null;
	public SeqLocator locator=null;
	
	public SeqLocatorMatchedAligns(Collection<SeqAlignment> a, SeqLocator l) throws SeqLMENotReplicatesException{
		aligns=a;
		locator = l;
		expts = new ArrayList<SeqExpt>();
		
		Set<SeqExpt> tmpExpts = new TreeSet<SeqExpt>();
		if(a!=null && !a.isEmpty()){
			for(SeqAlignment align : aligns)
				expts.add(align.getExpt());
			
			expts.addAll(tmpExpts);
			for(SeqExpt expt: expts){
				if(!expts.get(0).isReplicateOf(expt))
					throw new SeqLMENotReplicatesException("SeqLocatorMatchedAligns: input experiments are not replicates: "+locator.toString());
			}
		}
	}
	
	public int hashCode() { 
        return locator.hashCode();
    }
	public String toString(){return locator.toString();}
	
	public int compareTo(SeqLocatorMatchedAligns o) {
		return locator.compareTo(o.locator);
	}
	public boolean equals(Object o) {
		if(!(o instanceof SeqLocatorMatchedAligns)) { return false; }
		SeqLocatorMatchedAligns le =(SeqLocatorMatchedAligns)o;
		if(!locator.getExptName().equals(le.locator.getExptName())) { return false; }
        if(locator.getReplicates().size() != le.locator.getReplicates().size()) { return false; }
        for(String rep : locator.getReplicates()) { 
            if(!le.locator.getReplicates().contains(rep)) { return false; }
        }
        if(!locator.getAlignName().equals(le.locator.getAlignName())) { return false; }
        
        for(SeqAlignment a : aligns)
        	if(!le.aligns.contains(a))
        		return false;
        
        return true;
	}
	
	public class SeqLMENotReplicatesException extends Exception{
		public SeqLMENotReplicatesException(String message) {
	        super(message);
	    }
	}
}
