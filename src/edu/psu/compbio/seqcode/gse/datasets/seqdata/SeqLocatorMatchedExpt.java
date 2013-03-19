package edu.psu.compbio.seqcode.gse.datasets.seqdata;

import java.util.List;


/**
 * Simple pairing of a SeqLocator with the SeqExpt common to all replicates
 * @author mahony
 *
 */
public class SeqLocatorMatchedExpt  implements Comparable<SeqLocatorMatchedExpt>{

	public List<SeqExpt> expts=null;
	public SeqLocator locator=null;
	
	public SeqLocatorMatchedExpt(List<SeqExpt> e, SeqLocator l){
		expts = e;
		locator = l;
	}
	
	public int hashCode() { 
        return locator.hashCode();
    }
	public String toString(){return locator.toString();}
	
	public int compareTo(SeqLocatorMatchedExpt o) {
		return locator.compareTo(o.locator);
	}
	public boolean equals(Object o) {
		if(!(o instanceof SeqLocatorMatchedExpt)) { return false; }
		SeqLocatorMatchedExpt le =(SeqLocatorMatchedExpt)o;
		if(!locator.getExptName().equals(le.locator.getExptName())) { return false; }
        if(locator.getReplicates().size() != le.locator.getReplicates().size()) { return false; }
        for(String rep : locator.getReplicates()) { 
            if(!le.locator.getReplicates().contains(rep)) { return false; }
        }
        if(!locator.getAlignName().equals(le.locator.getAlignName())) { return false; }
        return true;
	}
}
