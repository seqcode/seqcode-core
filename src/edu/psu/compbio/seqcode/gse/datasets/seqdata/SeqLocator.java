package edu.psu.compbio.seqcode.gse.datasets.seqdata;

import java.util.Collection;
import java.util.Set;
import java.util.TreeSet;
import java.util.Iterator;

/**
 * @author tdanford
 * @author mahony
 * 
 * SeqLocator describes a subset of SeqAlignments (rows from the 
 * seqalignment table) from a single SeqExpt, which satisfy the following conditions:
 *  
 * (1) all SeqAlignments share the same name (i.e. are from the same aligner)
 * (2) if the locator's 'reps' field is non-empty, then only those replicates whose
 *     names appear in that set are indicated.
 *     
 * Created as ChipSeqLocator on Jan 11, 2008
 */
public class SeqLocator implements Comparable<SeqLocator> {

    private String exptName;
    private String alignName;
    // The reps field has a particular semantics -- if it contains any entries, then 
    // the locator designates *only* those replicates (of the named experiment) 
    // that have the given alignment name too.
    // On the other hand, if "reps" is empty, then the locator designates *all* 
    // available replicates for which the given alignment name is valid.
    private Set<String> reps; 
    
    public SeqLocator(String ename, String rname, String aname) {
        exptName = ename;
        alignName = aname;
        reps = new TreeSet<String>();
        reps.add(rname);
    }

    public SeqLocator(String ename, Collection<String> rnames, String aname) {
        exptName = ename;
        alignName = aname;
        reps = new TreeSet<String>(rnames);
    }
    
    public boolean isSubset(SeqLocator loc) { 
        if(!exptName.equals(loc.exptName)) { return false; }
        if(!alignName.equals(loc.alignName)) { return false; }
        for(String rep : reps) { 
            if(!loc.reps.contains(rep)) { return false; }
        }
        return true;
    }
    
    public String getExptName() { return exptName; }
    public String getAlignName() { return alignName; }
    public Collection<String> getReplicates() { return reps; }    
    
    public String getReplicateString() {
        if(reps.isEmpty()) { return "all"; }
        StringBuilder sb = new StringBuilder();
        for(String rep : reps) { 
            if(sb.length() > 0) { sb.append(","); }
            sb.append(rep);
        }
        return sb.toString();
    }
    
    public int hashCode() { 
        int code = 17;
        code += exptName.hashCode(); code *= 37;
        code += alignName.hashCode(); code *= 37;
        
        for(String rep : reps) { code += rep.hashCode(); code *= 37; }
        return code;
    }
    
    public String toString() { 
    	return String.format("%s [rep: %s, align: %s]", exptName, getReplicateString(), alignName);
    }

    public int compareTo(SeqLocator other) {
        int c = exptName.compareTo(other.exptName);
        if (c == 0) {
            c = alignName.compareTo(other.alignName);
            if (c == 0) {
                Iterator<String> mine = reps.iterator();
                Iterator<String> others = other.reps.iterator();
                while (mine.hasNext() && others.hasNext()) {
                    c = mine.next().compareTo(others.next());
                    if (c != 0) {
                        break;
                    }
                }
                if (c == 0) {
                    if (mine.hasNext()) {
                        c = 1;
                    } else {
                        c = -1;
                    }
                }
            }
        }
        return c;
    }
    
    public boolean equals(Object o) { 
        if(!(o instanceof SeqLocator)) { return false; }
        SeqLocator loc =(SeqLocator)o;
        if(!exptName.equals(loc.exptName)) { return false; }
        if(reps.size() != loc.reps.size()) { return false; }
        for(String rep : reps) { 
            if(!loc.reps.contains(rep)) { return false; }
        }
        if(!alignName.equals(loc.alignName)) { return false; }
        return true;
    }
}
