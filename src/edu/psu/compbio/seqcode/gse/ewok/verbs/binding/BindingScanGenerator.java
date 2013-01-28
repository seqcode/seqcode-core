/*
 * Created on Dec 5, 2006
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.psu.compbio.seqcode.gse.ewok.verbs.binding;

import edu.psu.compbio.seqcode.gse.datasets.binding.*;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.*;
import edu.psu.compbio.seqcode.gse.utils.iterators.EmptyIterator;

import java.sql.SQLException;
import java.util.*;

public class BindingScanGenerator implements Expander<Region,BindingEvent> {

    private BindingScanLoader loader;
    private LinkedList<BindingScan> scans;
    
    public BindingScanGenerator(BindingScanLoader l, BindingScan s) { 
        loader = l;
        scans = new LinkedList<BindingScan>(); scans.add(s);
    }

    public BindingScanGenerator(BindingScanLoader l) { 
        loader = l;
        scans = new LinkedList<BindingScan>(); 
    }
    
    public void addBindingScan(BindingScan bs) { scans.addLast(bs); }

    public Iterator<BindingEvent> execute(Region a) {
        LinkedList<BindingEvent> evts = new LinkedList<BindingEvent>();
        
        for(BindingScan scan : scans) { 
            try {
                Collection<BindingEvent> events = loader.loadEvents(scan, a);
                evts.addAll(events);
            } catch (SQLException e) {
                e.printStackTrace();
            }
        }
        
        return evts.iterator();
    }
}
