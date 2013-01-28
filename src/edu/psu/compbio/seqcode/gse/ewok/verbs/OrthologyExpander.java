/*
 * Created on Apr 1, 2006
 */
package edu.psu.compbio.seqcode.gse.ewok.verbs;

import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.orthology.*;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;

/**
 * @author tdanford
 */
public class OrthologyExpander implements Expander<String,String> {
    
    private Genome firstGenome, secondGenome;
    private OrthologyLoader loader;
    private OrthologyMapping mapping;

    public OrthologyExpander(OrthologyLoader l, OrthologyMapping m, Genome first, Genome second) {
        loader = l;
        firstGenome = first;
        secondGenome = second;
        mapping = m;
    }

    /* (non-Javadoc)
     * @see edu.psu.compbio.seqcode.gse.ewok.verbs.Expander#execute(null)
     */
    public Iterator<String> execute(String input) {
        Collection<OrthologyPair> pairs = loader.getFirstNameGenomePairs(input, firstGenome);
        LinkedList<String> lst = new LinkedList<String>();
        
        for(OrthologyPair op : pairs) { 
            if(op.getGenome2().equals(secondGenome)) { 
                lst.addLast(op.getName2());
            }
        }
        
        return lst.iterator();
    }

}
