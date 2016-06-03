/*
 * Created on Nov 28, 2006
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package org.seqcode.gse.gsebricks.verbs.location;

import java.util.Iterator;

import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Gene;
import org.seqcode.gse.gsebricks.verbs.Expander;


public class RefGeneByNameGenerator implements Expander<String,Gene> { 

    private Genome genome;
    private RefGeneGenerator gen;
    
    public RefGeneByNameGenerator(Genome g) { 
        genome = g;
        gen = new RefGeneGenerator(g);
    }

    public Iterator<Gene> execute(String a) {
        return gen.byName(a);
    }
}
