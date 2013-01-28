/*
 * Created on Mar 15, 2006
 */
package edu.psu.compbio.seqcode.gse.ewok.verbs.classification;

import edu.psu.compbio.seqcode.gse.ewok.verbs.Mapper;

/**
 * @author tdanford
 */
public interface Classifier<X> extends Mapper<X,Integer> {
    public int getNumClasses();
}
