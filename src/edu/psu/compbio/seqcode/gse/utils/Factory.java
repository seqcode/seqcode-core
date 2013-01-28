/*
 * Created on Aug 24, 2005
 */
package edu.psu.compbio.seqcode.gse.utils;

/**
 * @author tdanford
 */
public interface Factory<X> {
    public X createObject();
}
