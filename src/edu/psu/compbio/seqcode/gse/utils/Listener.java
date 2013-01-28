/*
 * Created on Aug 22, 2005
 */
package edu.psu.compbio.seqcode.gse.utils;

/**
 * @author tdanford
 */
public interface Listener<Event> {
    public void eventRegistered(Event e);
}
