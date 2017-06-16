/*
 * Created on Aug 22, 2005
 */
package org.seqcode.gseutils;

/**
 * @author tdanford
 */
public interface Listener<Event> {
	public void eventRegistered(Event e);
}
