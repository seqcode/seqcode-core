/*
 * Author: tdanford
 * Date: Jan 19, 2009
 */
package org.seqcode.utils.models;

import java.util.Iterator;

public interface Timer {

	public void addTiming(Timing t);
	public void addTimings(Iterator<Timing> ts);
}
