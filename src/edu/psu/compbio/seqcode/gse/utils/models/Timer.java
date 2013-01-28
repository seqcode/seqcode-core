/*
 * Author: tdanford
 * Date: Jan 19, 2009
 */
package edu.psu.compbio.seqcode.gse.utils.models;

import java.util.Iterator;

public interface Timer {

	public void addTiming(Timing t);
	public void addTimings(Iterator<Timing> ts);
}
