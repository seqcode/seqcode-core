/*
 * Author: tdanford
 * Date: May 27, 2008
 */
package edu.psu.compbio.seqcode.gse.viz.components;

public interface SelectionListener<X> {
	public void selectionMade(SelectionEvent<X> evt);
}
