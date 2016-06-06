/*
 * Author: tdanford
 * Date: May 27, 2008
 */
package org.seqcode.viz.components;

public interface SelectionListener<X> {
	public void selectionMade(SelectionEvent<X> evt);
}
