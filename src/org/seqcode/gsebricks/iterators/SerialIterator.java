/*
 * Author: tdanford
 * Date: Apr 4, 2009
 */
package org.seqcode.gsebricks.iterators;

import java.util.Iterator;

import org.seqcode.gseutils.ArrayUtils;


public class SerialIterator<X> implements Iterator<X> {
	
	private Iterator<Iterator<X>> itritr;
	private Iterator<X> itr;
	
	public SerialIterator(Iterator<X>... itrs) { 
		this(ArrayUtils.asIterator(itrs)); 
	}
	
	public SerialIterator(Iterator<Iterator<X>> iis) { 
		itritr = iis;
		itr = null;
		while((itr == null || !itr.hasNext()) && itritr.hasNext()) { 
			itr = itritr.next();
		}
	}

	public boolean hasNext() {
		return itr != null && itr.hasNext();
	}

	public X next() {
		X val = itr.next();
		while(!itr.hasNext() && itritr.hasNext()) {
			itr = itritr.next();
		}
		return val;
	}

	public void remove() {
		itr.remove();
	}

}
