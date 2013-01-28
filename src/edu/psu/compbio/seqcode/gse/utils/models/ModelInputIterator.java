/*
 * Author: tdanford
 * Date: Dec 19, 2008
 */
package edu.psu.compbio.seqcode.gse.utils.models;

import java.util.*;

import edu.psu.compbio.seqcode.gse.utils.Closeable;

public class ModelInputIterator<T extends Model> implements Iterator<T>, Closeable {
	
	private ModelInput<T> input;
	private T next;
	
	public ModelInputIterator(ModelInput<T> inp) { 
		input = inp;
		next = input.readModel();
	}
	
	public void close() { 
		input.close();
		input = null;
		next = null;
	}
	
	public boolean isClosed() { 
		return input ==null;
	}

	public boolean hasNext() {
		return next != null;
	}

	public T next() {
		T val = next;
		next = input.readModel();
		return val;
	}

	public void remove() {
		throw new UnsupportedOperationException();
	}
}
