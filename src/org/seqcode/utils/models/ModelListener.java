/*
 * Author: tdanford
 * Date: Jan 26, 2009
 */
package org.seqcode.utils.models;

public interface ModelListener<T extends Model> {
	public void modelChanged(T model);
}
