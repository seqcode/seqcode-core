/*
 * Author: tdanford
 * Date: Jan 26, 2009
 */
package org.seqcode.gseutils.models;

public interface ModelListener<T extends Model> {
	public void modelChanged(T model);
}
