/*
 * Author: tdanford
 * Date: Jan 26, 2009
 */
package edu.psu.compbio.seqcode.gse.utils.models;

public interface ModelListener<T extends Model> {
	public void modelChanged(T model);
}
