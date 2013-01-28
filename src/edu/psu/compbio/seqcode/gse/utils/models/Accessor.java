/*
 * Author: tdanford
 * Date: Mar 21, 2009
 */
package edu.psu.compbio.seqcode.gse.utils.models;

public interface Accessor<T extends Model> {
	public Object get(T object);
	public void set(T object, Object val);
	public Class getType();
	public String getName();
	public String getBaseName();
}
