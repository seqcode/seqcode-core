/*
 * Author: tdanford
 * Date: Sep 29, 2008
 */
package org.seqcode.gseutils.models;

public interface Modeleable {
	public Class getModelClass();

	public Model asModel();
}
