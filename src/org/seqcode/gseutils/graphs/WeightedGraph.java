/*
 * Author: tdanford
 * Date: Nov 15, 2008
 */
package org.seqcode.gseutils.graphs;

public interface WeightedGraph extends Graph {
	public Double weight(String node);
	public Double weight(String n1, String n2);
	
}
