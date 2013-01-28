package edu.psu.compbio.seqcode.gse.utils.graphs;

public interface CycleChecker {
	public boolean containsCycle();
	public boolean checkCycleWithEdge(String v1, String v2);
}