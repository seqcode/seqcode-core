package edu.psu.compbio.seqcode.gse.utils.graphs;

import java.util.*;
import java.io.*;

public interface Graph { 
	public Set<String> getVertices();
	public Set<String> getNeighbors(String vertex);
	public boolean isNeighbor(String v, String n);
}

