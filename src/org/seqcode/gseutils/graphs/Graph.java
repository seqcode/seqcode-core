package org.seqcode.gseutils.graphs;

import java.util.*;
import java.io.*;

public interface Graph {
	public Set<String> getVertices();

	public Set<String> getNeighbors(String vertex);

	public boolean isNeighbor(String v, String n);
}
