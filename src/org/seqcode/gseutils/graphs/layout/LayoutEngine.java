package org.seqcode.gseutils.graphs.layout;


import java.io.*;
import java.util.*;

import org.seqcode.gseutils.graphs.*;

import java.awt.Point;

public interface LayoutEngine<G extends Graph> { 
	public GraphLayout layoutGraph(G graph);
	public Set<String> getParameters();
	public void setParameter(String key, String value);
}
