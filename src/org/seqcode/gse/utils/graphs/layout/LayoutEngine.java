package org.seqcode.gse.utils.graphs.layout;


import java.io.*;
import java.util.*;
import java.awt.Point;

import org.seqcode.gse.utils.graphs.*;

public interface LayoutEngine<G extends Graph> { 
	public GraphLayout layoutGraph(G graph);
	public Set<String> getParameters();
	public void setParameter(String key, String value);
}
