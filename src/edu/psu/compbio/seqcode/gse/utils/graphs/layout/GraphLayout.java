package edu.psu.compbio.seqcode.gse.utils.graphs.layout;

import java.util.*;
import java.awt.*;

import edu.psu.compbio.seqcode.gse.utils.graphs.*;

public interface GraphLayout<G extends Graph> { 
	public G getGraph();

	public void displayGraph(Graphics2D g2, Rectangle bounds);
	public void displayNode(String node, Graphics2D g2, Rectangle bounds);
	public void displayEdge(String head, String tail, Graphics2D g2, Rectangle bounds);
}

