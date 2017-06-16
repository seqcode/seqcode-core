package org.seqcode.gseutils.graphs.layout;

import java.util.*;

import org.seqcode.gseutils.graphs.*;

import java.awt.*;

public interface GraphLayout<G extends Graph> {
	public G getGraph();

	public void displayGraph(Graphics2D g2, Rectangle bounds);

	public void displayNode(String node, Graphics2D g2, Rectangle bounds);

	public void displayEdge(String head, String tail, Graphics2D g2, Rectangle bounds);
}
