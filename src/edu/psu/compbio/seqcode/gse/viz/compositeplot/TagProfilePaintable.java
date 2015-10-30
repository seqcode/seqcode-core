package edu.psu.compbio.seqcode.gse.viz.compositeplot;

import java.awt.Color;
import java.awt.Graphics;

import edu.psu.compbio.seqcode.gse.viz.paintable.AbstractPaintable;

public class TagProfilePaintable extends AbstractPaintable{

	private TagProfile profile;
	private Color watsonCol=Color.blue;
	private Color crickCol=Color.blue;
	private String style="line";  //Line or histogram
	private boolean drawPointMarkers=true;
	
	private int[] xs, ys;
	
	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {}

}
