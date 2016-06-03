package org.seqcode.gse.viz.metaprofile;

import java.awt.*;

import org.seqcode.gse.viz.paintable.*;


/**
 * ProfilePaintable wraps a Profile object, and paints it in a window region in the simplest possible way 
 * (connected lines between points).
 * 
 * Profile style can be one of two options: line or histogram
 * 
 * @author tdanford
 * Date: Aug 12, 2008
 */
public class ProfilePaintable extends AbstractPaintable {
	
	private PaintableScale scale;
	private BinningParameters params;
	private Profile profile;
	private Color col=Color.blue;
	private String style="line";
	private boolean drawPointMarkers=true;
	
	private int[] xs, ys;
	
	public ProfilePaintable(PaintableScale sc, Profile p) { 
		scale = sc;
		profile = p;
		params = profile.getBinningParameters();
		xs = new int[params.getNumBins()];
		ys = new int[params.getNumBins()];
	}
	
	public Color getColor(){return col;}
	public void setColor(Color c){col=c;}
	public void setDrawPointMarkers(boolean drawMarkers){drawPointMarkers = drawMarkers;}
	/**
	 * Set style of profile
	 * @param s (line or histogram)
	 */
	public void setStyle(String s){style=s;}
	
	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		int w = x2-x1, h = y2-y1;
		Graphics2D g2 = (Graphics2D)g;
		Stroke oldStroke = g2.getStroke();
		g2.setStroke(new BasicStroke((float)2.5));
		
		int binPix = w / (params.getNumBins()+1);
		
		// Make sure that the profile doesn't change out from underneath us... 
		synchronized(profile) { 
			if(profile.min()<scale.getMin())
				scale.setScale(profile.min(), scale.getMax());
			if(profile.max()>scale.getMax())
				scale.setScale(scale.getMin(), profile.max());
			
			for(int i = 0; i < params.getNumBins(); i++) { 
				int x = x1 + (i+1)*binPix;
				double yf= scale.fractionalOffset(profile.value(i));
				int y = y2 - (int)Math.round(yf * (double)h);
				xs[i] = x; ys[i] = y;
			}
		
			g2.setColor(col);
			
			if(style.equals("histogram")){
				//Histogram style
				for(int i = 0; i < xs.length; i++) {
					g2.fillRect(xs[i],ys[i], binPix, h-ys[i]);
				}
			}else{
				//Line Graph style
				for(int i = 1; i < xs.length; i++) { 
					g2.drawLine(xs[i-1], ys[i-1], xs[i], ys[i]);
				}
				if(drawPointMarkers){
					int rad = 2;
					int diam = rad*2;
					for(int i = 0; i < xs.length; i++) {
						g2.setColor(Color.white);
						g2.fillOval(xs[i]-rad, ys[i]-rad, diam, diam);
						g2.setColor(col);
						g2.drawOval(xs[i]-rad, ys[i]-rad, diam, diam);
					}
				}
			}
		}
		g2.setStroke(oldStroke);
	}

}
