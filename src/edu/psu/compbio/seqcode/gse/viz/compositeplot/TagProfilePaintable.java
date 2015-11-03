package edu.psu.compbio.seqcode.gse.viz.compositeplot;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.Stroke;

import edu.psu.compbio.seqcode.gse.viz.colors.Coloring;
import edu.psu.compbio.seqcode.gse.viz.paintable.AbstractPaintable;

public class TagProfilePaintable extends AbstractPaintable{

	private TagProfile profile;
	private Color watsonCol=Color.blue;
	private Color crickCol=Color.blue;
	private boolean isStranded;
	private boolean filledColumns=false;
	private double ymax;
	private boolean drawAxis=true;
	private boolean drawPointMarkers=false;
	private boolean transparent=false;
	private int hmargin= 50, wmargin=30;
	private int leftLimit, rightLimit;

	public TagProfilePaintable(TagProfile profile){
		this.profile = profile;
		leftLimit = profile.getLeftRelCoord();
		rightLimit = profile.getRightRelCoord();
		isStranded = profile.isStranded();
	}
	
	public void setWatsonColor(Color c){watsonCol= transparent ? Coloring.clearer(c) : c;}
	public void setCrickColor(Color c){crickCol= transparent ? Coloring.clearer(c) : c;}
	public void setFilledColumns(boolean b){filledColumns=b;}
	public void setHmargin(int hm){hmargin = hm;}
	public void setWmargin(int wm){wmargin = wm;}
	public void setYmax(double y){ymax = y;}
	public double getYmax(){return ymax;}
	public void autoYmax(boolean auto){
		if(auto){
			double[] w = profile.getWatsonTags();
			double[] c = profile.getCrickTags();
			ymax = -Double.MAX_VALUE;
			for(int i=0; i<w.length; i++){
				if(w[i]>ymax)
					ymax=w[i];
				if(c[i]>ymax)
					ymax=c[i];
			}
		}
	}
	public void setProfileLeftLimit(int relCoord){
		if(relCoord>profile.getLeftRelCoord())
			leftLimit=relCoord;
	}
	public void setProfileRightLimit(int relCoord){
		if(relCoord<profile.getRightRelCoord())
			rightLimit=relCoord;
	}
	
	public void setTransparent(boolean t){
		if(!transparent && t){
			watsonCol = Coloring.clearer(watsonCol);
			crickCol = Coloring.clearer(crickCol);
		}transparent = t;
	}

	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		Graphics2D g2 = (Graphics2D)g;
		g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, 
				RenderingHints.VALUE_ANTIALIAS_ON);
		Stroke oldStroke = g2.getStroke();
		g2.setStroke(new BasicStroke((float)2.5));

		int w = x2-x1, h = y2-y1;
		int tY1 = hmargin, tY2 = h-hmargin;
		int tX1 = wmargin, tX2 = w-wmargin;
		int tH = h-(hmargin*2);
		int tW = w-(wmargin/2);
		g2.setColor(Color.white);
		if(!transparent)
			g2.fillRect(x1, y1, w, y2-y1);
		int[] xs = new int[profileWidth()];
		int[] yws = new int[profileWidth()];
		int[] ycs = new int[profileWidth()];
		
		//Split the area into two sections: watson & crick
		int wPlotH = profile.isStranded() ? (int)(tH/2) : tH;
		int wPlotY2 =tY1+wPlotH;
		int cPlotH = tH-wPlotH;
		int cPlotY1 = wPlotY2;
		int binPix = (int)Math.round((double)tW / ((double)profileWidth()));
		int centerX = tX1 + -leftLimit*binPix;
		
		// Make sure that the profile doesn't change out from underneath us... 
		synchronized(profile) { 
			int i=0;
			for(int coord = leftLimit; coord <= rightLimit; coord++) { 
				int x = tX1 + i*binPix;
				if(coord==0)
					centerX = x;
				//watson
				double ywf= fractionalOffsetY(profile.getWatsonTagsAtRelCoord(coord));
				int yw = wPlotY2 - (int)Math.round(ywf * (double)wPlotH);
				//crick
				if(isStranded){
					double ycf= fractionalOffsetY(profile.getCrickTagsAtRelCoord(coord));
					int yc = tY2 - (int)Math.round((1-ycf) * (double)cPlotH);
					ycs[i] = yc;
				}
				xs[i] = x; yws[i] = yw;
				i++;
			}
		
			//Draw the points
			if(filledColumns){
				//Histogram style
				for(i = 0; i < xs.length; i++) {
					g2.setColor(watsonCol);
					g2.fillRect(xs[i],yws[i], binPix, wPlotY2-yws[i]+1);
					if(isStranded){
						g2.setColor(crickCol);
						g2.fillRect(xs[i],cPlotY1, binPix, ycs[i]-cPlotY1+1);
					}
				}
			}else{
				//Line Graph style
				for(i = 1; i < xs.length; i++) {
					g2.setColor(watsonCol);
					g2.drawLine(xs[i-1], yws[i-1], xs[i], yws[i]);
					if(isStranded){
						g2.setColor(crickCol);
						g2.drawLine(xs[i-1], ycs[i-1], xs[i], ycs[i]);
					}
				}
				if(drawPointMarkers){
					int rad = 2;
					int diam = rad*2;
					for(i = 0; i < xs.length; i++) {
						g2.setColor(Color.white);
						g2.fillOval(xs[i]-rad, yws[i]-rad, diam, diam);
						g2.setColor(watsonCol);
						g2.drawOval(xs[i]-rad, yws[i]-rad, diam, diam);
						if(isStranded){
							g2.setColor(Color.white);
							g2.fillOval(xs[i]-rad, ycs[i]-rad, diam, diam);
							g2.setColor(crickCol);
							g2.drawOval(xs[i]-rad, ycs[i]-rad, diam, diam);
						}
					}
				}
			}
		}
		
		if(drawAxis){
			g2.setColor(Color.black);
    		g2.setStroke(new BasicStroke(2.0f));
    		g2.drawLine(centerX, tY1, centerX, tY2); //Y-axis
    		g2.drawLine(tX1, wPlotY2, tX2, wPlotY2); //X-axis
    		g2.setFont(new Font("Ariel", Font.PLAIN, 14));
    		FontMetrics metrics = g2.getFontMetrics();
    		//g2.drawString(String.format("%f",ymax), tX1+1, tY1+(metrics.getHeight()));
    		//if(isStranded)
    		//	g2.drawString(String.format("%f",ymax), tX1+1, tY2);
    		
        	//g2.drawLine(wmargin+center, tStart+trackHeight+8, wmargin+center, tStart+trackHeight);	// 0 tick on x-axis
        	//g2.drawString(String.format("%d",leftLimit), tX1+1, tY1+(metrics.getHeight()));
        	//		leftLimit, wmargin-(metrics.stringWidth(startCoord.toString())/2), tStart+trackHeight+metrics.getHeight()+10);
        	//g2.drawString(rightLimit, w-wmargin-(metrics.stringWidth(endCoord.toString())/2), tStart+trackHeight+metrics.getHeight()+10);
		}
		g2.setStroke(oldStroke);
	}

	
	private double fractionalOffsetY(double val){
		return (val)/(ymax);
	}
	private int profileWidth(){
		return rightLimit-leftLimit+1;
	}
}
