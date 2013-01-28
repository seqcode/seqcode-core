package edu.psu.compbio.seqcode.projects.shaun;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;

import edu.psu.compbio.seqcode.gse.utils.RealValuedHistogram;
import edu.psu.compbio.seqcode.gse.viz.paintable.AbstractPaintable;

/* A Paintable to handle RealValuedHistograms
 * Probably a redundant class, but I can't find any other version. 
 * Written by: Shaun (July 18th 2008)
 */

public class HistogramPaintable extends AbstractPaintable{
	
	private int screenSizeX, screenSizeY;
	private final int topBorder=50, bottomBorder=100;
	private final int leftBorder=100, rightBorder=30;
	private int topBound, bottomBound, leftBound, rightBound;
	private int numXTicks = 11, numYTicks=8;
	private double maxY, minY;
	private double maxX, minX;
	private boolean zeroY=false;
	private String title;
	private RealValuedHistogram histo;
	private Color barCol = Color.green; 
	
	//Constructors
	public HistogramPaintable(RealValuedHistogram h){this(h, "", Color.green);}
	public HistogramPaintable(RealValuedHistogram h, String titleText){this(h, titleText, Color.green);}
	public HistogramPaintable(RealValuedHistogram h, String titleText, Color barColor){
		histo = h;
		title=titleText;
		findLimits();
		barCol = barColor;
	}
	
	public void paintItem (Graphics g, int x1, int y1, int x2, int y2){
		Graphics2D g2d = (Graphics2D)g;
		FontMetrics metrics = g2d.getFontMetrics();
		screenSizeX = x2-x1;
		screenSizeY = y2-y1;
		g2d.setColor(Color.white);
		g2d.fillRect(0, 0, screenSizeX, screenSizeY);
		topBound = topBorder;
		bottomBound = screenSizeY-bottomBorder;
		leftBound = leftBorder;
		rightBound = screenSizeX-rightBorder - ((screenSizeX-rightBorder-leftBorder)%histo.getNumBins());
		
		//Draw the title
		if(title.length()>0){
			g2d.setColor(Color.black);
			g2d.setFont(new Font("Ariel", Font.BOLD, 20));
			metrics = g2d.getFontMetrics();
			g2d.drawString(title, leftBound+((rightBound-leftBound)/2)-metrics.stringWidth(title)/2, topBound-10);
		}
		
		//Draw the histogram bars
		g2d.setColor(barCol);
		g2d.setStroke(new BasicStroke(1.0f));
		int maxBarLen = (bottomBound-topBound);
		for(int b=0; b<histo.getNumBins(); b++){
			double currHeight = histo.getBin(b);
			int currBarLen = (int)(maxBarLen * ((currHeight-minY)/(maxY-minY)));
			int currYPos = bottomBound-currBarLen;
			int binWidth = (rightBound-leftBound)/histo.getNumBins();
			int currXPos = leftBound + (b*binWidth);
			
			g2d.fillRect(currXPos, currYPos, binWidth, currBarLen);
			//System.out.println(b+":\t"+currHeight+"\t"+currBarLen);
		}
		
		//Draw the axes & ticks
		g2d.setColor(Color.black);
		g2d.setStroke(new BasicStroke(2.0f));
		g2d.setFont(new Font("Ariel", Font.PLAIN, 14));
		metrics = g2d.getFontMetrics();
		//y-axis
		g2d.drawLine(leftBound, topBound, leftBound, bottomBound);
		double [] yTicks = new double[numYTicks];
		int [] yTickPos = new int[numYTicks];
		for(int y=0; y<numYTicks; y++){
			yTicks[y] = minY + ((double)y*(maxY-minY)/((double)numYTicks-1));
			yTickPos[y] = bottomBound + (y*(topBound-bottomBound)/(numYTicks-1));
			g2d.drawLine(leftBound-3, yTickPos[y], leftBound+3, yTickPos[y]);
			int textPos =0;
			if(yTicks[y]>=1){
				textPos = leftBound-5-metrics.stringWidth(String.format("%.0f", new Double(yTicks[y])));
				g2d.drawString(String.format("%.0f", new Double(yTicks[y])), textPos, yTickPos[y]);
			}else{
				textPos = leftBound-5-metrics.stringWidth(String.format("%f", new Double(yTicks[y])));
				g2d.drawString(String.format("%f", new Double(yTicks[y])), textPos, yTickPos[y]);
			}
		}
		//zero line if necessary
		if(!zeroY && minY<0){
			g2d.setColor(Color.gray);
			int zeroPos = (int)(bottomBound + (((-1*minY)/(maxY-minY))*(double)(topBound-bottomBound)));
			g2d.drawLine(leftBound, zeroPos, rightBound, zeroPos);
			g2d.setColor(Color.black);
		}
		//x-axis
		g2d.drawLine(leftBound, bottomBound, rightBound, bottomBound);
		double [] xTicks = new double[numXTicks];
		int [] xTickPos = new int[numXTicks];
		for(int x=0; x<numXTicks; x++){
			xTicks[x] = Math.floor(minX + (x*(maxX-minX)/(numXTicks-1)));
			xTickPos[x] = leftBound + (x*(rightBound-leftBound)/(numXTicks-1));
			g2d.drawLine(xTickPos[x], bottomBound+3, xTickPos[x], bottomBound-3);
			g2d.drawString(String.format("%.0f", new Double(xTicks[x])), xTickPos[x], bottomBound+20);
		}//Extra tick for zero
		int zbin = histo.getBinContainingVal(0.0);
		int zTickPos = leftBound + (zbin*(rightBound-leftBound)/(histo.getNumBins()));
		g2d.drawLine(zTickPos, bottomBound+3, zTickPos, bottomBound-3);
		g2d.drawString("0", zTickPos, bottomBound+20);
	}
	
	public void zeroTheY(boolean z){zeroY=z; if(zeroY){minY=0;}}
	
	private void findLimits(){
		maxY = histo.getMaxBinCount(); 
		minY = zeroY ? 0 : histo.getMinBinCount();
		minX = histo.getHistoStart();
		maxX = histo.getHistoStop(); 
	}
}
