package edu.psu.compbio.seqcode.projects.shaun.viz;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.util.ArrayList;

import edu.psu.compbio.seqcode.gse.viz.paintable.AbstractPaintable;

public class ColorBarPlotPaintable extends AbstractPaintable{
	
	private int screenSizeX, screenSizeY;
	private final int topBorder=100, bottomBorder=10;
	private int rows, cols;
	private BarPlotData barplot;
		
	public ColorBarPlotPaintable(BarPlotData bpd){
		barplot = bpd;
		rows=bpd.getRows();
		cols=bpd.getCols();		
	}
	
	public int getYborders(){
		return(topBorder+bottomBorder);
	}
	
	public void paintItem (Graphics g, int x1, int y1, int x2, int y2){
		Graphics2D g2d = (Graphics2D)g;
		screenSizeX = x2-x1;
		screenSizeY = y2-y1;
		g2d.setColor(Color.white);
		g2d.fillRect(0, 0, screenSizeX, screenSizeY);
		//g2d.setBackground(Color.white);
		ArrayList<Integer> colPos = new ArrayList<Integer>();
		
		//Define sizings
		int rectHeight = (screenSizeY-(topBorder+bottomBorder));
		int subSection = screenSizeX/(cols+1);
		int rectSpacing = (int)(subSection*0.2);
		int rectWidth = (screenSizeX-((cols+1)*rectSpacing))/cols;
		int cp = rectSpacing;
		for(int c=0; c<cols; c++){
			colPos.add(new Integer(cp));
			cp+=rectWidth+rectSpacing;
		}
		int thick= Math.max(1, rectHeight/rows);
		//int thick=1;
		g2d.setStroke(new BasicStroke(thick));
		
		//First column is a special coloring case; go through twice for positive & negative
		Color maxColor = barplot.colMaxColors.get(0);
		Color minColor = Color.black;
		double maxVal=Math.log(barplot.colMaxVals.get(0)), minVal=0;
		//double maxVal=barplot.colMaxVals.get(0), minVal=0;
		//double maxVal=colMaxVals.get(0), minVal=colMinVals.get(0);
		for(int r=0; r<rows; r++){
			double val = Math.log(barplot.data.get(0).get(r));
			//double val = barplot.data.get(0).get(r);
			if(val>=minVal){
				double sVal = (val-minVal)/(maxVal-minVal);
				int red = (int)(maxColor.getRed() * sVal + minColor.getRed() * (1 - sVal));
		        int green = (int)(maxColor.getGreen() * sVal + minColor.getGreen() * (1 - sVal));
		        int blue = (int)(maxColor.getBlue() *sVal + minColor.getBlue() * (1 - sVal));
		        Color currColor = new Color(red, green, blue);
				g2d.setColor(currColor);			
				int currY= topBorder+(rectHeight*r)/rows;
				g2d.drawLine(colPos.get(0), currY, colPos.get(0)+rectWidth, currY);
			}
		}
		maxColor = Color.black;
		minColor = barplot.colMinColors.get(0);
		maxVal=0; minVal=barplot.colMinVals.get(0);
		for(int r=0; r<rows; r++){
			double val = barplot.data.get(0).get(r);
			if(val<maxVal){
				double sVal = (val-minVal)/(maxVal-minVal);
				int red = (int)(maxColor.getRed() * sVal + minColor.getRed() * (1 - sVal));
		        int green = (int)(maxColor.getGreen() * sVal + minColor.getGreen() * (1 - sVal));
		        int blue = (int)(maxColor.getBlue() *sVal + minColor.getBlue() * (1 - sVal));
		        Color currColor = new Color(red, green, blue);
				g2d.setColor(currColor);			
				int currY= topBorder+(rectHeight*r)/rows;
				g2d.drawLine(colPos.get(0), currY, colPos.get(0)+rectWidth, currY);
			}
		}
		
		//Draw other indicators
		for(int c=1; c<cols; c++){
			maxVal=barplot.colMaxVals.get(c); minVal = 0; 
			//minVal=barplot.colMinVals.get(c);
			//Random Color
			if(c==1)
				maxColor = Color.red;
			else if(c==2)
				maxColor = Color.BLUE;
			else
				maxColor = barplot.colMaxColors.get(c);
			minColor = barplot.colMinColors.get(c);
			for(int r=0; r<rows; r++){
				double val = barplot.data.get(c).get(r);
				double sVal = (val-minVal)/(maxVal-minVal);
				int red = (int)(maxColor.getRed() * sVal + minColor.getRed() * (1 - sVal));
		        int green = (int)(maxColor.getGreen() * sVal + minColor.getGreen() * (1 - sVal));
		        int blue = (int)(maxColor.getBlue() *sVal + minColor.getBlue() * (1 - sVal));
		        Color currColor = new Color(red, green, blue);
		        g2d.setColor(currColor);
				if(val>minVal){
					int currY= topBorder+(rectHeight*r)/rows;
					g2d.drawLine(colPos.get(c), currY, colPos.get(c)+rectWidth, currY);
				}				
			}
		}
		
		//Set up bounding rectangles
		g2d.setColor(Color.black);
		for(int c=0; c<cols; c++){
			g2d.drawRect(colPos.get(c), topBorder, rectWidth, rectHeight);
		}
		//Draw Labels
		g2d.setFont(new Font("Serif", Font.BOLD, 60));
		for(int c=0; c<cols; c++){
			g2d.drawString(barplot.colNames.get(c), colPos.get(c), g2d.getFont().getSize()+topBorder/4);
		}		
	}
}
