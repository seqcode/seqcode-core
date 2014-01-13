package edu.psu.compbio.seqcode.projects.akshay.bayesments.bayesnet;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.Stroke;
import java.awt.geom.Line2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;



public class EMplotter {
	
	public static void main(String[] args) throws IOException{
		int w=500;
		int h=400+50+50;
		int[] mu = {40,30,60,380,0};
		int[] x = {10,20,30,40,50};
		int componentThickness =2;
		BufferedImage im = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
		Graphics g = im.getGraphics();
	    Graphics2D g2 = (Graphics2D)g;
	    g2.setRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
	    g2.setColor(Color.white);
	    g2.fillRect(0, 0, w, h);
	    Stroke defaultStroke = g2.getStroke();
	    Stroke componentStroke = new BasicStroke(componentThickness);
	    g2.setColor(Color.black);
	    g2.drawLine(50, 50+400, 450, 50+400);
	    g2.drawLine(50, 50, 50, 50+400);
	    
	    //g2.setStroke(componentStroke);
	    g2.setStroke(defaultStroke);
	    for(int i=0; i<4; i++){
	    	g2.setColor(Color.black.darker());
	    	g2.draw(new Line2D.Double(50+x[i], h-50-mu[i], 50+x[i+1], h-50-mu[i+1]));
	    	g2.setFont(new Font("Ariel", Font.PLAIN, 6));
	    	g2.drawString(Integer.toString(x[i]), 50+x[i], h-50+15);
	    }
	
	    File f = new File("/Users/akshaykakumanu/Desktop/temp");
	    
	    ImageIO.write(im, "png", f );
	    
		}
}
