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

import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.Config;



public class EMplotter {
	
	protected double[][][] trainMUc;
	protected double[][][] trainMUf;
	protected double[][][] trainSIGMAc;
	protected double[][][] trainSIGMAf;
	protected double[][] trainPIj;
	protected double[][][] trainBjk;
	protected int w; // plot width
	protected int h; // plot height
	protected int wmargin; //wifth margin
	protected int hmargin; //height margon
	protected String out_base; // out put base string name
	
	public EMplotter(Config config, double[][][] trainMUc, double[][][] trainMUf, double[][][] trainSIGMAc, double[][][] trainSIGMAf, double[][] trainPIj, double[][][] trainBjk) {
		this.trainMUc = trainMUc;
		this.trainMUf = trainMUf;
		this.trainSIGMAc = trainSIGMAc;
		this.trainSIGMAf = trainSIGMAf;
		this.trainPIj = trainPIj;
		this.trainBjk = trainBjk;
		
		
	}
	
	
	public void plot(){
		BufferedImage im = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
		Graphics g = im.getGraphics();
		Graphics2D g2 = (Graphics2D)g;
		g2.setRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
		g2.setColor(Color.white);
		g2.fillRect(0, 0, w, h);
		
		int trackheight = h-2*hmargin;
		int trackwidth = w-2*wmargin;
		
		// Drawing Margins for the graph
		g2.setColor(Color.black);
		g2.drawLine(wmargin, hmargin+trackheight, w-wmargin,hmargin+trackheight ); //x-axis
		g2.drawLine(wmargin, hmargin+trackheight, wmargin, hmargin); //y- axis
		
		
		
		 
	}
	
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
	
	    int[] temp = new int[4];
	    System.out.println(temp[2]);
	    System.out.println("is null"); 
	    //File f = new File("/Users/akshaykakumanu/Desktop/temp");
	    
	    //ImageIO.write(im, "png", f );
	    
		}
}
