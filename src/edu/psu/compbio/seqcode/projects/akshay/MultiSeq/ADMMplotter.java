package edu.psu.compbio.seqcode.projects.akshay.MultiSeq;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

public class ADMMplotter {
	
	protected int w = 800; // plot width
	protected int h = 800; // plot height
	protected int wmargin = 80; //width margin
	protected int hmargin = 60; //height margin
	protected String out_base; // out put base string name
	
	/**
	 * Plot mehtod
	 * @param Xax
	 * @param Yax
	 * @param Xlab
	 * @param Ylab
	 * @param fp
	 */
	public void plot(double[][] Xax, double[][] Yax, String Xlab, String Ylab, File fp){
		try{
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
			
			double scaley = trackheight/getMax(Yax);
			double scalex = trackwidth/getMax(Xax);
			
			for(int j=0; j<Xax[0].length; j++){
				for (int i=0; i<Xax.length; i++){
					g2.drawLine((int)(Xax[i][j]*scalex)+wmargin, h-hmargin-(int)(Yax[i][j]*scaley), (int)(Xax[i][j+1]*scalex)+wmargin, h-hmargin-(int)(Yax[i][j+1]*scaley));
				}
			}
			
			double legendy = trackheight/5;
			double legendx = trackwidth/5;
			
			//Y-axis legend
			g2.setFont(new Font("Ariel", Font.PLAIN, 15));
			for(int i=0; i<5; i++){
				g2.drawString(Integer.toString((int)getMax(Yax)*(i+1)/5), (int)(wmargin-20),h-hmargin-(int)((i+1)*legendy));
			}
			
			//X-axis legend
			for(int i=0; i<5; i++){
				g2.drawString(Integer.toString((int)getMax(Xax)*(i+1)/5), (int)(wmargin+(i+1)*legendx),h-hmargin+20);
			}
			
			ImageIO.write(im, "png", fp );
		}
		catch(IOException e){
			e.printStackTrace();
		}
		
		 
	}
	
	/**
	 * Get the maximum value in a given 2-d array
	 * @param vec
	 * @return
	 */
	public double getMax(double[][] vec){
		double ret=-100.0;
		for(int i=0; i<vec.length; i++){
			for(int j=0; j<vec[0].length; j++){
				if(vec[i][j]> ret){
					ret = vec[i][j];
				}
			}
		}
		return ret;
	}
	

}
