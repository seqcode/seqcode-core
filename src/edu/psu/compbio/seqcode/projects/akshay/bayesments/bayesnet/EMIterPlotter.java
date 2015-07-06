package edu.psu.compbio.seqcode.projects.akshay.bayesments.bayesnet;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.BayesmentsConfig;

public class EMIterPlotter {
	
	protected int w; // plot width
	protected int h; // plot height
	protected int wmargin; //width margin
	protected int hmargin; //height margin
	protected String out_base; // out put base string name
	protected int itrs; // total number of EM iterations
	//2-d arrays of the values to be plotted
	// rows: The number of parameters to be plotted on this graphic/plot
	// columns: The total number of EM iterations
	protected double[][] Xaxes;
	protected double[][] Yaxes;
	protected String parameter_name;
	protected File image_file;
	
	public EMIterPlotter(BayesmentsConfig config, double[][] Xaxes, double[][] Yaxes, String param_name) {
		this.w = config.W;
		this.h = config.H;
		this.wmargin = config.W_MARGIN;
		this.hmargin = config.H_MARGIN;
		this.out_base = config.getOutBase();
		this.itrs = config.getNumItrs();
		this.Xaxes = Xaxes;
		this.Yaxes = Yaxes;
		this.parameter_name = param_name;
		
		File f = config.getOutputImagesDir();
		if(!f.exists()){f.mkdirs();}
		
		this.image_file = new File(f.getAbsoluteFile()+File.separator+out_base+"_"+this.parameter_name+".png");
	}
	
	
	private double getMax(double[][] vec){
		double ret=-100.0;
		for(int i=0; i<vec.length; i++){
			for(int j=0; j<itrs+1; j++){
				if(vec[i][j]> ret){
					ret = vec[i][j];
				}
			}
		}
		return ret;
	}
	
	
	public void plot(){
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
			
			double scaley = trackheight/getMax(this.Yaxes);
			double scalex = trackwidth/getMax(this.Xaxes);
			
			for(int i=0; i<this.Xaxes.length; i++){
				for (int j=0; j<itrs; j++){
					g2.drawLine((int)(this.Xaxes[i][j]*scalex)+wmargin, h-hmargin-(int)(this.Yaxes[i][j]*scaley), (int)(this.Xaxes[i][j+1]*scalex)+wmargin, h-hmargin-(int)(this.Yaxes[i][j+1]*scaley));
				}
			}
			
			double legendy = trackheight/5;
			double legendx = trackwidth/5;
			
			//Y-axis legend
			g2.setFont(new Font("Ariel", Font.PLAIN, 15));
			for(int i=0; i<5; i++){
				g2.drawString(Integer.toString((int)getMax(this.Yaxes)*(i+1)/5), (int)(wmargin-20),h-hmargin-(int)((i+1)*legendy));
			}
			
			//X-axis legend
			for(int i=0; i<5; i++){
				g2.drawString(Integer.toString((int)getMax(this.Xaxes)*(i+1)/5), (int)(wmargin+(i+1)*legendx),h-hmargin+20);
			}
			
			ImageIO.write(im, "png", this.image_file );
		}
		catch(IOException e){
			e.printStackTrace();
		}
	}

}
