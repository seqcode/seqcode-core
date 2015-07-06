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

import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.BayesmentsConfig;

/**
 * Plots all the parameters in the Bayesian network over all rouds of traning iterations
 * @author akshaykakumanu
 *
 */

public class EMplotter {
	
	protected double[][][] trainMUc;
	protected double[][][] trainMUf;
	protected double[][][] trainSIGMAc;
	protected double[][][] trainSIGMAf;
	protected double[][] trainPIj;
	protected double[][][] trainBjk;
	protected int w; // plot width
	protected int h; // plot height
	protected int wmargin; //width margin
	protected int hmargin; //height margin
	protected String out_base; // out put base string name
	protected int itrs;
	
	/**
	 * Constructor method
	 * @param config
	 * @param trainMUc
	 * @param trainMUf
	 * @param trainSIGMAc
	 * @param trainSIGMAf
	 * @param trainPIj
	 * @param trainBjk
	 * @param C
	 * @param F
	 */
	public EMplotter(BayesmentsConfig config, double[][][] trainMUc, double[][][] trainMUf, double[][][] trainSIGMAc, double[][][] trainSIGMAf, double[][] trainPIj, double[][][] trainBjk, int C, int F, int nChromStates, int nFacStates)  {
		this.trainMUc = trainMUc;
		this.trainMUf = trainMUf;
		this.trainSIGMAc = trainSIGMAc;
		this.trainSIGMAf = trainSIGMAf;
		this.trainPIj = trainPIj;
		this.trainBjk = trainBjk;
		this.itrs = config.getNumItrs();
		int numChromStates = nChromStates;
		int numFacBindingStates = nFacStates;
		this.w = config.W;
		this.h = config.H;
		this.wmargin = config.W_MARGIN;
		this.hmargin = config.H_MARGIN;
		this.out_base = config.getOutBase();
	
		
		double[][] Xaxes;
		double[][] Yaxes;
		
		//Create Image directory if it does not already exist
		File f = config.getOutputImagesDir();
		if(!f.exists()){f.mkdirs();}
		
		
		// Draw PIj graphs
		//Initialzie the Xaxis and Yaxis vectors
		Xaxes = new double[numChromStates][itrs+1];  // plus 1 for initial random parameters
		Yaxes = new double[numChromStates][itrs+1];

		for(int j=0; j<numChromStates; j++){
			for(int itr =0; itr<itrs+1; itr++){
				Xaxes[j][itr] = itr;
				Yaxes[j][itr] = trainPIj[itr][j];
			}
		}
		File fp = new File(f.getAbsoluteFile()+File.separator+out_base+"_PI-c.png");
		plot(Xaxes,Yaxes,"number of interations","PI-C values",fp);
		
		//Draw MUc graphs
		
		//Re-Initialize the Xaxis and Yaxis vectors
		Xaxes = null; 
		Yaxes=null;
		Xaxes = new double[C*numChromStates][itrs+1];
		Yaxes = new double[C*numChromStates][itrs+1];
		int count=0;
		for(int j=0; j<numChromStates; j++){
			for(int c=0; c<C; c++){
				for(int itr=0; itr<itrs+1; itr++){
					Xaxes[count][itr] = itr;
					Yaxes[count][itr] = trainMUc[itr][j][c];
				}
				count++;
			}
		}
		fp = new File(f.getAbsolutePath()+File.separator+out_base+"_MU-c.png");
		plot(Xaxes,Yaxes,"number of interations","MU-C values",fp);
		
		//Draw MUf graphs
		
		//Re-Initialize the Xaxis and Yaxis vectors
		Xaxes = null; 
		Yaxes=null;
		Xaxes = new double[F*numFacBindingStates][itrs+1];
		Yaxes = new double[C*numFacBindingStates][itrs+1];
		
		count=0;
		for(int k=0; k<numFacBindingStates; k++){
			for(int nf=0; nf<F; nf++){
				for(int itr=0; itr<itrs+1; itr++){
					Xaxes[count][itr] = itr;
					Yaxes[count][itr] = trainMUf[itr][k][nf];
				}
				count++;
			}
		}
		fp = new File(f.getAbsolutePath()+File.separator+out_base+"_MU-f.png");
		plot(Xaxes,Yaxes,"number of interations","MU-F values",fp);
		
		//Draw SIGMAc graphs
		Xaxes = null; 
		Yaxes=null;
		Xaxes = new double[C*numChromStates][itrs+1];
		Yaxes = new double[C*numChromStates][itrs+1];
		count=0;
		for(int j=0; j<numChromStates; j++){
			for(int c=0; c<C; c++){
				for(int itr=0; itr<itrs+1; itr++){
					Xaxes[count][itr] = itr;
					Yaxes[count][itr] = trainSIGMAc[itr][j][c];
				}
				count++;
			}
		}
		fp = new File(f.getAbsolutePath()+File.separator+out_base+"_SIGMA-c.png");
		plot(Xaxes,Yaxes,"number of interations","SIGMA-C values",fp);
		
		//Draw SIGMAf graphs
		Xaxes = null; 
		Yaxes=null;
		Xaxes = new double[F*numFacBindingStates][itrs+1];
		Yaxes = new double[C*numFacBindingStates][itrs+1];
		
		count=0;
		for(int k=0; k<numFacBindingStates; k++){
			for(int nf=0; nf<F; nf++){
				for(int itr=0; itr<itrs+1; itr++){
					Xaxes[count][itr] = itr;
					Yaxes[count][itr] = trainSIGMAf[itr][k][nf];
				}
				count++;
			}
		}
		fp = new File(f.getAbsolutePath()+File.separator+out_base+"_SIGMA-f.png");
		plot(Xaxes,Yaxes,"number of interations","SIGMA-F values",fp);
		
		//Draw Bjk graps
		Xaxes = null; 
		Yaxes=null;
		Xaxes = new double[numFacBindingStates*numChromStates][itrs+1];
		Yaxes = new double[numFacBindingStates*numChromStates][itrs+1];
		
		count=0;
		for(int j=0; j<numChromStates; j++){
			for(int k=0; k<numFacBindingStates; k++){
				for(int itr=0; itr<itrs+1; itr++){
					Xaxes[count][itr] = itr;
					Yaxes[count][itr] = trainBjk[itr][j][k];
				}
				count++;
			}
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
			for(int j=0; j<itrs+1; j++){
				if(vec[i][j]> ret){
					ret = vec[i][j];
				}
			}
		}
		return ret;
	}
	
	/**
	 * Plot mehtod
	 * @param Xax
	 * @param Yax
	 * @param Xlab
	 * @param Ylab
	 * @param fp
	 */
	private void plot(double[][] Xax, double[][] Yax, String Xlab, String Ylab, File fp){
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
			
			for(int i=0; i<Xax.length; i++){
				for (int j=0; j<itrs; j++){
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
	
}
