package edu.psu.compbio.seqcode.projects.shaun.viz;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Polygon;
import java.awt.geom.AffineTransform;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.JFrame;

import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.viz.paintable.AbstractPaintable;
import edu.psu.compbio.seqcode.gse.viz.paintable.PaintableFrame;
import edu.psu.compbio.seqcode.projects.shaun.MultidataSpatialPaintable;

public class HoxGenePainter extends AbstractPaintable{

	private MultidataSpatialPaintable painter;
	private static PaintableFrame plotter;
	private static int deftScreenSizeX=800, deftScreenSizeY=1200;
	private int screenSizeX, screenSizeY;
	HashMap<String, ArrayList<Double>> expression = new HashMap<String, ArrayList<Double>>();
	private final int geneBoxHeight=34, geneBoxWidth=96;
	private final int geneArcWidth=6, geneArcHeight=6;
	private final int expBoxHeight=30, expBoxWidth = 15;
	private final int clusterSpacing = 120;
	private final int colorbarHeight=20, colorbarWidth=200;
	private final double maxExp=7.0, minExp=-5.0;
	private Color expMaxColor = Color.yellow;
	private Color expMidColor = Color.black;
	private Color expMinColor = Color.blue;
	private Color geneColor = Color.gray;
	private final int topBorder=50, bottomBorder=50;
	private final int leftBorder=25, rightBorder=25;
	private int topBound, bottomBound, leftBound, rightBound, baseLine, midLine;
	private static ArrayList<String> timeLabels= new ArrayList<String>();
	private static int numExpr =6;
	private static int Hoxes[][]={
		{1,3,4,8,9,10,11,12,13},
		{4,5,6,8,9,10,11,12,13},
		{1,2,3,4,5,6,7,8,9,13},
		{1,2,3,4,5,6,7,9,10,11,13}		
	};
	
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("data")) { 
            System.err.println("Usage:\n " +
                               "HoxPainter " +
                               "--data <file name> ");
            return;
        }
        String dfile = ap.getKeyValue("data");
        HashMap<String, ArrayList<Double>> expr = loadFile(dfile);
			
        //Paint the picture
        HoxGenePainter painter = new HoxGenePainter(expr);
		plotter = new PaintableFrame("Hox Cartoon", painter);
		plotter.setSize(deftScreenSizeX, deftScreenSizeY);
		plotter.setVisible(true);
		plotter.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
	}
	
	public HoxGenePainter(HashMap<String, ArrayList<Double>> x){
		expression = x;
	}
	
	public void paintItem (Graphics g, int x1, int y1, int x2, int y2){
		Graphics2D g2d = (Graphics2D)g;
		FontMetrics metrics = g2d.getFontMetrics();
		screenSizeX = x2-x1;
		screenSizeY = y2-y1;
		g2d.setColor(Color.white);
		g2d.fillRect(0, 0, screenSizeX, screenSizeY);
		topBound = topBorder+140;
		bottomBound = screenSizeY-bottomBorder;
		leftBound = leftBorder+30;
		rightBound = screenSizeX-rightBorder;
		baseLine = (topBound+bottomBound)/2;
		midLine = (leftBound+rightBound)/2;
	
		
		//Draw the cluster lines & background gene boxes
		for(int i=0; i<4; i++){
			String currHox="Hoxa";
			if(i==3){currHox = "Hoxa";}
			if(i==2){currHox = "Hoxb";}
			if(i==1){currHox = "Hoxc";}
			if(i==0){currHox = "Hoxd";}
			//Background lines
			g2d.setColor(geneColor);
			g2d.setStroke(new BasicStroke(5.0f));
			int xPos = (midLine-(3*clusterSpacing/2))+(i*clusterSpacing);
			int yEnd = topBound+(geneBoxHeight*13)-2;
			g2d.drawLine(xPos, topBound+2, xPos, yEnd);
			//Background rounded boxes
			for(int j=1; j<=13; j++){
				boolean found =false;
				for(int z=0; z<Hoxes[i].length; z++){
					if(Hoxes[i][z]==j)
						found=true;
				}
				if(found){
					int boxX = xPos-(geneBoxWidth/2);
					int boxY = topBound+geneBoxHeight*(j-1);
					g2d.setColor(geneColor);
					g2d.setStroke(new BasicStroke(2.0f));
					g2d.fillRoundRect(boxX, boxY, geneBoxWidth, geneBoxHeight, geneArcWidth, geneArcHeight);
					g2d.setColor(Color.darkGray);
					g2d.drawRoundRect(boxX, boxY, geneBoxWidth, geneBoxHeight, geneArcWidth, geneArcHeight);
					
					//Default expression boxes
					int eY = boxY+(geneBoxHeight-expBoxHeight)/2;
					for(int e=0; e<numExpr; e++){
						int eX = xPos - ((numExpr*expBoxWidth)/2) +(e*expBoxWidth);
						g2d.setColor(Color.lightGray);
						g2d.fillRect(eX, eY,expBoxWidth, expBoxHeight);
						g2d.setColor(Color.white);
						g2d.setStroke(new BasicStroke(1.0f));
						g2d.drawRect(eX, eY,expBoxWidth, expBoxHeight);
					}
					
					//Find appropriate expression values
					String currHoxGene = currHox+j;
					ArrayList<Double> evals = expression.get(currHoxGene);
					if(evals != null){
						for(int e=0; e<numExpr; e++){
							Double v = evals.get(e);
							Color currCol = expColor(v);
							int eX = xPos - ((numExpr*expBoxWidth)/2) +(e*expBoxWidth);
							g2d.setColor(currCol);
							g2d.fillRect(eX, eY,expBoxWidth, expBoxHeight);
							g2d.setColor(Color.white);
							g2d.setStroke(new BasicStroke(1.0f));
							g2d.drawRect(eX, eY,expBoxWidth, expBoxHeight);
						}
					}
				}
			}
		}
		
		//Paralog numbers
		g2d.setColor(Color.black);
		g2d.setFont(new Font("Ariel", Font.BOLD, 20));
		metrics = g2d.getFontMetrics();
		int xPos = (midLine-(5*clusterSpacing/2));
		g2d.drawString("Paralog", xPos-(metrics.stringWidth("Paralog")/2), topBound-(geneBoxHeight/2));
		for(int j=1; j<=13; j++){
			g2d.setFont(new Font("Ariel", Font.PLAIN, 20));
			//Vert g2d.drawString(String.format("%d",j), xPos, topBound+(geneBoxHeight*(j-1))+(geneBoxHeight/2));
			AffineTransform oldtrans = g2d.getTransform();
	        AffineTransform newtrans = new AffineTransform();
	        newtrans.translate(xPos, topBound+(geneBoxHeight*(j-1))+(geneBoxHeight/2));
	        newtrans.rotate(Math.toRadians(90));
	        g2d.setTransform(newtrans);
	        g2d.drawString(String.format("%d",j), 0,0);
	        g2d.setTransform(oldtrans);
		}
		//Big Hox Titles
		g2d.setColor(Color.black);
		g2d.setFont(new Font("Ariel", Font.BOLD, 26));
		metrics = g2d.getFontMetrics();
		int yPos = topBound-90;
		for(int i=0; i<4; i++){
			xPos = (midLine-(3*clusterSpacing/2))+(i*clusterSpacing);
			String text="HoxA";
			if(i==3){text = "HoxA";}
			if(i==2){text = "HoxB";}
			if(i==1){text = "HoxC";}
			if(i==0){text = "HoxD";}
			g2d.drawString(text, xPos-(metrics.stringWidth(text)/2), yPos);
		}	
		//Colorbar
		drawExpColorBar(g2d, midLine-(colorbarWidth/2), topBound-160);
		//Time Labels
		for(int i=0; i<4; i++){
			xPos = (midLine-(3*clusterSpacing/2))+(i*clusterSpacing);
			g2d.setColor(Color.black);
			AffineTransform oldtrans = g2d.getTransform();
	        AffineTransform newtrans = new AffineTransform();
	        newtrans.translate(xPos, topBound-10);
	        //Vert
	        //newtrans.rotate(Math.toRadians(-90));
	        newtrans.rotate(Math.toRadians(90));
	        g2d.setTransform(newtrans);
	        g2d.setFont(new Font("Ariel", Font.PLAIN, 16));
	        metrics = g2d.getFontMetrics();
	        g2d.drawString("Day",-1*metrics.stringWidth("Day")/2,-1*(geneBoxWidth/2));
	        g2d.setFont(new Font("Ariel", Font.PLAIN, 16));
	        metrics = g2d.getFontMetrics();
	        for(int e=0; e<numExpr; e++){
	        	//Vert			int etY = (e*expBoxWidth)-((numExpr*expBoxWidth)/2)+(expBoxWidth/2)+(metrics.getHeight()/2);
	        	int etY = ((numExpr-e-1)*expBoxWidth)-((numExpr*expBoxWidth)/2)+(expBoxWidth/2)+(metrics.getHeight()/2);
				String text = timeLabels.get(e);
				g2d.drawString(text,0,etY);
			}
	        g2d.setTransform(oldtrans);
		}
		//Anterior/posterior arrow
		xPos = (midLine+(5*clusterSpacing/2));
		g2d.setColor(Color.black);
		g2d.setStroke(new BasicStroke(5.0f));
		g2d.drawLine(xPos, topBound+((geneBoxHeight*13)/4), xPos, topBound+((geneBoxHeight*13)/4*3));
		drawArrow(g2d, xPos, topBound+((geneBoxHeight*13)/4)+10, xPos, topBound+((geneBoxHeight*13)/4), 5.0f);
		g2d.setFont(new Font("Ariel", Font.PLAIN, 20));
		metrics = g2d.getFontMetrics();
		AffineTransform oldtrans = g2d.getTransform();
        AffineTransform newtrans = new AffineTransform();
        newtrans.translate(xPos, topBound+((geneBoxHeight*13)/4));
        //Vert  newtrans.rotate(Math.toRadians(-90));
        newtrans.rotate(Math.toRadians(90));
        
        g2d.setTransform(newtrans);
        //Vert g2d.drawString("3' (anterior)",10,metrics.getHeight()/2);
        //Vert g2d.drawString("(posterior) 5'",0-((geneBoxHeight*13)/2)-10-metrics.stringWidth("(posterior) 5'"),metrics.getHeight()/2);
        g2d.drawString("(anterior) 3'",10-metrics.stringWidth("(anterior) 3'"),metrics.getHeight()/2);
        g2d.drawString("5' (posterior)",((geneBoxHeight*13)/2)+10,metrics.getHeight()/2);
        g2d.setTransform(oldtrans);
	}
	
	private Color expColor(double v){
		Color c;
		if(v>0){
			Color maxColor = expMaxColor;
			Color minColor = expMidColor;
			
			double sVal = v>maxExp ? 1 : (v)/(maxExp);
			int red = (int)(maxColor.getRed() * sVal + minColor.getRed() * (1 - sVal));
		    int green = (int)(maxColor.getGreen() * sVal + minColor.getGreen() * (1 - sVal));
		    int blue = (int)(maxColor.getBlue() *sVal + minColor.getBlue() * (1 - sVal));
		    c = new Color(red, green, blue);
		}else{
			Color maxColor = expMidColor;
			Color minColor = expMinColor;
			double sVal = v<minExp ? 1 : (minExp-v)/(minExp);
			int red = (int)(maxColor.getRed() * sVal + minColor.getRed() * (1 - sVal));
	        int green = (int)(maxColor.getGreen() * sVal + minColor.getGreen() * (1 - sVal));
	        int blue = (int)(maxColor.getBlue() *sVal + minColor.getBlue() * (1 - sVal));
	        c = new Color(red, green, blue);				
		}
		return(c);
	}
	private void drawExpColorBar(Graphics2D g2d, int x, int y){
		//Draw colors 
		GradientPaint colorbar = new GradientPaint(x, y, expMinColor, x+colorbarWidth/2, y, expMidColor, false);
		g2d.setPaint(colorbar);
		g2d.fillRect(x, y, colorbarWidth/2, colorbarHeight);
		colorbar = new GradientPaint(x+colorbarWidth/2, y, expMidColor, x+colorbarWidth, y, expMaxColor, false);
		g2d.setPaint(colorbar);
		g2d.fillRect(x+(colorbarWidth/2), y, colorbarWidth/2, colorbarHeight);
		
		//Draw border
		g2d.setPaint(Color.black);
		g2d.setColor(Color.black);
		g2d.setStroke(new BasicStroke(3.0f));
		g2d.drawRect(x, y, colorbarWidth, colorbarHeight);
		
		//Legend
		g2d.setFont(new Font("Ariel", Font.PLAIN, 12));
		FontMetrics metrics = g2d.getFontMetrics();
		int textY = y+colorbarHeight+ (metrics.getHeight());
		g2d.drawString("0", x+(colorbarWidth/2)-(metrics.stringWidth("0")/2), textY);
		g2d.drawString(String.format("%.1f",minExp), x-(metrics.stringWidth(String.format(".1f",minExp))/2), textY);
		g2d.drawString(String.format("%.1f",maxExp), x+colorbarWidth-(metrics.stringWidth(String.format(".1f",maxExp))/2), textY);
		
		//Title
		g2d.setFont(new Font("Ariel", Font.ITALIC, 12));
		metrics = g2d.getFontMetrics();
		g2d.drawString("log-foldchange", x+(colorbarWidth/2)-(metrics.stringWidth("log-foldchange")/2), y- (metrics.getHeight())/2);
	}
	private static HashMap<String, ArrayList<Double>> loadFile(String inF){
		HashMap<String, ArrayList<Double>> ex = new HashMap<String, ArrayList<Double>>();
		try{
			File aFile = new File(inF);
			if(aFile.isFile()){
				BufferedReader reader;
				reader = new BufferedReader(new FileReader(aFile));
				String line;
				//Labels
				line= reader.readLine();
				String [] tokens = line.split("\\s+");
				for(int i=1; i<tokens.length; i++){
					timeLabels.add(tokens[i]);
				}numExpr = timeLabels.size();
				//Expression
				while((line= reader.readLine())!=null){
					tokens = line.split("\\s+");
					String name = tokens[0];
					ArrayList<Double> vals = new ArrayList<Double>();
					for(int i=1; i<tokens.length; i++){
						vals.add(new Double(tokens[i]));
					}
					ex.put(name, vals);
				}
				reader.close();
			}
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return(ex);	
	}
	
	public static void drawArrow(Graphics2D g2d, int xCenter, int yCenter, int x, int y, float stroke) {
		double aDir=Math.atan2(xCenter-x,yCenter-y);
		g2d.drawLine(x,y,xCenter,yCenter);
		g2d.setStroke(new BasicStroke(stroke));					// make the arrow head solid even if dash pattern has been specified
		Polygon tmpPoly=new Polygon();
		int i1=12+(int)(stroke*2);
		int i2=6+(int)stroke;							// make the arrow head the same size regardless of the length length
		tmpPoly.addPoint(x,y);							// arrow tip
		tmpPoly.addPoint(x+xCor(i1,aDir+.5),y+yCor(i1,aDir+.5));
		tmpPoly.addPoint(x+xCor(i2,aDir),y+yCor(i2,aDir));
		tmpPoly.addPoint(x+xCor(i1,aDir-.5),y+yCor(i1,aDir-.5));
		tmpPoly.addPoint(x,y);							// arrow tip
		g2d.drawPolygon(tmpPoly);
		g2d.fillPolygon(tmpPoly);						// remove this line to leave arrow head unpainted
   	}
    private static int yCor(int len, double dir) {return (int)(len * Math.cos(dir));}
	private static int xCor(int len, double dir) {return (int)(len * Math.sin(dir));}
}
