package edu.psu.compbio.seqcode.gse.viz.metagenes.swing;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.event.ActionEvent;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;

import javax.imageio.ImageIO;
import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JFileChooser;
import javax.swing.JPanel;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

import edu.psu.compbio.seqcode.gse.viz.metagenes.BinningParameters;
import edu.psu.compbio.seqcode.gse.viz.metagenes.Profile;
import edu.psu.compbio.seqcode.gse.viz.metagenes.ProfilePaintable;
import edu.psu.compbio.seqcode.gse.viz.paintable.PaintableChangedEvent;
import edu.psu.compbio.seqcode.gse.viz.paintable.PaintableChangedListener;
import edu.psu.compbio.seqcode.gse.viz.paintable.PaintableScale;

public class ProfilePanel extends JPanel implements PaintableChangedListener {
	
	private Profile profile;
	private ProfilePaintable profilePainter;
	private PaintableScale scale;
	private int fontSize=12;
	private int border=20;
	private int lineHeight=20, lineWidth=4;
	private Color fontColor=Color.black;
	private Color peakColor=Color.blue;
	private String style = "Line";
	private boolean transparent = false;
	
	public ProfilePanel(Profile p, PaintableScale sc) { 
		profile = p;
		scale = sc;
		profilePainter = new ProfilePaintable(scale, profile);
		profilePainter.addPaintableChangedListener(this);
		
		setPreferredSize(new Dimension(500, 300));
	}
	
	public void updateFontSize(int size) {
		fontSize = size;
		repaint();
	}

	public void updateColor(Color c) {
		peakColor=c;
		repaint();
	}

	public void setStyle(String s) {
		style = s; 
		repaint();
	}
	public void setTransparent(boolean c){
		transparent = c;
	}
	
	public Action createSaveImageAction() { 
	    return new AbstractAction("Save Meta-Point Image...") { 
            /**
             * Comment for <code>serialVersionUID</code>
             */
            private static final long serialVersionUID = 1L;

            public void actionPerformed(ActionEvent e) { 
                String pwdName = System.getProperty("user.dir");
                JFileChooser chooser;
                if(pwdName != null) { 
                    chooser = new JFileChooser(new File(pwdName));
                } else {
                    chooser = new JFileChooser();
                }
                
                int v = 
                    chooser.showSaveDialog(null);
                if(v == JFileChooser.APPROVE_OPTION) { 
                    File f = chooser.getSelectedFile();
                    try {
                        saveImage(f, getWidth(), getHeight(), true);
                        //System.out.println("Saved Image [" + sImageWidth + " by " + sImageHeight +  "]");
                    } catch(IOException ie) {
                        ie.printStackTrace(System.err);
                    }
                }
                
            }
        };
	}
	public void saveImage(File f, int w, int h, boolean raster) 
    throws IOException { 
		if(raster){
			if(transparent)
				this.setOpaque(false);
			this.setSize(new Dimension(w, h));
			repaint();
	        BufferedImage im = 
	            new BufferedImage(w, h, BufferedImage.TYPE_INT_ARGB);
	        Graphics2D graphics = im.createGraphics();
	        graphics.setRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
	        this.print(graphics);
	        graphics.dispose();
	        ImageIO.write(im, "png", f);
		}else{
	        DOMImplementation domImpl =
	            GenericDOMImplementation.getDOMImplementation();
	        // Create an instance of org.w3c.dom.Document
	        Document document = domImpl.createDocument(null, "svg", null);
	        // Create an instance of the SVG Generator
	        SVGGraphics2D svgGenerator = new SVGGraphics2D(document);
	        svgGenerator.setSVGCanvasSize(new Dimension(w,h));
	        // Ask the test to render into the SVG Graphics2D implementation
	        if(!transparent){
	        	svgGenerator.setColor(Color.white);        
	        	svgGenerator.fillRect(0,0,w,h);
	        }
	        this.paintComponent(svgGenerator);
	
	        // Finally, stream out SVG to the standard output using UTF-8
	        // character to byte encoding
	        boolean useCSS = true; // we want to use CSS style attribute
	        FileOutputStream outStream = new FileOutputStream(f);
	        Writer out = new OutputStreamWriter(outStream, "UTF-8");
	        svgGenerator.stream(out, useCSS);
	        outStream.flush();
	        outStream.close();
		}
	}
	
	public Color getPeakColor(){return peakColor;}
	
	protected void paintComponent(Graphics g) {
		int w = getWidth(), h = getHeight();
		super.paintComponent(g);
		Graphics2D g2 = (Graphics2D)g;
		if(!transparent){
			g2.setColor(Color.white);
			g2.fillRect(0, 0, w, h);
		}
		profilePainter.setColor(peakColor);
		profilePainter.setStyle(style);
		profilePainter.paintItem(g, border, 0, w, h-border);
		
		int binPix = w / (profile.length()+1);
		
		//Paint labels & Y-axis stuff 
		FontMetrics metrics;
		g2.setFont(new Font("Arial", Font.PLAIN, fontSize));
		metrics = g2.getFontMetrics();
		g2.setColor(Color.black);
		if(profile.max()<10 && profile.max()>1)
			g2.drawString(String.format("%.2f", scale.getMax()), border/2, fontSize);
		else if(profile.max()>=1 || profile.max()==0)
			g2.drawString(String.format("%.0f", scale.getMax()), border/2, fontSize);
		else
			g2.drawString(String.format("%.2e", scale.getMax()), border/2, fontSize);
		if(profile.min()==0 || profile.min()<=-1)
			g2.drawString(String.format("%.0f", scale.getMin()), border/2, h-border-1);
		else
			g2.drawString(String.format("%.2e", scale.getMin()), border/2, h-border-1);
		String counter=String.format("%d datapoints", profile.getNumProfiles());
		g2.drawString(counter, w-border-metrics.stringWidth(counter), fontSize);
		//X-axis
		int xaxispos=h-border;
		if(profile.min()<0){
			double frac =scale.getMax()/(scale.getMax()-scale.getMin());
			xaxispos = (int)((double)((h-border-1))*frac);
			g2.drawString("0", border/2, xaxispos);
		}
		g2.setColor(Color.DARK_GRAY);
		g2.drawLine(border, h-border, (binPix*profile.length())+border, h-border);
		BinningParameters bps = profile.getBinningParameters();
		String minVal = String.format("%d", -1*(bps.getWindowSize()/2));
		String maxVal = String.format("%d", (bps.getWindowSize()/2));
		g2.drawString(maxVal, border+(binPix*profile.length())-metrics.stringWidth(maxVal), h);
		g2.drawString(minVal, border, h);
		
		//Title 
		
		//Draw marker line
		g2.setColor(Color.black);
		if(profile.isStranded()){
			g2.setStroke(new BasicStroke((float)lineWidth));	
			int [] a=new int [7];
			int [] b=new int [7];
			arrangeArrow(a, b, lineHeight, border+(binPix*(profile.length()/2)), border+(binPix*(profile.length()/2))+150, h-border);
            g2.drawPolyline(a, b, 7);
		}else{
            g2.fillRect(border+(binPix*(profile.length()/2))-lineWidth/2, h-lineHeight-(border/2), lineWidth, lineHeight);
		}
	}

	public void paintableChanged(PaintableChangedEvent pce) {
		repaint();
	}
	
	private void arrangeArrow(int[] a, int[] b, int height, int gx1, int gx2, int my) { 
        double arrowHt =0.1 *height;
        double arrowWd = 2;
        int a1, a2, a3;
        int startX = gx1;
        a1 = startX; 
        a2 = (int) Math.round(startX + (arrowWd * 6)); 
        a3 = (int) Math.round(startX + (arrowWd * 10)); 
        
        a[0] = a1; a[1] = a1;
        a[2] = a2; a[3] = a2;
        a[4] = a3; a[5] = a2; a[6] = a2;

        int b1 = (int) Math.round(my);
        int b2 = (int) Math.round(my - (arrowHt * 13));
        int b3 = (int) Math.round(my - (arrowHt * 10));
        int b4 = (int) Math.round(my - (arrowHt * 16));
        
        b[0] = b1; b[1] = b2;
        b[2] = b2; b[3] = b3; b[4] = b2;
        b[5] = b4; b[6] = b2;
    }
}
