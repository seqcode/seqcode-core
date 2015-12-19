package edu.psu.compbio.seqcode.gse.datasets.motifs;

import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.font.*;


public class WeightMatrixPainter {
	private static Color AColor = Color.RED;
	private static Color CColor = Color.BLUE;
	private static Color GColor = Color.ORANGE;
	private static Color TColor = Color.GREEN;
	private static Color GapColor = Color.WHITE;
	private static Color NColor = Color.GRAY;
	public final static int X_MARGIN = 10;
	public final static int Y_MARGIN = 2;
	public final static int YLABEL_SIZE = 12;
	public void paint(WeightMatrix wm, Graphics g, int x1, int y1, int x2, int y2) {this.paint(wm, g, x1, y1, x2, y2, null, false);}
	public void paint(WeightMatrix wm, Graphics g, int x1, int y1, int x2, int y2, String name) {this.paint(wm, g, x1, y1, x2, y2, name, false);}
	/* paints a representation of a weight matrix wherein the height of the letters roughly indicates
       their relative probability at each base.  The image is painted in g in the
       bounding box defined by the upper left corner x1,y1 and lower right corner x2,y2 */
    public void paint(WeightMatrix wm, Graphics g, int x1, int y1, int x2, int y2, String name, boolean drawAxes) {
        wm.toLogOdds();

        String label = wm.toString();
        if(name!=null)
        	label = name;
        
        int w = x2 - x1;
        int h = y2 - y1;
        Font labelFont = new Font("Arial",Font.PLAIN,YLABEL_SIZE);
        g.setFont(labelFont);
        FontMetrics fontmetrics = g.getFontMetrics();
        LineMetrics linemetrics = fontmetrics.getLineMetrics(label,g);
        g.setColor(Color.BLACK);
        g.drawString(label,x1+X_MARGIN + w/2 - fontmetrics.charsWidth(label.toCharArray(),0,label.length()) / 2,y2-Y_MARGIN);
        int labelHeight = fontmetrics.getHeight() + Y_MARGIN;
        int posAxisHeight = drawAxes ? fontmetrics.getHeight() + Y_MARGIN : 0;
        int boxHeight = h - labelHeight - posAxisHeight;

        if(drawAxes){
        	int block = (w-X_MARGIN)/wm.length();
        	int xpos = x1+X_MARGIN+(block/2);
        	for (int pos = 0; pos < wm.length(); pos++) {
        		String pstr = new Integer((pos+1)).toString();
        		g.drawString(pstr,xpos - fontmetrics.stringWidth(pstr)/2,y2-labelHeight-Y_MARGIN);
        		xpos+=block;
        	}
        	g.setColor(Color.black);
        	int ypos = y2 - labelHeight-posAxisHeight;
        	g.drawLine(x1+X_MARGIN, ypos, x2, ypos);
        	g.drawLine(x1+X_MARGIN, y1, x1+X_MARGIN, ypos);
        	g.drawString("2",x1,y1+fontmetrics.getAscent());
        	g.drawString("1",x1,y1+(boxHeight/2)+fontmetrics.getAscent());
        }
        Font baseFont = new Font("Arial",Font.BOLD, boxHeight );
        int pixelsPerLetter = (w-X_MARGIN*2) / wm.length();
        double xfactor = ((float)pixelsPerLetter) / baseFont.getSize() * 1.3;
        //        System.err.println("Xfactor is " + xfactor + " and base sizeis " + baseFont.getSize());
        double vals[] = new double[4];
        for (int pos = 0; pos < wm.length(); pos++) {
            Character[] letters = WeightMatrix.getLetterOrder(wm,pos);
            double total = 0;
            int ypos = y2 - labelHeight-posAxisHeight;
            for (int j = 3; j >= 0; j--) {
                vals[j] = Math.exp(wm.matrix[pos][letters[j]]);
                total += vals[j];
            }
            
            double bits = 0;
            for (int j = 3; j >= 0; j--) {
                vals[j] = vals[j] / total;
                bits -= vals[j] * (Math.log(vals[j]) / Math.log(2));
            }
            double totalscale = (2.0 - bits) / 2.0;
            for (int j = 3; j >= 0; j--) {
                if (letters[j] == 'A') {
                    g.setColor(AColor);
                } else if (letters[j] == 'C') {
                    g.setColor(CColor);
                } else if (letters[j] == 'G') {
                    g.setColor(GColor);
                }  else if (letters[j] == 'T') {
                    g.setColor(TColor);
                }        
                double val = vals[j];
                AffineTransform transform = new AffineTransform(xfactor,0,0,val * totalscale,0,0);
                Font thisFont = baseFont.deriveFont(transform);
                int letterHeight = (int) (thisFont.getSize() * val * totalscale);        
                if (letterHeight>1){
	                g.setFont(thisFont);
	                int offset = 0;
	                /*if (letters[j] == 'G') 
	                    offset = -(int)(pixelsPerLetter*0.05);
	                else if (letters[j] == 'T')
	                    offset = (int)(pixelsPerLetter*0.05);
	                    */
	                g.drawString(letters[j].toString(),x1+X_MARGIN +offset+ pos * pixelsPerLetter,ypos);
                }else if (letterHeight==1){
                	g.fillRect(x1+X_MARGIN+ pos*pixelsPerLetter, ypos, (int)(pixelsPerLetter*0.9), letterHeight);
                }
                ypos -= letterHeight;
            }                
        }
    }

}
