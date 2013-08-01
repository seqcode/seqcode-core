/*
 * Created on Mar 28, 2006
 */
package edu.psu.compbio.seqcode.gse.seqview.paintable;

import java.awt.*;
import java.text.*;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;

/**
 * @author tdanford
 */
public class HashMarkPaintable extends RegionPaintable {
    
    private static NumberFormat nf;
    private PaintableProperties props;
    private boolean labelAbove;
    static {
        nf = NumberFormat.getInstance();
    }    

    public HashMarkPaintable () {
        super();
        props = new PaintableProperties();
        labelAbove = true;
    }

    public PaintableProperties getProperties() {
        return props;
    }

    public void setLabelAbove (boolean b) {labelAbove = b;}

    public boolean canPaint() {
        return true;
    }

    public void setRegion(Region r) {
        super.setRegion(r);
        setCanPaint(true);
        setWantsPaint(true);
        notifyListeners();
    }

    public int getMaxVertSpace() {return 25;}
    public int getMinVertSpace() {return 25;}
    
    /* (non-Javadoc)
     * @see edu.psu.compbio.seqcode.gse.viz.paintable.Paintable#paintItem(java.awt.Graphics, int, int, int, int)
     */
    public void paintItem(Graphics2D g, int x1, int y1, int x2, int y2) {
        int rs = getRegion().getStart(), re = getRegion().getEnd();
        int rw = re - rs + 1;
        int h = y2 - y1;
        int w = x2 - x1;
        
        double log10 = Math.log(10.0);
        double logScale = Math.log((double)rw) / log10;
        Font oldFont = g.getFont();
        g.setFont(new Font("Arial",Font.PLAIN,(int)(h * .4)));
        FontMetrics fontmetrics = g.getFontMetrics();
        char[] nines = {'9','9','9','9','9','9','9','9','9','9','9','9','9'};


        int tickSize = Math.max(1, pow10((int)Math.floor(logScale) - 2));
        int tickDigits = (int)Math.ceil(Math.log(re / tickSize) / log10) + 1;
        int spaceUsed = fontmetrics.charsWidth(nines,0,tickDigits);
        int spaceAvail = (int) (.5 * w / (rw / ((double)tickSize)));

        while(spaceUsed > spaceAvail) {
            tickSize *= 10;
            tickDigits = (int)Math.ceil(Math.log(re / tickSize) / log10) + 1;
            spaceUsed = fontmetrics.charsWidth(nines,0,tickDigits);
            spaceAvail = (int) (.7 * w / (rw / ((double)tickSize)));
        }

        tickSize = Math.max(tickSize, 1);
        
        int hh = h/2;
        int my, ly;
        if (labelAbove) {
            my = y1 + hh;
            ly = my - 2;
        } else {
            my = y1;
            ly = y2;
        }
        g.setColor(Color.white);
        g.fillRect(x1, my, w, hh);
        g.setColor(Color.DARK_GRAY);
        g.drawRect(x1, my, w, hh);
        
        int bpStart = tickSize * (int)Math.ceil((double)rs / (double)tickSize);
        double scale = (double)w / (double)rw;
        int tickPixWidth = (int)Math.round((double)tickSize * scale);
        
        boolean paintBlack = true;
        String axisLabel = getRegion().getChrom() + " x " + String.valueOf(nf.format(tickSize));
        if (!axisLabel.matches("^chr.*")) {
            axisLabel = "chr" + axisLabel;
        }
        g.drawString(axisLabel, x1 + 2, ly);

        for(int bp = bpStart; bp <= re; bp += tickSize) { 
            int xOffset = x1 + (int)Math.round((double)(bp-rs) * scale);
            
            if(paintBlack) {
            	g.setColor(Color.black);
                g.fillRect(xOffset, my, tickPixWidth, hh);
            }
            
            g.setColor(paintBlack ? Color.white : Color.black);
            g.drawString(String.valueOf(bp / tickSize), xOffset, my + hh - 1);

            paintBlack = !paintBlack;
        }
        g.setFont(oldFont);
    }
    
    private int pow10(int exp) { 
        int value = 1;
        for(int i = 0; i < exp; i++) { 
            value *= 10;
        }
        return value;
    }
}

