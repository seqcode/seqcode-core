package org.seqcode.projects.seqview.paintable;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.*;

import org.seqcode.genome.location.Gene;
import org.seqcode.genome.location.Region;
import org.seqcode.projects.seqview.model.RegionExpanderModel;
import org.seqcode.utils.*;
import org.seqcode.viz.DynamicAttribute;


public class GenePainter extends RegionPaintable {
    
    private RegionExpanderModel<Gene> model;
    private NonOverlappingLayout<Gene> layout;
    private DynamicAttribute attrib;
    private double htRat, wdRat;
    private Vector<Gene> genes;    
    private GeneProperties props;

    public GenePainter(RegionExpanderModel<Gene> model) {
        super();
        layout = new NonOverlappingLayout<Gene>();
        this.model = model;
        attrib = DynamicAttribute.getGlobalAttributes();
        htRat = .03;
        wdRat = .03;
        model.addEventListener(this);
        props = new GeneProperties();
        initLabels();
    }

    public GeneProperties getProperties () {return props;}
    public void setProperties(GeneProperties p) {props = p;}
    
    public void clickedOnItem(ActionEvent e) {
        System.err.println(e);
        String geneName = e.getActionCommand();
    }
    
    public void cleanup() { 
        super.cleanup();
        model.removeEventListener(this);
    }

    public void removeEventListener(Listener<EventObject> l) {
        super.removeEventListener(l);
        if (!hasListeners()) {
            model.removeEventListener(this);
        }
    }

    public synchronized void eventRegistered(EventObject e) {
        
        if (e.getSource() == model &&
            model.isReady()) {            
        	System.out.println(this.getLabel());
            setCanPaint(true);
            setWantsPaint(true);
            genes = null;
            setLayoutGenes();
            notifyListeners();
        }
    }    

    public int getMaxVertSpace() { 
        int numTracks = layout.getNumTracks();
        return Math.min(Math.max(40,numTracks * 12),120);
    }
    public int getMinVertSpace() { 
        int numTracks = layout.getNumTracks();
        return Math.min(Math.max(40,numTracks * 12),120);
    }

    private void setLayoutGenes() {
        if (canPaint() && genes == null) {
            Iterator<Gene> itr = model.getResults();
            genes = new Vector<Gene>();
            while(itr.hasNext()) { genes.add(itr.next()); }
            layout.setRegions(genes); 
        }
    }
    
    public int getNumTracks() { return layout.getNumTracks(); }

    public void paintItem(Graphics2D g, 
                          int x1, int y1, 
                          int x2, int y2) {
        
        if (!canPaint()) {
            return;
        }
        int w = x2 - x1, h = y2 - y1;
        int my = y1 + (h / 2);
 
        int numTracks = layout.getNumTracks();
        
        int trackHeight = h;

        if(numTracks > 1) { 
            trackHeight = (int)Math.floor((double)h / (double)numTracks);
        }
        
        int geneHeight = Math.max(2, (int)Math.floor((double)trackHeight * 0.80));
        int halfHeight = trackHeight / 2;
        int halfGeneHeight = geneHeight / 2;

        for(int track=0; track<numTracks; track++) {
            g.setColor(Color.black);
            my = y1 + (halfHeight + (track * trackHeight));
            g.drawLine(x1, my, x2, my);
        }
        
        Region region = model.getRegion();
        int rs = region.getStart(), re = region.getEnd();
        int rw = re - rs + 1;

        Font oldFont = g.getFont();
        g.setFont(attrib.getRegionLabelFont(w,h));
        FontMetrics fontmetrics = g.getFontMetrics();
        //--------------------------------------------------
        
        double xScale = (double)w / (double)rw;

        clearLabels();
        boolean drewAnything = false;
        
        for(Gene gene : genes) { 
            int gs = gene.getStart(), ge = gene.getEnd();
            int track = 0;
            if(!layout.hasTrack(gene)) { 
                System.err.println("No track assigned to gene: " + gene.getName());
            } else { 
                track = layout.getTrack(gene);
            }
            
            char strand = gene.getStrand();
            
            my = y1 + (halfHeight + (track * trackHeight));
            
            int gx1 = x1 + (int)Math.round((double)(gs - rs) * xScale);
            int gx2 = x1 + (int)Math.round((double)(ge - rs) * xScale);
            int gxw = gx2 - gx1;
            int depth = Math.min(halfGeneHeight, gxw/3);

            int gleft = Math.max(gx1, x1), gright = Math.min(gx2, x2);
            int gtop = my - (halfGeneHeight/2), gbottom = gtop + (geneHeight/2);

            //g.setColor(Color.white);
            //g.fillRect(gx1, my - (halfGeneHeight / 2), gxw, geneHeight / 2);
            g.setColor(Color.GRAY);
            int rectwidth = gright - gleft;
            int rectheight = gbottom - gtop;
            g.fillRect(gleft, gtop, rectwidth, gbottom - gtop);                      
            drewAnything = true;

            g.setColor(Color.black);
            //g.drawRect(gx1, my - (halfGeneHeight / 2), gxw, geneHeight / 2);
            g.drawLine(gleft, gtop, gright, gtop);
            g.drawLine(gleft, gbottom, gright, gbottom);
            if(gx1 >= x1) { g.drawLine(gx1, gtop, gx1, gbottom); }
            if(gx2 <= x2) { g.drawLine(gx2, gtop, gx2, gbottom); }
            
            // Somewhat prettier gene-name output
            //g.drawString(gene.getName(), gx1 + (gxw/3), my);
            int nx = Math.max(x1 + 3, gx1 + 3);  // gotta do the Math.max(), to make sure the name is on the screen.
            int ny = my + (halfGeneHeight / 2) - 2;
            int fontsize = g.getFont().getSize();
            
            ArrayList<String> aliases = new ArrayList<String>();
            aliases.add(gene.getID());
            aliases.add(gene.getName());
            aliases.addAll(gene.getAliases());
            String first = aliases.remove(0);
            String todraw = null;
            addLabel(gleft,gtop,rectwidth,gbottom-gtop,first);
            if (first.length() * fontsize < rectwidth) {
                todraw = first;
            }
           
            for (String s : aliases) {
                String newtodraw = todraw + ", " + s;
                if (props.DrawAllGeneNames && todraw != null && 
                    fontmetrics.charsWidth((newtodraw).toCharArray(),0,newtodraw.length()) < rectwidth) {
                    todraw = newtodraw;
                } 
                addLabel(gleft,gtop,rectwidth,gbottom-gtop,s);
            }
            if (todraw != null && props.DrawGeneNames) {
                g.drawString(todraw, nx, ny);
            }
            
            // arrows go here
            
            //double arrowHt = htRat * h;
            //double arrowWd = wdRat * w;
            double arrowHt = htRat * geneHeight;
            double arrowWd = wdRat * gxw;
            
            if(arrowWd > arrowHt) { arrowWd = arrowHt; }
            
            int a[];
            // forward arrow
            if(strand == '+') {
            	// original arrows
                //g.drawLine(gx1, my-halfHeight, gx1 + depth, my);
                //g.drawLine(gx1 + depth, my, gx1, my+halfHeight);
            	int startX = gx1 - 5;
            	//
                int a1 = startX; //(int) Math.round(gx1 / (gx1+depth));
                //int a2 = (int) Math.round(startX + (arrowWd * 6)); //(int) Math.round((gx1 + 25) / (gx1+depth));
                //int a3 = (int) Math.round(startX + (arrowWd * 8)); //(int) Math.round((gx1 + 35) / (gx1+depth));

                int a2 = (int) Math.round(startX + (arrowWd * 8)); 
                int a3 = (int) Math.round(startX + (arrowWd * 12)); 
        		
                int[] t = {a1, a1, a2, a2, a3, a2, a2};
                a = t;
            }
            // backward arrow
            else { 
            	// original arrows
                //g.drawLine(gx2, my-halfHeight, gx2-depth, my);
                //g.drawLine(gx2-depth, my, gx2, my+halfHeight);
            	int startX = gx2 + 5;
                int a1 = startX; //(int) Math.round(gx1 / (gx1+depth));
                //int a2 = (int) Math.round(startX - (arrowWd * 6)); //(int) Math.round((gx1 + 25) / (gx1+depth));
                //int a3 = (int) Math.round(startX - (arrowWd * 8)); //(int) Math.round((gx1 + 35) / (gx1+depth));

                int a2 = (int) Math.round(startX - (arrowWd * 8)); 
                int a3 = (int) Math.round(startX - (arrowWd * 12)); 

                int[] t = {a1, a1, a2, a2, a3, a2, a2};
                a = t;
            }
            
            /*
              int b1 = (int) Math.round(my);
              int b2 = (int) Math.round(my - (arrowHt * 13));
              int b3 = (int) Math.round(my - (arrowHt * 10));
              int b4 = (int) Math.round(my - (arrowHt * 16));
            */

            int b1 = (int) Math.round(my);
            int b2 = (int) Math.round(my - (arrowHt * 13));
            int b3 = (int) Math.round(my - (arrowHt * 10));
            int b4 = (int) Math.round(my - (arrowHt * 16));
    		
            int[] b = {b1, b2, b2, b3, b2, b4, b2};
            g.drawPolyline(a, b, 7);
        }

        if (drewAnything && props.DrawTrackLabel) {
            g.setColor(Color.BLACK);
            g.setFont(attrib.getLargeLabelFont(w,h));
            g.drawString(getLabel(),x1,y2);
        }
        g.setFont(oldFont);
        
    }

}

