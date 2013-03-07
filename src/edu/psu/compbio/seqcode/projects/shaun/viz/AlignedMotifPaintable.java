package edu.psu.compbio.seqcode.projects.shaun.viz;

import java.sql.SQLException;
import java.util.*;
import java.awt.*;



import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.chipseq.*;
import edu.psu.compbio.seqcode.gse.utils.ArrayUtils;
import edu.psu.compbio.seqcode.gse.utils.Closeable;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.viz.colors.Coloring;
import edu.psu.compbio.seqcode.gse.viz.paintable.*;

public class AlignedMotifPaintable extends AbstractPaintable {
	
	private Region region;
	private LinkedList<Pair<String,StrandedRegion>> motifs;
	private ThinOverlapPaintable upper, lower;
	
	private Color watson, crick, watsonLine, crickLine;
	private Color rareColDR5 = new Color(0,153,0);
	private Color rareColDR2 = new Color(153,153,255);
	private Color rareColDR2Line, rareColDR5Line;
	private boolean reverseIt=false;
	private int motifWidth=5;
	
	public AlignedMotifPaintable(Region reg, Collection<Pair<String,StrandedRegion>> mtfs, ThinOverlapPaintable up, ThinOverlapPaintable down, boolean reverse, int motifWidth) { 
		region = reg;
		motifs = new LinkedList<Pair<String,StrandedRegion>>(mtfs);
		upper = up;
		lower = down;
		watson = Color.green;
		crick = Color.blue;
		reverseIt = reverse;
		this.motifWidth = motifWidth;
		int[] rgba = Coloring.rgba(watson);
		rgba[3] = 25;
		watsonLine = Coloring.asColor(rgba);
		rgba = Coloring.rgba(crick);
		rgba[3] = 25;
		crickLine = Coloring.asColor(rgba);
		rgba = Coloring.rgba(rareColDR2);
		rgba[3] = 25;
		rareColDR2Line = Coloring.asColor(rgba);
		rgba = Coloring.rgba(rareColDR5);
		rgba[3] = 25;
		rareColDR5Line = Coloring.asColor(rgba);	
	}

	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		int h = y2-y1, w = x2-x1;
		
		Graphics2D g2 = (Graphics2D)g;
		
		int[] upperProfile = upper.getDataProfile();
		int[] lowerProfile = lower.getDataProfile();

		/*Watson/Crick differences
		for(Pair<String,StrandedRegion> mpair : motifs) { 
			StrandedRegion motif = mpair.getLast();
			int midbp = motif.getMidpoint().getLocation();
			int x = x1 + getX(midbp, w);
			
			int uy = upperProfile[x-x1];
			int ly = lowerProfile[x-x1];
			
			g.setColor(motif.getStrand() == '+' ? watsonLine : crickLine);
			g.drawLine(x, uy+1, x, ly-1);
			
			g.setColor(motif.getStrand() == '+' ? watson : crick);
			drawMotif(g2, motif.getStrand() == '+', x, y1, y2);
		}*/
		//Motif type differences
		HashMap<String, Color> motifsProcessed= new HashMap<String, Color>();
		for(Pair<String, StrandedRegion> mpair : motifs){
			if(!motifsProcessed.containsKey(mpair.car())){
        		if(mpair.car().equals("DR2"))
        			motifsProcessed.put(mpair.car(), rareColDR2);
        		else if(mpair.car().equals("DR5"))
        			motifsProcessed.put(mpair.car(), rareColDR5);
        		else
        			motifsProcessed.put(mpair.car(), new Color((int)(Math.random()*256),(int)(Math.random()*256),(int)(Math.random()*256)));
        	}
        	Color currCol = motifsProcessed.get(mpair.car());
        	StrandedRegion motif = mpair.getLast();
			int midbp = motif.getMidpoint().getLocation();
			int x = x1 + getX(midbp, w);
			int uy = upperProfile[x-x1];
			int ly = lowerProfile[x-x1];
			if(mpair.car().equals("DR2"))
				g.setColor(rareColDR2Line);
			if(mpair.car().equals("DR5"))
				g.setColor(rareColDR5Line);
			g.drawLine(x, uy+1, x, ly-1);
			
			g.setColor(currCol);
			drawMotif(g2, motif.getStrand() == '+', x, y1, y2);
        	       	
        }
	}
	
	private void drawMotif(Graphics2D g2, boolean strand, int x, int y1, int y2) { 
		int[] xs = new int[3];
		int[] ys = new int[3];
		
		int h = y2-y1;
		int mid = y1 + (h/2);
		
		if((strand && !reverseIt) || (!strand && reverseIt)){
    		int [] xPoints = {x, x, x+motifWidth};
    		int [] yPoints = {y1, y2, (y1+y2)/2};
    		g2.fillPolygon(xPoints, yPoints, 3);
    	}else{
    		int [] xPoints = {x, x, x-motifWidth};
    		int [] yPoints = {y1, y2, (y1+y2)/2};
    		g2.fillPolygon(xPoints, yPoints, 3);
    	}
	}
	
	private int getX(int bp, int w) { 
		double f = (double)(bp-region.getStart()) / (double)region.getWidth();
		if(reverseIt){
			f=1.0-f;
		}
		return (int)Math.round(f * (double)w);
	}

}
