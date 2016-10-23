/*
 * Author: tdanford
 * Date: Sep 24, 2008
 */
package org.seqcode.viz.genomicplot;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Vector;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.QuadCurve2D;

import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.deepseq.StrandedPair;
import org.seqcode.deepseq.experiments.Sample;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.ScoredRegion;
import org.seqcode.genome.location.ScoredStrandedRegion;
import org.seqcode.gsebricks.verbs.*;
import org.seqcode.gseutils.Closeable;
import org.seqcode.gseutils.Pair;
import org.seqcode.viz.colors.Coloring;
import org.seqcode.viz.paintable.*;


public class ThinOverlapPaintable extends AbstractPaintable {
	
	private boolean axis =true;
	private boolean filledColumns=false;
	private Region region;
	private Vector<Region> highlighted;
	private List<Pair<Point,Point>> inters = new ArrayList<Pair<Point, Point>>();
	
	private Color bgColor, highlightColor, loopColor, interColor;
	private int bgThick, highlightThick;
	private boolean reverse;
	private int maxOverlap;
	private boolean drawPairedCurves=false;
	
	private Sample sample; 
	private int threePrimeExt;
	private RunningOverlapSum overlaps;
	private List<ScoredRegion> pairedRegs=null;
	private int[] dataProfile;  // y pixel values of where the data line is drawn.
	private boolean[][] bgPoints;
	private boolean[][] hlPoints;
	
	
	public ThinOverlapPaintable(Region basereg, Collection<Region> hls, Sample s, int threePrimeExt) { 
		this(basereg, hls, s, threePrimeExt, false);
	}
	public ThinOverlapPaintable(Region basereg,Collection<Region> hls, Sample s, int threePrimeExt, boolean drawPairs) { 
		this(basereg, hls, s, threePrimeExt, drawPairs, null);
		
	}
	public ThinOverlapPaintable(Region basereg, Collection<Region> hls, Sample s, int threePrimeExt, boolean drawPairs, List<Pair<Point,Point>> inters) { 
		region = basereg;
		this.inters = inters;
		sample = s;
		this.threePrimeExt = threePrimeExt;
		drawPairedCurves=drawPairs;
		if(hls != null && hls.size()>0)
			highlighted = new Vector<Region>(hls);
		else
			highlighted = new Vector<Region>();
		bgColor = Coloring.clearer(Color.gray);
		highlightColor = Coloring.clearer(Color.red);
		maxOverlap = 1;
		bgThick = 1; 
		highlightThick = 2;
		overlaps = new RunningOverlapSum(region.getGenome(), region.getChrom());
		reverse = false;
		
		loadData();
		
	}
	
	public void setBgColor(Color c){bgColor = c;}
	public void setLoopColor(Color c){loopColor = c;}
	public void setInterColor(Color c){interColor = c;}
	public void setHighlightColor(Color c){highlightColor = c;}
	public void setBgThick(int t){bgThick=t;}
	public void setHighlightThick(int t){highlightThick=t;}
	public void setFilledColumns(boolean b){filledColumns=b;}
	
	public void setReverse(boolean r) { 
		reverse = r;
		dispatchChangedEvent();
	}
	
	public void loadData() { 
		System.out.println("Loading data...");
		
		//Load hits in region, then convert to ScoredStrandedRegions
		List<StrandedBaseCount> hits = sample.getBases(region);
		float count = 0;
		for(StrandedBaseCount sbc : hits){ 
			ScoredStrandedRegion hit = sbc.expandToScoredStrandedRegion(region.getGenome(), region.getChrom(), 0, threePrimeExt);
			overlaps.addRegion(hit, (float)hit.getScore());
			count += hit.getScore();
		}
		System.out.println(String.format("\tLoaded %.2f hits.", count));
		maxOverlap = overlaps.getMaxOverlap();
		System.out.println("Max: "+maxOverlap);
		
		//Load pairs, convert to ScoredRegions
		if(drawPairedCurves){
			List<StrandedPair> pairs = sample.getPairs(region);
			if(pairs!=null){
				pairedRegs = new ArrayList<ScoredRegion>();
				for(StrandedPair sp : pairs){
					ScoredRegion pair = sp.toContiguousRegion();
					if(pair !=null)
						pairedRegs.add(pair);
				}
			}
		}
	}
	
	
	public void setMaxOverlap(int m){maxOverlap=m;}
	public void synchronizeMaxOverlaps(ThinOverlapPaintable p) { 
		int maxo = Math.max(maxOverlap, p.maxOverlap);
		maxOverlap = maxo;
		p.maxOverlap = maxo;
		dispatchChangedEvent();
		p.dispatchChangedEvent();
	}

	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		Graphics2D g2 = (Graphics2D)g;
		g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, 
				RenderingHints.VALUE_ANTIALIAS_ON);
		
		int w = x2-x1;
		dataProfile = new int[w+2];
		g2.setColor(Color.white);
		g2.fillRect(x1, y1, w, y2-y1);
		int[][] changes = overlaps.getChangePoints();
		
		//Split the area into two sections: reads & pairs. The latter is usually zero area.
		int readh = drawPairedCurves ? (int)((y2-y1)*0.7) : (y2-y1);
		int ready2 = y1+readh;
		int pairh = (y2-y1)-readh;
		int pairy1 = ready2;

		bgPoints = new boolean [w+1][readh+1];
		hlPoints = new boolean [w+1][readh+1];
		for(int a=0; a<=w; a++)
			for(int b=0; b<=readh; b++){
				bgPoints[a][b]=false;
				hlPoints[a][b]=false;
			}
		int runningOverlap = 0;
		int px = -1, py = -1, pbp = -1;
		
		//Map out the points that must be drawn
		for(int i = 0; i < changes.length; i++) {
			int bp = changes[i][0];
			runningOverlap += changes[i][1];
			int x = Math.min(Math.max(getX(bp, w), 0), w);	
			double yf = Math.min((double)runningOverlap / (double)maxOverlap, 1.0);
			int y = (int)Math.round(yf * (double)readh);
			
			if(px == -1) {
				if(reverse)
					px=w;
				else
					px = 0; 
				py = 0; pbp = bp;
			}
						
			fillPoints(px, py, x, py, isHighlighted(pbp));
			fillPoints(x, py, x, y, isHighlighted(bp));
						
			for(int k = px; k < x; k++) 
				setDataProfile(k, py);
			setDataProfile(x, y);
			
			if(i == changes.length-1) {	// last line 
				fillPoints(x, y, w, y, isHighlighted(bp));
				for(int k = x; k < x2; k++)
					setDataProfile(k, y);
			}
			px = x; py = y; pbp = bp;
		}

		//Fill in the circles...
		for(int a=0; a<=w; a++){
			int diam, rad;
			int x=a+x1;
			
			if(filledColumns){
				//Columns
				int maxY=-1;
				for(int b=0; b<=readh; b++){
					if(hlPoints[a][b] || bgPoints[a][b])
						maxY=b;
				}
				if(maxY!=-1){
					int y = ready2-maxY;
					if(hlPoints[a][maxY]){
						diam = highlightThick;
						rad = Math.max(1, diam/2);
						g2.setColor(highlightColor);
						g2.fillRect(x-rad, y-rad, diam, maxY);
					}else if(bgPoints[a][maxY]){
						diam = bgThick;
						rad = Math.max(1, diam/2);
						g2.setColor(bgColor);
						g2.fillRect(x-rad, y-rad, diam, maxY);
					}
				}
			}else{
				//Dots
				for(int b=0; b<readh; b++){
					int y = ready2-b;
					if(hlPoints[a][b]){
						diam = highlightThick;
						rad = Math.max(1, diam/2);
						g2.setColor(highlightColor);
						g2.fillOval(x-rad, y-rad, diam, diam);
					}else if(bgPoints[a][b]){
						diam = bgThick;
						rad = Math.max(1, diam/2);
						g2.setColor(bgColor);
						g2.fillOval(x-rad, y-rad, diam, diam);
					}
				}
			}
		}
		
		//Paired-end curves
		if(drawPairedCurves){
			for(ScoredRegion pair : pairedRegs){
				if(region.contains(pair)){
					int xA = x1+Math.min(Math.max(getX(pair.getStart(), w), 0), w);
					int xB = x1+Math.min(Math.max(getX(pair.getEnd(), w), 0), w);
					int xMid = (xA+xB)/2;
					int yMid = (int)(((double)pair.getWidth()/(double)region.getWidth()) * pairh*2) + pairy1;
					
					g2.setColor(loopColor);
		    		g2.setStroke(new BasicStroke(1.0f));
		    		QuadCurve2D loop = new QuadCurve2D.Float(xA, pairy1, xMid, yMid, xB, pairy1);
		    		g2.draw(loop);
				}
			}
			for(Pair<Point,Point> inter : inters){
				if(region.contains(inter.car()) && region.contains(inter.cdr())){
					int xA = x1+Math.min(Math.max(getX(inter.car().getLocation(), w), 0), w);
					int xB = x1+Math.min(Math.max(getX(inter.cdr().getLocation(), w), 0), w);
					int xMid = (xA+xB)/2;
					int yMid = (int)(((double)(inter.cdr().getLocation()-inter.car().getLocation())/(double)region.getWidth()) * pairh*2) + pairy1;
					
					g2.setColor(interColor);
		    		g2.setStroke(new BasicStroke(1.0f));
		    		QuadCurve2D loop = new QuadCurve2D.Float(xA, pairy1, xMid, yMid, xB, pairy1);
		    		g2.draw(loop);
				}
			}
		}
		
		if(axis){
			g2.setColor(Color.black);
    		g2.setStroke(new BasicStroke(1.0f));
    		g2.drawLine(x1, y1, x1,ready2);
    		g2.setFont(new Font("Ariel", Font.PLAIN, 16));
    		FontMetrics metrics = g2.getFontMetrics();
    		g2.drawString(String.format("%d",maxOverlap), x1+1, y1+(metrics.getHeight()/2));
		}
	}
	
	public int[] getDataProfile() { return dataProfile; } 
	
	private void setDataProfile(int x, int y) { 
		if(x >= 0 && x < dataProfile.length) {
			dataProfile[x] = Math.max(dataProfile[x], y);
		} else { 
			//System.err.println(String.format("%d not in %d profile", x, dataProfile.length));
		}
	}
	
	private boolean isHighlighted(int bp) { 
		for(Region r : highlighted) { 
			if(r.getStart() <= bp && r.getEnd() >= bp) { 
				return true;
			}
		}
		return false;
	}
	
	private int getX(int bp, int pixwidth) {
		double f = (double)(bp-region.getStart()) / (double)region.getWidth();
		if(reverse) { 
			f = 1.0 - f;
		}
		return (int)Math.round(f * (double)pixwidth);
	}
	
	private void fillPoints(int p1x, int p1y, int p2x, int p2y, boolean highlighted){
		int dx = Math.abs(p2x-p1x), dy = Math.abs(p2y-p1y);
		int divs = Math.max(dx, dy);
		
		for(int i = 0; i <= divs; i++) { 
			double f = (double)i / (double)divs;
			double nf = 1.0-f;
			int x = (int)Math.round(f*(double)p1x + nf*(double)p2x);
			int y = (int)Math.round(f*(double)p1y + nf*(double)p2y);
			//System.out.println(x+"\t"+y);
			
			if(highlighted)
				hlPoints[x][y]=true;
			else
				bgPoints[x][y]=true;
		}
	}
	
	private void strokeLine(Graphics2D g2, 
				java.awt.Point p1, java.awt.Point p2, 
				int thick1, int thick2, 
				Color c1, Color c2) { 
		
		int dx = Math.abs(p2.x-p1.x), dy = Math.abs(p2.y-p1.y);
		int divs = Math.max(dx, dy);
		int[] c1array = Coloring.rgba(c1);
		int[] c2array = Coloring.rgba(c2);
		int[] carray = new int[4];
		
		for(int i = 0; i <= divs; i++) { 
			double f = (double)i / (double)divs;
			double nf = 1.0-f;
			int x = (int)Math.round(f*(double)p1.x + nf*(double)p2.x);
			int y = (int)Math.round(f*(double)p1.y + nf*(double)p2.y);
			int diam = (int)Math.round(f*(double)thick1 + nf*(double)thick2);
			int rad = Math.max(1, diam/2);
			
			for(int j =0; j < carray.length; j++) { 
				carray[j] = (int)Math.round(f*(double)c1array[j] + nf*(double)c2array[j]);
			}
			Color c = Coloring.asColor(carray);
			g2.setColor(c);
			g2.fillOval(x-rad, y-rad, diam, diam);
		}
	}
	
	/**
	 * A one-shot (i.e., one use) expander, that automatically closes its inner expander
	 * when it's been called the first time.  Only used in the constructor, above.
	 * 
	 * @author tdanford
	 *
	 * @param <R>
	 */
	private static class ClosingRegionExpanderWrapper<R extends Region> 
		implements Expander<Region,Region> {  
		
		private Expander<Region,R> exp;
		private Mapper<R,R> mapper;
		
		public ClosingRegionExpanderWrapper(Expander<Region,R> e) { 
			exp = e;
			mapper = null;
		}

		public ClosingRegionExpanderWrapper(Expander<Region,R> e, Mapper<R,R> m) { 
			exp = e;
			mapper = m;
		}

		public Iterator<Region> execute(Region a) {
			Iterator<R> itr = exp.execute(a);
			if(mapper != null) { 
				itr = new MapperIterator<R,R>(mapper,itr);
			}
			Iterator<Region> regs = new MapperIterator<R,Region>(new CastingMapper<R,Region>(), itr);
			LinkedList<Region> regList = new LinkedList<Region>();
			while(regs.hasNext()) { regList.addLast(regs.next()); }
			
			if(exp instanceof Closeable) { 
				((Closeable)exp).close();
			}
			System.out.println(String.format("ClosingRegionExpanderWrapper: %d", regList.size()));
			
			return regList.iterator();
		}
	}
}


