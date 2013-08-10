package edu.psu.compbio.seqcode.gse.seqview.paintable;

import java.io.File;
import java.awt.*;
import java.awt.geom.QuadCurve2D;
import java.util.*;
import java.util.List;

import edu.psu.compbio.seqcode.gse.projects.readdb.PairedHit;
import edu.psu.compbio.seqcode.gse.seqview.model.InteractionArcModel;
import edu.psu.compbio.seqcode.gse.seqview.model.SeqHistogramModel;
import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.utils.probability.NormalDistribution;
import edu.psu.compbio.seqcode.gse.utils.stats.StatUtil;
import edu.psu.compbio.seqcode.gse.viz.DynamicAttribute;


public class SeqHistogramPainter extends RegionPaintable {

	private SeqHistogramModel histomodel;
	private InteractionArcModel arcmodel;
	private DynamicAttribute attrib;
	protected static List configurationFields = null;
	private SeqHistogramProperties props;
	private double[] gaussian;	//Gaussian kernel for density estimation
	private int kernelVar = 0;
	private boolean alignPaired = false;

	public SeqHistogramPainter(SeqHistogramModel hmodel, InteractionArcModel pmodel) {
		super();
		histomodel = hmodel;
		arcmodel = pmodel;
		props = new SeqHistogramProperties();
		histomodel.addEventListener(this);
		if(arcmodel!=null)
			arcmodel.addEventListener(this);
		attrib = DynamicAttribute.getGlobalAttributes();
		alignPaired = arcmodel!=null;
		if(alignPaired){
			props.Stranded=false;
			props.DrawPairedCurves=true;
			props.MaxReadCount = 100;
		}
	}
	public SeqHistogramProperties getProperties() {return props;}
	public void setProperties(SeqHistogramProperties p) {props = p;}
	public void savePropsInDir(File dir) {
		super.savePropsInDir(dir);
		saveModelPropsInDir(dir,histomodel);
	}
	public void loadPropsInDir(File dir) {
		super.loadPropsInDir(dir);
		loadModelPropsInDir(dir,histomodel);
	}    
	public void cleanup() { 
		super.cleanup();
		histomodel.removeEventListener(this);
		if(arcmodel!=null){
			arcmodel.removeEventListener(this);
		}
	}
	public boolean canPaint() {
		return histomodel.isReady() && (arcmodel==null || arcmodel.isReady());
	}
	public synchronized void eventRegistered(EventObject e) {        
		if (e.getSource() == histomodel || (arcmodel!=null && e.getSource() == arcmodel )) {
			if(histomodel.isReady() && (arcmodel==null || arcmodel.isReady())){
				setCanPaint(true);
				setWantsPaint(true);
				notifyListeners();
			}
		}
	}
	//pre-calculate and store the Guassian kernel prob., for efficiency
	private void initGaussianKernel(int var){
		kernelVar = var;
		gaussian = new double[250]; 
		NormalDistribution gaussianDist = new NormalDistribution(0, var*var);
		for (int i=0;i<gaussian.length;i++)
			gaussian[i]=gaussianDist.calcProbability((double)i);
	}
	// convert a weight histogram to a gaussian kernel density profile
	private Map<Integer,Float> convertKernelDensity(Map<Integer,Float> data){
		Map<Integer,Float> results = new TreeMap<Integer,Float>();
		int min= Integer.MAX_VALUE;
		int max= Integer.MIN_VALUE;
		for (int pos: data.keySet()){
			if (min>pos)
				min = pos;
			if (max<pos)
				max= pos;
		}
		// get all reads, convert to basepair-resolution density
		double[] profile = new double[max-min+1+100];	// add 50bp padding to the ends
		for (int pos: data.keySet()){
			profile[pos-min+50]=(double)data.get(pos);
		}

		double[] densities = StatUtil.symmetricKernelSmoother(profile, gaussian);
		// set density values back to the positions (only at certain resolution)
		// so that we can paint efficiently
		int step = 1;
		if (max-min>1024)
			step = (max-min)/1024;
		for (int i=min-50; i<=max+50; i+=step){
			results.put(i, (float)densities[i-min+50]);
		}
		return results;
	}
	public void removeEventListener(Listener<EventObject> l) {
		super.removeEventListener(l);
		if (!hasListeners()) {
			histomodel.removeEventListener(this);
			if(arcmodel!=null){
				arcmodel.removeEventListener(this);System.out.println("ListenerRemoved");
			}
		}
	}
	public void paintItem(Graphics2D g, 
			int x1, int y1, 
			int x2, int y2) {

		g.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_GASP);
		g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
		
		if (!canPaint()) {
			return;
		}
		boolean stranded = getProperties().Stranded;
		int trackWidth = x2 - x1;
		int trackHeight = y2-y1;
		boolean autoUpdateBins = getProperties().BinAutoUpdate;
		boolean logscale =false;

		Map<Integer,Float> plus = histomodel.getPlus(), minus = histomodel.getMinus();
		
		double maxobs = 0;
		for (int i : plus.keySet()) {
			if (plus.get(i) > maxobs) {
				maxobs = plus.get(i);
			}
		}
		for (int i : minus.keySet()) {
			if (minus.get(i) > maxobs) {
				maxobs = minus.get(i);
			}
		}
		int maxhits = props.MaxReadCount;
		if (maxhits < 0) {maxhits = (int)maxobs;}
		int regionStart = histomodel.getRegion().getStart();
		int regionEnd = histomodel.getRegion().getEnd();
		
		//Set y coordinates for drawing
		int readh = (props.DrawPairedCurves && alignPaired) ? (int)((y2-y1)*0.7) : (y2-y1);
		int ready2 = y1+readh;
		int pairh = (y2-y1)-readh;
		int pairy1 = ready2;
		int midpoint = (y1 + ready2) / 2;

		Stroke oldStroke = g.getStroke();
		int linewidth = getProperties().getLineWidth();
		if (linewidth < 0) {
			linewidth = trackWidth / (histomodel.getRegion().getWidth() / histomodel.getProperties().BinWidth);
		}
		if (linewidth < 1) {
			linewidth = 1;
		}
		int actualBinWidth = histomodel.getProperties().BinWidth;
		if(autoUpdateBins){
			if (trackWidth / linewidth < histomodel.getRegion().getWidth() / histomodel.getProperties().BinWidth) {
				actualBinWidth = histomodel.getRegion().getWidth() / (trackWidth / linewidth);
				combineBins(plus, actualBinWidth);
				combineBins(minus, actualBinWidth);
			}
		}
		int binPixels = trackWidth/(histomodel.getRegion().getWidth() / histomodel.getProperties().BinWidth);
		if(binPixels<1)
			binPixels=1;

		g.setStroke(new BasicStroke((float)linewidth));

		// Draw Gaussian density probabilities (if viewing <10Kbp)
		int gk_var = histomodel.getProperties().GaussianKernelVariance; 
		if (gk_var!=0 && regionEnd-regionStart<=getProperties().getDrawGaussianMaxWindow()){
			if (gaussian==null){
				this.initGaussianKernel(gk_var);
			}
			else if (gk_var!=kernelVar){
				this.initGaussianKernel(gk_var);
			}

			Map<Integer, Float> density_p = convertKernelDensity(plus);
			Map<Integer, Float> density_m = convertKernelDensity(minus);
			// find max prob
			float max= Integer.MIN_VALUE;
			for (float prob: density_p.values())
				if (max<prob)
					max= prob;
			for (float prob: density_m.values())
				if (max<prob)
					max= prob;
			double scaling = (maxobs/max)*2;

			if (stranded) {
				g.setColor(getProperties().getMinusColor());
				int x_prev = -1;
				int y_prev = -1;
				for (int pos : density_p.keySet()) {
					double val = density_p.get(pos)*scaling;
					int xpix = getXPos(pos, regionStart, regionEnd, x1, x2);
					int ypix = getYPos(val, 0, maxhits, y1, midpoint, logscale);
					if (x_prev!=-1)
						g.drawLine(x_prev, y_prev, xpix, ypix);
					x_prev = xpix;
					y_prev = ypix;
				}
				g.setColor(getProperties().getPlusColor());
				x_prev = -1;
				y_prev = -1;
				for (int pos : density_m.keySet()) {
					double val = density_m.get(pos)*scaling;
					int xpix = getXPos(pos, regionStart, regionEnd, x1, x2);
					int ypix = midpoint + (ready2 - getYPos(val, 0, maxhits, midpoint, ready2, logscale));
					if (x_prev!=-1)
						g.drawLine(x_prev, y_prev, xpix, ypix);
					x_prev = xpix;
					y_prev = ypix;                    
				}
			}  	
		}
		//Draw the read density
		if (stranded) {
			HashMap<Integer,Double> plotXVals = new HashMap<Integer,Double>();
			//Plus strand : first screen out overlapping rects (take max), then plot
			g.setColor(getProperties().getMinusColor());
			for (int pos : plus.keySet()) {
				double val = plus.get(pos);
				int xpix = getXPos(pos, regionStart, regionEnd, x1, x2);
				if(!plotXVals.containsKey(xpix) || plotXVals.get(xpix)<val){
					plotXVals.put(xpix,val);
				}
			}
			for(int xpix : plotXVals.keySet()){
				double val = plotXVals.get(xpix); 
				int ypix = getYPos(val, 0, maxhits, y1, midpoint, logscale);
				g.fillRect(xpix, ypix, binPixels, midpoint-ypix);
			}
			plotXVals.clear();
			//Minus strand
			g.setColor(getProperties().getPlusColor());
			for (int pos : minus.keySet()) {
				double val = minus.get(pos);
				int xpix = getXPos(pos, regionStart, regionEnd, x1, x2);
				if(!plotXVals.containsKey(xpix) || plotXVals.get(xpix)<val){
					plotXVals.put(xpix,val);
				}
			}for(int xpix : plotXVals.keySet()){
				double val = plotXVals.get(xpix); 
				int ypix = midpoint + (ready2 - getYPos(val, 0, maxhits, midpoint, ready2, logscale));
				g.fillRect(xpix, midpoint, binPixels, ypix-midpoint);
			}
			//Line & trimmings
			g.setColor(Color.darkGray);
			g.drawLine(x1, midpoint, x2, midpoint);
			g.setFont(attrib.getPointLabelFont(trackWidth,trackHeight));
			int fh = g.getFont().getSize();
			int step = Math.max(1,(int)Math.round(maxhits / 1)); //y-axis scale
			for (int i = step; i <= Math.ceil(maxhits); i += step) {
				int ypos = getYPos(i, 0, maxhits, y1+fh, midpoint, logscale);
				g.drawString(Integer.toString(i),5,ypos);
				ypos = midpoint + (ready2 - getYPos(i, 0, maxhits, midpoint, ready2, logscale));
				g.drawString(Integer.toString(i),5,ypos);
			}
			g.setColor(Color.black);
		} else {
			//Plot density : first screen out overlapping rects (take max), then plot
			HashMap<Integer,Double> plotXVals = new HashMap<Integer,Double>();
			g.setColor(getProperties().getUnstrandedColor());
			for (int pos : plus.keySet()) {
				if (true/*!model.getProperties().ShowInteractionHistogram || !pvals.containsKey(pos)*/) {
					double val = plus.get(pos);
					if (minus.containsKey(pos)) {
						val += minus.get(pos);
					}
					int xpix = getXPos(pos, regionStart, regionEnd, x1, x2);
					if(!plotXVals.containsKey(xpix) || plotXVals.get(xpix)<val){
						plotXVals.put(xpix,val);
					}
				}
			}
			for (int pos : minus.keySet()) {
				if (plus.containsKey(pos)) {
					continue;
				}
				double val = minus.get(pos);
				int xpix = getXPos(pos, regionStart, regionEnd, x1, x2);
				if(!plotXVals.containsKey(xpix) || plotXVals.get(xpix)<val){
					plotXVals.put(xpix,val);
				}
			}
			for(int xpix : plotXVals.keySet()){
				double val = plotXVals.get(xpix); 
				int ypix = getYPos(val, 0, maxhits, y1, ready2, logscale);
				g.fillRect(xpix, ypix, binPixels, ready2-ypix);
			}
			
			for(int xpix : plotXVals.keySet()){
				double val = plotXVals.get(xpix); 
				int ypix = getYPos(val, 0, maxhits, y1, ready2, logscale);
				g.fillRect(xpix, ypix, binPixels, ready2-ypix);
			}

			//Line & trimmings
			g.setColor(Color.darkGray);
			g.drawLine(x1, ready2, x2, ready2);
			g.setFont(attrib.getPointLabelFont(trackWidth,trackHeight));
			int fh = g.getFont().getSize();
			int step = Math.max(1,(int)Math.round(maxhits / 1)); //y-axis scale labels
			for (int i = step; i <= Math.ceil(maxhits); i += step) {
				int ypos = getYPos(i, 0, maxhits, y1+fh, ready2, logscale);
				g.drawString(Integer.toString(i),5,ypos);
			}
			g.setColor(Color.black);
		}

		//Draw pairing arcs
		if(props.DrawPairedCurves && alignPaired){
			List<PairedHit> pairs = arcmodel.getResults();
			for(PairedHit pair : pairs){
				if(pair.leftChrom==pair.rightChrom && 
						pair.leftPos>=arcmodel.getRegion().getStart() && pair.leftPos<=arcmodel.getRegion().getEnd() &&
						pair.rightPos>=arcmodel.getRegion().getStart() && pair.rightPos<=arcmodel.getRegion().getEnd() ){
					
					int lPos = pair.lesserPos();
					int rPos = pair.greaterPos();
					if(pair.pairCode==1){
						g.setColor(getProperties().getMateArcColor());
					}else{
						//Junction-mapped reads are a different color,
						//and it looks better if they go from the end of one block to the start of the other 
						g.setColor(getProperties().getSplitReadArcColor());
						
						if(pair.lesserStrand())
							lPos = pair.lesserPos()+pair.lesserLength()-1;
						if(!pair.greaterStrand())
							rPos = pair.greaterPos()-pair.greaterLength()+1;
					}
		    		int xA =getXPos(lPos, regionStart, regionEnd, x1, x2);
					int xB =getXPos(rPos, regionStart, regionEnd, x1, x2);
					double pwidth = rPos-lPos;
					int xMid = (xA+xB)/2;
					int yMid = (int)(((double)pwidth/(double)(regionEnd-regionStart)) * pairh*2) + pairy1;
					
					g.setStroke(new BasicStroke(1.0f));
		    		QuadCurve2D loop = new QuadCurve2D.Float(xA, pairy1, xMid, yMid, xB, pairy1);
		    		g.draw(loop);
				}
			}
		}
		
		
		//Draw labels
		g.setStroke(oldStroke);
		if (getProperties().DrawTrackLabel) {
			g.setFont(attrib.getRegionLabelFont(trackWidth,trackHeight));
			g.setColor(Color.black);
			g.drawString(getLabel(),x1 + g.getFont().getSize()*2,y1 + g.getFont().getSize());
		}if(getProperties().DrawBinSize) {
			g.setFont(attrib.getPointLabelFont(trackWidth,trackHeight));
			g.setColor(getProperties().getUnstrandedColor());
			String binString = new String("(bin size: "+actualBinWidth+"bp)");
			g.drawString(binString,x2 - g.getFontMetrics().stringWidth(binString) - g.getFont().getSize()*2,y1 + g.getFont().getSize());
		}
	}

	/**
	 * If the binWidth was too small then there aren't enough pixels to show
	 * each bin and so they'll be drawn overlapping.  This method
	 * resizes the bins into bins of width binWidth
	 * by combining skinny bins that map to the same fat bin
	 */
	private void combineBins(Map<Integer,Float> map, int binWidth) {
		if (map == null) {
			throw new NullPointerException("null map");
		}
		java.util.List<Integer> keys = new ArrayList<Integer>();
		keys.addAll(map.keySet());
		for (int p : keys) {
			if (p % binWidth == 0) {
				continue;
			} else {
				int newbin = (p / binWidth) * binWidth;
				if (map.containsKey(newbin)) {
					map.put(newbin, map.get(newbin) + map.get(p));
				} else {
					map.put(newbin, map.get(p));
				}
				map.remove(p);
			}
		}
	}

}