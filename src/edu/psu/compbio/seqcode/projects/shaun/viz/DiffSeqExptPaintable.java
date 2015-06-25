/*
 * Author: tdanford
 * Date: Sep 24, 2008
 */
package edu.psu.compbio.seqcode.projects.shaun.viz;

import java.util.*;
import java.util.List;
import java.awt.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.*;
import edu.psu.compbio.seqcode.gse.projects.gps.DeepSeqExpt;
import edu.psu.compbio.seqcode.gse.projects.gps.ReadHit;
import edu.psu.compbio.seqcode.gse.utils.Closeable;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.viz.paintable.*;

public class DiffSeqExptPaintable extends AbstractPaintable {
	
	public static void main(String[] args) { 
		try {
			Genome g = Genome.findGenome("mm8");

			Region r = new Region(g, "6", 52000000, 52200000);
			
			List<SeqLocator> expts = new ArrayList<SeqLocator>();
			List<SeqLocator> bases = new ArrayList<SeqLocator>();
			Set<String> reps = new TreeSet<String>();
			expts.add(new SeqLocator("PPG Day2+8h RAR HBG3",reps,  "bowtie_unique"));
			bases.add(new SeqLocator("PPG Day2 RAR HBG3", reps, "bowtie_unique"));
            
			DiffSeqExptPaintable p = new DiffSeqExptPaintable(g, r, expts, bases, 2000, 1000, 1);

			PaintableFrame pf = new PaintableFrame("Test", p);
			
		} catch (NotFoundException e) {
			e.printStackTrace();
		}
	}

	private Genome gen;
	private boolean axis =true;
	private boolean filledColumns=true;
	private Region region;
	private Vector<Region> highlighted;
	private int binWidth, binStep;
	private double scalingFactor;
	private Color negColor, posColor;
	private boolean reverse=false;
	protected DeepSeqExpt signal;
	protected DeepSeqExpt control;
	private int maxLogFold=8;
	private int minLogFold=-8;
	private int sigPerBaseMax=2;
	private int ctrlPerBaseMax=2;
	private boolean needleFiltering=true;
	private double LOG2 = Math.log(2);
	private double[] diff;
	private int[] dataProfile;  // y pixel values of where the data line is drawn.

	//Load data from ReadDB sources
	public DiffSeqExptPaintable(Genome gen, Region basereg, List<SeqLocator> exptLocs, List<SeqLocator> baseLocs, int win, int step, double scaling) { 
		this.gen = gen;
		region = basereg;
		signal = new DeepSeqExpt(gen, exptLocs, "readdb", 1);
		control = new DeepSeqExpt(gen, baseLocs, "readdb", 1);
		binWidth = win;
		binStep = step;
		scalingFactor = scaling;
		negColor = Color.blue;
		posColor = Color.yellow;
		reverse = false;
		
		loadData();
	}
	//Load data from preformatted file
	public DiffSeqExptPaintable(Genome gen, Region basereg, File data, int win, int step, double scaling) { 
		this.gen = gen;
		region = basereg;
		signal = null;
		control = null;
		binWidth = win;
		binStep = step;
		scalingFactor = scaling;
		negColor = Color.blue;
		posColor = Color.yellow;
		reverse = false;
		loadDataFromFile(data);
	}
	
	public void setNegColor(Color c){negColor = c;}
	public void setPosColor(Color c){posColor = c;}
	public void setSigPerBaseMax(int pb){sigPerBaseMax = pb;}
	public void setCtrlPerBaseMax(int pb){ctrlPerBaseMax = pb;}
	public void setReverse(boolean r) { 
		reverse = r;
		dispatchChangedEvent();
	}
	public void setMaxLogFold(int m){maxLogFold=m;}
	public void setMinLogFold(int m){minLogFold=m;}
	public void synchronizeMaxLogFold(DiffSeqExptPaintable p) { 
		int maxo = Math.max(maxLogFold, p.maxLogFold);
		maxLogFold = maxo;
		p.maxLogFold = maxo;
		int mino = Math.max(minLogFold, p.minLogFold);
		minLogFold = mino;
		p.minLogFold = mino;
		dispatchChangedEvent();
		p.dispatchChangedEvent();
	}

	public void loadData() { 
		System.err.println("Loading data in "+region);
		List<ReadHit> AHits = signal.loadHits(region);
		List<ReadHit> BHits = control.loadHits(region);
		
		double [] ABinnedReads = makeHitLandscape(AHits, region, sigPerBaseMax, '.');
		double [] BBinnedReads = makeHitLandscape(BHits, region, ctrlPerBaseMax, '.');
		int numBins = (int)(region.getWidth()/binStep);
		diff = new double[numBins+1];
		for(int x=0; x<=numBins; x++){
			if(ABinnedReads[x]==0 && BBinnedReads[x]==0)
				diff[x] = 0;
			else if(BBinnedReads[x]==0)
				diff[x] = Math.min(Math.log(ABinnedReads[x]), Math.log(maxLogFold)/LOG2);
			else if(ABinnedReads[x]==0)
				diff[x] = -1* Math.min(Math.log(BBinnedReads[x]*scalingFactor), -Math.log(minLogFold)/LOG2);
			else
				diff[x] = Math.log(ABinnedReads[x]/(BBinnedReads[x]*scalingFactor))/LOG2;
		}
	}
	
	public void loadDataFromFile(File data) { 
		//Load data
		HashMap<String, Double> dataHash = new HashMap<String,Double>();
		try {			
			if(!data.isFile()){System.err.println("Invalid file name: "+data.getAbsolutePath());System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(data));
			String line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            if(!words[0].equals("Position"))
	            	dataHash.put(words[0], new Double(words[1]));
	        }
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
        
		//format
		int numBins = (int)(region.getWidth()/binStep);
		diff = new double[numBins+1];
		String chr = region.getChrom();
		int start = region.getStart();
		for(int x=0; x<=numBins; x++){
			int curr =start +(x*binStep); 
			String currPos = chr+":"+curr;  
			if(dataHash.containsKey(currPos))
				diff[x] = dataHash.get(currPos);
			else
				diff[x] = 0;
		}
	}
	
	//Makes integer arrays corresponding to the read landscape over the current region
	private double [] makeHitLandscape(List<ReadHit> hits, Region currReg, int perBaseMax, char strand){
		int numBins = (int)(currReg.getWidth()/binStep);
		int [] counts = new int[currReg.getWidth()+1];
		double [] landscape = new double[numBins+1];
		for(int i=0; i<=numBins; i++){ landscape[i]=0;}
		for(int i=0; i<=currReg.getWidth(); i++){counts[i]=0;}
		int ext = binWidth/2;
		for(ReadHit r : hits){
			if(strand=='.' || r.getStrand()==strand){
				int rpos = r.getFivePrime(); 
				int offset=inBounds(rpos-currReg.getStart(),0,currReg.getWidth());
				counts[offset]++;//small issue here... counts will not be corrected for scalingFactor in control
				if(!needleFiltering ||  (counts[offset] <= perBaseMax)){
					int binstart = inBounds((int)((double)(offset-ext)/binStep), 0, numBins);
					int binend = inBounds((int)((double)(offset+ext)/binStep), 0, numBins);
					for(int i=binstart; i<=binend; i++){
						landscape[i]+=r.getWeight();
					}
				}
			}
		}
		return(landscape);
	}
	private final int inBounds(int x, int min, int max){
		if(x<min){return min;}
		if(x>max){return max;}
		return x;
	}

	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		Graphics2D g2 = (Graphics2D)g;
		g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, 
				RenderingHints.VALUE_ANTIALIAS_ON);
		
		int w = x2-x1;
		dataProfile = new int[w+2];
		g2.setColor(Color.white);
		g2.fillRect(x1, y1, w, y2-y1);
		
		//Split the area into two sections: positive & negative.
		int posDiffH = (y2-y1)/ ((maxLogFold-minLogFold)/maxLogFold);
		int posDiffY1 = y1;
		int posDiffY2 = y1+posDiffH;
		int negDiffH = (y2-y1)-posDiffH;
		int negDiffY1 = posDiffY2;
		int negDiffY2 = y2;

		//Map out the points that must be drawn
		int numBins = (int)(region.getWidth()/binStep);
		for(int a=0; a<=numBins; a++){
			int bp = region.getStart() + (a*binStep)+(binWidth/2);
			int currX1 = x1+Math.min(Math.max(getX(bp-(binStep/2), w), 0), w);
			int currX2 = x1+Math.min(Math.max(getX(bp+(binStep/2), w), 0), w);
			
			
			double foldDiff = diff[a];
			if(foldDiff > (double)maxLogFold)
				foldDiff = (double)maxLogFold;
			else if(foldDiff < (double)minLogFold)
				foldDiff = (double)minLogFold;
			double yf = foldDiff<0 ? Math.abs(foldDiff / -(double)minLogFold) : Math.abs(foldDiff / (double)maxLogFold);
			int currY1 =foldDiff<0 ? negDiffY1 : posDiffY2 - (int)Math.round(yf * (double)posDiffH);
			int currY2 = foldDiff<0 ?  negDiffY1 + (int)Math.round(yf * (double)negDiffH) : posDiffY2;
			
			Color currColor = foldDiff<0 ? negColor : posColor;
			
			//Columns
			g2.setColor(currColor);
			int x = Math.min(currX1,  currX2);
			int y = Math.min(currY1,  currY2);
			int width = Math.abs(currX2-currX1);
			int height = Math.abs(currY2-currY1);
			g2.fillRect(x, y, width, height);
		}
		
		
		if(axis){
			g2.setColor(Color.black);
    		g2.setStroke(new BasicStroke(1.0f));
    		g2.drawLine(x1, posDiffY1, x1,negDiffY2);
    		g2.setFont(new Font("Ariel", Font.PLAIN, 16));
    		FontMetrics metrics = g2.getFontMetrics();
    		g2.drawString(String.format("%d",maxLogFold), x1+1, y1+(metrics.getHeight()/2));
    		g2.drawString(String.format("%d",minLogFold), x1+1, y2+(metrics.getHeight()/2));
		}
	}
	
	public int[] getDataProfile() { return dataProfile; } 
	
	private void setDataProfile(int x, int y) { 
		if(x >= 0 && x < dataProfile.length) {
			dataProfile[x] = Math.max(dataProfile[x], y);
		} else { 
			System.err.println(String.format("%d not in %d profile", x, dataProfile.length));
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


