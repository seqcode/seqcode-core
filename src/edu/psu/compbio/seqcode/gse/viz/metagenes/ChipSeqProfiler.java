package edu.psu.compbio.seqcode.gse.viz.metagenes;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.LinkedList;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedPoint;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqExptHandler;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqHit;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.verbs.chipseq.SeqExpander;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.RealValuedHistogram;

public class ChipSeqProfiler implements PointProfiler<Point, Profile>{

	private Genome gen;
	private BinningParameters params=null;
	private int readLen=26, readExt=174;
	private double totalHits=0, backTotalHits=0;
	private ArrayList<SeqLocator> locs;
	private ArrayList<SeqLocator> backlocs;
	private ArrayList<SeqExptHandler> handles;
	private ArrayList<SeqExptHandler> ctrl_handles;
	private boolean zScoring=false;

	public ChipSeqProfiler(BinningParameters bp, Genome g, ArrayList<SeqLocator> experiments, ArrayList<SeqLocator> controls, int readLength, int readExtension, boolean z){
		this(bp, g, experiments, controls, readLength, readExtension);
		zScoring=z;
	}
	public ChipSeqProfiler(BinningParameters bp, Genome g, ArrayList<SeqLocator> experiments, ArrayList<SeqLocator> controls){
		this(bp, g, experiments, controls, 26, 174);
	}
	public ChipSeqProfiler(BinningParameters bp, Genome g, ArrayList<SeqLocator> experiments, ArrayList<SeqLocator> controls, int readLength, int readExtension){
		gen=g;
		params=bp; 
		readLen =readLength;
		readExt=readExtension;
		totalHits=0; backTotalHits=0;
		locs = experiments;
		backlocs = controls;
		
		handles = new ArrayList<SeqExptHandler>();
		ctrl_handles = new ArrayList<SeqExptHandler>();
		try {
			//Load experiments
            for(SeqLocator l : locs){
				System.err.print(String.format("%s\t", l.getExptName()));
				SeqExptHandler curr = new SeqExptHandler(gen, l);
				curr.setReadLength(readLen);
                curr.setReadExtension(readExt);
				handles.add(curr);
				totalHits += curr.getHitCount();
			}System.err.print(String.format("%.0f reads loaded\n", totalHits));
			//Load controls
			if(backlocs!=null){
	            for(SeqLocator l : backlocs){
					System.err.print(String.format("%s\t", l.getExptName()));
					SeqExptHandler curr = new SeqExptHandler(gen, l);
					curr.setReadLength(readLen);
	                curr.setReadExtension(readExt);
					ctrl_handles.add(curr);
					backTotalHits += curr.getHitCount();
				}System.err.print(String.format("%.0f reads loaded\n", backTotalHits));
			}
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public BinningParameters getBinningParameters() {
		return params;
	}

	public Profile execute(Point a) {
		double[] array = new double[params.getNumBins()];
		for(int i = 0; i < array.length; i++) { array[i] = 0; }
		
		int window = params.getWindowSize();
		int left = window/2;
		int right = window-left-1;
		
		boolean strand = (a instanceof StrandedPoint) ? 
				((StrandedPoint)a).getStrand() == '+' : true;
		
		int start = Math.max(0, a.getLocation()-left);
		int end = Math.min(a.getLocation()+right, a.getGenome().getChromLength(a.getChrom())-1);
		Region query = new Region(gen, a.getChrom(), start, end);
		Region extQuery = new Region(a.getGenome(), a.getChrom(), start-(readLen+readExt)>0 ? start-(readLen+readExt) : 1, end+(readLen+readExt) < a.getGenome().getChromLength(a.getChrom()) ? end+(readLen+readExt) : a.getGenome().getChromLength(a.getChrom()) );
		
		LinkedList<StrandedRegion> hits = new LinkedList<StrandedRegion>();
		LinkedList<StrandedRegion> backhits = new LinkedList<StrandedRegion>();
		for(SeqExptHandler e: handles){
			hits.addAll(e.loadExtendedHits(extQuery));
		}
		if(backTotalHits>0){
			for(SeqExptHandler e: ctrl_handles){
				backhits.addAll(e.loadExtendedHits(query));
			}
		}
		
		if(zScoring){
			int winWidth = params.getBinSize();
			int [] ipHitCounts = makeHitLandscape(hits, query, winWidth);
            int [] backHitCounts = null;
            if (backTotalHits>0) {
                backHitCounts = makeHitLandscape(backhits, query, winWidth);
            }
			
			for(int i=query.getStart(); i<query.getEnd()-winWidth; i+=winWidth){
				Region currWin = new Region(gen, query.getChrom(), i, (i+winWidth-1));
				
				int binid = (int)Math.max(0, ((double)(currWin.getStart()-query.getStart())/winWidth));
				double ipWinHits=(double)ipHitCounts[binid];
				double backWinHits= backTotalHits>0 ? ((double)backHitCounts[binid]) : 0;
				
				if(ipWinHits>0){
					double Z0 = binomialSampleEquality(ipWinHits, backWinHits, totalHits, backTotalHits);
					int startOffset = Math.max(0, currWin.getStart()-query.getStart());
					int endOffset = Math.max(0, Math.min(query.getEnd(), currWin.getEnd()-query.getStart()));
					
					int centerOffset = (startOffset+endOffset)/2;
					
					if(!strand) { 
						int tmp = window-centerOffset;
						centerOffset=tmp;
					}
					
					int bin = params.findBin(centerOffset);
					addToArray(bin, bin, array, Z0);
				}
			}
		}else{
			//Experiment hits
			for(StrandedRegion r : hits){
				if(r.overlaps(query)){
					int startOffset = Math.max(0, r.getStart()-query.getStart());
					int endOffset = Math.max(0, Math.min(query.getEnd(), r.getEnd()-query.getStart()));
					
					if(!strand) { 
						int tmpEnd = window-startOffset;
						int tmpStart = window-endOffset;
						startOffset = tmpStart;
						endOffset = tmpEnd;
					}
					
					int startbin = params.findBin(startOffset);
					int endbin = params.findBin(endOffset);
					
					addToArray(startbin, endbin, array, 1/totalHits);
				}
			}
			//Control hits
			if(backTotalHits>0){
				for(StrandedRegion r : backhits){
					if(r.overlaps(query)){
						int startOffset = Math.max(0, r.getStart()-query.getStart());
						int endOffset = Math.max(0, Math.min(query.getEnd(), r.getEnd()-query.getStart()));
						
						if(!strand) { 
							int tmpEnd = window-startOffset;
							int tmpStart = window-endOffset;
							startOffset = tmpStart;
							endOffset = tmpEnd;
						}
						
						int startbin = params.findBin(startOffset);
						int endbin = params.findBin(endOffset);
						
						addToArray(startbin, endbin, array, -1*(1/backTotalHits));
					}
				}
			}
		}
		return new PointProfile(a, params, array, (a instanceof StrandedPoint));
	}

	private void addToArray(int i, int j, double[] array, double value) { 
		for(int k = i; k <= j; k++) { 
			array[k] += value;
		}
	}
	
	/* Binomial test for differences between two population proportions */
	private double binomialSampleEquality(double X1, double X2, double n1, double n2){
		double P1 = X1/n1;
		double P2 = X2/n2;
		double P = (X1+X2)/(n1+n2);
		double Z0 = (P1-P2)/(Math.sqrt(P*(1-P)*((1/n1)+(1/n2))));
		return(Z0);
	}
	
	private int [] makeHitLandscape(LinkedList<StrandedRegion> hits, Region currReg, int winWidth){
		int numBins = (int)(currReg.getWidth()/winWidth);
		int [] land = new int[numBins+1];
		for(int i=0; i<=numBins; i++){land[i]=0;}
		for(StrandedRegion r : hits){
			if(r.overlaps(currReg)){
				int binstart = (int)Math.max(0, ((double)((r.getStart()-currReg.getStart())/winWidth)-(Math.floor(winWidth/winWidth)-1)));
				int binend = (int)(Math.min((double)(r.getEnd()-currReg.getStart()), (double)currReg.getWidth())/winWidth);
				for(int i=binstart; i<=binend; i++){
					land[i]++; 
				}
			}
		}
		return(land);
	}
	
	//No cleanup
	public void cleanup(){
		for(SeqExptHandler h : handles)
			h.close();
	}
}
