package org.seqcode.viz.metaprofile;

import java.util.*;

import org.seqcode.data.seqdata.*;
import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.Sample;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.genome.sequence.SequenceUtils;


/**
 * Stranded5PrimeProfiler profiles the occurrence of sequencing reads around points
 *  
 * @author: tdanford
 * Date: Aug 19, 2008
 */

public class Stranded5PrimeProfiler implements PointProfiler<Point,PointProfile> {
	
	private Genome genome;
	private ExperimentManager manager=null;
	private BinningParameters params;
	private char strand;
	private char base='.'; //Only plot tags that have this base at baseRelPosition
	private int baseRelPosition=0;  //Only plot tags that have this base at baseRelPosition
	private SequenceGenerator seqgen=null;
	private int fivePrimeShift = 0;
	private boolean cpmNorm = false;
	private boolean signal = true; // Profile signal
	
	
	public Stranded5PrimeProfiler(GenomeConfig genConfig, BinningParameters ps, ExperimentManager man, char strand, int fivePrimeShift, char base, int baseRelPosition, boolean profileSignal, boolean countsPerMillion) {
		genome = genConfig.getGenome();
		manager = man;
		params = ps;
		this.strand = strand;
		this.base = base;
		this.baseRelPosition = baseRelPosition;
		this.fivePrimeShift=fivePrimeShift;
		this.signal=profileSignal;
		this.cpmNorm=countsPerMillion;
		
		if(base != '.'){
			seqgen = genConfig.getSequenceGenerator();
		}
	}
	
	
	
	public BinningParameters getBinningParameters() {
		return params;
	}

	public PointProfile execute(Point a) {
		int window = params.getWindowSize();
		int upstream = window/2;
		int downstream = window-upstream-1;
		char pointStrand = '+';
		
		if(a instanceof StrandedPoint)
			pointStrand = ((StrandedPoint)a).getStrand();
		boolean wantPosStrandReads = this.strand=='+';
		if(pointStrand == '-')
			wantPosStrandReads = !wantPosStrandReads;
		char wantedStrand = wantPosStrandReads?'+':'-';
				
		int start = pointStrand == '+' ?  Math.max(1, a.getLocation()-upstream) : Math.max(1, a.getLocation()-downstream);
		int end = pointStrand == '+' ?  Math.min(a.getLocation()+downstream, a.getGenome().getChromLength(a.getChrom())) : Math.min(a.getLocation()+upstream, a.getGenome().getChromLength(a.getChrom()));
		Region query = new Region(a.getGenome(), a.getChrom(), start, end);
		int ext = 200;
		Region extQuery = new Region(a.getGenome(), a.getChrom(), start-ext>0 ? start-ext : 1, end+ext < a.getGenome().getChromLength(a.getChrom()) ? end+ext : a.getGenome().getChromLength(a.getChrom()) );
		
		//Get sequence if required
		char [] seq=null;
		if(seqgen!=null)
			seq = seqgen.execute(extQuery).toCharArray();
		
		//Populate the tag density
		double[] array = new double[params.getNumBins()];
		for(int i = 0; i < array.length; i++) { array[i] = 0; }
		for(ControlledExperiment expt : manager.getReplicates()){
			Sample currSample = signal ? expt.getSignal() : expt.getControl();	
			double normFactor = cpmNorm ? 1000000/currSample.getHitCount() : 1.0;
			List<StrandedBaseCount> sbcs = currSample.getBases(extQuery);
			for(StrandedBaseCount sbc : sbcs){
				if(base=='.' || (seq!=null && getBaseAtPosition(sbc, baseRelPosition, extQuery, seq)==base)){
					SeqHit hit = new SeqHit(genome, a.getChrom(), sbc);
					if (this.strand=='.' || hit.getStrand()==wantedStrand){  //only count one strand
						if (start<=hit.getFivePrime() && end>=hit.getFivePrime()){
							int hitPos = hit.getStrand()=='+' ? hit.getFivePrime()+fivePrimeShift : hit.getFivePrime()-fivePrimeShift; 
							int hit5Prime = hitPos-start;
							if(pointStrand=='-')
								hit5Prime = end-hitPos;
							array[params.findBin(hit5Prime)]+=(hit.getWeight()*normFactor);
						}				
					}
				}
			}
		}
		
		return new PointProfile(a, params, array, (a instanceof StrandedPoint));
	}
	
	private char getBaseAtPosition(StrandedBaseCount a, int position, Region queryReg, char[] seq){
		char b = '.';
		int wantedPos = a.getStrand()=='+' ? 
				a.getCoordinate()+position : 
				a.getCoordinate()-position;
		if(wantedPos>=queryReg.getStart() && wantedPos<queryReg.getEnd()){
			b = a.getStrand()=='+' ? seq[wantedPos-queryReg.getStart()] : SequenceUtils.complementChar(seq[wantedPos-queryReg.getStart()]);
		}
		return b;
	}
	
	public void cleanup(){
	}
}
