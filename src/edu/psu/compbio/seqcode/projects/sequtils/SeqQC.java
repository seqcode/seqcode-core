package edu.psu.compbio.seqcode.projects.sequtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ChromosomeGenerator;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.multigps.framework.Config;
import edu.psu.compbio.seqcode.projects.multigps.framework.StrandedBaseCount;

public class SeqQC {

	Config config;
	Genome gen;
	ExperimentManager manager;
	float[] startcounts =null;
	float densityWindow = 500;
	
	/**
	 * Constructor
	 * @param c
	 * @param man
	 */
	public SeqQC(Config c, ExperimentManager man){
		config = c;
		gen = config.getGenome();
		manager = man;
	}
	
	/**
	 * estimateLibrarySize
	 */
	public void estimateLibrarySize(){
	
		//If we have multiple experiments, process one at a time
		for(ControlledExperiment expt : manager.getExperimentSet().getReplicates()){
			
			// 1: Find the density of covered bases surrounding each read
			List<DensityCountPair> densities = new ArrayList<DensityCountPair>();
			Iterator<Region> chroms = new ChromosomeGenerator().execute(gen);
			while (chroms.hasNext()) {
				Region currentRegion = chroms.next();
				//Split the job up into large chunks
	            for(int x=currentRegion.getStart(); x<=currentRegion.getEnd(); x+=config.MAXSECTION){
	                int y = x+config.MAXSECTION; 
	                if(y>currentRegion.getEnd()){y=currentRegion.getEnd();}
	                Region currSubRegion = new Region(gen, currentRegion.getChrom(), x, y);
				
					List<StrandedBaseCount> ipHits = expt.getSignal().getUnstrandedBases(currSubRegion);
					makeHitStartArray(ipHits, currSubRegion, '+');
					float posCounts[] = startcounts.clone();
					makeHitStartArray(ipHits, currSubRegion, '-');
					float negCounts[] = startcounts.clone();
		            
					int halfDWin = (int)densityWindow/2;
					for(int i=currSubRegion.getStart()+halfDWin; i<currSubRegion.getEnd()-halfDWin; i++){  //Note that this means a tiny fraction of reads that are located on the large block boundaries will be ignored
						if(posCounts[i]>0){ //Treat fragments on opposite strands separately
							float dens =0; 
							for(int j=i-halfDWin; j<i+halfDWin; j++){
								if(posCounts[j]>0 && j!=i)
									dens++;
								if(negCounts[j]>0)
									dens++;
							}
							dens /= (densityWindow*2);
							densities.add(new DensityCountPair(dens, posCounts[i]));
						}
						if(negCounts[i]>0){ //Treat fragments on opposite strands separately
							float dens =0; 
							for(int j=i-halfDWin; j<i+halfDWin; j++){
								if(posCounts[j]>0)	
									dens++;
								if(negCounts[j]>0 && j!=i)
									dens++;
							}
							dens /= (densityWindow*2);
							densities.add(new DensityCountPair(dens, negCounts[i]));
						}
					}
	            }
			}
			Collections.sort(densities); //Sort the density pairs in increasing order
			
			//2: Generate a read count per base distribution for the lowest density sites. 
		}
	}
	
	//Makes integer arrays corresponding to the read starts over the current region
	protected void makeHitStartArray(List<StrandedBaseCount> hits, Region currReg, char strand){
		startcounts = new float[currReg.getWidth()+1];
		for(int i=0; i<=currReg.getWidth(); i++){startcounts[i]=0;}
		for(StrandedBaseCount r : hits){
			if(strand=='.' || r.getStrand()==strand){
				int offset=inBounds(r.getCoordinate()-currReg.getStart(),0,currReg.getWidth());
				startcounts[offset]+=r.getCount();
			}
		}
	}
	protected final int inBounds(int x, int min, int max){
		if(x<min){return min;}
		if(x>max){return max;}
		return x;
	}
	
	/**

	 */
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
		if(args.length==0 || ap.hasKey("h") || !ap.hasKey("emp")){
			System.err.println("SeqQC:\n" +
					"");
		}else{
			Config con = new Config(args, false);
			ExperimentManager man = new ExperimentManager(con);
			
			
		}
	}
	
	
	private class DensityCountPair implements Comparable<DensityCountPair>{
		private float density;
		private float count;
		
		public DensityCountPair(float d, float c){
			density = d;
			count =c;
		}
		
		public void setCount(float count) {
			this.count = count;
		}
		public void setDensity(float density) {
			this.density = density;
		}
		public float getCount() {
			return count;
		}
		public float getDensity() {
			return density;
		}
		
		// sort according to density
		public int compareTo(DensityCountPair dc) {
			if(density < dc.density)
				return -1;
			else if (density > dc.density)
				return 1;
			else return 0;
		}

	}
}
