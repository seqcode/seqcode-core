package org.seqcode.motifs.scores2motifs;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.RepeatMaskedRegion;
import org.seqcode.genome.location.ScoredPoint;
import org.seqcode.genome.sequence.SequenceUtils;
import org.seqcode.gseutils.Pair;
import org.seqcode.motifs.MemeER;


public class Scores2Motifs {
	
	protected S2MConfig s2mConfig;
	protected List<Hill> hills;
	protected String[] randomSequences = new String[MemeER.MOTIF_FINDING_NEGSEQ];
	
	
	/** List of all discriminative motifs for all the labels*/ 
	protected List<WeightMatrix> discrimMotifs = new ArrayList<WeightMatrix>();
	
	/** The ROC value of the identified motifs that indicates their significance */
	protected HashMap<String, Double> discrimMotifsRocs = new HashMap<String,Double>();
	
	

	/**
	 * Constructor
	 * @param h2mcon
	 */
	public Scores2Motifs(S2MConfig h2mcon) {
		s2mConfig = h2mcon;
	}
	
	/**
	 * execute: Run the method
	 * @throws IOException
	 */
	public void execute() throws IOException{
		makeMemeDirs();
		
		//Load a set of random regions
		setRandomRegs();
		
		//Define the hills according to one of the input options
		setHills();
		
		//Cluster hills
		
		
		//Run motif-finding
		
		
	}

	
	/**
	 * setHills:
	 * 	There is currently only one way to establish hills; by reading in a list of scored binding locations from a file.
	 * 	However, in the future we will want to extend this to determining hills from a scored landscape, or even from a model + sequence
	 * @return
	 */
	protected List<Hill> setHills(){
		hills = new ArrayList<Hill>();
		
		if(s2mConfig.getPeaksFile() != null){
			//Load the scored points from a file
			List<ScoredPoint> peaks = processScoredPoints();

			for(ScoredPoint sp : peaks){
				//Expand to region
				Region searchWin = sp.expand(s2mConfig.getMemeSearchWin());
				//Load sequences for each point
				String seq = s2mConfig.getSeqGen().execute(searchWin);
				//Define the Hill objects
				hills.add(new Hill(searchWin, sp.getScore(), seq));
			}
		}
		
		//Add k-mer profiles to hills
		
		
		//sort by score
		
		
		
		return hills;
	}
	
	
	
	
	
	/**
	 * Read the scored points from a file and optionally screen repeats. 
	 * @return
	 */
	protected List<ScoredPoint> processScoredPoints(){
		if(s2mConfig.getPeaksFile() != null){
			// Loading peaks
			List<ScoredPoint> sp = RegionFileUtilities.loadScoredPointsFromFile(s2mConfig.getGenome(), s2mConfig.getPeaksFile());
			
			//Filter repeats if appropriate
			if(s2mConfig.getScreenReps()){
				if(sp.size() > 0){

					Iterator<ScoredPoint> peakItr = sp.iterator();
					while(peakItr.hasNext()){
						ScoredPoint currPeak = peakItr.next();
						Region currReg = currPeak.expand(s2mConfig.getRepMaskWin());
						double repLen = 0;
						Iterator<RepeatMaskedRegion> repItr = s2mConfig.getRepMask().execute(currReg);
						while (repItr.hasNext()) {
							RepeatMaskedRegion currRep = repItr.next();
							if (currRep.overlaps(currReg)) {
								repLen += (double) currRep.getOverlapSize(currReg);
							}
						}
						if (repLen / (double) currReg.getWidth() > s2mConfig.getRepPropLimit()){
							peakItr.remove();
						}
					}
				}
			}
			return sp;
		}
		return null;
	}
	
	// Generate random sequences later need for meme analysis to assess motifs
	protected void setRandomRegs(){
		List<Region> randomRegions = MemeER.randomRegionPick(s2mConfig.getGenome(), null, MemeER.MOTIF_FINDING_NEGSEQ,s2mConfig.getMemeSearchWin());
		for(int i=0; i<randomRegions.size(); i++){
			randomSequences[i] = s2mConfig.getSeqGen().execute(randomRegions.get(i));
		}
	}

	
	
	//Make MEME directories
	protected void makeMemeDirs(){
		File memeDir = new File(s2mConfig.getOutDir().getAbsoluteFile()+File.separator+"meme");
				memeDir.mkdirs();
	}

}
