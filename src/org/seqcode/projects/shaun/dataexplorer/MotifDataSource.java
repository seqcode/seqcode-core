package org.seqcode.projects.shaun.dataexplorer;

import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.gse.datasets.motifs.MarkovBackgroundModel;
import org.seqcode.gse.datasets.motifs.WeightMatrix;
import org.seqcode.gse.gsebricks.verbs.motifs.WeightMatrixScoreProfile;
import org.seqcode.gse.gsebricks.verbs.motifs.WeightMatrixScorer;
import org.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import org.seqcode.projects.shaun.MarkovMotifThresholdFinder;
import org.seqcode.projects.shaun.Score2Sp;

public class MotifDataSource extends DataSource{

	private WeightMatrix motif;
	private WeightMatrixScorer scorer;
	private SequenceGenerator seqgen;
	private Score2Sp spCalc;
	
	public MotifDataSource(WeightMatrix wm, MarkovBackgroundModel back, double threshold, double weight) {
		super(wm.getName(), threshold, weight);
		motif=wm;
		
		scorer = new WeightMatrixScorer(motif);
		seqgen = new SequenceGenerator();
		MarkovMotifThresholdFinder mmtf = new MarkovMotifThresholdFinder(wm,back);
		spCalc = mmtf.getMotifROC();
	}

	public double genData(Region r) {
		String seq = seqgen.execute(r);
		WeightMatrixScoreProfile profiler = scorer.execute(seq);
		//return(profiler.getMaxScore());
		return(spCalc.getSp(profiler.getMaxScore()));
	}

	public double[] genDataVector(Region r, int binSize) {
		int numBins = (r.getWidth()/binSize)+1;
		double [] vec = new double[numBins];
		for(int i=0; i<numBins; i++){
			String seq = seqgen.execute(new Region(r.getGenome(), r.getChrom(), r.getStart()+(i*binSize), r.getStart()+((i+1)*binSize)-1+motif.length()));
			WeightMatrixScoreProfile profiler = scorer.execute(seq);
			//return(profiler.getMaxScore());
			vec[i]=spCalc.getSp(profiler.getMaxScore());
		}
		if(r instanceof StrandedRegion)
			if(((StrandedRegion) r).getStrand()=='-')
				return reverseVec(vec);
		return(vec);
	}
	public void cleanup(){}
}
