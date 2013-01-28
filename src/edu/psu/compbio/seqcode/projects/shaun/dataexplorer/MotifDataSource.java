package edu.psu.compbio.seqcode.projects.shaun.dataexplorer;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.motifs.MarkovBackgroundModel;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.ewok.verbs.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.motifs.WeightMatrixScoreProfile;
import edu.psu.compbio.seqcode.gse.ewok.verbs.motifs.WeightMatrixScorer;
import edu.psu.compbio.seqcode.projects.shaun.MarkovMotifThresholdFinder;
import edu.psu.compbio.seqcode.projects.shaun.Score2Sp;

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
