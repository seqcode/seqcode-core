package edu.psu.compbio.seqcode.projects.shaun.dataexplorer;

import java.util.Iterator;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.ScoredRegion;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ScoredRegionGenerator;
import edu.psu.compbio.seqcode.gse.tools.sequence.StrandedRegionsToFasta;

public class ConsDataSource extends DataSource{

		private Genome gen;
		private ScoredRegionGenerator scorer;
		
		public ConsDataSource(Genome gen, String tableName, double threshold, double weight){
			super(tableName, threshold, weight);
			
			scorer = new ScoredRegionGenerator(gen, tableName);
		}
		
		public double genData(Region r) {
			double consAvg=0;
			
			double cSum=0;
			Iterator<ScoredRegion> iter = scorer.execute(r);
			while(iter.hasNext()){
				ScoredRegion sr = (ScoredRegion) iter.next();
				//System.out.println(r.getLocationString()+"\t"+sr.getLocationString()+"\t"+sr.getOverlapSize(r)+"\t"+sr.getScore());
				cSum+=sr.getScore()*(sr.getOverlapSize(r));
			}
			consAvg=cSum/(double)r.getWidth();
			return consAvg;
		}

		public double[] genDataVector(Region r, int binSize) {
			int numBins = (r.getWidth()/binSize)+1;
			double [] vec = new double[numBins];
			for(int i=0; i<numBins; i++){
				Region s = new Region(r.getGenome(), r.getChrom(), r.getStart()+(i*binSize), r.getStart()+((i+1)*binSize)-1);
				double cSum=0;
				Iterator<ScoredRegion> iter = scorer.execute(s);
				while(iter.hasNext()){
					ScoredRegion sr = (ScoredRegion) iter.next();
					cSum+=sr.getScore()*(sr.getOverlapSize(s));
				}
				vec[i]=cSum/(double)r.getWidth();
			}
			if(r instanceof StrandedRegion)
				if(((StrandedRegion) r).getStrand()=='-')
					return reverseVec(vec);
			return(vec);
		}
		public void cleanup(){}
}
