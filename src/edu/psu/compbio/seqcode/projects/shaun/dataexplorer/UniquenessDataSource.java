package edu.psu.compbio.seqcode.projects.shaun.dataexplorer;

import java.util.Iterator;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.ScoredRegion;
import edu.psu.compbio.seqcode.genome.location.StrandedRegion;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.ScoredRegionGenerator;

public class UniquenessDataSource extends DataSource{

		private Genome gen;
		private ScoredRegionGenerator scorer;
		
		public UniquenessDataSource(Genome gen, String tableName, double threshold, double weight){
			super(tableName, threshold, weight);
			
			scorer = new ScoredRegionGenerator(gen, tableName);
		}
		
		public double genData(Region r) {
			double uProp=0;
			double uCount=0;
			Iterator<ScoredRegion> iter = scorer.execute(r);
			while(iter.hasNext()){
				ScoredRegion sr = (ScoredRegion) iter.next();
				//System.out.println(r.getLocationString()+"\t"+sr.getLocationString()+"\t"+sr.getOverlapSize(r)+"\t"+sr.getScore());
				if(sr.getScore()==1)
					uCount+=sr.getOverlapSize(r);
			}
			uProp=uCount/(double)r.getWidth();
			return uProp;
		}

		public double[] genDataVector(Region r, int binSize) {
			int numBins = (r.getWidth()/binSize)+1;
			double [] vec = new double[numBins];
			for(int i=0; i<numBins; i++){
				Region s = new Region(r.getGenome(), r.getChrom(), r.getStart()+(i*binSize), r.getStart()+((i+1)*binSize)-1);
				double uCount=0;
				Iterator<ScoredRegion> iter = scorer.execute(s);
				while(iter.hasNext()){
					ScoredRegion sr = (ScoredRegion) iter.next();
					if(sr.getScore()==1)
						uCount+=sr.getOverlapSize(r);
				}
				vec[i]=uCount/(double)s.getWidth();
			}
			if(r instanceof StrandedRegion)
				if(((StrandedRegion) r).getStrand()=='-')
					return reverseVec(vec);
			return(vec);
		}
		public void cleanup(){}
}
