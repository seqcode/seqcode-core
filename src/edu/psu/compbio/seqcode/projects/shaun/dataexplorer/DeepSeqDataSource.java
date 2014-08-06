package edu.psu.compbio.seqcode.projects.shaun.dataexplorer;

import java.util.List;

import cern.jet.random.Binomial;
import cern.jet.random.engine.DRand;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.StrandedRegion;
import edu.psu.compbio.seqcode.gse.projects.gps.DeepSeqExpt;
import edu.psu.compbio.seqcode.gse.projects.gps.ExtReadHit;
import edu.psu.compbio.seqcode.gse.utils.Pair;

public class DeepSeqDataSource extends DataSource{

	private Pair<DeepSeqExpt, DeepSeqExpt> experiment;
	private double totalIP=0, totalCtrl=0;
	private double scale;
	private boolean controlled = true; //if true, return -log(P-value), if false, return IP weight counts 
	
	public DeepSeqDataSource(DeepSeqExpt expt, DeepSeqExpt ctrl, String exptName, double threshold, double weight){
		super(exptName, threshold, weight);
		experiment = new Pair<DeepSeqExpt, DeepSeqExpt>(expt, ctrl);
		totalIP=expt.getWeightTotal();
		totalCtrl=ctrl==null ? 0 : ctrl.getWeightTotal();
		controlled = ctrl==null ? false : true;
		scale = totalCtrl/totalIP;
		
		System.out.println(exptName+"\t"+totalIP+" signal\t"+totalCtrl+" control");
	}
	
	
	public double genData(Region r) {
		if(controlled){
			double currSig = getIP().sumWeights(r);
			double currBack = getCtrl().sumWeights(r);
			//System.out.println(currSig+"\t"+currBack);
			double backScale=currBack/scale;
			double score = currSig>0 ? -1*Math.log(binomialPValue(backScale, backScale+currSig)) : 0;
			return score;
		}else{
			double currSig = getIP().sumWeights(r);
			//return (currSig/(r.getWidth()*totalIP));			
			return (currSig);
		}
	}
	public double[] genDataVector(Region r, int binSize) {
		int numBins = (r.getWidth()/binSize)+1;
		double [] vec = new double[numBins];
		for(int i=0; i<numBins; i++){vec[i]=0;}
		
		if(controlled){
			double [] dvec = new double[numBins];
			double [] cvec = new double[numBins];
			for(int i=0; i<numBins; i++){dvec[i]=0; cvec[i]=0;}

			List<ExtReadHit> ipHits = getIP().loadExtHits(r);
			List<ExtReadHit> ctrlHits = getCtrl().loadExtHits(r);
			
			for(ExtReadHit h : ipHits){
				int startBin=Math.max(0, (h.getStart()-r.getStart())/binSize);
				int endBin=Math.max(0, Math.min(numBins-1, (h.getEnd()-r.getStart())/binSize));
				for(int x=startBin; x<=endBin; x++){dvec[x]++;}
			}for(ExtReadHit h : ctrlHits){
				int startBin=Math.max(0, (h.getStart()-r.getStart())/binSize);
				int endBin=Math.max(0, Math.min(numBins-1, (h.getEnd()-r.getStart())/binSize));
				for(int x=startBin; x<=endBin; x++){cvec[x]++;}
			}
			
			for(int i=0; i<numBins; i++){
				double backScale=cvec[i]/scale;
				vec[i] = dvec[i]>0 ? -1*Math.log(binomialPValue(backScale, backScale+dvec[i])) : 0;
			}
		}else{
			List<ExtReadHit> ipHits = getIP().loadExtHits(r);
			for(ExtReadHit h : ipHits){
				int startBin=Math.max(0, (h.getStart()-r.getStart())/binSize);
				int endBin=Math.max(0, Math.min(numBins-1, (h.getEnd()-r.getStart())/binSize));
				for(int x=startBin; x<=endBin; x++){vec[x]++;}
			}
		}
		if(r instanceof StrandedRegion)
			if(((StrandedRegion) r).getStrand()=='-')
				return reverseVec(vec);
		return(vec);
	}

	private DeepSeqExpt getIP(){return experiment.car();}
	private DeepSeqExpt getCtrl(){return experiment.cdr();}
	public void cleanup(){
		if(getIP()!=null)
			getIP().closeLoaders();
		if(getCtrl()!=null)
			getCtrl().closeLoaders();
	}
	
	//Binomial CDF assuming scaled control. Uses COLT binomial test
	// k=scaled control, n=scaled control+signal
	protected double binomialPValue(double k, double n){
		double pval=1;
		Binomial b = new Binomial((int)Math.ceil(n), 0.5, new DRand());
		pval = b.cdf((int) Math.ceil(k));
		return(pval);		
	}
}
