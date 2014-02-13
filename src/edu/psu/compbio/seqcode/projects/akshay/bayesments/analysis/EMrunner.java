package edu.psu.compbio.seqcode.projects.akshay.bayesments.analysis;

import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.bayesnet.EMtrain;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.bayesnet.MAPassignment;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.features.GenomicLocations;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.features.Sequences;
import edu.psu.compbio.seqcode.projects.akshay.bayesments.framework.Config;

public class EMrunner {
	
	protected EMtrain model;
	protected boolean onlyChrom;
	protected int num_chrom_itrs;
	protected int num_seq_itrs;
	protected int total_itrs;
	protected GenomicLocations chromdata;
	protected Sequences seqdata;
	protected boolean finishedTraining;
	protected ExperimentManager manager;
	protected List<WeightMatrix> globalListOfMotifs;
	
	protected Config config;
	
	public EMrunner(Config conf, GenomicLocations chromdata, ExperimentManager manager) {
		config = conf;
		this.chromdata = chromdata;
		this.onlyChrom = conf.runOnlyChrom();
		this.num_chrom_itrs = conf.getNumChromIters();
		this.num_seq_itrs = conf.getNumSeqIters();
		this.total_itrs = conf.getNumItrs();
		this.manager = manager;
	}
	
	public void trainModel(){
		//Initializing the EM train class
		// Since we are just initializing the model, it is not in seq state
		this.model = new EMtrain(this.config, this.chromdata,this.manager);
		// Running the EM only on chromatin features
		model.runEM(num_chrom_itrs);
		if(!this.onlyChrom){ // Now we have to initialze the seq data and features
			// Before doing anything lets fetch DNA sequences at the peak pair locations
			this.seqdata = new Sequences(config,chromdata.getLocations());
			// Firsly, get the top 200 (or all, whichever is maximum) sites for all chromatin state assignments
				// To do that , lets us first do a Map assignment
			MAPassignment onlyChromAddignment = new MAPassignment(model,config);
			double[][] assignment = onlyChromAddignment.getMapAssignments();
				// Now get the top 200 seqs and add them to alltopseqs
			List<String> alltopseqs = new ArrayList<String>();
			for(int j=0; j<config.getNumChrmStates(); j++){
				List<String> seqs = new ArrayList<String>();
				for(int i=0; i<chromdata.getNumTrainingExamples();i++){
					if(j==assignment[i][0]){
						seqs.add(seqdata.getIthSeq(i));
					}
					if(seqs.size()==200){
						break;
					}
				}
				for(String s : seqs){
					alltopseqs.add(s);
				}
			}
			//Secondly, run meme on these top seqs and store them in globalListOfMotifs
			this.globalListOfMotifs = this.runMeme(alltopseqs);
			//Thirdly, scan these motifs over the peak- pair locations and convert them into chip like signals
			this.seqdata.setMotifs(globalListOfMotifs);
			this.seqdata.setXc();
			//Finaly, run the remaining EM iterations with thses new seqence features
			this.model.setSeqMode(seqdata, seqdata.getXs());
			model.runEM(num_seq_itrs);
		}
		
		
		
	}
	
	private List<WeightMatrix> runMeme(List<String> seqs){
		List<WeightMatrix> ret;
		MEMERunner meme = new MEMERunner(config);
		Pair<List<WeightMatrix>, List<WeightMatrix>> wm =meme.execute(seqs, "results", false);
		ret = wm.car();
		return ret;
		
	}
	
	//Accessors
	public double[] getPIj(){return model.getPIj();}
	public double[][] getBjk(){return model.getBjk();}
	public double[][] getMUc(){return model.getMUc();}
	public double[][] getMUf(){return model.getMUf();}
	public double[][] getMUs(){return model.getMUs();}
	public double[][] getSIGMAc(){return model.getSIGMAc();}
	public double[][] getSIGMAf(){return model.getSIGMAf();}
	public double[][] getSIGMAs(){return model.getSIGMAs();}
	public EMtrain getModel(){return this.model;}
	
}
