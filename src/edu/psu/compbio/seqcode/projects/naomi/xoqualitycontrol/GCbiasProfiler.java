package edu.psu.compbio.seqcode.projects.naomi.xoqualitycontrol;

import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;

import edu.psu.compbio.seqcode.deepseq.StrandedBaseCount;
import edu.psu.compbio.seqcode.deepseq.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.deepseq.experiments.Sample;
import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.ChromosomeGenerator;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;

public class GCbiasProfiler {
	
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected SequenceGenerator<Region> seqgen;
	
	protected double[][] GCvsFragmentRate; 
	
	public GCbiasProfiler (GenomeConfig gcon, ExptConfig econ, SequenceGenerator<Region> seqgenerator){	
		gconfig = gcon;
		econfig = econ;
		seqgen = seqgenerator;
	}
	
	public void getSinglePositionModel(int a, int l){
		
		ExperimentManager manager = new ExperimentManager(econfig);
		Genome genome = gconfig.getGenome();
		seqgen = new SequenceGenerator<Region>(genome);
		
		Iterator<Region> chroms = new ChromosomeGenerator<Genome>().execute(genome);
		double chromNum = (double) genome.getChromList().size();
		
		GCvsFragmentRate = new double[l+1][2];
		for (int i =0; i<l+1; i++){
			GCvsFragmentRate[i][0] = (double) i/(double) l;
			GCvsFragmentRate[i][1] = 0;
		}
				
		while (chroms.hasNext()) {
	
			Region currChrom = chroms.next();
			int currchromSize = currChrom.getWidth()+1;
			
			Map<Sample, List<StrandedBaseCount>> sampleCountsMap = new HashMap<Sample, List<StrandedBaseCount>>();			
			for (ControlledExperiment rep: manager.getReplicates()){
				sampleCountsMap.put(rep.getSignal(), rep.getSignal().getBases(currChrom));
//				sampleCountsMap.put(rep.getControl(), rep.getControl().getBases(currChrom));			
			}
			
			//first dimention contains hit counts, second dimention contains strand
			Map<Sample, double[][]> strandedSampleCounts = new HashMap<Sample, double[][]>();
			
			for (ControlledExperiment rep: manager.getReplicates()){
				List<StrandedBaseCount> currentCounts = sampleCountsMap.get(rep.getSignal());
				double[][] strandedCounts = new double[currchromSize][2];
				for (StrandedBaseCount hits :currentCounts){
					if (hits.getStrand()=='+')
						strandedCounts[hits.getCoordinate()][0]+=hits.getCount();
					else
						strandedCounts[hits.getCoordinate()][1]+=hits.getCount();
				}
				strandedSampleCounts.put(rep.getSignal(),strandedCounts);
				currentCounts=null;
			}			
			
			double[][] gcSinglePositionModel = new double[l+1][3];
			for (int i = 0; i <l+1; i++){
				for (int j = 0 ; j< 3; j++)
					gcSinglePositionModel[i][j] = 0;
			}
						
			for (ControlledExperiment rep: manager.getReplicates()){
				
				List<StrandedBaseCount> controlCounts = rep.getControl().getBases(currChrom);
				double[][] strandedCounts = strandedSampleCounts.get(rep.getSignal());
				
				// sampling 1/10 of hits in controlCounts
				for (int randomIteration = 0; randomIteration < controlCounts.size()/10; randomIteration ++){
				
					Random randomizer = new Random();
					StrandedBaseCount randomBase = controlCounts.get(randomizer.nextInt(controlCounts.size()));
					Region reg = null;
					if (randomBase.getStrand()=='+' && randomBase.getCoordinate()+a+l<currchromSize)
						reg = new Region(genome,currChrom.getChrom(),randomBase.getCoordinate()+a,randomBase.getCoordinate()+a+l);
					else if (randomBase.getStrand()=='-' && randomBase.getCoordinate()-a-l>=0)
						reg = new Region(genome,currChrom.getChrom(),randomBase.getCoordinate()-a-l, randomBase.getCoordinate()-a);
					if (reg !=null){
						String seq = seqgen.execute(reg);
						seq =seq.toUpperCase();
						int gcScore = 0;
						for(int i=0; i<seq.length()-1; i++){
							if (seq.charAt(i) == 'C'|| seq.charAt(i) == 'G')
								gcScore++;
						}
						
						//calculate Posi(N gc)
						gcSinglePositionModel[gcScore][0]++;
					
						//calculate Frag(F gc)
						if (randomBase.getStrand()=='+')
							gcSinglePositionModel[gcScore][1] += strandedCounts[randomBase.getCoordinate()][0];
						else
							gcSinglePositionModel[gcScore][1] += strandedCounts[randomBase.getCoordinate()][1];
					}
				}//end of random iteration from controlCounts

				//calculate Rate (lambda gc)
				for (int i = 0; i<l+1;i++){
					if (gcSinglePositionModel[i][1] != 0)
						gcSinglePositionModel[i][2] = (double) gcSinglePositionModel[i][1]/ (double) gcSinglePositionModel[i][0];
				}				
			}//end of sample hits iteration
			for (int i = 0; i<l+1;i++){
				GCvsFragmentRate[i][1] += gcSinglePositionModel[i][2]/chromNum;
			}
		}// end of chromosome iteration
		manager.close();
	}
	
	public void printGCvsFragmentRate(){
		
		for (int i = 0; i<GCvsFragmentRate.length; i++){
			System.out.println(GCvsFragmentRate[i][0]+" : "+GCvsFragmentRate[i][1]);
		}		
	}
	
	public static void main(String[] args) {
		
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig  econf = new ExptConfig(gconf.getGenome(), args);
		SequenceGenerator<Region> sequencegen = new SequenceGenerator<Region>(gconf.getGenome());
		ArgParser ap = new ArgParser(args);
		
		int posStart = 0;
		if (ap.hasKey("a"))
			posStart = Args.parseInteger(args,"a",3);
		
		int length = 180;
		if (ap.hasKey("l"))
			length = Args.parseInteger(args,"l",3);
		
		GCbiasProfiler gfProfile = new GCbiasProfiler(gconf, econf, sequencegen); 
		gfProfile.getSinglePositionModel(posStart,length);
		gfProfile.printGCvsFragmentRate();	
		
	}
}
