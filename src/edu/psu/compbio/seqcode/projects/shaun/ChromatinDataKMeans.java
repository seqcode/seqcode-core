package edu.psu.compbio.seqcode.projects.shaun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.Vector;

import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.machinelearning.clustering.ClusterRepresentative;
import edu.psu.compbio.seqcode.machinelearning.clustering.PairwiseElementMetric;
import edu.psu.compbio.seqcode.machinelearning.clustering.kmeans.KMeansClustering;
import edu.psu.compbio.seqcode.machinelearning.clustering.vectorcluster.DefaultVectorClusterElement;
import edu.psu.compbio.seqcode.machinelearning.clustering.vectorcluster.EuclideanDistance;
import edu.psu.compbio.seqcode.machinelearning.clustering.vectorcluster.Mean;
import edu.psu.compbio.seqcode.machinelearning.clustering.vectorcluster.VectorClusterElement;

public class ChromatinDataKMeans {

	private int K=2;
	private String inFileName=null;
	private int DIM=1;
	private int TIMES = 1;
	private ArrayList<ArrayList<DefaultVectorClusterElement>> allData = new ArrayList<ArrayList<DefaultVectorClusterElement>>();
	private ArrayList<VectorClusterElement> concatData = new ArrayList<VectorClusterElement>();
	private ArrayList<VectorClusterElement> trainData = new ArrayList<VectorClusterElement>();
	private int TRAIN_REPEATS=10;
	private int MAX_ITER = 100;
	private int MAX_TRAIN_DATA=100000;
	
	public ChromatinDataKMeans(String file, int K, int TIMES, int DIM){
		this.inFileName=file;
		this.K = K;
		this.TIMES = TIMES;
		this.DIM = DIM;	
		for(int t=0; t<TIMES; t++){
			allData.add(new ArrayList<DefaultVectorClusterElement>());
		}
	}
	
	public void cluster(){
		//Initialize
		PairwiseElementMetric<VectorClusterElement> metric = new EuclideanDistance<VectorClusterElement>();
		ClusterRepresentative<VectorClusterElement> crep = new Mean();
		Vector<VectorClusterElement> bestClusterMeans=null;
		Double SSD = Double.MAX_VALUE;
		
		for(int i=0; i<=TRAIN_REPEATS; i++){
			
			//For now, we will cluster datapoints from all times together. 
			trainData = new ArrayList<VectorClusterElement>();
			Random generator = new Random();
			Collections.shuffle(concatData, generator);
			for(int s=0; s<MAX_TRAIN_DATA && s<concatData.size(); s++){
				trainData.add(concatData.get(s));
			}
			System.out.println("Training round: "+i+" training set = "+trainData.size());
			
			//Random starts
			Vector<VectorClusterElement> starts = new Vector<VectorClusterElement>();
			for(int s=0; s<K; s++){
				int r = generator.nextInt(trainData.size());
				starts.add(trainData.get(r));
			}
			
					
			//Initialize clustering
			KMeansClustering<VectorClusterElement> kmc = new KMeansClustering<VectorClusterElement>(metric, crep, starts);
			kmc.setIterations(MAX_ITER);
			
			//Cluster!
			kmc.clusterElements(trainData, 0.01);			
			
			Double currSSD = kmc.sumOfSquaredDistance();
			if(currSSD<SSD){
				SSD = currSSD;
				bestClusterMeans = kmc.getClusterMeans();
			}System.out.println("SSD: "+currSSD);
		}
		//Print the clusters
		int ccount=0;
		for(VectorClusterElement c : bestClusterMeans){
			System.out.print("Cluster"+ccount);
			for(int d=0; d<c.dimension(); d++)
				System.out.print("\t"+c.getValue(d));
			System.out.println("");
			ccount++;
		}
	}
	public void loadData(){
		try{
		    File dFile = new File(inFileName);
		    if(!dFile.isFile()){System.err.println("Invalid data file name");System.exit(1);}
		    BufferedReader reader = new BufferedReader(new FileReader(dFile));
		    String line;// = reader.readLine(); //Ignore first line
		    while ((line = reader.readLine()) != null) {
		    	line = line.trim();
		    	String[] words = line.split("\\s+");
		    	
		    	if(!words[0].startsWith("#") && words.length==(TIMES+(TIMES*DIM))){
		    		int currTime = 0;
		    		//	Ignore first T fields (fillers)
		    		for(int i=TIMES; i<words.length; i+=DIM){
		    			Double [] currDat = new Double[DIM];
		    			for(int d=0; d<DIM; d++){
		    				currDat[d] = new Double(words[d+i]);
		    			}DefaultVectorClusterElement elem = new DefaultVectorClusterElement(currDat);
		    			allData.get(currTime).add(elem);
		    			currTime++;
		    		}
		    	}
		    }
		    for(int d=0; d<TIMES; d++){
		    	concatData.addAll(allData.get(d));
		    }
		} catch (FileNotFoundException e) {
		    // TODO Auto-generated catch block
		    e.printStackTrace();
		} catch (IOException e) {
		    // TODO Auto-generated catch block
		    e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
		String inFileName = Args.parseString(args,"in",null);
		int k = ap.hasKey("k") ? new Integer(ap.getKeyValue("k")).intValue():-1;
		int t = ap.hasKey("t") ? new Integer(ap.getKeyValue("t")).intValue():-1;
		int d = ap.hasKey("dim") ? new Integer(ap.getKeyValue("dim")).intValue():-1;
		
		if(inFileName !=null || k!=-1 || t!=-1 || d!=-1){
			ChromatinDataKMeans clusterer = new ChromatinDataKMeans(inFileName, k,t,d);
			clusterer.loadData();
			clusterer.cluster();
		}else{
			System.err.println("Options:");
			System.err.println("--in: infile");
			System.err.println("--k: cluster number");
			System.err.println("--t: number of times");
			System.err.println("--dim: dimensionality of data");
		}
	}
}

