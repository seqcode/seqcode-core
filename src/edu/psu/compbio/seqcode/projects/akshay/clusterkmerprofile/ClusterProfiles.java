package edu.psu.compbio.seqcode.projects.akshay.clusterkmerprofile;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.machinelearning.clustering.Cluster;
import edu.psu.compbio.seqcode.machinelearning.clustering.ClusterRepresentative;
import edu.psu.compbio.seqcode.machinelearning.clustering.ClusteringMethod;
import edu.psu.compbio.seqcode.machinelearning.clustering.kmeans.KMeansClustering;

public class ClusterProfiles {
	private KmerProfileEucDistComparator comparator;
	private ClusteringMethod<int[]> method;
	private KmerProfileAvgDistRep rep;
	private int K;
	private ArrayList<int[]> profiles;
	
	
	//Settors
	public void setProfiles(ArrayList<int[]> data){profiles=data;}
	public void setNumClusters(int nc){K=nc;}
	
	//public static void main(String args[]){
		
		
	//	int K = Args.parseInteger(args, "numClusters", 3);
	//	int 
		
		
		
	//}
	
	public void cluster(){
		Collection<Cluster<int[]>> clusters = method.clusterElements(profiles);
		
	}
	
	
	public void drawClustersHeatmap(Collection<Cluster<int[]>> clus, String PicFilename){
		
		
		
	}
	
	
	
	
	public ClusterProfiles(int itrs, int k, ArrayList<int[]> pfls) {
		setProfiles(pfls);
		setNumClusters(k);
		comparator = new KmerProfileEucDistComparator();
		rep = new KmerProfileAvgDistRep(comparator);
		Random generator = new Random();
		
		List<int[]> starts = new ArrayList<int[]>();
		for(int s=0; s<K; s++){
			int r = generator.nextInt(profiles.size());
			starts.add(profiles.get(r));
		}
		
		method = new KMeansClustering<int[]>(comparator,rep,starts);
		((KMeansClustering)method).setIterations(itrs);
	}

}
