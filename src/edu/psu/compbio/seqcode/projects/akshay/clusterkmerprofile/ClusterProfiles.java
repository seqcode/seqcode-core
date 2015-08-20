package edu.psu.compbio.seqcode.projects.akshay.clusterkmerprofile;

import java.awt.Color;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.Vector;

import org.tc33.jheatchart.HeatChart;

import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;
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
	// Out tag name
	private String outtag="out";
	// The name of the directory that holds all the results
	// Assumes that the directory has already been created
	private File outbase;
	private HashMap<Integer,String> indToLocation;
	private HashMap<Integer,Double> intToLogitScore;
	private int kmerLen=4;
	
	// Heatmap options
	public final int W_MARGIN=80;
	public final int H_MARGIN=60;
	public final int W=800;
	public final int H=800;
	
	
	
	//Settors
	public void setProfiles(ArrayList<int[]> data){profiles=data;}
	public void setNumClusters(int nc){K=nc;}
	public void setOuttag(String otag){outtag=otag;}
	public void setOutDir(File odir){outbase = odir;}
	public void setProfileInds(HashMap<Integer,String> plfsinds){indToLocation=plfsinds;}
	public void setProfileScores(HashMap<Integer,Double> plfscore){intToLogitScore = plfscore;}
	public void setKmerModLen(int k){kmerLen = k;}
	
	
	/**
	 * The method that should be executed after initiating the class object
	 * @throws IOException 
	 */
	public void execute() throws IOException{
		Collection<Cluster<int[]>> clusters = method.clusterElements(profiles);
		Vector<int[]> clustermeans = ((KMeansClustering)method).getClusterMeans();
		
		//Print the clusters
		writeClusters(clustermeans);
		
		//Plot the clusters
		Mappable orderedClusters = reorderKmerProfileMaps(clusters);
		drawClusterHeatmap(orderedClusters);
	}
	
	
	
	// Slave methods
	private void writeClusters(Vector<int[]> clusMeans) throws IOException{
		File clusout = new File(outbase.getAbsolutePath()+File.separator+outtag+"_clusterAssignment.list");
		FileWriter ow = new FileWriter(clusout);
		BufferedWriter bw = new BufferedWriter(ow);
		for(int p=0; p<profiles.size(); p++){
			int memebership = getClusterAssignment(profiles.get(p),clusMeans);
			bw.write(indToLocation.get(p)+"\t"+Integer.toString(memebership)+"\t"+Double.toString(intToLogitScore.get(p))+"\n");
		}
		bw.close();
	}
	
	private int getClusterAssignment(int[] pfl, Vector<int[]> clusMeans){
		int minCluster = -1;
        double minDist = 0.0;
        for(int i = 0; i < clusMeans.size(); i++) { 
            double clustDist = comparator.evaluate(pfl, clusMeans.get(i));
            if(minCluster == -1 || clustDist < minDist) { 
                minDist = clustDist;
                minCluster = i;
            }
        }
        return minCluster;
	}
	
	
	private Mappable reorderKmerProfileMaps(Collection<Cluster<int[]>> clus){
		Mappable ret = null; 
		
		//Mappable features
		double[][] matrix;
		String[] rnames;
		String[] cnames;
		
		// Which colums to retain while drawing the heatmap
		boolean[] keepCol = new boolean[profiles.get(0).length];
		for(int i=0; i<keepCol.length; i++){
			keepCol[i] = false;
		}
		
		//Which clusters to these columns belong to...
		ArrayList<Integer> colCluster = new ArrayList<Integer>();
		for(int i=0; i<colCluster.size(); i++){
			colCluster.add(K+1);
		}
		
		
		int clusID=1;
		for(Cluster<int[]> c : clus){
			for(int i=0; i<keepCol.length; i++){
				for(int[] elems : c.getElements()){
					if(elems[i] > 0){
						keepCol[i] = true;
						if(clusID < colCluster.get(i)){
							colCluster.set(i, clusID);
						}
					}
				}
			}
			clusID++;
		}
		
		// Now reorder 
		ArrayIndexComparator comp = new ArrayIndexComparator(colCluster);
		Integer[] indexes = comp.createIndexArray();
		Arrays.sort(indexes, comp);
		
		int sparseLenght = 0;
		for(int i=0;i<indexes.length; i++){
			if(keepCol[i])
				sparseLenght++;
			else
				break;
		}
		
		matrix = new double[profiles.size()][sparseLenght];
		rnames = new String[profiles.size()];
		cnames = new String[sparseLenght];
		
		//fill the cnames
		for(int j=0; j<sparseLenght; j++){
			cnames[j] = RegionFileUtilities.int2seq(indexes[j], kmerLen);
		}
		
		int rowInd = 0;
		for(Cluster<int[]> c : clus){
			for(int[] elems : c.getElements()){
				for(int j=0; j<sparseLenght; j++){
					matrix[rowInd][j] = elems[indexes[j]];
				}
				rnames[rowInd] = indToLocation.get(rowInd);
				rowInd++;
			}
		}
		
		ret = new Mappable(matrix, rnames, cnames);
		return ret;
	}
	
	
	
	public void drawClusterHeatmap(Mappable plotMat) throws IOException{
		double[][] matrix = plotMat.matrix;
		HeatChart map = new HeatChart(matrix);
		map.setHighValueColour(new Color(10));
		map.setLowValueColour(new Color(20));
		map.setChartMargin(100);
		map.setAxisLabelsFont(new Font("Ariel",Font.PLAIN,55));
		map.setXValues(plotMat.rownmanes);
		map.setYValues(plotMat.colnames);
		File f = new File(outbase.getAbsolutePath()+File.separator+outtag+"_clusters.png");
		map.saveToFile(f);
	}
	
	
	
	/**
	 * Constructor that sets up the k-means object
	 * @param itrs
	 * @param k
	 * @param pfls
	 * @param otag
	 */
	public ClusterProfiles(int itrs, int k, ArrayList<int[]> pfls, HashMap<Integer,String> pflsIndsMap, int kmerL, HashMap<Integer,Double> pflscores, String otag, File odir) {
		setProfiles(pfls);
		setNumClusters(k);
		setOuttag(otag);
		setOutDir(odir);
		setProfileScores(pflscores);
		setKmerModLen(kmerL);
		
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
	
	public class Mappable{
		public double[][] matrix;
		public String[] rownmanes;
		public String[] colnames;
		public Mappable(double[][] m, String[] rnames, String[] cnames) {
			matrix = m;
			rownmanes = rnames;
			colnames = cnames;
		}
	}
	
	public class ArrayIndexComparator implements Comparator<Integer>{
		ArrayList<Integer> list;
		
		public ArrayIndexComparator(ArrayList<Integer> ls) {
			list = ls;
		}
		
		public Integer[] createIndexArray(){
			Integer[] indexes = new Integer[list.size()];
			for(int i=0; i<indexes.length; i++){
				indexes[i] = i;
			}
			return indexes;
		}

		@Override
		public int compare(Integer o1, Integer o2) {
			return list.get(o1).compareTo(list.get(o2));
		}
		
	}

}
