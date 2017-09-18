package org.seqcode.ml.clustering.kmedoids;

import java.util.Collection;
import java.util.Vector;

import org.seqcode.ml.clustering.Cluster;
import org.seqcode.ml.clustering.ClusterRepresentative;
import org.seqcode.ml.clustering.ClusteringMethod;
import org.seqcode.ml.clustering.DefaultCluster;
import org.seqcode.ml.clustering.PairwiseElementMetric;
import org.seqcode.ml.clustering.vectorcluster.VectorClusterElement;
import java.util.Iterator;
import java.util.Set;


/** Modified from MMeansClustering
 * @author Naomi Yamada
 */
public class KMedoidsClustering<X> implements ClusteringMethod<X> {
	
    private PairwiseElementMetric<X> metric;
    private Vector<X> startMedoids;
    private int numClusters;
    private int iterations;
    private Vector<X> elmts;
    private Vector<DefaultCluster<X>> clusters;
    private Vector<X> clusterMedoids;
	
    public KMedoidsClustering(PairwiseElementMetric<X> m, Collection<X> starts) { 
        metric = m;
        numClusters = starts.size(); 
        clusters = new Vector<DefaultCluster<X>>();
        for(int c=0; c<numClusters; c++){clusters.add(new DefaultCluster<X>());}
        clusterMedoids = new Vector<X>(starts);
        startMedoids = new Vector<X>(starts);
        iterations = 10;
        elmts = new Vector<X>();
    }
	   
    public void setIterations(int i) { iterations = i; }

    public Collection<Cluster<X>> clusterElements(Collection<X> e) {return(clusterElements(e, 0));}
    public Collection<Cluster<X>> clusterElements(Collection<X> e, double convergenceDifference) {
        init(e);
        boolean converged=false;
        for(int i = 0; i < iterations && !converged; i++) {
        	Vector<X> oldClusterMedoids = (Vector<X>) clusterMedoids.clone();
        	
        	//K-means
        	assignToClusters();
        	getClusterMedoids();
            
            //Check convergence
        	if (oldClusterMedoids.equals(clusterMedoids))
        		converged=true;
        }
        return new Vector<Cluster<X>>(clusters);
    }
	
    private void assignToClusters() { 
        for(int i = 0; i < numClusters; i++) { clusters.get(i).clear(); }
        for(int k = 0; k < elmts.size(); k++) { 
            X e = elmts.get(k);
            int minCluster = -1;
            double minDist = 0.0;
            for(int i = 0; i < numClusters; i++) { 
                double clustDist = metric.evaluate(e, clusterMedoids.get(i));
                if(minCluster == -1 || clustDist < minDist) { 
                    minDist = clustDist;
                    minCluster = i;
                }
            }
            clusters.get(minCluster).addElement(e);
        }
    }
	
    public Vector<X> getClusterMedoids() { 
    	for (int i=0; i < numClusters; i++){
    		X newMedoid = null;
			double clustMinDist = 0.0;
			for (X m : clusters.get(i).getElements()){ // Candidate medoid
				double clustTotalDist =0.0;
				for (X e : clusters.get(i).getElements())
					if ( m != e)
						clustTotalDist += metric.evaluate(m, e);
				if (newMedoid == null || clustTotalDist < clustMinDist){
					clustMinDist = clustTotalDist;
					newMedoid = m;
				}
			}
			clusterMedoids.set(i, newMedoid);
    	}return(clusterMedoids);
    }
	
    private void init(Collection<X> e) {
        elmts = new Vector<X>(e);
        for(int i = 0; i < numClusters; i++) { 
            clusters.set(i, new DefaultCluster<X>());
            clusterMedoids.set(i, startMedoids.get(i));
        }
    }
    
    public double sumOfSquaredDistance(){
    	double totalDist =0;
    	for(int k = 0; k < elmts.size(); k++) { 
            X e = elmts.get(k);
            int minCluster = -1;
            double minDist = 0.0;
            for(int i = 0; i < numClusters; i++) { 
                double clustDist = metric.evaluate(e, clusterMedoids.get(i));
                if(minCluster == -1 || clustDist < minDist) { 
                    minDist = clustDist;
                    minCluster = i;
                }
            }totalDist += minDist*minDist;
        }
    	return(totalDist);
    }
    
    public double silhouette(){
    	double sh = 0.0;
    	for(int k=0; k< elmts.size(); k++){
    		X e = elmts.get(k);
    		double[] avgDist = new double[numClusters];
    		int minCluster = -1;
    		double minDist = 0.0;
    		for(int i = 0; i < numClusters; i++) { 
                double clustDist = metric.evaluate(e, clusterMedoids.get(i));
                if(minCluster == -1 || clustDist < minDist) { 
                    minDist = clustDist;
                    minCluster = i;
                }
                double currDist = 0.0;
                for(X ve : clusters.get(i).getElements()){
                	currDist = currDist + metric.evaluate(e, ve);
                }
                currDist = currDist/clusters.get(i).size();
                avgDist[i] = currDist;
            }
    		double neighborDistance = Double.MAX_VALUE;
    		for(int i=0; i<numClusters; i++){
    			if(i!=minCluster){ // Not its cluster
    				if(avgDist[i]<neighborDistance){
    					neighborDistance = avgDist[i];
    				}
    			}
    		}
    		
    		sh = sh + (neighborDistance - avgDist[minCluster])/Math.max(neighborDistance, avgDist[minCluster]);
    	}
    	sh = sh/elmts.size();
    	return sh;
    }
}

