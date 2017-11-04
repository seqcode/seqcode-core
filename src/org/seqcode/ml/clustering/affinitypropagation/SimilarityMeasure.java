package org.seqcode.ml.clustering.affinitypropagation;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import org.seqcode.ml.clustering.Clusterable;
import org.seqcode.ml.clustering.ClusterablePair;
import org.seqcode.ml.clustering.PairwiseElementMetric;


/**
 * 
 * @author reeder
 *
 */
public abstract class SimilarityMeasure<X extends Clusterable> implements PairwiseElementMetric<X> {

    protected boolean noiseAdded;
	public static final double NEGINF = -Double.MAX_VALUE;
	protected int[] assgn;
	protected int[] exidx;

    public abstract void addNoise();
    public abstract double get(int idx1, int idx2);
    public abstract int size();
    public abstract String getName(int idx);
    public abstract boolean exists(ClusterablePair p);
    public abstract double evaluate(ClusterablePair p);
    
    /**
     * 
     * @param assgn an array of the index of the exemplar for each point
     */
    public void putAssignments(int[] assgn) {
    	this.assgn = assgn;
    }
    
    /**
     * 
     * @param ex an array of the indices of the exemplars
     */
    public void putExemplars(int[] exidx) {
    	this.exidx = exidx;
    }
    
    public void printAssignments(PrintStream p) {
    	for (int i=0; i<this.size(); i++) {
    		p.println(getName(i)+": "+assgn[i]);
    	}
    }
    
    public void printExemplars(PrintStream p) {
    	p.println("Number of clusters: "+exidx.length);
    	for (int i=0; i<exidx.length; i++) {
    		p.println("Cluster "+i+": "+getName(exidx[i]));
    	}
    }
    
    public void printClusterCenterIndices(PrintStream p) {
    	for (int i=0; i<this.size(); i++) {
    		p.println(exidx[assgn[i]]+1);
    	}
    }

    public List<APExemplar> getExemplars(){
    	List<APExemplar> exs = new ArrayList<APExemplar>();
    	for (int i=0; i<exidx.length; i++) 
    		exs.add(new APExemplar(i));
    	return exs;
    }
    
    public List<APAssignment> getAssignments(){
    	List<APAssignment> assigns = new ArrayList<APAssignment>();
    	for (int i=0; i<this.size(); i++)
    		assigns.add(new APAssignment(i));
    	return assigns;
    }

    
    public class APExemplar{
    	public int index;
    	public String name;
    	public APExemplar(int i){
    		index = exidx[i];
    		name = getName(exidx[i]);
    	}
    }
    
    public class APAssignment{
    	public int index;
    	public String name;
    	public APExemplar exemplar;
    	public APAssignment(int i){
    		index = i;
    		name = getName(i);
    		exemplar = new APExemplar(assgn[i]);
    	}
    }
}