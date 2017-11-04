package org.seqcode.ml.clustering.affinitypropagation;

import java.util.HashMap;
import java.util.List;
import java.util.Vector;

import org.seqcode.ml.clustering.Clusterable;
import org.seqcode.ml.clustering.ClusterablePair;
import org.seqcode.ml.clustering.SimpleClusterable;


/**
 * @author mahony
 *
 */
public class MatrixSimilarityMeasure<X extends Clusterable> extends SimilarityMeasure<X> {

	Vector<Vector<Double>> simValues = new Vector<Vector<Double>>();
	Vector<Clusterable> objects = new Vector<Clusterable>();
	HashMap<ClusterablePair, Double> valuemap;
	double prefvalue;
	
	//names index matrix. If matrix is NxN, names should be N
	public MatrixSimilarityMeasure(List<String> names, double[][] matrix, double prefvalue) {
		this.prefvalue = prefvalue;
		valuemap = new HashMap<ClusterablePair, Double>();
		
		for(Integer i=0; i<names.size(); i++){
			SimpleClusterable object1 = new SimpleClusterable(names.get(i));
			if (!objects.contains(object1))
				objects.add(object1);
			
			for(Integer j=0; j<names.size(); j++){
				SimpleClusterable object2 = new SimpleClusterable(names.get(j));
				if (!objects.contains(object2))
					objects.add(object2);
				
				valuemap.put(new ClusterablePair(object1, object2), matrix[i][j]);
			}
		}
	}
	
	
	
	@Override
	public void addNoise() {
		noiseAdded = true;

	}

	@Override
	public double get(int idx0, int idx1) {
		System.out.println("get shouldn't be called");
		if (simValues.get(idx0).get(idx1)==null) {
			return NEGINF;
		} else {
			return simValues.get(idx0).get(idx1);
		}
	}

	@Override
	public int size() {
		return objects.size();
	}
	
	public String getName(int idx) {
		return objects.get(idx).name();
	}

	/*public double evaluate(Clusterable e1, Clusterable e2) {
		if (valuemap.containsKey(e1.name()+e2.name()));
		return 0;
	}*/

	public double evaluate(X e1, X e2) {
		if (e1.name().equals(e2.name())) {
			return prefvalue;
		} else if (valuemap.containsKey(e1.name()+"SEP"+e2.name())) {
			return valuemap.get(e1.name()+"SEP"+e2.name());
		} else if (valuemap.containsKey(e1.name()+"SEP"+e1.name())) {
			return valuemap.get(e1.name()+"SEP"+e2.name());
		} else {
			return NEGINF;
		}
	}
	
	public double evaluate(ClusterablePair p) {
		if (p.symmetric()) {
			return prefvalue;
		} else if (valuemap.containsKey(p)) {
			return valuemap.get(p);
		} else {
			return NEGINF;
		}
	}
	
	public boolean exists(X e1, X e2) {
		return (e1.name().equals(e2.name())) || (valuemap.containsKey(e1.name()+"SEP"+e2.name())) || 
			(valuemap.containsKey(e1.name()+"SEP"+e1.name()));
	}
	
	public boolean exists(ClusterablePair p) {
		return p.symmetric() || valuemap.containsKey(p);
	}
	
	public Vector<Clusterable> objects() {
		return objects;
	}
	
}
