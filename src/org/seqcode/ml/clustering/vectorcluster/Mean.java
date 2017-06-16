package org.seqcode.ml.clustering.vectorcluster;

import java.util.Set;

import org.seqcode.ml.clustering.Cluster;
import org.seqcode.ml.clustering.ClusterRepresentative;
import org.seqcode.ml.clustering.SingletonCluster;

public class Mean implements ClusterRepresentative<VectorClusterElement> {

	public Mean() {
	}

	public VectorClusterElement getRepresentative(Cluster<VectorClusterElement> c) {

		if (c instanceof SingletonCluster) {
			return (VectorClusterElement) ((SingletonCluster) c).getValue();
		}

		double[] array = null;
		Set<VectorClusterElement> elmts = c.getElements();
		if (elmts.size() == 0) {
			throw new IllegalArgumentException();
		}
		for (VectorClusterElement ce : elmts) {
			if (array == null) {
				array = new double[ce.dimension()];
				for (int i = 0; i < array.length; i++) {
					if (ce.isMissingValue(i)) {
						array[i] = 0.0;
					} else {
						array[i] = ce.getValue(i);
					}
				}
			} else {
				if (ce.dimension() != array.length) {
					throw new IllegalArgumentException(ce.toString());
				}

				for (int i = 0; i < array.length; i++) {
					if (!ce.isMissingValue(i)) {
						array[i] += ce.getValue(i);
					}
				}
			}
		}

		for (int i = 0; i < array.length; i++) {
			array[i] /= (double) elmts.size();
		}

		return new DefaultVectorClusterElement(array);
	}
}
