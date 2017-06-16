package org.seqcode.ml.clustering;

public class ClusterablePair {

	private Clusterable a, b;

	public ClusterablePair(Clusterable a, Clusterable b) {
		this.a = a;
		this.b = b;
	}

	public int hashCode() {
		return (a.name() + b.name()).hashCode();
	}

	public boolean symmetric() {
		return a.equals(b);
	}

	public boolean equals(Object o) {
		if (o instanceof ClusterablePair) {
			return a.equals(((ClusterablePair) o).a) && b.equals(((ClusterablePair) o).b);
		} else {
			return false;
		}
	}

}
