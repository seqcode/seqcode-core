package edu.psu.compbio.seqcode.gse.clustering.affinitypropagation;

import edu.psu.compbio.seqcode.gse.clustering.Clusterable;

public class Pair {

	private Clusterable a, b;
	
	public Pair(Clusterable a, Clusterable b) {
		this.a = a;
		this.b = b;
	}
	
	public int hashCode() {
		return (a.name()+b.name()).hashCode();
	}
	
	public boolean symmetric() {
		return a.equals(b);
	}
	
	 public boolean equals(Object o) {
		 if (o instanceof Pair) {
			 return a.equals(((Pair)o).a) && b.equals(((Pair)o).b);
		 } else {
			 return false;
		 }
	 }
	
}
