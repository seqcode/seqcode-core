package org.seqcode.projects.seqview.components;

import java.util.Comparator;

import org.seqcode.genome.location.Gene;

class GeneorRegionComparator implements Comparator<Object> {

	public int compare(Object a, Object b) {
		if (a instanceof Gene && b instanceof Gene) {
			Gene x = (Gene) a;
			Gene y = (Gene) b;
			String k1 = x.getName() != null ? x.getName() : x.toString();
			String k2 = y.getName() != null ? y.getName() : y.toString();
			return k1.compareTo(k2);
		} else {
			if (a == null) {
				return 1;
			} else {
				if (b == null) {
					return -1;
				} else {
					return a.toString().compareTo(b.toString());
				}
			}
		}
	}

}
