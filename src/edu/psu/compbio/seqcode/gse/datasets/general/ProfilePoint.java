package edu.psu.compbio.seqcode.gse.datasets.general;

import edu.psu.compbio.seqcode.gse.datasets.species.Genome;

public class ProfilePoint extends Point {
	
	double profile;

	public ProfilePoint(Genome g, String c, int position, double profile) {
		super(g, c, position);
		this.profile = profile;
	}
	
	public double getProfile() {
		return profile;
	}

}
