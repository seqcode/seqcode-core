package org.seqcode.projects.seqview.model;

import org.seqcode.genome.location.Region;

public class RegionMultipleExpanderModel<OUT> extends MultipleExpanderModel<Region, OUT> implements RegionModel {
	private Region region;

	public RegionMultipleExpanderModel() {
	}

	public void setRegion(Region r) throws NullPointerException {
		if (r == null) {
			throw new NullPointerException("Region can't be null");
		}
		region = r;
		setInput(r);
	}

	public void resetRegion(Region r) throws NullPointerException {
		if (r == null) {
			throw new NullPointerException("Region can't be null");
		}
		region = r;
		setInput(r);
	}

	public Region getRegion() {
		return region;
	}
}
