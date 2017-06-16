package org.seqcode.projects.seqview.components;

import org.seqcode.genome.location.Region;

public interface RegionList {
	public void addRegion(Region r);

	public int regionListSize();

	public Region regionAt(int i);
}
