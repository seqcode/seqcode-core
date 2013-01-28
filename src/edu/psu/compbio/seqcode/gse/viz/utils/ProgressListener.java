package edu.psu.compbio.seqcode.gse.viz.utils;

public interface ProgressListener {

	public void registerEventMax(String key, int max);
	public void progressMade(ProgressEvent e);
}
