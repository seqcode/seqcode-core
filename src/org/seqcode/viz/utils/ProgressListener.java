package org.seqcode.viz.utils;

public interface ProgressListener {

	public void registerEventMax(String key, int max);

	public void progressMade(ProgressEvent e);
}
