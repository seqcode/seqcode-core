package org.seqcode.viz.metaprofile;

/**
 * Profile: meta-profiles are made from these.
 * 
 * @author: tdanford Date: Aug 12, 2008
 */
public interface Profile {
	public String getName();

	public double value(int i);

	public int length();

	public double max();

	public double min();

	public int getNumProfiles();

	public void setStranded(boolean s);

	public boolean isStranded();

	public BinningParameters getBinningParameters();

}