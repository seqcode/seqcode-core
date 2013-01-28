package edu.psu.compbio.seqcode.gse.viz.preferences;

public interface PreferencesListener {

	public void preferencesUpdated(PreferencesEvent evt);
	public void preferencesUpdateCanceled(PreferencesEvent evt);
}
