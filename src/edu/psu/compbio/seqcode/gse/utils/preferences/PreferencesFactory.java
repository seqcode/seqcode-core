/*
 * Created on Aug 27, 2005
 */
package edu.psu.compbio.seqcode.gse.utils.preferences;

import edu.psu.compbio.seqcode.gse.utils.Factory;

/**
 * @author tdanford
 */
public interface PreferencesFactory<X> extends Factory<Preferences<X>> {
    public Preferences<X> createObject(X v);
}
