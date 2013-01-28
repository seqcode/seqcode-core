/*
 * Created on Aug 24, 2005
 */
package edu.psu.compbio.seqcode.gse.utils.preferences;

import edu.psu.compbio.seqcode.gse.utils.Factory;
import edu.psu.compbio.seqcode.gse.utils.Saveable;

import java.io.*;

/**
 * @author tdanford
 */
public interface Preferences<X> extends Factory<X>, Saveable {
    public String getName();
    public PreferencesPanel createPanel();
    public void saveFromPanel(PreferencesPanel pp);
}
