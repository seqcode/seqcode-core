/*
 * Created on Apr 2, 2006
 */
package edu.psu.compbio.seqcode.gse.ewok.utils;

import java.sql.SQLException;
import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.binding.BindingExtent;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.locators.*;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.binding.PeakCaller;
import edu.psu.compbio.seqcode.gse.ewok.verbs.binding.RegionProber;
import edu.psu.compbio.seqcode.gse.utils.database.UnknownRoleException;

/**
 * @author tdanford
 */
public interface EwokBinding<Locator extends ExptLocator, P extends Probe> {
    public RegionProber<P> getProber(Locator loc);
    public Expander<Region,BindingExtent> getPeakCaller(Locator loc, double[] params);
    public Collection<Locator> getAllLocators() throws SQLException, UnknownRoleException;
    public EwokBase getBase();
}