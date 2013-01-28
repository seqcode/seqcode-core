/*
 * Created on Oct 17, 2005
 */
package edu.psu.compbio.seqcode.gse.datasets.alignments;

import java.sql.SQLException;

import edu.psu.compbio.seqcode.gse.utils.database.UnknownRoleException;

/**
 * @author tdanford
 */
public interface PairwiseDistributionLearner {
    public int getWidth();
    public void addDistribution(String other);
    public PairwiseLetterDistribution[] getDistribution(String other);
    public void learnFromAlignment(Alignment a, int ungappedStart)
            throws SQLException, UnknownRoleException;
}