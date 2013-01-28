package edu.psu.compbio.seqcode.gse.tools.chipchip;

import java.sql.*;
import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.chipchip.*;
import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.tools.chipchip.NormalizeExperiment;
import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.utils.database.*;

/**
 *  MixedChannels creates a new experiment where
 *   ip' = a * ip + b * wce;
 *   wce' = c * ip + d * wce;
 * This is useful if you want to artifically reduce signal levels to see what
 * it does to your analysis.
 */


public class MixedChannels extends NormalizeExperiment {

    private int oldexptid, newexptid;
    private Organism org;
    private Genome genome;
    private double a, b, c, d;


    public MixedChannels (Experiment oldexpt,
                          Experiment newexpt,
                          double a, double b, double c, double d) throws SQLException, NotFoundException {
        oldexptid = oldexpt.getDBID();
        newexptid = newexpt.getDBID();
        this.a = a;
        this.b = b;
        this.c = c;
        this.d = d;
    }

    public void doNorm() throws SQLException {        
        java.sql.Connection cxn = DatabaseFactory.getConnection("chipchip");
        Statement stmt = cxn.createStatement();
        String sql = String.format("insert into data(experiment, probe, channelone, channeltwo, channelratio) " +
                                   "(select %d as experiemnt, probe, channelone * %f + channeltwo * %f, channelone * %f + channeltwo * %f, channelratio from data where experiment = %d)",
                                   newexptid,
                                   a,b,c,d,
                                   oldexptid);        
        stmt.execute(sql);
        sql = String.format("update data set ratio = channelone / channeltwo, mor = channelone / channeltwo where experiment = %d",newexptid);
        stmt.execute(sql);
    }

}