package edu.psu.compbio.seqcode.gse.tools.chipchip;

import java.sql.*;
import java.util.*;
import java.io.*;

import edu.psu.compbio.seqcode.gse.datasets.chipchip.*;
import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.species.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.*;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.utils.database.*;

/** 
 * Driver class and super class for producing new versions of a previously loaded experiment
 * with some normalization performed.
 *
 * Usage:
 * java edu.psu.compbio.seqcode.gse.tools.chipchip.NormalizeExperiment --species "species;genome" \
 *   --old "oldname;oldver;oldrep" --new "newname;newver;newrep" ...
 * Then you specify the type of normalization (one of the following):
 * --crosstalk --iptowce .1 --wcetoip .02
 * --control "controlname;controlversion;controlreplicate"
 * --mixed "a;b;c;d"
 * for crosstalk normalization, control experiment normalization (subtracts from ratios
 * by ratios in control expt).  Mixed generates a new experiment in which
 * the channels are a mix of the old channels:
 *   ip' = a * ip + b * wce;
 *   wce' = c * ip + d * wce;
 * 
 *
 */
public abstract class NormalizeExperiment {

    public static void main(String args[]) throws Exception {

        String on = null, ov = null, or = null, 
            nn = null, nv = null, nr = null, genomever = null;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--new")) {
                String pieces[] = args[++i].split(";");
                nn = pieces[0];
                nv = pieces[1];
                nr = pieces[2];
            }
            if (args[i].equals("--old")) {
                String pieces[] = args[++i].split(";");
                on = pieces[0];
                ov = pieces[1];
                or = pieces[2];
            }
        }
        ChipChipMetadataLoader loader = new ChipChipMetadataLoader();
        MetadataLoader general = new MetadataLoader();
        Experiment oldexpt = loader.loadExperiment(on,ov,or);
        FragDist fd = loader.loadFragDist(oldexpt.getFragDist());
        Pair<Organism,Genome> pair = Args.parseGenome(args);
        Organism org = pair.car();
        Genome genome = pair.cdr();

        Experiment newexpt = loader.getExperiment(nn,nv,nr,
                                                  org.getName(),
                                                  fd.getName(),
                                                  fd.getVersion(),
                                                  general.loadExptTarget(oldexpt.getFactorOne()).getName(),
                                                  general.loadExptTarget(oldexpt.getFactorTwo()).getName(),
                                                  general.loadCellLine(oldexpt.getCellsOne()).getName(),
                                                  general.loadCellLine(oldexpt.getCellsTwo()).getName(),
                                                  general.loadExptCondition(oldexpt.getConditionOne()).getName(),
                                                  general.loadExptCondition(oldexpt.getConditionTwo()).getName(),
                                                  true);
        int oldexptid = oldexpt.getDBID();
        int newexptid = newexpt.getDBID();
        System.err.println("old id is " + oldexptid + " new id is " +  newexptid);
        ResultSet rs;
        java.sql.Connection cxn = DatabaseFactory.getConnection("chipchip");
        PreparedStatement getToGenome = cxn.prepareStatement("select genome from exptToGenome where experiment = ?");
        PreparedStatement addToGenome = cxn.prepareStatement("insert into exptToGenome(genome,experiment) values(?,?)");
        getToGenome.setInt(1,oldexptid);
        rs = getToGenome.executeQuery();
        while (rs.next()) {
            addToGenome.setInt(1,rs.getInt(1));
            addToGenome.setInt(2,newexptid);
            addToGenome.execute();
        }
        rs.close();
        getToGenome.close();
        addToGenome.close();
        loader.mapExptToGenome(newexptid,genome.getDBID());
        NormalizeExperiment normalizer = parseNormalize(oldexpt, newexpt, genome, args);
        
        normalizer.doNorm();
    }
    public static NormalizeExperiment parseNormalize(Experiment oldexpt,
                                                     Experiment newexpt,
                                                     Genome genome,
                                                     String args[]) throws SQLException, NotFoundException, IOException {
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--crosstalk")) {
                double ipToWCE = Args.parseDouble(args,"iptowce",.1);
                double wceToIP = Args.parseDouble(args,"wcetoip",.1);
                CrossTalkNormalization ctn = new CrossTalkNormalization(oldexpt,
                                                                        newexpt,
                                                                        ipToWCE,
                                                                        wceToIP);
                return ctn;
            }          
            if (args[i].equals("--control")) {
                String control = Args.parseString(args,"control",null);
                String pieces[] = control.split(";");
                return new ControlExptNormalization(oldexpt,
                                                    newexpt,
                                                    pieces[0], pieces[1], pieces[2],
                                                    Args.parseFlags(args).contains("ratio"),
                                                    !Args.parseFlags(args).contains("subtractive"),
                                                    !Args.parseFlags(args).contains("usecontrolip"),
                                                    Args.parseFlags(args).contains("usecontrolboth"));
            }
            if (args[i].equals("--mixed")) {
                String pieces[] = args[++i].split(";");
                return new MixedChannels(oldexpt,
                                         newexpt,
                                         Double.parseDouble(pieces[0]),
                                         Double.parseDouble(pieces[1]),
                                         Double.parseDouble(pieces[2]),
                                         Double.parseDouble(pieces[3]));
            }

        }
        throw new RuntimeException("You didn't specify any normalizer");
    }

    public abstract void doNorm() throws SQLException ;
}
                                