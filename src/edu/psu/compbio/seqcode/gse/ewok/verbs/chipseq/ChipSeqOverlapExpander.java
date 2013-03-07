package edu.psu.compbio.seqcode.gse.ewok.verbs.chipseq;

import java.sql.*;
import java.util.*;
import java.io.IOException;

import edu.psu.compbio.seqcode.gse.datasets.general.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Mapper;
import edu.psu.compbio.seqcode.gse.utils.Closeable;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.database.*;

public class ChipSeqOverlapExpander implements Closeable, Mapper<StrandedRegion, RunningOverlapSum> {

    private SeqDataLoader loader;
	private LinkedList<SeqAlignment> alignments;
    private SeqLocator locator;
    private java.sql.Connection cxn;
    private PreparedStatement stmt;
    private int extension;
    private Genome lastGenome;

    public ChipSeqOverlapExpander(SeqLocator loc, int extension)  throws SQLException, IOException { 
        this.extension = extension;
        loader = new SeqDataLoader();
        alignments = null;
        stmt = null;
        locator = loc;
    }   
    private void getAligns(Genome genome) throws SQLException {
        if (alignments != null && genome.equals(lastGenome)) {
            return;
        }
        alignments = new LinkedList<SeqAlignment>();
        try {
            alignments.addAll(locator.loadAlignments(loader, genome));
        } catch (SQLException e) {
            e.printStackTrace(System.err);
        } catch (NotFoundException e) {
            e.printStackTrace();
        }
        cxn = DatabaseFactory.getConnection("seqdata");
        StringBuffer alignIDs = new StringBuffer();
        if (alignments.size() == 1) {
            alignIDs.append("alignment = " + alignments.get(0).getDBID());
        } else {

            alignIDs.append("alignment in (");
            for (int i = 0; i < alignments.size(); i++) {
                if (i == 0) {
                    alignIDs.append(alignments.get(i).getDBID());
                } else {
                    alignIDs.append("," + alignments.get(i).getDBID());
                }
            }
        }
        stmt = cxn.prepareStatement("select startpos, stoppos from chipseqhits where " + alignIDs.toString() +
                                    " and chromosome = ? and startpos > ? and stoppos < ? and strand = ?");
        stmt.setFetchSize(1000);
    }

	public RunningOverlapSum execute(StrandedRegion a) {
		try {
            Genome g = a.getGenome();
            getAligns(g);
            stmt.setInt(1, g.getChromID(a.getChrom()));
            stmt.setInt(2, a.getStart());
            stmt.setInt(3, a.getEnd());
            stmt.setString(4, a.getStrand() == '+' ? "+" : a.getStrand() == '-' ? "-" : " ");
            RunningOverlapSum sum = new RunningOverlapSum(g, a.getChrom());
            ResultSet rs = stmt.executeQuery();
            if (a.getStrand() == '+') {
                while (rs.next()) {
                    sum.addInterval(rs.getInt(1), rs.getInt(2) + extension);
                }
            } else if (a.getStrand() == '-') {
                while (rs.next()) {
                    sum.addInterval(rs.getInt(1) - extension, rs.getInt(2));
                }
            } 
            rs.close();
            return sum;
		} catch (SQLException e) {
			e.printStackTrace();
            throw new DatabaseException(e.toString(), e);
		}
	}
    public void close() {
		if(loader != null) { 
			loader.close();
            loader = null;
            if (alignments != null) { alignments.clear();}
            if (stmt != null) {
                try {
                    stmt.close();
                } catch (SQLException e) {
                    e.printStackTrace();
                } finally {
                    stmt = null;
                }
            }
            DatabaseFactory.freeConnection(cxn);
            cxn = null;
        }
	}
    public boolean isClosed() {
        return loader == null;
    }
}