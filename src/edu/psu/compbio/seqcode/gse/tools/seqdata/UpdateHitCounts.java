package edu.psu.compbio.seqcode.gse.tools.seqdata;

import java.io.IOException;
import java.sql.PreparedStatement;
import java.sql.SQLException;

import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqAlignment;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqDataLoader;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.database.DatabaseFactory;

/**
 * Updates the hit counts and weights stored in the db.
 * 
 * @author mahony
 *
 */
public class UpdateHitCounts {
	public static void main(String args[]) throws SQLException, NotFoundException, IOException {
		ArgParser ap = new ArgParser(args);
        if(ap.hasKey("id")) { 
	        java.sql.Connection cxn = DatabaseFactory.getConnection("seqdata");
	        cxn.setAutoCommit(false);
	        Integer id = Args.parseInteger(args,"id", -1);
	        
	        SeqDataLoader loader = new SeqDataLoader();
	        SeqAlignment align = loader.loadAlignment(id);
	        
	        Integer oldsinglecount = align.getNumHits();
	        Integer oldpaircount = align.getNumPairs();
	        Float oldsingleweight = (float)align.getTotalWeight();
	        Float oldpairweight = (float)align.getTotalPairWeight();
	        
	        Integer singlecount = Args.parseInteger(args,"singlecount", oldsinglecount);
	        Integer paircount = Args.parseInteger(args,"paircount", oldpaircount);
	        Float singleweight = Args.parseFloat(args,"singleweight", oldsingleweight);
	        Float pairweight = Args.parseFloat(args,"pairweight", oldpairweight);
	        
	        
	        PreparedStatement update = SeqAlignment.createUpdateHitsAndWeights(cxn);
	        System.err.println("Updating counts for alignment: "+id+" ("+align.getName()+")");
	        System.err.println("\tnumhits="+singlecount);
	        System.err.println("\ttotalweight="+singleweight);
	        System.err.println("\tnumpairs="+paircount);
	        System.err.println("\ttotalpairweight="+pairweight);
	        update.setInt(1, singlecount);
	        update.setFloat(2, singleweight);
	        update.setInt(3, paircount);
	        update.setFloat(4, pairweight);
	        update.setInt(5, id);
	        update.execute();
	        update.close();
	        
	        loader.close();
	        cxn.commit();
	        cxn.close();
        }else{
        	System.err.println("UpdateHitCounts:\n" +
        			"\t--id <ReadDB ID>\n" +
        			"\t--singlecount <hitcount, single>\n" +
        			"\t--singleweight <hitweight, single>\n" +
        			"\t--paircount <hitcount, pairs>\n" +
        			"\t--pairweight <hitweight, pairs>\n" +
        			"");
        }
	}
}
