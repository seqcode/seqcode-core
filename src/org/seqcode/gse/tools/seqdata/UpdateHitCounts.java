package org.seqcode.gse.tools.seqdata;

import java.io.IOException;
import java.sql.SQLException;

import org.seqcode.gse.datasets.seqdata.SeqAlignment;
import org.seqcode.gse.datasets.seqdata.SeqDataLoader;
import org.seqcode.gse.datasets.seqdata.SeqDataModifier;
import org.seqcode.gse.tools.utils.Args;
import org.seqcode.gse.utils.ArgParser;
import org.seqcode.gse.utils.NotFoundException;
import org.seqcode.gse.utils.database.DatabaseConnectionManager;


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
	        java.sql.Connection cxn = DatabaseConnectionManager.getConnection("seqdata");
	        cxn.setAutoCommit(true);
	        Integer id = Args.parseInteger(args,"id", -1);
	        
	        SeqDataLoader loader = new SeqDataLoader();
	        SeqDataModifier modifier = new SeqDataModifier(loader);
	        SeqAlignment align = loader.loadAlignment(id);
	        
	        Integer oldsinglecount = align.getNumHits();
	        Integer oldpaircount = align.getNumPairs();
	        Float oldsingleweight = (float)align.getTotalWeight();
	        Float oldpairweight = (float)align.getTotalPairWeight();
	        
	        Integer singlecount = Args.parseInteger(args,"singlecount", oldsinglecount);
	        Integer singletype2count = Args.parseInteger(args,"singletype2count", oldsinglecount);
	        Integer paircount = Args.parseInteger(args,"paircount", oldpaircount);
	        Float singleweight = Args.parseFloat(args,"singleweight", oldsingleweight);
	        Float singletype2weight = Args.parseFloat(args,"singletype2weight", oldsingleweight);
	        Float pairweight = Args.parseFloat(args,"pairweight", oldpairweight);
	        
	        modifier.updateSeqAlignmentHitCounts(align, singlecount, singleweight, singletype2count, singletype2weight, paircount, pairweight);
		     
	        loader.close();
	        cxn.close();
        }else{
        	System.err.println("UpdateHitCounts:\n" +
        			"\t--id <ReadDB ID>\n" +
        			"\t--singlecount <hitcount, single>\n" +
        			"\t--singleweight <hitweight, single>\n" +
        			"\t--singletype2count <hitcount, single type2>\n" +
        			"\t--singletype2weight <hitweight, single type2>\n" +
        			"\t--paircount <hitcount, pairs>\n" +
        			"\t--pairweight <hitweight, pairs>\n" +
        			"");
        }
	}
}
