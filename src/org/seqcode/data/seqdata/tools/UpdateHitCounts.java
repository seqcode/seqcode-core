package org.seqcode.data.seqdata.tools;

import java.io.IOException;
import java.sql.SQLException;

import org.seqcode.data.connections.DatabaseConnectionManager;
import org.seqcode.data.seqdata.SeqAlignment;
import org.seqcode.data.seqdata.SeqDataLoader;
import org.seqcode.data.seqdata.SeqDataModifier;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;
import org.seqcode.gseutils.NotFoundException;

/**
 * Updates the hit counts and weights stored in the db.
 * 
 * @author mahony
 *
 */
public class UpdateHitCounts {
	public static void main(String args[]) throws SQLException, NotFoundException, IOException {
		ArgParser ap = new ArgParser(args);
		if (ap.hasKey("id")) {
			Integer id = Args.parseInteger(args, "id", -1);

			SeqDataLoader loader = new SeqDataLoader();
			SeqDataModifier modifier = new SeqDataModifier(loader);
			SeqAlignment align = loader.loadAlignment(id);

			Integer oldsinglecount = align.getNumHits();
			Integer oldpaircount = align.getNumPairs();
			Float oldsingleweight = (float) align.getTotalWeight();
			Float oldpairweight = (float) align.getTotalPairWeight();

			Integer singlecount = Args.parseInteger(args, "singlecount", oldsinglecount);
			Integer singletype2count = Args.parseInteger(args, "singletype2count", oldsinglecount);
			Integer paircount = Args.parseInteger(args, "paircount", oldpaircount);
			Float singleweight = Args.parseFloat(args, "singleweight", oldsingleweight);
			Float singletype2weight = Args.parseFloat(args, "singletype2weight", oldsingleweight);
			Float pairweight = Args.parseFloat(args, "pairweight", oldpairweight);

			modifier.updateSeqAlignmentHitCounts(align, singlecount, singleweight, singletype2count, singletype2weight,
					paircount, pairweight);

			loader.close();
		} else {
			System.err.println("UpdateHitCounts:\n" + "\t--id <ReadDB ID>\n" + "\t--singlecount <hitcount, single>\n"
					+ "\t--singleweight <hitweight, single>\n" + "\t--singletype2count <hitcount, single type2>\n"
					+ "\t--singletype2weight <hitweight, single type2>\n" + "\t--paircount <hitcount, pairs>\n"
					+ "\t--pairweight <hitweight, pairs>\n" + "");
		}
	}
}
