/*
 * Created on Oct 3, 2007
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package edu.psu.compbio.seqcode.gse.ewok.verbs;

import java.io.IOException;
import java.sql.SQLException;
import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.chipseq.*;
import edu.psu.compbio.seqcode.gse.datasets.general.*;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.verbs.chipseq.ChipSeqExpander;
import edu.psu.compbio.seqcode.gse.utils.Closeable;

public class H2AZMidpointExpander implements Expander<Region,Point>, Closeable {
    
    private ChipSeqExpander exp;
    
    public H2AZMidpointExpander(String exptName, String alignName) throws SQLException { 
        ChipSeqLocator loc = new ChipSeqLocator(exptName, alignName);
        try {
			exp = new ChipSeqExpander(loc);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }

    public Iterator<Point> execute(Region r) {
        Iterator<ChipSeqHit> hits = exp.execute(r);
        LinkedList<Point> pts = new LinkedList<Point>();
        Genome g = r.getGenome();
        
        while(hits.hasNext()) { 
            ChipSeqHit hit = hits.next();
            int mdp = 0;
            if(hit.getStrand() == '+') { 
                mdp = hit.getStart() + 73;
            } else { 
                mdp = hit.getEnd() - 73;
            }
            pts.addLast(new Point(g, r.getChrom(), mdp));
        }
        return pts.iterator();
    }

    public void close() {
        exp.close();
        exp = null;
    }

    public boolean isClosed() {
        return exp==null;
    }

}
