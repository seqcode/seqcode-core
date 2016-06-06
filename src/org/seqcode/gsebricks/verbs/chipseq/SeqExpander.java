package org.seqcode.gsebricks.verbs.chipseq;

import java.sql.SQLException;
import java.util.*;
import java.io.*;

import org.seqcode.data.seqdata.*;
import org.seqcode.genome.Genome;
import org.seqcode.genome.location.Region;
import org.seqcode.gsebricks.verbs.Expander;
import org.seqcode.utils.Closeable;
import org.seqcode.utils.NotFoundException;

/**
 * SeqExpander: legacy class for interacting with SeqAlignments
 * @author mahony
 * 
 */
public class SeqExpander implements Expander<Region, SeqHit>, Closeable {

    private SeqDataLoader loader;
    private Genome lastGenome;
    private LinkedList<SeqAlignment> alignments;
    private SeqLocator locator;
    private boolean loadR2;


    public SeqExpander(SeqLocator loc, boolean loadR2) throws SQLException, IOException {
        loader = new SeqDataLoader();
        locator = loc;
        alignments = null;
        lastGenome = null;
        this.loadR2 = loadR2;
    }
    private void getAligns(Genome genome) throws SQLException {
        if (alignments != null && genome.equals(lastGenome)) {
            return;
        }
        alignments = new LinkedList<SeqAlignment>();
        lastGenome = genome;
        try {
			alignments.addAll(loader.loadAlignments(locator, genome));
		} catch (NotFoundException e) {
			e.printStackTrace();
		}
        
    }


    public Iterator<SeqHit> execute(Region a) {
    	return this.getHits(a);
    }
    
    /**
     * Get single hits
     * @param a
     * @return
     */
    public Iterator<SeqHit> getHits(Region a) {
        try {
            getAligns(a.getGenome());
            Collection<SeqHit> hits = loader.loadByRegion(alignments, a, loadR2);
            return hits.iterator();
        }
        catch (Exception e) {
            e.printStackTrace();
            return new LinkedList<SeqHit>().iterator();
        }
    }

    /**
     * Get pairs
     * @param a
     * @return
     */
    public Iterator<SeqHitPair> getPairs(Region a) {
        try {
            getAligns(a.getGenome());
            Collection<SeqHitPair> pairs = loader.loadPairsByRegion(alignments, a);
            return pairs.iterator();
        }
        catch (Exception e) {
            e.printStackTrace();
            return new LinkedList<SeqHitPair>().iterator();
        }
    }
    
    /**
     * Get the single hit count
     * @param a
     * @return
     */
    public int getHitCount(Region a) {
        int hits = 0;
        try {
            getAligns(a.getGenome());
            hits = loader.countByRegion(alignments, a, loadR2);
            return hits;
        }
        catch (Exception e) {
            e.printStackTrace();
            return 0;
        }
    }

    
    public Collection<Genome> alignedGenomes() {
        LinkedList<Genome> genomes = new LinkedList<Genome>();
        if (alignments != null) {
            for (SeqAlignment align : alignments) {
                genomes.add(align.getGenome());
            }
        }
        return genomes;
    }


    public void close() {
        loader.close();
        loader = null;
        if (alignments != null) {
            alignments.clear();
        }
    }


    public boolean isClosed() {
        return loader == null;
    }
}
