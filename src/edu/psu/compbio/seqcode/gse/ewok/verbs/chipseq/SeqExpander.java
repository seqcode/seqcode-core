package edu.psu.compbio.seqcode.gse.ewok.verbs.chipseq;

import java.sql.SQLException;
import java.util.*;
import java.io.*;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Expander;
import edu.psu.compbio.seqcode.gse.utils.Closeable;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

public class SeqExpander implements Expander<Region, SeqHit>, Closeable {

    private SeqDataLoader loader;
    private Genome lastGenome;
    private LinkedList<SeqAlignment> alignments;
    private SeqLocator locator;
    private boolean closeLoader;


    public SeqExpander(SeqLocator loc) throws SQLException, IOException {
        loader = new SeqDataLoader();
        closeLoader = true;
        locator = loc;
        alignments = null;
        lastGenome = null;
    }
    private void getAligns(Genome genome) throws SQLException {
        if (alignments != null && genome.equals(lastGenome)) {
            return;
        }
        alignments = new LinkedList<SeqAlignment>();
        lastGenome = genome;
        try {
            alignments.addAll(locator.loadAlignments(loader, genome));
        } catch (SQLException e) {
            e.printStackTrace(System.err);
        } catch (NotFoundException e) {
            e.printStackTrace();
        }
    }

    public SeqExpander(SeqDataLoader l, SeqAlignment a, boolean closeLoader) {
        loader = l;
        alignments = new LinkedList<SeqAlignment>();
        alignments.add(a);
        this.closeLoader = closeLoader;
    }


    public Iterator<SeqHit> execute(Region a) {
        try {
            getAligns(a.getGenome());
            Collection<SeqHit> hits = loader.loadByRegion(alignments, a);
            return hits.iterator();
        }
        catch (Exception e) {
            e.printStackTrace();
            return new LinkedList<SeqHit>().iterator();
        }
    }


    public int getHitCount(Region a) {
        int hits = 0;
        try {
            getAligns(a.getGenome());
            hits = loader.countByRegion(alignments, a);
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
        if (closeLoader) {
            loader.close();
        }
        loader = null;
        if (alignments != null) {
            alignments.clear();
        }
    }


    public boolean isClosed() {
        return loader == null;
    }
}
