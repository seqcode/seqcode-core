package edu.psu.compbio.seqcode.gse.seqview.model;

import java.io.IOException;
import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.projects.readdb.Client;
import edu.psu.compbio.seqcode.gse.projects.readdb.ClientException;
import edu.psu.compbio.seqcode.gse.projects.readdb.PairedHit;
import edu.psu.compbio.seqcode.gse.projects.readdb.PairedHitLeftComparator;

public class SeqPairedEndModel extends SeqViewModel implements RegionModel, Runnable {

    private Client client;
    private Set<SeqAlignment> alignments;
    private Set<String> ids;
    private Region region;
    private boolean newinput;
    private List<PairedHit> results, otherchrom;
    private Comparator<PairedHit> comparator;
    private SeqPairedEndModelProperties props;
    
    public SeqPairedEndModel (Collection<SeqAlignment> alignments) throws IOException, ClientException{
        client = new Client();
        comparator = new PairedHitLeftComparator();
        this.alignments = new HashSet<SeqAlignment>();
        this.alignments.addAll(alignments);
        ids = new HashSet<String>();
        for (SeqAlignment a : alignments) {
            ids.add(Integer.toString(a.getDBID()));
        }
        results = null;
        otherchrom = null;
        props = new SeqPairedEndModelProperties();
    }
    public SeqPairedEndModelProperties getProperties() {return props;}

    public void clearValues() {
        results = null;
        otherchrom = null;
    }
    public Region getRegion() {return region;}
    public void setRegion(Region r) {
        if (newinput == false) {
            if (!r.equals(region)) {
                region = r;
                newinput = true;
            } else {
                notifyListeners();
            }
        }
    }
    public void resetRegion(Region r) {
        if (newinput == false) {
        	region = r;
        	newinput = true;
        }
    }
    public List<PairedHit> getResults () {return results;}
    public List<PairedHit> getOtherChromResults() {return otherchrom;}
    public boolean connectionOpen(){return client.connectionAlive();}
    public void reconnect(){client.reConnect();}
    
    public boolean isReady() {return !newinput;}
    public List<PairedHit> dedup(List<PairedHit> hits) {
        ArrayList<PairedHit> deduped = new ArrayList<PairedHit>();
        deduped.add(hits.get(0));
        for (int i = 1; i < hits.size(); i++) {
            PairedHit a = hits.get(i);
            PairedHit b = deduped.get(deduped.size() - 1);
            if (a.leftPos != b.leftPos ||
                a.rightPos != b.rightPos ||
                a.leftStrand != b.leftStrand ||
                a.rightStrand != b.rightStrand) {
                deduped.add(a);
                            }
        }
        return deduped;
    }
    public synchronized void run() {
        while(keepRunning()) {
            try {
                if (!newinput) {
                    wait();
                }
            } catch (InterruptedException ex) { }
            if (newinput) {
                try {
                    results = new ArrayList<PairedHit>();
                    otherchrom = new ArrayList<PairedHit>();
                    double mindist = getProperties().MinimumDistance;
                    if (mindist < 1) {
                        mindist = mindist * region.getWidth();
                    }
                    for (String alignid : ids) {
                        List<PairedHit> r = client.getPairedHits(alignid,
                                                                 region.getGenome().getChromID(region.getChrom()),
                                                                 true,
                                                                 region.getStart(),
                                                                 region.getEnd(),
                                                                 null,
                                                                 null);
                        for (PairedHit h : r) {
                            if (h.leftChrom == h.rightChrom) { 
                                if (h.rightPos >= region.getStart() &&
                                    h.rightPos <= region.getEnd() && 
                                    Math.abs(h.leftPos - h.rightPos) > mindist) {
                                    results.add(h);
                                }
                            } else {
                                otherchrom.add(h);
                            }
                        }
                    }
                    Collections.sort(results, comparator);
                    Collections.sort(otherchrom, comparator);
                    
                    if (getProperties().DeDuplicateByPosition && results.size() > 0) {
                        results = dedup(results);
                        if (otherchrom.size() > 0) {
                        	otherchrom = dedup(otherchrom);
                        }
                    }
                    if (getProperties().LeftAlwaysLesser) {
                        for (PairedHit h : results) {
                            if (h.leftPos > h.rightPos) {
                                h.flipSides();
                            }
                        }
                        Collections.sort(results, comparator);
                    }
                } catch (Exception ex) {
                    ex.printStackTrace();
                    // assign empty output.  This is useful because Client
                    // throws an exception for non-existant chromosomes, such
                    // as those for which there were no alignment results
                    results = new ArrayList<PairedHit>();
                }
                newinput = false;
                notifyListeners();
            }
        }
        client.close();
    }
}
