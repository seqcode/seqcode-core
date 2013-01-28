package edu.psu.compbio.seqcode.gse.ewok.verbs;

import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.*;
import edu.psu.compbio.seqcode.gse.datasets.binding.BindingEvent;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.ChipChipMLE;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.ewok.*;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;
import edu.psu.compbio.seqcode.gse.utils.*;

public class MLEBindingGenerator implements Expander<Region,BindingEvent> {

    private ChipChipMLE mle;
    private double sizethresh, probthresh;

    public MLEBindingGenerator (ChipChipMLE mle, double probthresh, double sizethresh) {
        this.mle = mle;
        this.sizethresh = sizethresh;
        this.probthresh = probthresh;
    }

    public Iterator<BindingEvent> execute(Region r) {
        ArrayList results = new ArrayList<BindingEvent>();
        try {
            long t1 = System.currentTimeMillis();
            mle.window(r.getChrom(),
                       r.getStart(),
                       r.getEnd(),
                       sizethresh, probthresh);
            long t2 = System.currentTimeMillis();
            //            System.err.println("MLEBindingGenerator window() took " + (t2 - t1));
            int count = mle.getCount();
            for (int i = 0; i < count; i++) {
                results.add(new BindingEvent(r.getGenome(),
                                             r.getChrom(),
                                             mle.getPos(i),
                                             mle.getPos(i),
                                             mle.getSize(i),
                                             mle.getConf(i),
                                             "MLE"));
            }
        } catch (NotFoundException ex) {
            throw new RuntimeException("Couldn't window()" +ex,ex);
        }
        return results.iterator();
    }
    public void setPeaks(boolean peaks) {}
}
