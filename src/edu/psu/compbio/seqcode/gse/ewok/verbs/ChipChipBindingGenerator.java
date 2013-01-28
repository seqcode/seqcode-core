package edu.psu.compbio.seqcode.gse.ewok.verbs;

import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.binding.BindingExtent;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.ChipChipData;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.ewok.*;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;
import edu.psu.compbio.seqcode.gse.utils.*;

public class ChipChipBindingGenerator implements Expander<Region,BindingExtent>{
    private ChipChipData data;
    private double t1, t2;
    private boolean peaks;

    /* thresholds are usually conf, size, so here, thresh2 is the center
       size and thresh1 is the outer size */
    public ChipChipBindingGenerator(ChipChipData data, double thresh1, double thresh2, boolean peaks) {
        this.data = data;
        this.t1 = thresh1;
        this.t2 = thresh2;
        this.peaks = peaks;
    }
    public void setPeaks(boolean peaks) {
        this.peaks = peaks;
    }
    
    public Iterator<BindingExtent> execute(Region r) {
        ArrayList results = new ArrayList<BindingExtent>();
        //        System.err.println("Getting " + data.getName() + " binding for " + r);
        try {
            data.window(r.getChrom(),
                        r.getStart(),
                        r.getEnd());
            double avgs[] = new double[data.getCount()];
            for (int i = 0; i < data.getCount(); i++) {
                int numrepl = data.getReplicates(i);
                double avg = 0; int count = 0, j = 0;
                for (j = 0; j < numrepl; j++) {
                    avg += data.getRatio(i,j);
                }
                avg = avg / numrepl;
                avgs[i] = avg;
            }
            int lastboundpos = -1;
            double lsize = 0, rsize = 0;
            for (int i = 1; i < data.getCount() - 1; i++) {
                boolean leftok = false, rightok = false;
                if (avgs[i] >= t2) {
                    //                    System.err.println("   " + data.getName() + "  " + i + " : " + avgs[i-1] + "," + avgs[i] +"," + avgs[i+1]);
                    if ((i == 0) || (Math.abs(data.getPos(i) - data.getPos(i-1)) > 500)) {
                        leftok = false;
                    } else {
                        leftok = avgs[i] > avgs[i-1];
                        lsize = avgs[i-1];
                    }
                    if ((i == data.getCount() - 1) || (Math.abs(data.getPos(i) - data.getPos(i+1)) > 500)) {
                        rightok = false;
                    } else {
                        rightok = avgs[i] > avgs[i+1];
                        rsize = avgs[i+1];
                    }
                    double outer = Math.max(lsize,rsize);
                    if ((!peaks || (leftok && rightok)) && 
                        (data.getPos(i) != lastboundpos) &&
                        (outer >= t1)) {
                        results.add(new BindingExtent(r.getGenome(),r.getChrom(),
                                                     data.getPos(i),data.getPos(i),
                                                     avgs[i],outer,
                                                     "Raw",
                                                      data.getPos(i-1),data.getPos(i+1)));
                    }
                }
            }
        } catch (NotFoundException ex) {
            throw new RuntimeException("Couldn't window()" +ex,ex);
        }
        return results.iterator();
    }

}
