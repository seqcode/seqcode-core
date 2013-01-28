package edu.psu.compbio.seqcode.gse.conservation;

import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.*;
import edu.psu.compbio.seqcode.gse.datasets.binding.BindingExtent;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.ChipChipBayes;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.ewok.*;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Expander;
import edu.psu.compbio.seqcode.gse.utils.*;

public class CustomBayesBindingGenerator implements Expander<Region,BindingExtent> {

    private ChipChipBayes data;
    private double probThresh, strengthThresh;

    public CustomBayesBindingGenerator(ChipChipBayes data, double pt, double st) { 
        this.data = data;
        probThresh = pt;
        strengthThresh = st;
    }
    
    public CustomBayesBindingGenerator(ChipChipBayes data) { 
        this.data = data;
        probThresh = 0.05;
        strengthThresh = 1.0;        
    }
    
    public boolean meetsCriteria(ChipChipBayes data, int i) { 
        return data.getPosterior(i) >= probThresh && data.getStrength(i) >= strengthThresh;
    }
    
    public Iterator <BindingExtent> execute (Region r) {
        LinkedList<BindingExtent> results = new LinkedList<BindingExtent>();
        try {
            data.window(r.getChrom(),
                        r.getStart(),
                        r.getEnd());
            int start = -1;
            int end = -1;
            double maxProb = 0.0, maxRatio = 0.0;

            for (int i = 1; i < data.getCount()-1; i++) {
                if(meetsCriteria(data, i)) { 
                    if(start == -1) { 
                        start = data.getPos(i-1);
                        maxProb = data.getPosterior(i);
                        maxRatio = data.getStrength(i);
                    } else { 
                        maxProb = Math.max(maxProb, data.getPosterior(i));
                        maxRatio = Math.max(maxRatio, data.getStrength(i));
                    }
                    end = data.getPos(i+1);
                } else { 
                    if(start != -1) { 
                        results.add(new BindingExtent(r.getGenome(),
                                r.getChrom(),
                                start, end,
                                maxRatio, maxProb, "Bayes", start, end));
                    }
                    start = -1;
                    maxProb = maxRatio = 0.0;
                }
            }
            
            if(start != -1) { 
                results.add(new BindingExtent(r.getGenome(),
                        r.getChrom(),
                        start, end,
                        maxRatio, maxProb, "Bayes", start, end));                
            }
            
        } catch (NotFoundException ex) {
            throw new RuntimeException("Couldn't window()" +ex,ex);
        }
        return results.iterator();        
    }
}
