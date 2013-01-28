/*
 * Created on Mar 3, 2006
 */
package edu.psu.compbio.seqcode.gse.ewok.verbs.binding;

import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.binding.BindingEvent;
import edu.psu.compbio.seqcode.gse.datasets.binding.ProbedBindingEvent;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.Probe;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Distiller;

/**
 * @author tdanford
 */
public class BayesPeakFinder implements PeakFinder<Probe> {
    
    private int strIndex, postIndex;
    private double postThresh, strThresh;
    private int maxDist;
    private String key;
    private Genome genome;

    private LinkedList<Probe> currentProbes;
    private String chrom;
    private double W, L;
    private double maxStr, maxPost;
    
    public BayesPeakFinder(String exptName, double pt, double st) {
        key = exptName;
        maxDist = 60;
        strIndex = 0;
        postIndex = 1;
        postThresh = pt;
        strThresh = st;

		currentProbes = new LinkedList<Probe>();
		chrom = "null_chrom";
		W = L = 0.0;
		maxStr = maxPost = 0.0;
    }
    
    public double getPostThreshold() { return postThresh; }
    public double getStrThreshold() { return strThresh; }
    
    public void setPostThreshold(double p) { postThresh = p; }
    public void setStrThreshold(double s) { strThresh = s; }
    
    public void reset() { 
        currentProbes.clear();
        chrom = "null_chrom";
        W = L = 0.0;
        maxStr = maxPost = 0.0;
    }
    
    /* (non-Javadoc)
     * @see edu.psu.compbio.seqcode.gse.ewok.verbs.Filter#execute(null)
     */ 
    public BindingEvent execute(Probe p) {
        double post = p.getValue(key)[postIndex], str = p.getValue(key)[strIndex];
        BindingEvent ret = null;
		//System.out.println("[ " + p.toString() + " ] ");
        if(!chrom.equals(p.getChrom()) || 
                post < postThresh || 
                str < strThresh ||
                (!currentProbes.isEmpty() &&
                        Math.abs(p.getLocation() - currentProbes.getLast().getLocation()) > maxDist)) {  

			ret = getCurrent();
            currentProbes.clear();
            chrom = p.getChrom();
            
        } else { 
            if(currentProbes.isEmpty()) { 
                chrom = p.getChrom();
                maxStr = str; 
                maxPost = post;
                W = str * post;
                L = str * post * p.getLocation();                
            } else { 
                W += (str * post);
                L += (str * post * p.getLocation());
                if(str > maxStr && post > maxPost) { 
                    maxStr = str; 
                    maxPost = post;
                }
            }
            currentProbes.addLast(p);
        }
        return ret;
    }

    /* (non-Javadoc)
     * @see edu.psu.compbio.seqcode.gse.ewok.verbs.Distiller#getCurrent()
     */
    public BindingEvent getCurrent() {
        int pos = (int)Math.round(W / L);
		if(currentProbes.isEmpty()) { return null; }
        return new ProbedBindingEvent(genome, chrom, pos, pos, maxStr, maxPost, key, currentProbes);
    }

    public static class PosteriorTunable implements TunablePeakFinder<Probe> {
        
        private BayesPeakFinder inner;
        
        public PosteriorTunable(BayesPeakFinder bpf) { 
            inner = bpf;
        }

        /* (non-Javadoc)
         * @see edu.psu.compbio.seqcode.gse.ewok.verbs.binding.TunablePeakFinder#getMaxParameter()
         */
        public double getMaxParameter() {
            return Double.MAX_VALUE;
        }

        /* (non-Javadoc)
         * @see edu.psu.compbio.seqcode.gse.ewok.verbs.binding.TunablePeakFinder#getMinParameter()
         */
        public double getMinParameter() {
            return 0.0;
        }

        /* (non-Javadoc)
         * @see edu.psu.compbio.seqcode.gse.ewok.verbs.binding.TunablePeakFinder#getCurrentParameter()
         */
        public double getCurrentParameter() {
            return inner.getPostThreshold();
        }

        /* (non-Javadoc)
         * @see edu.psu.compbio.seqcode.gse.ewok.verbs.binding.TunablePeakFinder#setParameter(double)
         */
        public void setParameter(double p) {
            inner.setPostThreshold(p);
        }

        /* (non-Javadoc)
         * @see edu.psu.compbio.seqcode.gse.ewok.verbs.Distiller#execute(java.lang.Object)
         */
        public BindingEvent execute(Probe a) {
            return inner.execute(a);
        }

        /* (non-Javadoc)
         * @see edu.psu.compbio.seqcode.gse.ewok.verbs.Distiller#getCurrent()
         */
        public BindingEvent getCurrent() {
            return inner.getCurrent();
        }

        public void reset() {
            inner.reset();
        }
    }
    
    public static class StrengthTunable implements TunablePeakFinder<Probe> {
        
        private BayesPeakFinder inner;
        
        public StrengthTunable(BayesPeakFinder bpf) { 
            inner = bpf;
        }

        /* (non-Javadoc)
         * @see edu.psu.compbio.seqcode.gse.ewok.verbs.binding.TunablePeakFinder#getMaxParameter()
         */
        public double getMaxParameter() {
            return Double.MAX_VALUE;
        }

        /* (non-Javadoc)
         * @see edu.psu.compbio.seqcode.gse.ewok.verbs.binding.TunablePeakFinder#getMinParameter()
         */
        public double getMinParameter() {
            return 0.0;
        }

        /* (non-Javadoc)
         * @see edu.psu.compbio.seqcode.gse.ewok.verbs.binding.TunablePeakFinder#getCurrentParameter()
         */
        public double getCurrentParameter() {
            return inner.getStrThreshold();
        }

        /* (non-Javadoc)
         * @see edu.psu.compbio.seqcode.gse.ewok.verbs.binding.TunablePeakFinder#setParameter(double)
         */
        public void setParameter(double p) {
            inner.setStrThreshold(p);
        }

        /* (non-Javadoc)
         * @see edu.psu.compbio.seqcode.gse.ewok.verbs.Distiller#execute(java.lang.Object)
         */
        public BindingEvent execute(Probe a) {
            return inner.execute(a);
        }

        /* (non-Javadoc)
         * @see edu.psu.compbio.seqcode.gse.ewok.verbs.Distiller#getCurrent()
         */
        public BindingEvent getCurrent() {
            return inner.getCurrent();
        }

        public void reset() {
            inner.reset();
        }
    }
}
