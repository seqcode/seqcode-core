package org.seqcode.gse.gsebricks.verbs.location;

import java.util.*;

import org.seqcode.genome.location.Region;
import org.seqcode.gse.gsebricks.*;
import org.seqcode.gse.gsebricks.types.*;
import org.seqcode.gse.gsebricks.verbs.Mapper;
import org.seqcode.gse.gsebricks.verbs.MapperIterator;


public class RegionExpander<X extends Region> 
    implements Mapper<X,Region>, SelfDescribingVerb, DefaultConstantsParameterized {

    private int minWidth;
    
    public RegionExpander() { 
        minWidth = 0;
    }

    public RegionExpander(int minw) {
        minWidth = minw;
    }

    public Iterator<Region> execute(Iterator<X> r) {
        return new MapperIterator<X,Region>(this,r);
    }
    
    public Region execute(X r) {
        if(r.getWidth() >= minWidth) { 
            return r;
        } else {
            int diff = minWidth - r.getWidth();
            int startdelta = diff/2;
            int enddelta = diff-startdelta;
            return r.expand(startdelta, enddelta);
        }
    }
    
    /*
     * Methods for the implementation of SelfDescribingVerb
     */

    public EchoType getInputClass() {
        return new ClassType(Region.class);
    }

    public EchoType getOutputClass() {
        return new ClassType(Region.class);
    }
    
    private static final EchoType[] pclasses = { new ClassType(Integer.class) };
    private static final String[] pnames = { "MinWidth" };
    
    private static final EchoType[] inputClasses = { new ClassType(Region.class) };
    private static final String[] inputNames = { "Regions" };
    
    private static final SelfDescribingConstant[] defConsts = { new ValueWrapper(1000) };
    private static final String[] defConstNames = { "MinWidth" };
    
    public String[] getInputNames() { 
        return inputNames; 
    }
    
    public EchoType[] getInputClasses() { 
        return inputClasses; 
    }

    public EchoType[] getParameterClasses() {
        return pclasses;
    }

    public String[] getParameterNames() {
        return pnames;
    }

    public void init(Map<String, Object> params) {
        int minw = (Integer)params.get(pnames[0]);
        minWidth = minw;
    }

    public String[] defaultConstantNames() {
        return defConstNames;
    }

    public SelfDescribingConstant[] defaultConstants() {
        return defConsts;
    }
}

