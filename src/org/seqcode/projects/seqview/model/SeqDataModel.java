/*
 * Created on Jan 11, 2008
 *
 * TODO 
 * 
 * To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package org.seqcode.projects.seqview.model;

import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.ArrayList;
import java.util.TreeMap;
import java.util.Collection;

import org.seqcode.data.readdb.*;
import org.seqcode.data.seqdata.SeqAlignment;
import org.seqcode.data.seqdata.SeqHit;
import org.seqcode.genome.location.Region;
import org.seqcode.gsebricks.verbs.WeightedRunningOverlapSum;



/**
 * SeqDataModel:
 * SeqView data model for SeqHits.
 * If an extension is specified, then each return hit will
 * be extended by that many bp.
 */
public class SeqDataModel extends SeqViewModel implements RegionModel, Runnable {
    
    private Client client;
    private Collection<String> alignids;
    private Collection<SeqAlignment> alignments;
    private SeqAlignment align;

    private WeightedRunningOverlapSum totalSum;
    private WeightedRunningOverlapSum watsonSum;
    private WeightedRunningOverlapSum crickSum;
    private int extension;
    private int shift;
    private Region region;
    private boolean newinput, reloadInput, doSums, doHits;
    private ArrayList<SeqHit> results;
    private SeqDataProperties props;

    public SeqDataModel(Client c, Collection<SeqAlignment> alignments) throws IOException, ClientException {
        client = c;
        extension = 0;
        totalSum = null;
        watsonSum = null;
        crickSum = null;
        shift=0;
        newinput = false;
        reloadInput = false;
        doSums = true;
        doHits = true;
        alignids = new ArrayList<String>();
        this.alignments = alignments;
        align = null;
        for (SeqAlignment a : alignments) {
            alignids.add(Integer.toString(a.getDBID()));
            if (align == null) {
                align = a;
            }
            if(!client.exists(Integer.toString(a.getDBID()))){
            	System.err.println("SeqHistogramModel: Error: "+a.getExpt().getName()+";"+a.getExpt().getReplicate()+";"+a.getName()+"\tRDBID:"+a.getDBID()+" does not exist in ReadDB.");
            	dataError=true;
            }
        }

        results = new ArrayList<SeqHit>();
        props = new SeqDataProperties();
    }
    public SeqDataProperties getProperties() {return props;}

    public WeightedRunningOverlapSum getTotalRunningOverlap() { return totalSum; }
    public WeightedRunningOverlapSum getWatsonRunningOverlap() { return watsonSum; }
    public WeightedRunningOverlapSum getCrickRunningOverlap() { return crickSum; }
    
    public double getTotalMaxOverlap() { return totalSum.getMaxOverlap(); }
    public double getWatsonMaxOverlap() { return watsonSum.getMaxOverlap(); }
    public double getCrickMaxOverlap() { return crickSum.getMaxOverlap(); }
    
    public void setExtensionAndShift(int extension, int shift) { 
        if (this.extension == extension && this.shift == shift) {
            return;
        }
    	this.extension = extension;
        this.shift = shift;
    }
    
    protected void clearValues() {
        Region r = getRegion();
        totalSum = new WeightedRunningOverlapSum(r.getGenome(), r.getChrom());
        watsonSum = new WeightedRunningOverlapSum(r.getGenome(), r.getChrom());
        crickSum = new WeightedRunningOverlapSum(r.getGenome(), r.getChrom());
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
    
    public boolean isReady() {return !newinput;}
    public synchronized void run() {
        while(keepRunning()) {
            try {
                if (!newinput) {
                    wait();
                }
            } catch (InterruptedException ex) {

            }
            if (newinput) {
                try {
                    setExtensionAndShift(getProperties().ExtendRead,getProperties().ShiftRead);
                    if (doSums) {
                    	if(props.ShowType1Reads){
                    		mapToSum(totalSum, client.getWeightHistogram(alignids,
                    	                                             region.getGenome().getChromID(region.getChrom()),
                    	                                             false, //type1
                                                                     false, // paired
                                                                     extension,
                                                                     1, // binsize
                                                                     props.DeDuplicate,
                                                                     region.getStart(),
                                                                     region.getEnd(),
                                                                     null,
                                                                     null));
                    		mapToSum(watsonSum, client.getWeightHistogram(alignids,
                                                                      region.getGenome().getChromID(region.getChrom()),
                                                                      false, //type1
                                                                      false, // paired
                                                                      extension,
                                                                      1, // binsize
                                                                      props.DeDuplicate,
                                                                      region.getStart(),
                                                                      region.getEnd(),
                                                                      null,
                                                                      true));
                    		mapToSum(crickSum, client.getWeightHistogram(alignids,
                                                                      region.getGenome().getChromID(region.getChrom()),
                                                                      false, //type1
                                                                      false, // paired
                                                                      extension,
                                                                      1, // binsize
                                                                      props.DeDuplicate,
                                                                      region.getStart(),
                                                                      region.getEnd(),
                                                                      null,
                                                                      false));  
                    	}
                    	if(props.ShowType2Reads){
                    		mapToSum(totalSum, client.getWeightHistogram(alignids,
                    	                                             region.getGenome().getChromID(region.getChrom()),
                    	                                             true, //type2
                                                                     false, // paired
                                                                     extension,
                                                                     1, // binsize
                                                                     props.DeDuplicate,
                                                                     region.getStart(),
                                                                     region.getEnd(),
                                                                     null,
                                                                     null));
                    		mapToSum(watsonSum, client.getWeightHistogram(alignids,
                                                                      region.getGenome().getChromID(region.getChrom()),
                                                                      true, //type2
                                                                      false, // paired
                                                                      extension,
                                                                      1, // binsize
                                                                      props.DeDuplicate,
                                                                      region.getStart(),
                                                                      region.getEnd(),
                                                                      null,
                                                                      true));
                    		mapToSum(crickSum, client.getWeightHistogram(alignids,
                                                                      region.getGenome().getChromID(region.getChrom()),
                                                                      true, //type2
                                                                      false, // paired
                                                                      extension,
                                                                      1, // binsize
                                                                      props.DeDuplicate,
                                                                      region.getStart(),
                                                                      region.getEnd(),
                                                                      null,
                                                                      false));  
                    	}
                    }
                    if (doHits) {
                        results.clear();
                        for (String aid : alignids) {
                        	if(props.ShowType1Reads){
	                            for (SingleHit hit : client.getSingleHits(aid,
	                                                                      region.getGenome().getChromID(region.getChrom()),
	                                                                      false,
	                                                                      region.getStart(),
	                                                                      region.getEnd(),
	                                                                      null,
	                                                                      null)) {
	                                results.add(convert(hit));                                
	                            }
                        	}
                        	if(props.ShowType2Reads){
	                            for (SingleHit hit : client.getSingleHits(aid,
	                                                                      region.getGenome().getChromID(region.getChrom()),
	                                                                      true,
	                                                                      region.getStart(),
	                                                                      region.getEnd(),
	                                                                      null,
	                                                                      null)) {
	                                results.add(convert(hit));                                
	                            }
                        	}
                        }
                        if(props.DeDuplicate>0)
                        	results = deDuplicateSingleHits(results, props.DeDuplicate);
                    }
                } catch (Exception ex) {
                    ex.printStackTrace();
                }
                newinput = false;
                notifyListeners();
            }
        }
    }
    
	public void setDoSums(boolean b) {doSums = b;}
    public void setDoHits(boolean b) {doHits = b;}
    public Iterator<SeqHit> getResults() {
        return results.iterator();
    }
    private void mapToSum(WeightedRunningOverlapSum sum, TreeMap<Integer,Float> map) {
        for (Integer i : map.keySet()) {
            sum.addWeightedInterval(i, i, map.get(i));
        }
    }
    private SeqHit convert(SingleHit hit) {
        hit.length += extension;
        hit.pos += shift * (hit.strand ? 1 : -1);
        SeqHit out = new SeqHit(region.getGenome(),
                                        region.getChrom(),
                                        hit.strand ? hit.pos : hit.pos - hit.length + 1,
                                        hit.strand ? hit.pos + hit.length + 1 : hit.pos,
                                        hit.strand ? '+' : '-',
                                        hit.weight);
        return out;
    }
    //Hack to allow de-duplication of reads when drawing all reads
    private ArrayList<SeqHit> deDuplicateSingleHits(ArrayList<SeqHit> res, int deDup) {
		ArrayList<SeqHit> out = new ArrayList<SeqHit>();
		HashMap<String, Integer> hitCounts = new HashMap<String, Integer>();
		for(SeqHit h : res){
			if(hitCounts.containsKey(h.getLocationString())){
				int c = hitCounts.get(h.getLocationString());
				hitCounts.put(h.getLocationString(), c+1);
				if(c <= deDup){
					out.add(h);
				}
			}else{
				hitCounts.put(h.getLocationString(), 1);
				out.add(h);
			}
		}
		return out;
	}
    protected void doReload() {
    	synchronized(this) {
    		reloadInput = true;
			this.notifyAll();
		}
    }    
}
