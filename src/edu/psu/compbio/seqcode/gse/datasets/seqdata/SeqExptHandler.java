package edu.psu.compbio.seqcode.gse.datasets.seqdata;

import java.sql.SQLException;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import java.util.Vector;
import java.io.IOException;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqAlignment;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqExpt;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqHit;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqDataLoader;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

public class SeqExptHandler {
	
    public static final double defaultReadLength = 26.0, defaultReadExtension = 174.0;
    private double readLength = defaultReadLength;
    private double readExtension = defaultReadExtension;
    private double readShift = 0;
    private String exptName;
    private Genome currentGen = null;
    private double hitCount = 0;
    private double hitWeight = 0;
    private double totalSeq = 0;
    
    private SeqLocator loc = null;
    private SeqDataLoader loader;
    private LinkedList<SeqAlignment> alignments;


    public SeqExptHandler(Genome g, String exptName) throws NotFoundException, IOException {
        this(g, exptName, "");
    }


    public SeqExptHandler(Genome g, String exptName, String replicate) throws NotFoundException, IOException {
        currentGen = g;
        this.exptName = exptName;

        SeqDataLoader seqLoader = null;
        try {
            seqLoader = new SeqDataLoader();
        }
        catch (SQLException sqlex) {
            sqlex.printStackTrace();
        }

        if (seqLoader != null) {
            try {
                if (replicate.equals("")) {
                    List<SeqLocator> locs = new Vector<SeqLocator>();

                    List<SeqExpt> expts = new Vector<SeqExpt>();
                    expts.addAll(seqLoader.loadExperiments(exptName));

                    for (SeqExpt expt : expts) {
                        Collection<SeqAlignment> aligns;
                        aligns = seqLoader.loadAllAlignments(expt);

                        for (SeqAlignment currentAlign : aligns) {
                            if (currentAlign.getGenome().equals(g)) {
                				SeqLocator currentLoc = new SeqLocator(expt.getName(), 
                                                                               expt.getReplicate(), currentAlign.getName());
                                locs.add(currentLoc);
                                break;
                            }
                        }
                    }

                    List<SeqLocator> collapsedLocs = new Vector<SeqLocator>(this.collapseLocatorsByName(locs));
                    if (collapsedLocs.size() != 1) {
                        System.err.println(collapsedLocs.size() + " collapsed locators");
                        System.exit(0);
                    }
                    loc = collapsedLocs.get(0);
                }
                else {
                    SeqExpt expt = seqLoader.loadExperiment(exptName, replicate);

                    SeqAlignment align = null;
                    Collection<SeqAlignment> aligns;
                    aligns = seqLoader.loadAllAlignments(expt);

                    for (SeqAlignment currentAlign : aligns) {
                        if (currentAlign.getGenome().equals(g)) {
                            align = currentAlign;
                            break;
                        }
                    }

                    loc = new SeqLocator(expt.getName(), expt.getReplicate(), align.getName());
                }

                loader = new SeqDataLoader();
                alignments = new LinkedList<SeqAlignment>();

                if (loc.getReplicates().isEmpty()) {
                    Collection<SeqExpt> expts = loader.loadExperiments(loc.getExptName());
                    for (SeqExpt expt : expts) {
                        SeqAlignment alignment = loader.loadAlignment(expt, loc.getAlignName(), g);
                        if (alignment != null) {
                            alignments.add(alignment);
                        }
                    }
                } else {
                    for (String repName : loc.getReplicates()) {
                        SeqExpt expt = loader.loadExperiment(loc.getExptName(), repName);
                        SeqAlignment alignment = loader.loadAlignment(expt, loc.getAlignName(), g);
                        if (alignment != null) {
                            alignments.add(alignment);
                        }
                    }
                }
                countHits();
            }
            catch (SQLException e) {
                e.printStackTrace();
            }
            catch (NotFoundException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        }
    }


    public SeqExptHandler(Genome g, SeqLocator locator) throws NotFoundException, SQLException, IOException {
        currentGen = g;
        exptName = locator.getExptName();
        loader = new SeqDataLoader();
        alignments = new LinkedList<SeqAlignment>();
        if (locator.getAlignName() == null) {
            if(locator.getReplicates().isEmpty()) { //No alignment name, no replicate names
                Collection<SeqExpt> expts = loader.loadExperiments(locator.getExptName());
                List<SeqLocator> locs = new Vector<SeqLocator>();
                for (SeqExpt expt : expts) {
                    Collection<SeqAlignment> aligns;
                    aligns = loader.loadAllAlignments(expt);
                    for (SeqAlignment currentAlign : aligns) {
                        if (currentAlign.getGenome().equals(g)) {
            				SeqLocator currentLoc = new SeqLocator(expt.getName(), 
                                                                           expt.getReplicate(), currentAlign.getName());
                            locs.add(currentLoc);
                            alignments.add(currentAlign);
                            break;
                        }
                    }
                }
                List<SeqLocator> collapsedLocs = new Vector<SeqLocator>(this.collapseLocatorsByName(locs));
                if (collapsedLocs.size() != 1) {
                    System.err.println(collapsedLocs.size() + " collapsed locators");
                    System.exit(0);
                }
                locator = collapsedLocs.get(0);
            }
            else { // No alignment name, given replicate names
                for (String repName : locator.getReplicates()) {
                    SeqExpt expt = loader.loadExperiment(locator.getExptName(), repName);
                    SeqAlignment alignment = loader.loadAlignment(expt, locator.getAlignName(), g);
                    if (alignment != null) {
                        locator = new SeqLocator(locator.getExptName(), locator.getReplicates(), alignment.getName());
                        alignments.add(alignment);
                        break;
                    }
                }
            }
        }
        else {
            if (locator.getReplicates().isEmpty()) {// Given alignment name, no
                // replicate names
                Collection<SeqExpt> expts = loader.loadExperiments(locator.getExptName());
                List<SeqLocator> locs = new Vector<SeqLocator>();
                for (SeqExpt expt : expts) {
                    Collection<SeqAlignment> aligns;
                    aligns = loader.loadAllAlignments(expt);
                    for (SeqAlignment currentAlign : aligns) {
                        if (currentAlign.getGenome().equals(g) && currentAlign.getName().equals(locator.getAlignName())) {
            				SeqLocator currentLoc = new SeqLocator(expt.getName(), 
                                                                           expt.getReplicate(), currentAlign.getName());
                            locs.add(currentLoc);
                            alignments.add(currentAlign);
                            break;
                        }
                    }
                }
                List<SeqLocator> collapsedLocs = new Vector<SeqLocator>(this.collapseLocatorsByName(locs));
                if (collapsedLocs.size() != 1) {
                    System.err.println(collapsedLocs.size() + " collapsed locators");
                    System.exit(0);
                }
                locator = collapsedLocs.get(0);
            }
            else {
                for (String replicate : locator.getReplicates()) {// Given alignment name, 
                    // given replicate names
                    alignments.add(loader.loadAlignment(loader.loadExperiment(locator.getExptName(), 
                                                                              replicate), locator.getAlignName(), g));
                }
            }
        }
        loc = locator;
        //countHits();
    }


    public double getHitCount() {
        return (hitCount);
    }

    public double getHitWeight() {
        return (hitWeight);
    }

  
    public double getTotalSeq() {
        return (totalSeq);
    }


    public void setReadLength(double rl) {
        readLength = rl;
    }


    public void setReadExtension(double re) {
        readExtension = re;
    }


    public void setReadShift(double rs) {
        readShift = rs;
        readExtension = readShift * 2;
    }


    public double getReadLength() {
        return (readLength);
    }


    public double getExtendedReadLength() {
        return (readLength + readExtension);
    }


    public double getExtension() {
        return readExtension;
    }


    public Genome getGenome() {
        return currentGen;
    }


    public String getExptName() {
        return exptName;
    }


    public SeqLocator getLocator() {
        return loc;
    }


    public double countHits() {
        hitCount = 0;
        try {
            for (SeqAlignment alignment : alignments) {
                double currHits = (double) loader.countAllHits(alignment);
                hitCount += currHits;
            }
        }
        catch (IOException e) {
            e.printStackTrace();
            return 0;
        }
     
        return hitCount;
    }


    public double countHits(Region a) {
        double hits = 0;
        try {
            hits = (double) loader.countByRegion(alignments, a);
            return hits;
        }
        catch (IOException e) {
            e.printStackTrace();
            return 0;
        }
    }

    public double weighHits() {
		hitWeight = 0;
		try {
			for (SeqAlignment align : alignments) {
				double currWeight = loader.weighAllHits(align);
				hitWeight += currWeight;
			}
        }
        catch (IOException e) {
            e.printStackTrace();
            return 0;
        }
		return hitWeight;
	}


    public double weighHits(Region a) {
        double hits = 0;
        try {
            hits = (double) loader.weightByRegion(alignments, a);
            return hits;
        }
        catch (IOException e) {
            e.printStackTrace();
            return 0;
        }
    }


  
    public LinkedList<StrandedRegion> loadHits(Region a) {
        try {
            LinkedList<StrandedRegion> total = new LinkedList<StrandedRegion>();
            Collection<SeqHit> hits = loader.loadByRegion(alignments, a);
            for (SeqHit curr : hits) {
                total.add(hit2region(0, curr));
            }
            return total;
        }
        catch (IOException e) {
            e.printStackTrace();
            return null;
        }
    }


   public LinkedList<SeqHit> loadExtendedHits(Region a) {
     try {
       LinkedList<SeqHit> total = new LinkedList<SeqHit>();
       Collection<SeqHit> hits = loader.loadByRegion(alignments, a);
       for (SeqHit curr : hits) {
           total.add(curr.extendHit((int) readExtension));
       }
       return total;
     } catch (IOException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
		return new LinkedList<SeqHit>();
	}
   }


   public LinkedList<SeqHit> loadShiftedExtendedHits(Region a) {
     try {
       LinkedList<SeqHit> total = new LinkedList<SeqHit>();
       Collection<SeqHit> hits = loader.loadByRegion(alignments, a);
       for (SeqHit curr : hits) {
           total.add(curr.shiftExtendHit((int) readExtension, (int) readShift));
       }
       return total;
     } catch (IOException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
		return new LinkedList<SeqHit>();
	}
   }


    public int[] hitCountLandscape(Region a) {
        int[] hits = new int[a.getWidth() + 1];
        for (int i = 0; i <= a.getWidth(); i++) {
            hits[i] = 0;
        }

        try {
            Iterator<SeqHit> h = loader.loadByRegion(alignments, a).iterator();
            while (h.hasNext()) {
                SeqHit csh = h.next();
                int pos = csh.getStart() - a.getStart();
                if (pos < 0) {
                    pos = 0;
                }
                if (pos > a.getWidth())
                    pos = a.getWidth();
                hits[pos]++;
            }
        }
        catch (IOException e) {
            e.printStackTrace();
            return null;
        }
        return hits;
    }


    // Extension can be zero
    // Strand can be +/-/., with '.' for both strands
    public int[] hitDepthLandscape(Region a, int ext, char strand) {
        int[] hits = new int[a.getWidth() + 1];
        for (int i = 0; i <= a.getWidth(); i++)
            hits[i] = 0;
        try {
            Iterator<SeqHit> h = loader.loadByRegion(alignments, a).iterator();
            while (h.hasNext()) {
                SeqHit csh = h.next();
                int start, end;
                if (csh.getStrand() == '+') {
                    start = Math.max(0, csh.getStart() - a.getStart());
                    end = Math.min(a.getWidth(), (csh.getEnd() + ext) - a.getStart());
                }
                else {
                    start = Math.max(0, (csh.getStart() - ext) - a.getStart());
                    end = Math.min(a.getWidth(), csh.getEnd() - a.getStart());
                }
                if (strand == '.' || csh.getStrand() == strand)
                    for (int p = start; p <= end; p++)
                        hits[p]++;
            }
        }
        catch (IOException e) {
            e.printStackTrace();
            return null;
        }
        return hits;
    }

    //Bins the start positions of hits in a region
	//Strand can be +/-/., with '.' for both strands
	public int [] binHitStarts(Region a, int binSize, char strand) {
		int numBins = ((a.getWidth())/binSize)+1;
		int [] bins = new int [numBins];
		for(int i=0; i<numBins; i++)
			bins[i]=0;
		try {
			Iterator<SeqHit> h = loader.loadByRegion(alignments, a).iterator();
			while(h.hasNext()){
				SeqHit csh = h.next();
				int start, end;
				if(csh.getStrand()=='+'){
					start = Math.max(0, csh.getStart()-a.getStart());					
				}else{
					start =Math.min(a.getWidth(), csh.getEnd()-a.getStart());
				}
				if(strand == '.' || csh.getStrand() == strand)
					bins[start/binSize]++;
			}		
        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }
		return bins;
	}

    private Collection<SeqLocator> collapseLocatorsByName(Collection<SeqLocator> locs) {
        LinkedHashMap<String, Map<String, Set<String>>> map = new LinkedHashMap<String, Map<String, Set<String>>>();

        for (SeqLocator loc : locs) {
            String exptName = loc.getExptName();
            String alignName = loc.getAlignName();
            if (!map.containsKey(exptName)) {
                map.put(exptName, new LinkedHashMap<String, Set<String>>());
            }
            if (!map.get(exptName).containsKey(alignName)) {
                map.get(exptName).put(alignName, new TreeSet<String>());
            }
            map.get(exptName).get(alignName).addAll(loc.getReplicates());
        }

        LinkedList<SeqLocator> collapsed = new LinkedList<SeqLocator>();

        for (String exptName : map.keySet()) {
            for (String alignName : map.get(exptName).keySet()) {
                SeqLocator newloc = new SeqLocator(exptName, map.get(exptName).get(alignName), alignName);
                collapsed.add(newloc);
            }
        }

        return collapsed;
    }


    private StrandedRegion hit2region(int ext, SeqHit hit) {
        if (hit.getStrand() == '+') {
            return new StrandedRegion(hit.getGenome(), hit.getChrom(), hit.getStart(), hit.getEnd() + ext, hit.getStrand());
        }
        else {
            return new StrandedRegion(hit.getGenome(), hit.getChrom(), hit.getStart() - ext, hit.getEnd(), hit.getStrand());
        }
    }


    private StrandedRegion hit2region(int ext, SeqHit hit, int shift) {
        if (hit.getStrand() == '+') {
            return new StrandedRegion(hit.getGenome(), hit.getChrom(), hit.getStart() + shift - (ext / 2), hit.getEnd()
                                      + shift + (ext / 2), hit.getStrand());
        }
        else {
            return new StrandedRegion(hit.getGenome(), hit.getChrom(), hit.getStart() - shift - (ext / 2), hit.getEnd()
                                      - shift + (ext / 2), hit.getStrand());
        }
    }
    
    public void close(){
    	loader.close();
    }
}
