package edu.psu.compbio.seqcode.projects.shaun;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Vector;

import cern.jet.random.Exponential;
import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;
import cern.jet.random.engine.RandomEngine;

import edu.psu.compbio.seqcode.gse.datasets.binding.BindingEvent;
import edu.psu.compbio.seqcode.gse.datasets.general.NamedRegion;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.RepeatMaskedRegion;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqAlignment;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqExpt;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqHit;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqDataLoader;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.RunningOverlapSum;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.RepeatMaskedGenerator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.chipseq.SeqExpander;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

public class SeqDataHandler {
	
	private double readLength = 26.0;
	private double readExtension = 174;
	private String exptName;
	private Organism currentOrg=null;
	private Genome currentGen=null;
	private double hitCount = 0;
	private double totalSeq = 0;
	private boolean uniqueFilter=false;
	private double [] cdf;
	private int cdfThres=-1;
	private RandomEngine re;
	private Poisson P;
	private SeqExpander expander=null;
	private SeqLocator loc = null;
	
	private boolean hitsCounted=false;
	private boolean cdfCompiled=false;
	
	
	public SeqDataHandler(Organism o, Genome g, String exptName){
		this(o, g, exptName, "");
	}
	 
	public SeqDataHandler(Organism o, Genome g,  String exptName, String replicate){
		currentOrg=o;
		currentGen=g;
		this.exptName=exptName;
		
		re = new DRand();
		P = new Poisson(0, re);
		
		SeqDataLoader chipSeqLoader = null;
    	try {
    		chipSeqLoader = new SeqDataLoader();
    	} 
    	catch (SQLException sqlex) {
    		sqlex.printStackTrace();
    	} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
        if (chipSeqLoader != null) { 
        	try{        		
        		if (replicate.equals("")) {
        			List<SeqLocator> locs = new Vector<SeqLocator>();
        			
        			List<SeqExpt> expts = new Vector<SeqExpt>();
        			expts.addAll(chipSeqLoader.loadExperiments(exptName));
					
        			for (SeqExpt expt : expts) {
        				Collection<SeqAlignment> aligns;
						aligns = chipSeqLoader.loadAllAlignments(expt);
						
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
        			SeqExpt expt;
					expt = chipSeqLoader.loadExperiment(exptName, replicate);
					
        			SeqAlignment align = null;
            		Collection<SeqAlignment> aligns;
					aligns = chipSeqLoader.loadAllAlignments(expt);
					
            		for (SeqAlignment currentAlign : aligns) {
            			if (currentAlign.getGenome().equals(g)) { 
    						align = currentAlign;
    						break;
    					}
            		}
            		
            		loc = new SeqLocator(expt.getName(), expt.getReplicate(), align.getName());
        		}
        		
        		expander = new SeqExpander(loc);
	        	
	        } catch (SQLException e) {
				e.printStackTrace();
			} catch (NotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} 
        }        
	}
	
	public double getHitCount(){return(hitCount);}
	public double getTotalSeq(){return(totalSeq);}
	public void setReadLength(double rl){readLength=rl;}
	public void setReadExtension(double re){readExtension=re;}
	public double getExtendedReadLength(){return(readLength+readExtension);}
	public double getExtension(){return readExtension;}
	public Organism getOrg(){return currentOrg;}
	public Genome getGenome(){return currentGen;}
	public String getExptName(){return exptName;}
	public SeqExpander getExpander(){return expander;}
	public SeqLocator getLocator(){return loc;}
	public void useUniqueFilter(boolean unique){uniqueFilter=unique;}
	public boolean hitsAreCounted(){return hitsCounted;}
	public boolean cdfIsComiled(){return cdfCompiled;}
	
	public double countHits(){
		
		hitCount=0;
		totalSeq=0;
		//For later: if uniqueFilter is set, scan only over regions covered by unique k-mers + extensions
		Iterator<NamedRegion> chroms = new ChromRegionIterator(currentGen);
		while (chroms.hasNext()) {
			NamedRegion currentChrom = chroms.next();
			
			hitCount += expander.getHitCount(currentChrom);
			totalSeq += (double)currentChrom.getWidth();
		}		
		
		if(hitCount >0){
			P.setMean((readLength+readExtension)*(hitCount/totalSeq));
			//System.out.println("E: "+ (readLength+readExtension)*(hitCount/totalSeq)+"\tHits: "+hitCount+"\tTotal: "+totalSeq);
		}
		
		if(hitCount>0){hitsCounted=true;}
		return hitCount;
	}
	
	public void compileCDF(int maxThres) {
		cdfThres = maxThres;
		cdf = new double [maxThres+1];
		TreeMap<Integer, Double> hist = new TreeMap<Integer, Double>();
		Iterator<NamedRegion> chroms = new ChromRegionIterator(currentGen);
		int maxOverlap=0;
		
		if(!hitsCounted){
			countHits();
		}
		
		//Calculate the PDF/histogram
		while (chroms.hasNext()) {
			NamedRegion currentChrom = chroms.next();
			RunningOverlapSum totalSum = this.createRunningOverlapSum(currentChrom); 
			
			maxOverlap = Math.max(maxOverlap,totalSum.getMaxOverlap());
			
			int[][] changePoints = totalSum.getChangePoints();
			if (changePoints.length > 0) {
    			
				for (int i = 0; i < changePoints.length-1; i++) {
     				int curr = changePoints[i][1];
     				int dist = (changePoints[i+1][0]-changePoints[i][0]);
     				if (curr > 0) {
     					if (!hist.containsKey(curr)) {
         					hist.put(curr, 0.0);
         				}
         				hist.put(curr, hist.get(curr) + (double)dist);
     				}
     			}
			}
		}	
		double totalCount=0;
		for (int i = 0; i <= maxOverlap; i++) {
			double count = 0;
			if (hist.containsKey(i)) {
				count = hist.get(i);
			}
			totalCount+=count;
		}
		double zeros = totalSeq-totalCount;
		hist.put(0, zeros);
		
		//Generate the CDF
		double cumul=0;
		for (int i = 0; i <= maxThres; i++) {
			if (hist.containsKey(i)) {
				cumul +=hist.get(i).doubleValue();  
			}
			//System.out.println(i+"\t"+cumul+"\t"+hist.get(i));
			cdf[i]=cumul/totalSeq;
		}//System.out.println("");
	}
	
	public void printCDF(String filename){
		try {
			FileWriter fout = new FileWriter(filename);
			
			for(int i=0; i<cdfThres+1; i++){
				fout.write(i+"\t"+cdf[i]+"\n");
			}
			
			fout.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public void loadCDFFile(String filename){
		double [] x = new double[10000];
		try {
			File rFile = new File(filename);
			if(rFile.isFile()){
				BufferedReader reader = new BufferedReader(new FileReader(rFile));
				
				String line;
				while((line= reader.readLine())!=null){
					String [] tokens = line.split("[\\s*\\t\\r\\n\\f]");
					x[new Integer(tokens[0]).intValue()] = new Double(tokens[1]).doubleValue();
					cdfThres = new Integer(tokens[0]).intValue();
				}
				reader.close();
			}
			
			cdf = new double[cdfThres+1];
			for(int c=0; c<=cdfThres; c++){
				cdf[c]=x[c];
			}
			cdfCompiled=true;
			
			System.out.println("CDF loaded");
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public double [] selfProbLandscape(Region r){
		int rlen = r.getWidth();
		double [] probs = new double [rlen];
		for(int i=0; i<rlen; i++){
			probs[i]=0;
		}
		
		if(hitCount ==0 || cdfThres==-1){
			countHits();
			if(hitCount==0){//If still 0 after counting
				System.err.println("No hits found in genome");
			}else{
				if(cdfThres==-1){
					compileCDF(100);
				}
			}
		}if(!r.getGenome().equals(currentGen)){
			System.err.println("Genome versions not matching!");
		}else{
			//Load the observed overlapping read landscape
			int ws = Math.max(0, r.getStart() - (int)readExtension);
			int we = r.getEnd() + (int)readExtension;
			Region w = new Region(r.getGenome(), r.getChrom(), ws, we);			
			Iterator<SeqHit> hits = expander.execute(w);
			RunningOverlapSum summer = new RunningOverlapSum(r.getGenome(), r.getChrom());
			while(hits.hasNext()) { 
				SeqHit hit = hits.next();
				int ehs = hit.getStrand()=='+' ?
						hit.getStart() : 
						hit.getStart()-(int)readExtension ;
				int ehe = hit.getStrand() == '+' ? 
						hit.getEnd() + (int)readExtension : 
						hit.getEnd();
				Region ehit = new Region(hit.getGenome(), hit.getChrom(), ehs, ehe);
				summer.addRegion(ehit);
			}
			
			//iterate through region, calculating probabilities for each base
			for(int b=0; b<rlen; b++){
				int count = summer.countOverlapping(r.getStart()+b, r.getStart()+b);
				if(count<cdfThres){
					probs[b]=cdf[count];
				}else{
					probs[b]=cdf[cdfThres];
				}
			}
		}
		return(probs);
	}
		
	public double [] hitLandscape(Region r){
		int rlen = r.getWidth();
		double [] hitL = new double [rlen];
		for(int i=0; i<rlen; i++){
			hitL[i]=0;
		}
		
		if(hitCount ==0){
			countHits();
		}if(hitCount==0){//If still 0 after counting
			System.err.println("No hits found in genome");
		}else if(!r.getGenome().equals(currentGen)){
			System.err.println("Genome versions not matching!");
		}else{
			//Load the observed overlapping read landscape
			int ws = Math.max(0, r.getStart() - (int)readExtension);
			int we = r.getEnd() + (int)readExtension;
			Region w = new Region(r.getGenome(), r.getChrom(), ws, we);			
			Iterator<SeqHit> hits = expander.execute(w);
			RunningOverlapSum summer = new RunningOverlapSum(r.getGenome(), r.getChrom());
			while(hits.hasNext()) { 
				SeqHit hit = hits.next();
				int ehs = hit.getStrand()=='+' ?
						hit.getStart() : 
						hit.getStart()-(int)readExtension ;
				int ehe = hit.getStrand() == '+' ? 
						hit.getEnd() + (int)readExtension : 
						hit.getEnd();
				Region ehit = new Region(hit.getGenome(), hit.getChrom(), ehs, ehe);
				summer.addRegion(ehit);
			}
			
			//iterate through region, calculating hit counts for each base
			for(int b=0; b<rlen; b++){
				int count =summer.countOverlapping(r.getStart()+b, r.getStart()+b);
				hitL[b]=count;
			}
		}
		return(hitL);
	}
	
	public double [] poissonLandscape(Region r){
		int rlen = r.getWidth();
		double [] probs = new double [rlen];
		for(int i=0; i<rlen; i++){
			probs[i]=0;
		}
		
		if(hitCount ==0){
			countHits();
		}if(hitCount==0){//If still 0 after counting
			System.err.println("No hits found in genome");
		}else if(!r.getGenome().equals(currentGen)){
			System.err.println("Genome versions not matching!");
		}else{
			//Load the observed overlapping read landscape
			int ws = Math.max(0, r.getStart() - (int)readExtension);
			int we = r.getEnd() + (int)readExtension;
			Region w = new Region(r.getGenome(), r.getChrom(), ws, we);			
			Iterator<SeqHit> hits = expander.execute(w);
			RunningOverlapSum summer = new RunningOverlapSum(r.getGenome(), r.getChrom());
			while(hits.hasNext()) { 
				SeqHit hit = hits.next();
				int ehs = hit.getStrand()=='+' ?
						hit.getStart() : 
						hit.getStart()-(int)readExtension ;
				int ehe = hit.getStrand() == '+' ? 
						hit.getEnd() + (int)readExtension : 
						hit.getEnd();
				Region ehit = new Region(hit.getGenome(), hit.getChrom(), ehs, ehe);
				summer.addRegion(ehit);
			}
			
			//iterate through region, calculating Poisson probabilities for each base
			for(int b=0; b<rlen; b++){
				int count = summer.countOverlapping(r.getStart()+b, r.getStart()+b);
				probs[b]=P.cdf(count);
			}
		}
		return(probs);
	}
	
	
	private Collection<SeqLocator> collapseLocatorsByName(Collection<SeqLocator> locs) { 
        LinkedHashMap<String,Map<String,Set<String>>> map = 
            new LinkedHashMap<String,Map<String,Set<String>>>();
        
        for(SeqLocator loc : locs) { 
            String exptName = loc.getExptName();
            String alignName = loc.getAlignName();
            if(!map.containsKey(exptName)) { map.put(exptName, new LinkedHashMap<String,Set<String>>()); }
            if(!map.get(exptName).containsKey(alignName)) { map.get(exptName).put(alignName, new TreeSet<String>()); }
            map.get(exptName).get(alignName).addAll(loc.getReplicates());
        }
        
        LinkedList<SeqLocator> collapsed = new LinkedList<SeqLocator>();
        
        for(String exptName : map.keySet()) { 
            for(String alignName : map.get(exptName).keySet()) { 
                SeqLocator newloc = new SeqLocator(exptName, map.get(exptName).get(alignName), alignName);
                collapsed.add(newloc);
            }
        }
        
        return collapsed;
    }
	
	private RunningOverlapSum createRunningOverlapSum(Region r) {
		RunningOverlapSum totalSum = new RunningOverlapSum(r.getGenome(), r.getChrom());

		Iterator<SeqHit> iter = expander.execute(r);
		while (iter.hasNext()) {
			SeqHit hit = iter.next();
			Region extended = null;
			if(hit.getStrand() == '+') { 
				extended = new Region(hit.getGenome(), hit.getChrom(), hit.getStart(), hit.getEnd() + (int)readExtension);
			} else {
				extended = new Region(hit.getGenome(), hit.getChrom(), hit.getStart()-(int)readExtension, hit.getEnd());
			}
			totalSum.addRegion(extended);
		}
		return totalSum;
	}
	
	public void printHitHistogram() {
		Iterator<NamedRegion> chroms = new ChromRegionIterator(currentGen);
		TreeMap<Integer, Double> sparseHist = new TreeMap<Integer, Double>();
		int maxOverlap=0;
		while (chroms.hasNext()) {
			NamedRegion currentChrom = chroms.next();
			RunningOverlapSum totalSum = this.createRunningOverlapSum(currentChrom); 
			
			maxOverlap = Math.max(maxOverlap,totalSum.getMaxOverlap());
			
			int[][] changePoints = totalSum.getChangePoints();
			if (changePoints.length > 0) {
    			
     			for (int i = 0; i < changePoints.length-1; i++) {
     				int curr = changePoints[i][1];
     				int dist = (changePoints[i+1][0]-changePoints[i][0]);
     				if (curr > 0) {
     					if (!sparseHist.containsKey(curr)) {
         					sparseHist.put(curr, 0.0);
         				}
         				sparseHist.put(curr, sparseHist.get(curr) + (double)dist);
     				}
     			}
			}
		}	
		double totalCount=0;
		for (int i = 0; i <= maxOverlap; i++) {
			double count = 0;
			if (sparseHist.containsKey(i)) {
				count = sparseHist.get(i);
			}
			totalCount+=count;		
			System.out.println(i + "\t" + count);
		}
		double zeros = totalSeq-totalCount;
		System.out.println("0\t" + zeros);
	}
	
	public void printPoissonModel(){
		for(int i=0; i<100; i++){
			double exp = P.pdf(i)*totalSeq;
			System.out.println(i + "\t" + exp);
		}
	}
	
}

