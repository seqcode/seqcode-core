package edu.psu.compbio.seqcode.projects.shaun;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.PropertyResourceBundle;
import java.util.TreeMap;

import edu.psu.compbio.seqcode.gse.datasets.chipseq.ChipSeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.deepseq.DeepSeqExpt;
import edu.psu.compbio.seqcode.gse.deepseq.ReadHit;
import edu.psu.compbio.seqcode.gse.projects.readdb.PairedHit;
import edu.psu.compbio.seqcode.gse.projects.readdb.PairedHits;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;

public class TatsuReadExperiment {
	private boolean needlefiltering=true;
	public Pair<Organism,Genome> pair;
	public List<ChipSeqLocator> expts;
	public ChipSeqCountReads counter;
	public Genome gen;
	public DeepSeqExpt signal;
	
	public static void main(String args[]) throws NotFoundException, SQLException, FileNotFoundException, IOException{
		TatsuReadExperiment re = new TatsuReadExperiment(args);
		//re.getReads("12", 58641618, 58651618, 0, 10);
		re.getReads("1", 10371090, 10371100, 0, 1);
		//re.getAllRegions("1",100000,0);
		//re.setPaired();
		//re.getPairedReads();
	}
	
	public TatsuReadExperiment(String args[]) throws NotFoundException, SQLException, FileNotFoundException, IOException {
		for (String string : args) {
			System.out.println(string);
		}    
		boolean  metaPeak=false;
        pair = Args.parseGenome(args);
        expts = Args.parseChipSeq(args,"expt");
        if (expts.size() == 0) {
            System.err.println("Usage:\n " +
                               "ChipSeqCountReads " +
                               "--species <organism name;genome version> "+
                               "--expt <solexa expt> " );
            return;
        }
        gen =pair.cdr();
        signal = new DeepSeqExpt(gen, expts, "readdb", 32);
        //Initialize the peak finder
	}
	
	public double[] getAllRegions(String chrom,int by,int strand){
		int clen = gen.getChromLength(chrom);
		int start = 1;
		int end = start+clen;
		return getMetaReads(chrom,start,end,by,strand);
	}
	
	public double[] getMetaReads(String chrom,int start, int end, int by,int strand){
		int width = end-start;
		int nument = (int)Math.floor(width/by);
		double[] hits = new double[nument+1];
	/*	Region x;
		x=new Region(gen,chrom,start,end);
		char strandch;
		if(strand==0){
			strandch='-';
		}else{
			strandch='+';
		}
		TreeMap<Integer,Float> hist= signal.makeHistogram(x, by,strandch);
		int i = 0;
		for (Iterator iterator = hist.keySet().iterator(); iterator.hasNext();) {
			Integer ival = (Integer)iterator.next();
			Float fp = hist.get(ival);
			hits[(int)(Math.floor((ival.intValue()-start)/by))]+=fp.doubleValue();
			i++;
		}*/
		return hits;
	}
	public void getPairedReads(){
		Region r = new Region(gen, "10", 3002630, 3002850);
		List<ReadHit> hits = signal.loadHits(r);
		for(ReadHit h : hits)
			System.out.println("ReadHit= "+h.toString());
	}
	public double[] getReads(String chrom,int start,int end,int strandint,int pileup){
		System.out.println("searching..");
		System.out.println(chrom);
		System.out.println(start);
		System.out.println(end);
		System.out.println(pileup);
		double iphittot = signal.getHitCount();
		System.out.println(iphittot);
		Region x = new Region(gen,chrom,start,end);
		char strand='.';
		if(strandint==0){
			strand='-';
		}else{
			strand='+';
		}
		List<ReadHit> hits = signal.loadHits(x);
		for(ReadHit h : hits)
			System.out.println("ReadHit= "+h.toString());
		ArrayList<Double> hls=makeHitLandscape(hits,x,pileup,strand);
		 double[] target = new double[hls.size()];
		 double sum=0;
		 System.out.println("ListSize= "+hls.size());
		 for (int i = 0; i < target.length; i++) {
		    target[i] = hls.get(i).doubleValue();                // java 1.5+ style (outboxing)
		    sum=sum+hls.get(i).doubleValue();
		    System.out.println("Landscape\t"+i+"\t"+hls.get(i));
		 }
		 System.out.println(sum);
		 return target;
	}
	
	public void setPaired(){signal.setPairedEnd(true);}
	
	public double[] getAnchoredReads(String chrom,int start,int end,int pileup){
		System.out.println("searching..");
		System.out.println(chrom);
		System.out.println(start);
		System.out.println(end);
		System.out.println(pileup);
		double iphittot = signal.getHitCount();
		System.out.println(iphittot);
		Region x = new Region(gen,chrom,start,end);
		char strand='+';
		//ArrayList<PairedHit> hits = (ArrayList)signal.loadAnchoredPairs(x);
		//System.out.println(hits.size());
		double[] toret = null;//new double[hits.size()*2];
		//for(int i = 0; i<hits.size(); i++){
		//	PairedHit p = hits.get(i);
		//	toret[2*i]=p.leftPos;
		//	toret[2*i+1]=p.rightPos;			
		//}
		 return toret;
	}
	
	protected ArrayList<Double> makeHitLandscape(List<ReadHit> hits, Region currReg, int perBaseMax, char strand){
		int numBins = (int)(currReg.getWidth());
		int [] counts = new int[currReg.getWidth()+1];
		ArrayList<Double> startcounts = new ArrayList<Double>();
		for(int i=0; i<=numBins; i++){startcounts.add(new Double(0));}
		for(int i=0; i<=currReg.getWidth(); i++){counts[i]=0;}
		for(ReadHit r : hits){
			if(strand=='.' || r.getStrand()==strand){
				int offset = r.getFivePrime()-currReg.getStart();
				if(offset >=0 && offset <=numBins){
					double currWeight = (needlefiltering && counts[offset]+r.getWeight()>perBaseMax) ?
							Math.max(0, perBaseMax-counts[offset]) : 
							r.getWeight();
					counts[offset]+=r.getWeight();
					double countnow=startcounts.get(offset);
					startcounts.set(offset,countnow+currWeight);
				}
			}
		}
		return startcounts;
	}
	
		
	
	protected final int inBounds(int x, int min, int max){
		if(x<min){return min;}
		if(x>max){return max;}
		return x;
	}
}