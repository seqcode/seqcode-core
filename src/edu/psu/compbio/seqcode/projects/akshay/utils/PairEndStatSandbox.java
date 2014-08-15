package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Genome.ChromosomeInfo;
import edu.psu.compbio.seqcode.genome.location.NamedRegion;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqAlignment;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqDataLoader;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqHitPair;
import edu.psu.compbio.seqcode.genome.location.Region;

import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.chipseq.SeqExpander;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.projects.shaun.Utilities;

public class PairEndStatSandbox {
	private SeqLocator locator;
	private SeqExpander expander;
	private boolean cleanExpander;
	private Genome gen;
	
	//stats
	//public double[] fragment_size;
	
	
	public PairEndStatSandbox(SeqLocator sloc, Genome g) throws SQLException, IOException {
		this.locator = sloc;
		this.expander = new SeqExpander(sloc);
		this.cleanExpander = true;
		this.gen = g;
		
	}
	
	// Fetchers
	
	public Iterator<SeqHitPair> loadHitsbyChrom(NamedRegion chrom){
		Iterator<SeqHitPair> ret = this.expander.getPairs(chrom);
		return ret;
	}
	
	
	
	
	
	
	//Stats Calculators
	
	public void printFragmentsSizes(){
		double[] fragment_size = new double[500];
		
		ChromRegionIterator chrItr = new ChromRegionIterator(this.gen);
		
		while(chrItr.hasNext()){
			NamedRegion chrom = chrItr.next();
			if(!chrom.getName().contains("random")){
				Iterator<SeqHitPair> pairItr = this.loadHitsbyChrom(chrom);
				while(pairItr.hasNext()){
					SeqHitPair pair = pairItr.next();
					if(pair.getCode() == 1 && pair.getMidpoint() != null){
						int dist  = pair.getLeft().getFivePrime() > pair.getRight().getFivePrime() ? pair.getLeft().getFivePrime() - pair.getRight().getFivePrime() : pair.getRight().getFivePrime() - pair.getLeft().getFivePrime();
						if(dist <=500){
							fragment_size[dist-1]++;
						}
					}
				
				}
			}
		}
		
		double total = 0;
		for(int i=0; i<fragment_size.length; i++){
			total = total+ fragment_size[i];
		}
		for(int i=0; i<fragment_size.length; i++){
			fragment_size[i] = fragment_size[i]/total;
		}
		
		for(int i=0; i<fragment_size.length; i++){
			System.out.println(Integer.toString(i+1)+"\t"+Double.toString(fragment_size[i]));
		}
		
		
	}
	
	// make sure win is -1 if you pass a regions file instead of peaks file
	public void printWindowCounts(String peaks_file, int win, int minFragLen, int maxFragLen){
		List<Region> regs  = Utilities.loadRegionsFromPeakFile(gen, peaks_file, win);
		for(Region r : regs){
			int count = 0;
			Iterator<SeqHitPair> pairItr = this.expander.getPairs(r);
			while(pairItr.hasNext()){
				SeqHitPair pair = pairItr.next();
				if(pair.getCode() == 1 && pair.getMidpoint() != null){
					int fragLen = Math.abs(pair.getRight().getFivePrime() - pair.getLeft().getFivePrime())+1;
					if(fragLen >= minFragLen && fragLen < maxFragLen){
						count++;
					}
				}
			}
			System.out.println(r.getLocationString()+"\t"+Integer.toString(count));
		}
	}
	
	
	
	
	
	
	public void clean(){
		if(isNotClean()){
			this.expander.close();
			this.cleanExpander = false;
		}
	}


    public boolean isNotClean() {
        return this.cleanExpander;
    }
	
	
	public static void main(String[] args) throws NotFoundException, SQLException, IOException{
		ArgParser ap = new ArgParser(args);
		SeqLocator expt = Args.parseSeqExpt(args, "expt").get(0);
		Genome g = Args.parseGenome(args).cdr();
		PairEndStatSandbox analyzer = new PairEndStatSandbox(expt, g);
		
		if(ap.hasKey("fragSize")){
			analyzer.printFragmentsSizes();
			analyzer.clean();
		}else if(ap.hasKey("windowCount")){
			String peaks_file;
			int win;
			if(ap.hasKey("peaks")){
				peaks_file = ap.getKeyValue("peaks");
				win = Args.parseInteger(args, "win", -1);
			}else{
				peaks_file = ap.getKeyValue("region");
				win = -1;
			}
			
			int minFragLen = Args.parseInteger(args, "minFrag", 140);
			int maxFragLen = Args.parseInteger(args, "maxFrag", 200);
			
			analyzer.printWindowCounts(peaks_file, win, minFragLen, maxFragLen);
			analyzer.clean();
			
		}
		
	}
}
