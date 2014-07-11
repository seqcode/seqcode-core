package edu.psu.compbio.seqcode.projects.akshay.Utils;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqAlignment;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqDataLoader;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqHitPair;

import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.verbs.chipseq.SeqExpander;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

public class PairEndStatSandbox {
	
	// Current this class uses the old seqexpander methods instead of the new hitloaders. ..
	// Change to the new hitloader methods as soon as possible
	
	private SeqLocator locator;
	
	private LinkedList<SeqAlignment> alignments;
	private Genome lastGenome;
	private SeqDataLoader loader;
	private boolean closeLoader;
	private Genome gen;
	
	//stats
	
	public double[] fragment_size; 
	
	public PairEndStatSandbox(SeqLocator expt, Genome g) throws SQLException, IOException {
		this.locator = expt;
		this.loader = new SeqDataLoader();
		closeLoader = true;
		this.gen = g;
		try {
			this.getAligns(g);
		} catch (SQLException e) {
			e.printStackTrace();
		}	
	}
	private void getAligns(Genome genome) throws SQLException {
		if (alignments != null && genome.equals(lastGenome)) {
            return;
        }
        alignments = new LinkedList<SeqAlignment>();
        lastGenome = genome;
        try {
            alignments.addAll(locator.loadAlignments(loader, genome));
        } catch (SQLException e) {
            e.printStackTrace(System.err);
        } catch (NotFoundException e) {
            e.printStackTrace();
        }
    }
	
	
	
	// Fetchers
	
	public List<SeqHitPair> loadHitsbyChrom(String chr){
		List<SeqHitPair> ret = new ArrayList<SeqHitPair>();
		int chrID = gen.getChromID(chr);
		for(SeqAlignment a : this.alignments){
			try {
				List<SeqHitPair> temp = loader.loadPairsByChrom(a, chrID);
				ret.addAll(temp);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		return ret;
		
	}
	
	
	
	
	
	
	//Stats Calculators
	
	public void setFragmentsSizes(){
		this.fragment_size = new double[500];
		for(String chr : this.gen.getChromList()){
			if(!chr.contains("random")){
				List<SeqHitPair> tmp = this.loadHitsbyChrom(chr);
				for(SeqHitPair pair : tmp){
					if(pair.getCode() == 1 && pair.getMidpoint() != null){
						int dist  = pair.getLeft().getFivePrime() > pair.getRight().getFivePrime() ? pair.getLeft().getFivePrime() - pair.getRight().getFivePrime() : pair.getRight().getFivePrime() - pair.getLeft().getFivePrime();
						if(dist <=500){
							this.fragment_size[dist-1]++;
						}
					}
				
				}
			}
		}
		
		double total = 0;
		for(int i=0; i<this.fragment_size.length; i++){
			total = total+ this.fragment_size[i];
		}
		for(int i=0; i<this.fragment_size.length; i++){
			this.fragment_size[i] = this.fragment_size[i]/total;
		}
	}
	
	
	
	
	// getters
	
	public double[] getFragmentSizes(){return this.fragment_size;}
	
	
	
	public void close() {
        if (closeLoader) {
            loader.close();
        }
        loader = null;
        if (alignments != null) {
            alignments.clear();
        }
    }


    public boolean isClosed() {
        return loader == null;
    }
	
	
	public static void main(String[] args) throws NotFoundException, SQLException, IOException{
		SeqLocator expt = Args.parseSeqExpt(args, "expt").get(0);
		Genome g = Args.parseGenome(args).cdr();
		PairEndStatSandbox analyzer = new PairEndStatSandbox(expt, g);
		analyzer.setFragmentsSizes();
		double[] frags = analyzer.getFragmentSizes();
		analyzer.close();
		for(int i=0; i<frags.length; i++){
			System.out.println(Integer.toString(i+1)+"\t"+Double.toString(frags[i]));
		}
		
		
	}
}
