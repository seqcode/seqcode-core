package org.seqcode.projects.shaun.rnaseq;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;

import org.seqcode.gse.projects.gps.DeepSeqExpt;
import org.seqcode.gse.projects.gps.ReadHit;
import org.seqcode.projects.shaun.rnaseq.genemodels.GeneTUnit;
import org.seqcode.projects.shaun.rnaseq.genemodels.TUnit;


public class RNASeqExonicCounter extends RNASeqAnalyzer{

	public RNASeqExonicCounter(String[] args) {
		super(args);
	}

	public void execute(){
		Double totalHitCount = 0.0;
		for(DeepSeqExpt e : experiments){
			totalHitCount += e.getHitCount();
		}
		
		System.err.println("Counting reads in "+knownGenes.size()+" known genes");
		
		//Dumb iteration over all genes for now
		for(GeneTUnit gene : knownGenes){
			ArrayList<ReadHit> hits = new ArrayList<ReadHit>();
			
			//If we have multiple experiments, just pool all data
			for(DeepSeqExpt e : experiments){
				hits.addAll(e.loadHits(gene.getCoords()));
			}
			
			//System.err.println(gene.getName()+"\t"+hits.size());
			
			//Assign hits to genes if they overlap exons. Each hit added only once
			Iterator<TUnit> ite = gene.getExonIterator();
			while(ite.hasNext()){
				TUnit exon = ite.next();
				for(ReadHit h : hits){
					if(h.getFivePrime()>=exon.getCoords().getStart() && h.getFivePrime()<=exon.getCoords().getEnd())
						gene.addHits(h.getWeight()*h.getReadLength());
				}
			}
		}
		
		//Print counts
		try{
			FileWriter fw = new FileWriter(outputName);
			fw.write(String.format("Name\tID\tCoord\tType\tExonHits\tExonLengths\tFPKM\n"));
			for(GeneTUnit gene : knownGenes){
				double fpkm = ((gene.getHits()/readLength)/((double)gene.getExonLength()/1000.0))/(totalHitCount/1000000);
				fw.write(String.format("%s\t%s\t%s\t%s\t%.2f\t%d\t%.2f\n", gene.getName(),gene.getID(),gene.getCoords().getLocationString(),gene.getType(),gene.getHits(),gene.getExonLength(), fpkm));
			}
			fw.close();
			System.out.println("Results written to "+outputName);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//Close loaders
		this.cleanup();
	}
	
	public static void main(String[] args) {
		RNASeqExonicCounter counter = new RNASeqExonicCounter(args);
		counter.execute();
	}
	
	public void printError() {
		printError("RNASeqExonicCounter");
	}	
}
