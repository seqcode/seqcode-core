package edu.psu.compbio.seqcode.projects.shaun;

import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Iterator;

import edu.psu.compbio.seqcode.gse.datasets.binding.BindingEvent;
import edu.psu.compbio.seqcode.gse.datasets.general.NamedRegion;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Gene;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.RefGeneGenerator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.chipseq.ChipSeqBindingGenerator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.chipseq.SeqExpander;

public class SeqDataLikelihoodRatio {
	private Organism currentOrg;
	private Genome currentGen;
	private SeqDataHandler ipChannel;
	private SeqDataHandler backChannel;
	private double maxLikelihood;
	private int cdfMax = 100;
	private int minPeakWidth=50;
	private int geneWindowSize = 100000; 
	private ArrayList<Peak> peaks=null;
	
	public SeqDataLikelihoodRatio(SeqDataHandler ip, SeqDataHandler back){
		ipChannel=ip;
		backChannel=back;
		
		//Check that the genome versions match
		if(ipChannel.getOrg().equals(backChannel.getOrg()) && ipChannel.getGenome().equals(backChannel.getGenome())){
			currentOrg = ipChannel.getOrg();
			currentGen = ipChannel.getGenome();
		}else{
			System.err.println("Organism/genome mismatch");System.exit(0);
		}
		
		//Initialize the handlers if not already done
		System.out.println("Counting hits");
		double ipHits, backHits;
		if(!ipChannel.hitsAreCounted()){
			ipHits = ipChannel.countHits();
		}else{
			ipHits = ipChannel.getHitCount();
		}
		//System.out.println("Total IP hits:\t"+ipHits);
		if(!backChannel.hitsAreCounted()){
			backHits = backChannel.countHits();
		}else{
			backHits = backChannel.getHitCount();
		}			
		//System.out.println("Total background hits:\t"+backHits);
		
		if(!ipChannel.cdfIsComiled()){
			//System.out.println("Generating IP channel CDF");
			ipChannel.compileCDF(cdfMax);
		}
		if(!backChannel.cdfIsComiled()){
			//System.out.println("Generating background channel CDF");
			backChannel.compileCDF(cdfMax);
		}
	}
	
	public void calcPeaks(double likelihoodThres){
		
		maxLikelihood = likelihoodThres;
		
		int numPeaks=0;
		peaks = new ArrayList<Peak>();
		System.out.println("Calculating peaks");
		
		////////// This section only needed if taking the short-cut for peaks /////////////////////////////
		////////// Replace this section with the commented section below for the long way /////////////////
		////////// Even if we stick with this short-cut for peaks, we should still do the long-way to put the landscape into the db //////////
		/*ChipSeqLocator loc =ipChannel.getLocator();
		ChipSeqExpander csex;
		try {
			csex = new ChipSeqExpander(loc);
			ChipSeqBindingGenerator csbg = new ChipSeqBindingGenerator(csex, 174, 5, ipChannel.getExptName());
			Iterator<NamedRegion> chroms = new ChromRegionIterator(currentGen);
			while (chroms.hasNext()) {
				NamedRegion currentChrom = chroms.next();
				int regionStart = currentChrom.getStart();
				System.out.println(currentChrom.getChrom()+"\t"+currentChrom.getWidth());
				
				Iterator <BindingEvent> iter = csbg.execute(currentChrom);
				while(iter.hasNext()){
					BindingEvent r = iter.next();
					regionStart = r.getStart();


		//END OF SHORT-CUT CODE
		///////////////////////////////////////////////////////////////////////////////////////////////////
			*/		
			Iterator<NamedRegion> chroms = new ChromRegionIterator(currentGen);
			while (chroms.hasNext()) {
				NamedRegion currentChrom = chroms.next();
				int regionStart = currentChrom.getStart();
				System.out.println(currentChrom.getChrom()+"\t"+currentChrom.getWidth());
				
				//Iterate over regions here... 
				//at present this is done over 10Mbp regions
				//Later, it could iterate over unique regions only
				for(int x=0; x<currentChrom.getWidth()/10000000; x++){
					int rstart = x*10000000;
					int rend = Math.min((x+1)*10000000,currentChrom.getWidth()-1);
					Region r = new Region(currentChrom.getGenome(), currentChrom.getChrom(), rstart, rend);
					regionStart=rstart;
					
					double [] ipP = ipChannel.selfProbLandscape(r);
					double [] ipH = ipChannel.hitLandscape(r);
					double [] backH = backChannel.hitLandscape(r);
					double [] backP = backChannel.selfProbLandscape(r);
					
					//Scan through the region, adding new peaks as necessary
					Peak currPeak=new Peak();
					boolean recording=false;
					for(int i=0; i<r.getWidth(); i++){
						
						double llratio = Math.log(1-ipP[i])-Math.log(1-backP[i]);
						
						if(llratio<=maxLikelihood){
							if(!recording){
								currPeak.chr=currentChrom.getChrom();
								currPeak.start = i+regionStart;
								currPeak.end = i+regionStart;
								currPeak.maxIPHeight = ipH[i];
								currPeak.maxBackHeight = backH[i];
								currPeak.maxLikelihood = llratio;
								recording=true;
							}else{
								currPeak.end = i+regionStart;
								if(ipH[i]>currPeak.maxIPHeight){currPeak.maxIPHeight = ipH[i];}
								if(backH[i]>currPeak.maxBackHeight){currPeak.maxBackHeight = backH[i];}
								if(llratio<currPeak.maxLikelihood){currPeak.maxLikelihood = llratio;}					
							}						
						}else{
							if(recording){
								recording=false;
								if((currPeak.end-currPeak.start)>=minPeakWidth){
									numPeaks++;
									peaks.add(currPeak);
								}
								currPeak = new Peak();						
							}
						}
					}				
				}
			}		
	//	} catch (SQLException e) {
	//		// TODO Auto-generated catch block
	//		e.printStackTrace();
	//	}
	}
	
	public void printPeaks(String [] geneAnnots){
		boolean geneScan=false;
		if(geneAnnots!=null && geneAnnots.length>0){
			geneScan=true;
		}
		
		//Scan for closest genes if necessary
		if(geneScan){
			for(int g=0; g<geneAnnots.length; g++){
				RefGeneGenerator<Region> rgg = new RefGeneGenerator<Region>(currentGen, geneAnnots[g]);
				Iterator<NamedRegion> chroms = new ChromRegionIterator(currentGen);
				while (chroms.hasNext()) {
					NamedRegion currentChrom = chroms.next();
					Iterator<Gene> geneIter = rgg.execute(currentChrom);
					while (geneIter.hasNext()) {
						Gene currentGene = geneIter.next();
						int currentStart; 
						if(currentGene.getStrand() == '+')
							currentStart = currentGene.getStart();
						else
							currentStart = currentGene.getEnd();
						
						//Look through the peaks
						Iterator<Peak> itp = peaks.iterator();
						while(itp.hasNext()){
							Peak p = itp.next();
							int peakMid = (p.start+p.end)/2;
							if(p.chr.equals(currentChrom.getChrom())){
								int currentDist = (int)Math.abs(peakMid - currentStart);
								if(currentDist<geneWindowSize){
									if(currentDist < p.distToGene){
										p.nearestGene =currentGene;
										p.distToGene=currentDist;
									}
								}
							}
						}						
					}
				}
			}
		}
		//Print out the list of peaks
		Iterator<Peak> itp = peaks.iterator();
		while(itp.hasNext()){
			Peak p = itp.next();
			if(geneScan){
				if(p.distToGene==geneWindowSize+1){
					System.out.println(p.chr+":"+p.start+"-"+p.end+"\t"+p.maxIPHeight+"\t"+p.maxBackHeight+"\t"+p.maxLikelihood+"\tNO_GENE\t"+p.distToGene);
				}else{
					System.out.println(p.chr+":"+p.start+"-"+p.end+"\t"+p.maxIPHeight+"\t"+p.maxBackHeight+"\t"+p.maxLikelihood+"\t"+p.nearestGene.getName()+"\t"+p.distToGene);
				}
			}else{System.out.println(p.chr+":"+p.start+"-"+p.end+"\t"+p.maxIPHeight+"\t"+p.maxBackHeight+"\t"+p.maxLikelihood);}
		}
	
	}
	
	public double [] llratioRegion(Region r){
		double [] ipP = ipChannel.selfProbLandscape(r);
		double [] backP = backChannel.selfProbLandscape(r);
		double [] ll = new double[r.getWidth()];
		
		for(int i=0; i<r.getWidth(); i++){
			ll[i]= Math.log(1-ipP[i])-Math.log(1-backP[i]);
		}
		
		return(ll);
	}
	
	public class Peak{
		public String chr;
		public int start;
		public int end;
		public double maxIPHeight;
		public double maxBackHeight;
		public double maxLikelihood;
		public Gene nearestGene=null;
		public int distToGene=geneWindowSize+1;
	}
}
