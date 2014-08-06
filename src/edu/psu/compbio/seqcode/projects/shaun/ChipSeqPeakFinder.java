package edu.psu.compbio.seqcode.projects.shaun;

import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.sql.SQLException;

import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;
import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Organism;
import edu.psu.compbio.seqcode.genome.location.Gene;
import edu.psu.compbio.seqcode.genome.location.NamedRegion;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqExptHandler;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqHit;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.*;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.NamedGeneratorFactory;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.NamedStrandedGeneratorFactory;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.NamedTypedGeneratorFactory;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.RefGeneGenerator;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.RepeatMaskedGenerator;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;

/**
 * Finds peaks in Solexa experiments. 
 * This peak finder reads in IP experiments and control/background experiments and tries to find peaks in IP that 
 * are significantly enriched over background. The genome is split into overlapping bins and the numbers of (extended)
 * reads falling into each bin are counted. Two thresholds are employed for peak finding. One is that the IP channel 
 * bin must contain significantly more reads than expected from a completely random scattering of reads (i.e. Poisson
 * distributed) and the background channel must have less than the corresponding limit. The other threshold is set using 
 * the binomial test for differences between two population proportions; the IP bin must be significantly different than
 * the background bin. Neighboring bins are merged. Predicted peak information is printed to a file. Sequences surrounding
 * the peaks are printed to a separate file.   
 * 
 *  The peak finder can read from multiple different experiments for IP and background channels. Provide each experiment
 *  name using separate "--expt" or "--back" options. If specifying the replicate number, use a semi-colon as in 
 *  "expt;replicate". If no replicate name is provided, ALL replicates for a given experiment name are loaded. 
 *  
 *  Example usage:
 *  ChipSeqPeakFinder --species "Mus musculus" --genome mm8 --expt PPG_Solexa_RAR_8hr [--back PPG_Solexa_WCE_ES+2d [--back PPG_Solexa_WCE_2+1]] --outpeak rar_8hr.peaks --outseq rar_8hr_peaks.seq --seqwin 200   
 * 
 */
public class ChipSeqPeakFinder {
	protected double winWidth=100;
	protected double winStep = 25;
	protected double sigThres = 2.33; //p<0.01?
	protected double poissThres=-9; //10th power for Poission threshold
	protected double readLength, readExtension, readShift;
	protected boolean shiftTags=false;
	protected int seqwin = -1; //-1 equals all of the peak region sequence is printed
    protected boolean noBackground;
    protected boolean geneOverlap; // if true, then only return genes that overlap the bound region; don't do the closest gene thing
    protected int maxGeneDistance=50000;
	
	protected Genome gen;
	protected double genomeLen=0;
	protected double ipGenomePoissonThres=0;
	protected double ipBasePoissonThres=0;
	protected double backGenomePoissonThres=0;
	protected double backBasePoissonThres=0;
	protected ArrayList<Integer> poissonThresModels=new ArrayList<Integer>(); //-1 for entire genome, 0 for current region, length of window otherwise
	protected int numPoissonModels=0;
	protected double[] ipPoissonThresholds;
	protected double[] backPoissonThresholds;
	protected double [] ipStackedHitCounts;
	protected double[] backStackedHitCounts;
	protected ArrayList<SeqExptHandler> IPhandles;
	protected ArrayList<SeqExptHandler> backhandles;
	protected double iphittot;
	protected double backhittot;
	protected boolean towerfiltering=true, needlefiltering=false, peakTrimming=true;
	protected boolean repeatMaskerAnnots = false, peakCenter = false;
	protected boolean addStrandedness = false;
	protected int towerWindow=(int)winWidth*10;
	protected String outPeakName="peak.out", outSeqName="seq.out";
	public String [] annots = new String[]{
		"refGene"//,  "knownGene"//, "mgcGenes", "ensGene"
	};
	protected boolean genePeaksOnly=false;
    public Collection<String> namedRegions, namedStrandedRegions, namedTypedRegions;
    
	
	/* command-line driver */
	public static void main(String[] args) throws SQLException, NotFoundException {
		boolean  metaPeak=false;
		Pair<Organism,Genome> pair = Args.parseGenome(args);
		if(pair==null){printError();return;}
		List<SeqLocator> expts = Args.parseSeqExpt(args,"expt");
        List<SeqLocator> back = Args.parseSeqExpt(args,"back");
        if (expts.size() == 0) {
            printError();
            return;
        }
        double rLen = Args.parseDouble(args,"readlen",SeqExptHandler.defaultReadLength);
        double rExt = Args.parseDouble(args,"readextend",SeqExptHandler.defaultReadExtension);
        double rShift =0;
        if(Args.parseArgs(args).contains("shifttags")){
        	rShift = Args.parseInteger(args,"shifttags", 0);
        }
    
        //Initialize the peak finder
        ChipSeqPeakFinder finder = new ChipSeqPeakFinder(pair.cdr(), expts, back, rLen, rExt, rShift);
        
        //Load options for peak calling 
        if(Args.parseArgs(args).contains("dynpoisson")){finder.addPoissonModel(Args.parseInteger(args,"dynpoisson",5000));}
        finder.setPeakFileName(Args.parseString(args,"outpeak",finder.outPeakName));
        finder.setSeqFileName(Args.parseString(args,"outseq",finder.outSeqName));
        finder.setSeqwin(Args.parseInteger(args,"seqwin",finder.seqwin));
        finder.setBinWidth(Args.parseDouble(args,"binwidth",finder.winWidth));
        finder.setBinStep(Args.parseDouble(args,"binstep",finder.winStep));
        finder.setMinZ(Args.parseDouble(args,"minz",finder.sigThres));
        finder.setPoissonThres(Args.parseDouble(args,"poisson",finder.poissThres));
        finder.setTowerFilter(!Args.parseFlags(args).contains("allowtowers"));
        finder.setNeedleFilter(!Args.parseFlags(args).contains("allowneedles"));
        finder.setGeneOverlap(Args.parseFlags(args).contains("geneOverlap"));
        finder.setMaxGeneDistance(Args.parseInteger(args,"maxgenedist",finder.maxGeneDistance));
        finder.peakCenter = Args.parseFlags(args).contains("peakcenter");
        Collection<String> genes = Args.parseStrings(args,"genes");
        if (genes.size() > 0) {
            finder.annots = new String[genes.size()];
            int i = 0;
            Iterator<String> iter = genes.iterator();
            while (iter.hasNext()) {
                finder.annots[i++] = iter.next();
            }
        } 
        finder.namedRegions = Args.parseStrings(args,"namedregion");
        finder.namedStrandedRegions = Args.parseStrings(args,"namedstrandedregion");
        finder.namedTypedRegions = Args.parseStrings(args,"namedtypedregion");
        finder.repeatMaskerAnnots = Args.parseFlags(args).contains("repeatmasker");
        metaPeak = Args.parseArgs(args).contains("printpeakdistrib");
        finder.setGeneScanning(Args.parseFlags(args).contains("scangenesonly"));
      
        
        //Run the peak finder
        ArrayList<ChipSeqPeak> enriched = finder.execute();
        System.out.println("Printing");
		finder.printEnrichedPeaks(enriched);
		finder.printPeakSequences(enriched);
		if(metaPeak){
			System.out.println("Printing a metapeak");
			finder.printIPMetaPeak(enriched, 2000, Args.parseInteger(args, "printpeakdistrib", 200));
		}		
		
		System.out.println("Finished!");		
	}
	
	public ChipSeqPeakFinder(){
		
	}
	/* Constructor requires a species, genome version and lists of ip and background solexa experiments */
	public ChipSeqPeakFinder(Genome gen, List<SeqLocator> ips, List<SeqLocator> backs, double rLen, double rExt, double rShift) {
		System.out.println("Initializing the Peak Finder");
		this.gen = gen;
        readLength = rLen;
        readExtension = rExt;
        readShift = rShift;
        if(readShift > 0)
        	shiftTags=true;
        iphittot=0; backhittot=0;
                
        //Loading experiments
        loadExperiments(ips, backs);
	}
	
	protected void loadExperiments(List<SeqLocator> ips, List<SeqLocator> backs){
	  IPhandles = new ArrayList<SeqExptHandler>();
	  backhandles = new ArrayList<SeqExptHandler>();
	  try {
	    genomeLen = gen.getGenomeLength(); 

	    for(SeqLocator ip : ips){
	      System.out.print(String.format("%s\t", ip.getExptName()));
	      SeqExptHandler curr = new SeqExptHandler(gen, ip);
	      curr.setReadLength(readLength);
	      if(shiftTags){
	        curr.setReadExtension(readShift*2);
	        curr.setReadShift(readShift);
	      }
	      else {
	        curr.setReadExtension(readExtension);
	      }	      
	      IPhandles.add(curr);
	    }

	    for(SeqLocator back : backs){
	      System.out.print(String.format("%s\t", back.getExptName()));
	      SeqExptHandler curr = new SeqExptHandler(gen, back);
	      curr.setReadLength(readLength);
	      if(shiftTags){
	        curr.setReadExtension(readShift*2);
	        curr.setReadShift(readShift);
	      } 
	      else {
	        curr.setReadExtension(readExtension);
	      }
	      backhandles.add(curr);
	    }
	    noBackground = backhandles.isEmpty();

	  } 
	  catch (NotFoundException e) {
	    e.printStackTrace();
	  } 
	  catch (SQLException e) {
	    // TODO Auto-generated catch block
	    e.printStackTrace();
	  } catch (IOException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	}
	}
	
	
	/* Run the peak-finder (assumes all initialization took place successfully) */
	public ArrayList<ChipSeqPeak> execute() {
	//Initialize the iterator for the test regions/genes
    Iterator<NamedRegion> testRegions=null;
    if(!genePeaksOnly)
      testRegions = new ChromRegionIterator(gen);
    else{
      ArrayList<NamedRegion> geneList = new ArrayList<NamedRegion>(); 
      for(int g=0; g<annots.length; g++){
        RefGeneGenerator<Region> rgg = new RefGeneGenerator<Region>(gen, annots[g]);
        ChromRegionIterator chroms = new ChromRegionIterator(gen);
        while(chroms.hasNext()){
          NamedRegion c = chroms.next();
          Iterator<Gene> geneIter = rgg.execute(c);
                    while (geneIter.hasNext()) {
                        Gene gene = geneIter.next();   
                        geneList.add(new NamedRegion(gene.getGenome(),gene.getChrom(), gene.getStart(),gene.getEnd(),gene.getName()));
                    }
        }
      }
      testRegions = geneList.iterator();
    }
    
    //weigh hits (for test regions only??)
    weightHits();
	  
    //Initialize Poisson Models
    addPoissonModel(-1); //customize further        

		ipBasePoissonThres = getPoissonThreshold(Math.pow(10, poissThres), iphittot, 1, genomeLen, 0.8, 1, 1);
		ipGenomePoissonThres = getPoissonThreshold(Math.pow(10, poissThres), iphittot, IPhandles.get(0).getExtendedReadLength(), genomeLen, 0.8, winWidth, winStep);
        if (noBackground) {
            backBasePoissonThres = 1;
            backGenomePoissonThres = 1;
        }  else {
        	backBasePoissonThres = getPoissonThreshold(Math.pow(10, poissThres), backhittot, 1, genomeLen, 0.8, 1, 1);
        	backGenomePoissonThres = getPoissonThreshold(Math.pow(10, poissThres), backhittot, backhandles.get(0).getExtendedReadLength(), genomeLen, 0.8, winWidth, winStep);
        }        
		System.out.println("IP Genome-wide Poisson Threshold: "+ipGenomePoissonThres+"\nBackground Genome-wide Poisson Threshold: "+backGenomePoissonThres);
		System.out.println("IP Per Base Poisson Threshold: "+ipBasePoissonThres+"\nBackground Per Base Poisson Threshold: "+backBasePoissonThres);
		System.out.println(String.format("Genome Length: %.0f", genomeLen));
		


		//Call the peak finder
		System.out.println("Finding enriched regions");
		ArrayList<ChipSeqPeak> enriched = callEnrichedRegions(testRegions);
		System.out.println("Matching peaks with closest genes");
		if (peakCenter) {
		  for (ChipSeqPeak peak : enriched) {
		    int pos = (peak.coords.getStart() + peak.coords.getEnd()) / 2;
		    peak.coords = new Region(peak.coords.getGenome(),
		        peak.coords.getChrom(), pos, pos);
		  }
		}

		findClosestGenes(enriched, annots);
		for (String s : namedRegions) {
		  addRegionAnnotations(enriched, (new NamedGeneratorFactory(s)).getExpander(gen));
		}
		for (String s : namedStrandedRegions) {
		  addRegionAnnotations(enriched, (new NamedStrandedGeneratorFactory(s)).getExpander(gen));
		}
		for (String s : namedTypedRegions) {
		  addRegionAnnotations(enriched, (new NamedTypedGeneratorFactory(s)).getExpander(gen));
		}
		if (repeatMaskerAnnots) {
		  addRegionAnnotations(enriched, new RepeatMaskedGenerator(gen));
		}

		return(enriched);
	}
	
	
	private void weightHits() {
	  for(SeqExptHandler handler : IPhandles) {
      handler.weighHits();
      //iphittot += 20958530; 
      iphittot += handler.getHitWeight();
    }
    System.out.print(String.format("%.0f ip read weights loaded\n", iphittot));
    for(SeqExptHandler handler : backhandles) {
      System.out.print(String.format("%s\t", handler.getExptName()));
      handler.weighHits();
      //backhittot += 3072065;
      backhittot += handler.getHitWeight();
    }
    System.out.print(String.format("%.0f bg read weights loaded\n", backhittot));
	}
	
	private void weightHitsByRegion(List<Region> regions) {
	  for(SeqExptHandler handler : IPhandles) {
	    System.out.print(String.format("%s\t", handler.getExptName()));
	    for (Region reg : regions) {

	      //TODO
	      //	        handler.weighHits(reg);
	      //iphittot += handler.getHitWeight();
	    }
	  }
	  System.out.print(String.format("%.0f ip read weights loaded\n", iphittot));
	  for(SeqExptHandler handler : backhandles) {
	    System.out.print(String.format("%s\t", handler.getExptName()));
	    for (Region reg : regions) {
	      //TODO
	      //	        handler.weighHits(reg);
	      //	        backhittot += handler.getHitWeight();
	    }
	  }
	  System.out.print(String.format("%.0f bg read weights loaded\n", backhittot));
	}


	//Some option setters
	public void setSeqwin(int s){seqwin=s;}
	public void setBinWidth(double b){winWidth=b;}
	public void setBinStep(double b){winStep=b;}
	public void setMinZ(double m){sigThres=m;}
	public void setPoissonThres(double p){poissThres=p;}
	public void setPeakFileName(String f){outPeakName=f;}
	public void setSeqFileName(String f){outSeqName=f;}
	public void setTowerFilter(boolean f){towerfiltering=f;}
	public void setNeedleFilter(boolean f){needlefiltering=f;}
    public void setGeneOverlap(boolean b){geneOverlap=b;}
    public void setMaxGeneDistance(int i){maxGeneDistance=i;}
    public void setGeneScanning(boolean g){genePeaksOnly = g;}
    
    /*Refreshes the Poisson thresholds based on current environment */
    protected void refreshPoissonThresholds(Region currReg, int currOffset){
    	for(int p=0; p<numPoissonModels; p++){
    		if(poissonThresModels.get(p)==-1){//Genome-wide
    			if(ipGenomePoissonThres==0){
    				ipGenomePoissonThres = getPoissonThreshold(Math.pow(10, poissThres), iphittot, IPhandles.get(0).getExtendedReadLength(), genomeLen, 0.8, winWidth, winStep);
    			}ipPoissonThresholds[p]= ipGenomePoissonThres;
    			if(!noBackground && backGenomePoissonThres==0){
    				backGenomePoissonThres = getPoissonThreshold(Math.pow(10, poissThres), backhittot, backhandles.get(0).getExtendedReadLength(), genomeLen, 0.8, winWidth, winStep);
    			}backPoissonThresholds[p]= backGenomePoissonThres;
    		}else if(poissonThresModels.get(p)==0){//Current region
    			double totalIP=0;
    			for(int i=0; i<ipStackedHitCounts.length; i++)
    				totalIP+=ipStackedHitCounts[i];
    			double estReads = totalIP/(((readLength+readExtension)/winStep)+1);
    			ipPoissonThresholds[p]= (double)getPoissonThreshold(Math.pow(10, poissThres), estReads, IPhandles.get(0).getExtendedReadLength(), currReg.getWidth(), 1.0, winWidth, winStep);
    			if(!noBackground){
	    			double totalBack=0;
	    			for(int i=0; i<backStackedHitCounts.length; i++)
	    				totalBack+=backStackedHitCounts[i];
	    			estReads = totalBack/(((readLength+readExtension)/winStep)+1);
	    			backPoissonThresholds[p]= (double)getPoissonThreshold(Math.pow(10, poissThres), estReads, backhandles.get(0).getExtendedReadLength(), currReg.getWidth(), 1.0, winWidth, winStep);
    			}
    		}else{//Window around current position
    			int win = poissonThresModels.get(p);
    			double winIP=0;
    			int istart = currOffset-(win/2)<0 ? 0 :currOffset-(win/2);
    			int istop = currOffset+(win/2)>=currReg.getWidth() ? currReg.getWidth()-1 :currOffset+(win/2);
    			int winTrueSize=0;
    			for(int i=istart; i<=istop; i++){
    				int binid = (int)Math.max(0, ((double)(i)/winStep));
    				winIP+=ipStackedHitCounts[binid]; winTrueSize+=winStep;
    			}
    			double estReads = winIP/(((readLength+readExtension)/winStep)+1);
    			ipPoissonThresholds[p]= (double)getPoissonThreshold(Math.pow(10, poissThres), estReads, IPhandles.get(0).getExtendedReadLength(), winTrueSize, 1.0, winWidth, winStep);
    			if(!noBackground){
	    			double winBack=0;
	    			istop = currOffset+(win/2)>=currReg.getWidth() ? currReg.getWidth()-1 :currOffset+(win/2);
	    			winTrueSize=0;
	    			for(int i=istart; i<=istop; i++){
	    				int binid = (int)Math.max(0, ((double)(i)/winStep));
	    				winBack+=backStackedHitCounts[binid]; winTrueSize+=winStep;
	    			}
	    			estReads = winIP/(((readLength+readExtension)/winStep)+1);
	    			backPoissonThresholds[p]= (double)getPoissonThreshold(Math.pow(10, poissThres), estReads, backhandles.get(0).getExtendedReadLength(), winTrueSize, 1.0, winWidth, winStep);
    			}
    		}
    	}
    }
	
	/*This is the actual peak-finder; it returns an arraylist of peaks fulfilling the thresholds*/
	public ArrayList<ChipSeqPeak> callEnrichedRegions(Iterator<NamedRegion> testRegions){
		double basesDone=0, printStep=10000000,  numPrint=0;;
		
		ArrayList<ChipSeqPeak> allres = new ArrayList<ChipSeqPeak>();
		ChipSeqPeak lastHit=null;
		double ipTotHits = iphittot;
        // set to 1 if no background so the binomial test doesn't return NaN
		double backTotHits = noBackground ? 1 : backhittot;        
		
		while (testRegions.hasNext()) {
			NamedRegion currentRegion = testRegions.next();
			
			lastHit=null;
			//Split the job up into chunks of 100Mbp
			int chunkSize = 1000000;
			for(int x=currentRegion.getStart(); x<=currentRegion.getEnd(); x+=chunkSize){
				int y = x+chunkSize; 
				if(y>currentRegion.getEnd()){y=currentRegion.getEnd();}
				Region currSubRegion = new Region(gen, currentRegion.getChrom(), x, y);
				ArrayList<ChipSeqPeak> currres = new ArrayList<ChipSeqPeak>();
				
				LinkedList<SeqHit> ipHits = new LinkedList<SeqHit>();
				LinkedList<SeqHit> backHits = new LinkedList<SeqHit>();
				
				for(SeqExptHandler IP: IPhandles){
					if(shiftTags){
						ipHits.addAll(IP.loadShiftedExtendedHits(currSubRegion));
					}else
						ipHits.addAll(IP.loadExtendedHits(currSubRegion));
				}
				for(SeqExptHandler back: backhandles){
					if(shiftTags)
						backHits.addAll(back.loadShiftedExtendedHits(currSubRegion));
					else
						backHits.addAll(back.loadExtendedHits(currSubRegion));
				}
                
				//Sets ipStackedHitCounts, backStackedHitCounts, ipHitCounts, backHitCounts
				ipStackedHitCounts = makeHitLandscape(ipHits, currSubRegion, (int)ipBasePoissonThres);
				backStackedHitCounts = null;
                if (!noBackground) {
                	backStackedHitCounts = makeHitLandscape(backHits, currSubRegion, (int)backBasePoissonThres);
                }
				
                //Scan regions
				for(int i=currSubRegion.getStart(); i<currSubRegion.getEnd()-(int)winWidth; i+=(int)winStep){
					Region currWin = new Region(gen, currentRegion.getChrom(), i, (int)(i+winWidth-1));
					
					int binid = (int)Math.max(0, ((double)(currWin.getStart()-currSubRegion.getStart())/winStep));
					double ipWinHits=(double)ipStackedHitCounts[binid];
					double backWinHits= noBackground ? 0 : ((double)backStackedHitCounts[binid]);
					
					//Should this region be added as a hit?
					//first pass: genome thresholds
					if(ipWinHits >= ipGenomePoissonThres && backWinHits < backGenomePoissonThres){
						//Second pass: refresh all thresholds
						refreshPoissonThresholds(currSubRegion, i-x);
						if(ipWinHits >= (double)getMaxPoissonThreshold(ipPoissonThresholds) && backWinHits < (double)getMaxPoissonThreshold(backPoissonThresholds)){
						
							double Z0 = binomialSampleEquality(ipWinHits, backWinHits, ipTotHits, backTotHits);
							double over=-1;
							if(noBackground || Z0 > sigThres){
								if(backWinHits > 0){
									over = (ipWinHits/ipTotHits)/(backWinHits/backTotHits);
								}
								//Is this hit close to a previously added one?
								if(lastHit!=null && currWin.getStart()-lastHit.coords.getEnd()<=winWidth){
									Region merge = new Region(gen, lastHit.coords.getChrom(), lastHit.coords.getStart(), currWin.getEnd());
									lastHit.coords=merge;
									if(ipWinHits>lastHit.ipHits){lastHit.ipHits=ipWinHits;}
									if(backWinHits>lastHit.backHits){lastHit.backHits=backWinHits;}
									if(Z0>lastHit.Z){lastHit.Z=Z0;}
									if(over>lastHit.overrep){lastHit.overrep=over;}								
									//Update the peak
									lastHit.peak = findPeakMaxHit(ipHits, lastHit);
								}else{
									ChipSeqPeak hit = new ChipSeqPeak();
									hit.coords = currWin;
									hit.Z = Z0;
									hit.overrep = over;
									hit.ipHits = ipWinHits;
									hit.backHits = backWinHits;
									//Find the peak
									hit.peak = findPeakMaxHit(ipHits, hit);
									currres.add(hit);
									lastHit=hit;								
								}
							}
						}
					}					
					
					//Print out progress
					basesDone+=winStep;
					if(basesDone > numPrint*printStep){
						if(numPrint%10==0){System.out.print(String.format("(%.0f)", (numPrint*printStep)));}
						else{System.out.print(".");}
						if(numPrint%50==0 && numPrint!=0){System.out.print("\n");}
						numPrint++;
					}
				}
				//A filter for those darned towers and other nasty crap
				if(towerfiltering){
                    //System.err.println("Filtering towers");
					ArrayList<ChipSeqPeak> tmpres = filterTowers(currres, ipHits, backStackedHitCounts, currSubRegion.getStart());
					currres = tmpres;
				}
				if(addStrandedness){
					addStrands(currres, ipHits);
				}
				if(peakTrimming){
					trimPeaks(currres, ipHits);
				}
					allres.addAll(currres);
				}
		}System.out.print("\n");
		Collections.sort(allres);
		return(allres);
	}
	
	/* Filter the tower; protected because this method only makes sense from within the callEnrichedRegions method */
	protected ArrayList<ChipSeqPeak> filterTowers(ArrayList<ChipSeqPeak> curr, LinkedList<SeqHit> ipHits, double [] backHitCounts, int chromOffset){
		ArrayList<ChipSeqPeak> res = new ArrayList<ChipSeqPeak>();
		Iterator<ChipSeqPeak> itr = curr.iterator();
		while(itr.hasNext()){
			ChipSeqPeak peak = itr.next();
			boolean isTower=false; //innocent until proven guilty
			
			ArrayList<SeqHit> winHits = overlappingHits(ipHits, peak.coords);
			
			//Filter 1: remove things that are ridiculously one-stranded in the IP channel
			double strandprop=forwardReadProportion(winHits);
			if(strandprop>=0.85 || strandprop<=0.15){
				isTower=true;
			}
			
			//Filter 2: remove "peaks" that neighbor background towers
			Region currTowerWin = new Region(gen, peak.coords.getChrom(), peak.coords.getStart()-towerWindow, peak.coords.getEnd()+towerWindow);
			for(int i=currTowerWin.getStart(); i<currTowerWin.getEnd(); i+=(int)winStep){
				int binid = (int)Math.max(0, ((double)(i-chromOffset)/winStep));
				if(!noBackground && binid < backHitCounts.length){
					double backWinHits=(double)backHitCounts[binid];
					if(backWinHits>=backGenomePoissonThres){
						isTower=true;
					}
				}
			}
			
			//Filter 3: Another class of tower are the needles that stack up over a very tight space. This isn't what we'd expect a TF binding site to look like.
			double all=0, mid=0, num=0;
			for(SeqHit hit : winHits){ //get mid-point of all hit starts
				if(hit.getStrand()=='+'){
					all+=(double)hit.getStart();
				}else{
					all+=(double)hit.getEnd();
				}
				num++;
			}mid=all/num;
			//if most of the hits are within one read-length of the mid-point, this is a needle
			double needles=0;
			for(SeqHit hit : winHits){ //get mid-point of all hit starts
				if(hit.getStrand()=='+' && Math.abs((double)hit.getStart()-mid)<=readLength){needles++;}
				else if(hit.getStrand()=='-' && Math.abs((double)hit.getEnd()-mid)<=readLength){needles++;}
			}if(needles/num>=0.95){isTower=true;}
			
			
			if(!isTower){
				res.add(peak);
			}
		}
		return(res);
	}
	
	/** print out the peaks.  Fields are
     * - coordinates
     * - region size
     * - location string
     * {- strand if CLIP}
     * - distance from location string to start of region
     * - number of hits (IP)
     * - number of hits (background)
     * - Z (compares probability of read being in this region in IP and background)
     * - overrepresentation factor
     * - closest gene
     * - distance to gene
     * - extra annotations if queried
     */
	public void printEnrichedPeaks(ArrayList<ChipSeqPeak> enriched){
	  FileWriter fout = null;
	  try {
			fout = new FileWriter(outPeakName);
			
			Iterator<ChipSeqPeak> itr = enriched.iterator();
			while(itr.hasNext()){
				ChipSeqPeak peak = itr.next();
				String closestGene = (peak.nearestGene==null) ? "NONE":peak.nearestGene.getName();
                StringBuilder annots = new StringBuilder();
                if (peak.annotations != null) {
                    for (Region r : peak.annotations) {
                        annots.append(r.toString() + ",");
                    }
                }
                if (annots.length() > 1) {
                    annots.setLength(annots.length() - 1);
                }
				int peakToStart = peak.peak.getLocation() - peak.coords.getStart();
				if(addStrandedness)
					fout.write(peak.coords.getLocationString()+"\t"+peak.strand+"\t"+peak.coords.getWidth()+"\t"+peak.peak.getLocationString()+"\t"+peakToStart+"\t"+peak.ipHits+"\t"+peak.backHits+"\t"+peak.Z+"\t"+peak.overrep+"\t"+closestGene+"\t"+peak.distToGene+"\t" + annots.toString() + "\n");
				else
				fout.write(peak.coords.getLocationString()+"\t"+peak.coords.getWidth()+"\t"+peak.peak.getLocationString()+"\t"+peakToStart+"\t"+peak.ipHits+"\t"+peak.backHits+"\t"+peak.Z+"\t"+peak.overrep+"\t"+closestGene+"\t"+peak.distToGene+"\t" + annots.toString() + "\n");
			}			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		finally {
		  if (fout != null) {
		    try {
		      fout.close();
		    }
		    catch (IOException ioex) {
		      ioex.printStackTrace();
		    }
		  }
		}
	}
	/* print out sequences around the peak */
	public void printPeakSequences(ArrayList<ChipSeqPeak> enriched){
		FileWriter fout = null;
	  try {
			fout = new FileWriter(outSeqName);
			SequenceGenerator seqgen = new SequenceGenerator();
			Iterator<ChipSeqPeak> itr = enriched.iterator();
			while(itr.hasNext()){
				ChipSeqPeak peak = itr.next();
				String closestGene = (peak.nearestGene==null) ? "NONE":peak.nearestGene.getName();
				
				Region peakWin;
				if(seqwin == -1){
					peakWin = peak.coords;
				}else{
					int start = peak.peak.getLocation()-((int)(seqwin/2));
					if(start<1){start=1;}
					int end = peak.peak.getLocation()+((int)seqwin/2)-1;
					if(end>gen.getChromLength(peak.coords.getChrom())){end =gen.getChromLength(peak.coords.getChrom());} 
					peakWin = new Region(peak.peak.getGenome(), peak.peak.getChrom(), start, end);
				}
				
				String seqName = new String(">"+peakWin.getLocationString()+"\t"+peakWin.getWidth());
				String seq = seqgen.execute(peakWin);
				fout.write(seqName +"\n"+ seq+"\n");
			}
		} 
	  catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	  finally {
	    if (fout != null) {
	      try {
	        fout.close();
	      }
	      catch (IOException ioex) {
	        ioex.printStackTrace();
	      }	      	     
	    }
	  }
	}
	
	/* Binomial test for differences between two population proportions */
	protected double binomialSampleEquality(double X1, double X2, double n1, double n2){
		double P1 = X1/n1;
		double P2 = X2/n2;
		double P = (X1+X2)/(n1+n2);
		double Z0 = (P1-P2)/(Math.sqrt(P*(1-P)*((1/n1)+(1/n2))));
		return(Z0);
	}
	
	protected ArrayList<SeqHit> overlappingHits(LinkedList<SeqHit> hits, Region window){
		ArrayList<SeqHit> sub = new ArrayList<SeqHit>();
		for(SeqHit r : hits){
			if(window.overlaps(r)){sub.add(r);}			
		}
		return(sub);
	}
	protected double [] makeHitLandscape(LinkedList<SeqHit> hits, Region currReg, int perBaseMax){
		int numBins = (int)(currReg.getWidth()/winStep);
		double [] counts = new double[currReg.getWidth()+1];
		double [] land = new double[numBins+1];
		for(int i=0; i<=numBins; i++){land[i]=0;}
		for(int i=0; i<=currReg.getWidth(); i++){counts[i]=0;}
		for(SeqHit r : hits){
			int offset=r.getStart()-currReg.getStart()<0 ? 0 : r.getStart()-currReg.getStart();
			counts[offset] = counts[offset] + r.getWeight();
			if(!needlefiltering || (counts[offset] <= perBaseMax)){
				//int binstart = (int)Math.max(0, ((double)(offset/winStep)-(Math.floor(winWidth/winStep)-1)));
				//int binend = (int)(Math.min((double)(r.getEnd()-currReg.getStart()), (double)currReg.getWidth())/winStep);
				int binstart = (int)Math.max(0, ((double)(offset/winStep)));
				int binend = (int)(Math.min((double)(r.getEnd()-currReg.getStart())/winStep, numBins-1));
				for(int i=binstart; i<=binend; i++){
					land[i] = land[i] + r.getWeight();
				}
			}
		}
		return(land);
	}

    /** add annotations to the ChipSeqPeaks.  An annotation is a subclass of Region
     * that is near or associated with the peak
     */
    protected void addRegionAnnotations(ArrayList<ChipSeqPeak> enriched,
                                      Expander<Region, ? extends Region> expander) {
        System.err.println("Adding annots with " + expander);
        for (ChipSeqPeak peak : enriched) {
            Region query;
            if (geneOverlap) {
                query = peak.coords;
            } else {
                query = peak.coords.expand(maxGeneDistance, maxGeneDistance);
            }
            Iterator<? extends Region> iter = expander.execute(query);
            while (iter.hasNext()) {
                peak.addAnnotation(iter.next());
            }
        }
    }
	
	//Arguments: a list of enriched peaks and an array of annotation source names
	protected void findClosestGenes(ArrayList<ChipSeqPeak> enriched, String [] geneAnnots){
		if(geneAnnots!=null && geneAnnots.length>0){
			for(int g=0; g<geneAnnots.length; g++){
				RefGeneGenerator<Region> rgg = new RefGeneGenerator<Region>(gen, geneAnnots[g]);
                for (ChipSeqPeak peak : enriched) {
                    if (geneOverlap) {
                        Iterator<Gene> geneIter = rgg.execute(peak.coords);
                        while (geneIter.hasNext()) {
                            Gene gene = geneIter.next();
                            int overlap = gene.getOverlapSize(peak.coords);
                            if (peak.nearestGene == null || overlap > peak.distToGene) {
                                peak.nearestGene = gene;
                                peak.distToGene = overlap;
                            }                           
                        }
                        if (peak.nearestGene != null) {
                            peak.distToGene = peak.nearestGene.distance(peak.coords);
                        }
                    } else {
                    	peak.distToGene = maxGeneDistance;
                    	Region query = peak.coords.expand(maxGeneDistance, maxGeneDistance);
                        Iterator<Gene> geneIter = rgg.execute(query);
                        //int midpoint = (peak.coords.getStart() + peak.coords.getEnd()) / 2;
                        while (geneIter.hasNext()) {
                            Gene gene = geneIter.next();
                            //int distance = Math.abs(midpoint - gene.getFivePrime());
                            int distance = Math.abs(peak.peak.getLocation() - gene.getFivePrime());
                            if (distance < peak.distToGene) {
                                peak.nearestGene = gene;
                                peak.distToGene = distance;
                            }                            
                        }
                    }
                }
			}
		}
	}
	/* Add strandedness of peaks (useful for CLIP)*/
	protected void addStrands(ArrayList<ChipSeqPeak> curr, LinkedList<SeqHit> ipHits){
		Iterator<ChipSeqPeak> itr = curr.iterator();
		while(itr.hasNext()){
			ChipSeqPeak peak = itr.next();
			ArrayList<SeqHit> winHits = overlappingHits(ipHits, peak.coords);
			
			//Filter 1: remove things that are ridiculously one-stranded in the IP channel
			double strandprop=forwardReadProportion(winHits);
			if(strandprop>=0.5){
				peak.strand='+';
			}else{
				peak.strand='-';
			}
		}
	}
	/* Trim the peaks back to the coordinates of the first & last read in the peak*/
	protected void trimPeaks(ArrayList<ChipSeqPeak> curr, LinkedList<SeqHit> ipHits){
		Iterator<ChipSeqPeak> itr = curr.iterator();
		while(itr.hasNext()){
			ChipSeqPeak peak = itr.next();
			ArrayList<SeqHit> winHits = overlappingHits(ipHits, peak.coords);
			StrandedRegion min=winHits.get(0);
			StrandedRegion max=winHits.get(0);
			for(StrandedRegion sr : winHits){
				if(sr.getStart()<min.getStart()){min=sr;}
				if(sr.getEnd()>max.getEnd()){max=sr;}
			}
			int startOff = peak.coords.getStart()-min.getStart()<0 ? peak.coords.getStart()-min.getStart(): 0;
			int endOff = max.getEnd()-peak.coords.getEnd()<0 ? max.getEnd()-peak.coords.getEnd():0;
			peak.coords = peak.coords.expand(startOff, endOff);			
		}
	}
	/* Returns the proportion of reads in the region that are in the forward direction*/
	protected double forwardReadProportion(ArrayList<SeqHit> winHits){
		double totalHits =(double)winHits.size();
		double forwardHits =0;
		for(SeqHit s : winHits){
			if(s.getStrand()=='+'){
				forwardHits++;
			}			
		}
		return(forwardHits/totalHits);
	}
	
	/* Find the exact peak locations based on left-right balance of forward/reverse reads */
	protected Point findPeakLRBalance(LinkedList<SeqHit> ipHits, ChipSeqPeak reg){
		ArrayList<SeqHit> winHits = overlappingHits(ipHits, reg.coords);
		int [] forward = new int [reg.coords.getWidth()+1];
		int [] reverse = new int [reg.coords.getWidth()+1];
		for(int s=0; s<=reg.coords.getWidth(); s++){forward[s]=0; reverse[s]=0;}
		
		for(SeqHit r : winHits){
			int i;
			if(r.getStrand()=='+'){
				i=Math.max(0, r.getStart()-reg.coords.getStart());
				forward[i]++;
			}else{
				i=Math.min(r.getEnd()-reg.coords.getStart(),reg.coords.getWidth());
				reverse[i]++;
			}			
		}
		int minBal = 10000, minPos = -1;
		for(int s=0; s<=reg.coords.getWidth(); s++){
			int left=0, right=0;
			for(int l=0; l<=s; l++){left+=forward[l];}
			for(int r=reg.coords.getWidth(); r>s; r--){right+=reverse[r];}
			if(Math.abs(left-right)<minBal){
				minBal =Math.abs(left-right);
				minPos = s;
			}
		}
		Point p = new Point(gen, reg.coords.getChrom(), minPos+reg.coords.getStart());
		return(p);
	}
	
	/* Find the exact peak locations based on maximum overlapping read counts. Careful; make sure the reads are extended... */
	protected Point findPeakMaxHit(LinkedList<SeqHit> ipHits, ChipSeqPeak reg){
		ArrayList<SeqHit> winHits = overlappingHits(ipHits, reg.coords);
		int [] sum = new int [reg.coords.getWidth()+1];
		for(int s=0; s<sum.length; s++){sum[s]=0;}
		
		for(SeqHit r : winHits){
			int start = r.getStart()-reg.coords.getStart(); 
			int stop= r.getEnd()-reg.coords.getStart();
			
			for(int i=start; i<stop; i++){
				if(i>=0 && i<sum.length){
					sum[i]++;
				}
			}
		}
		int max = 0, maxPos = -1;
		for(int s=0; s<sum.length; s++){
			if(sum[s]>max){
				max= sum[s];
				maxPos=s;
			}
		}
		Point p = new Point(gen, reg.coords.getChrom(), maxPos+reg.coords.getStart());
		return(p);
	}
	
	//Set Poisson thresholds
	public int getPoissonThreshold(double threshold, double count, double hitLength, double seqLen, double mappable, double binWidth, double binOffset){
		int countThres=0;
		DRand re = new DRand();
		Poisson P = new Poisson(0, re);
		double pMean = (count*(hitLength/binOffset + binWidth/binOffset))/(seqLen*mappable/binOffset); 
		P.setMean(pMean);
		double l=1;
		for(int b=1; l>threshold; b++){
			l=1-P.cdf(b);
			countThres=b;
		}
		return(Math.max(1,countThres));
	}
	public int getMaxPoissonThreshold(double [] p){
		double max=p[0];
		for(int i=1; i<p.length; i++)
			if(p[i]>max)
				max=p[i];
		return((int)max);
	}
	public void addPoissonModel(int code){
		poissonThresModels.add(code);
		numPoissonModels++;
		ipPoissonThresholds=new double[numPoissonModels];
		backPoissonThresholds=new double[numPoissonModels];
		if(code==-1){//Genome-wide initialization
			if(ipGenomePoissonThres==0){
				ipGenomePoissonThres = getPoissonThreshold(Math.pow(10, poissThres), iphittot, IPhandles.get(0).getExtendedReadLength(), genomeLen, 0.8, winWidth, winStep);
			}ipPoissonThresholds[numPoissonModels-1]= ipGenomePoissonThres;
			if(!noBackground && backGenomePoissonThres==0){
				backGenomePoissonThres = getPoissonThreshold(Math.pow(10, poissThres), backhittot, backhandles.get(0).getExtendedReadLength(), genomeLen, 0.8, winWidth, winStep);
			}backPoissonThresholds[numPoissonModels-1]= backGenomePoissonThres;
		}
	}	
	//Find the Poisson thresholds empirically
	public void empiricalPoissonCalc(double threshold, double mappable){
		System.out.println("Calculating Poission Thresholds");
		
		double totalIPBins=0, totalIPBinCounts=0, totalBackBins=0, totalBackBinCounts=0;
			
		Iterator<NamedRegion> chroms = new ChromRegionIterator(gen);
		while (chroms.hasNext()) {
			NamedRegion currentChrom = chroms.next();
			System.out.println(currentChrom.getChrom());
			
			for(int x=currentChrom.getStart(); x<=currentChrom.getEnd(); x+=100000000){
				int y = x+100000000; 
				if(y>currentChrom.getEnd()){y=currentChrom.getEnd();}
				Region currRegion = new Region(gen, currentChrom.getChrom(), x, y);
				
				LinkedList<SeqHit> ipHits = new LinkedList<SeqHit>();
				LinkedList<SeqHit> backHits = new LinkedList<SeqHit>();
				for(SeqExptHandler IP: IPhandles){
					ipHits.addAll(IP.loadExtendedHits(currRegion));
				}
				for(SeqExptHandler back: backhandles){
					backHits.addAll(back.loadExtendedHits(currRegion));
				}
				double [] ipHitCounts=makeHitLandscape(ipHits, currRegion, (int)ipBasePoissonThres);
				double [] backHitCounts=makeHitLandscape(backHits, currRegion, (int)backBasePoissonThres);
				for(int i=currRegion.getStart(); i<currRegion.getEnd()-winWidth; i+=winStep){
					Region currWin = new Region(gen, currentChrom.getChrom(), i, (int)(i+winWidth-1));
					totalIPBins++;
					totalBackBins++;
					
					int binid = (int)Math.max(0, (double)(currWin.getStart()-currRegion.getStart())/winStep);
					totalIPBinCounts+=(double)ipHitCounts[binid];
					totalBackBinCounts+=(double)backHitCounts[binid];
				}
			}
		}
		
		int countThres=0;
		DRand re = new DRand();
		Poisson P = new Poisson(0, re);
		double pMean=totalIPBinCounts/totalIPBins; System.out.println("IPBins: "+totalIPBins+" IPCounts: "+totalIPBinCounts+" IPMean: "+pMean);
		P.setMean(pMean);
		double l=1;
		for(int b=1; l>threshold; b++){
			l=1-P.cdf(b);
			countThres=b;
		}
		ipGenomePoissonThres = countThres;
		pMean=totalBackBinCounts/totalBackBins; System.out.println("BackBins: "+totalBackBins+" BackCounts: "+totalBackBinCounts+" BackMean: "+pMean);
		P.setMean(pMean);
		l=1;
		for(int b=1; l>threshold; b++){
			l=1-P.cdf(b);
			countThres=b;
		}
		backGenomePoissonThres = countThres;
	}
	
	//Print the distribution of reads around the peaks
	//Used for calculating the fragmentation distribution
	public void printIPMetaPeak(ArrayList<ChipSeqPeak> peaks, int window, int topX){
		int [] forward = new int [window+1];
		int [] reverse = new int [window+1];
		for(int w=0; w<=window; w++){
			forward[w]=0; reverse[w]=0;
		}
		//Assuming the peaks are pre-sorted
		for(SeqExptHandler IP: IPhandles){
			int count =0;
			for(ChipSeqPeak p : peaks){
				if(count<topX){
					int start = p.peak.getLocation()-(window/2)<1 ? 1 : p.peak.getLocation()-(window/2);
					int end = p.peak.getLocation()+(window/2) >= p.peak.getGenome().getChromLength(p.peak.getChrom()) ? p.peak.getGenome().getChromLength(p.peak.getChrom()) : p.peak.getLocation()+(window/2);
					Region r = new Region(p.peak.getGenome(), p.peak.getChrom(), start, end);
					int [] currFor = IP.hitDepthLandscape(r, 0, '+');
					int [] currRev = IP.hitDepthLandscape(r, 0, '-');
					//System.out.println(p.peak.getLocationString());
					//Careful here... assuming that no region is less than window/2 closer to the chromo start
					for(int w=0; w<=window; w++){
						forward[w]+=currFor[w];
						reverse[w]+=currRev[w];
					}
					count++;
				}
			}
		}
		
		//Print the distribs & find maxes
		int maxFor=0, maxRev=0;
		int maxForPos=0, maxRevPos=window;
		System.out.println("Tag distributions around peaks\n");
		System.out.println("RelPos\tForward\tReverse");
		for(int w=0; w<=window; w++){
			int rel = w-(window/2);
			System.out.println(rel+"\t"+forward[w]+"\t"+reverse[w]);
			if(forward[w]>maxFor && rel<-30){maxFor=forward[w]; maxForPos=w;}
			if(reverse[w]>maxRev && rel>30){maxRev=reverse[w];maxRevPos=w;}			
		}int diff = (maxRevPos-maxForPos)/2;
		System.out.println("Forward-Reverse Shift:\t"+diff);
	}
	
	public static void printError(){
		System.err.println("Usage:\n " +
                "ChipSeqPeakFinder \n" +
                "--species <organism name;genome version> "+
                "--expt <solexa expt> " +
                "--back <background expt> \n"+
                "--outpeak <output file> "+
                "--outseq <output file> "+
                "--seqwin <sequence length> \n"+
                "--binwidth <width of bins> "+
                "--binstep <offset of bins> "+
                "--minz <minimum z score> \n"+
                "--poisson <Poisson threshold (10 to the power of)> " +
                "--dynpoisson <dynamic Poisson threshold (Xbp)> \n" +
                "--readlen <length> " +
                "--readextend <extension> \n"+
                "--allowtowers [don't filter towers] " +
                "--allowneedles [don't filter needles] " +
                "--printpeakdistrib <number of top peaks to print meta-peak with> \n" +
                "--genes <name of gene annotation to examine> " +
                "--scangenesonly <only look for peaks within genes (CLIP-seq specific)> " +
                "--maxgenedist <distance from gene> " +
                "--geneOverlap [only look at genes that overlap a peak] \n" +
                "--shifttags <distance to shift tags (overrides extension)> \n");
	}
}


