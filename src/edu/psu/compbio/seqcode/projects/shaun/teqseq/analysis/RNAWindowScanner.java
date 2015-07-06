package edu.psu.compbio.seqcode.projects.shaun.teqseq.analysis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Species;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.AlignBlock;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.AlignHit;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.GenomeLoader;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.GenomeSegmenter;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.ReadLoader;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.ReadUtils;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.SAMReadLoader;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.TUnit;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.TUnitGene;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.core.TUnitJunction;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.exptprops.ObsExpSeqBiasModel;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.exptprops.SeqBiasModel;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.geneloaders.AGene;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.geneloaders.AIsoform;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.geneloaders.ARegion;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.geneloaders.AnnotationFilter;
import edu.psu.compbio.seqcode.projects.shaun.teqseq.geneloaders.GTFAnnotationLoader;

/**
 * RNAWindowScanner: The aim of this class is to take a collection of experiments, scan each one independently using a sliding window 
 * to pick out regions and junctions of read enrichment, and combine all such regions found across all experiments (with gene annotation) 
 * into a single set of non-overlapping TUnits. 
 * The execute method also scans the experiment collection again using the discovered TUnits, counting enrichment in each.  
 *   
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class RNAWindowScanner {

	//Basics
	protected String outName;
	protected List<ReadLoader> loaders;
	protected Map<String, List<Region>> scanRegions;
	protected GenomeLoader gLoad;
	protected Collection<AGene> genes;
	protected Map<ReadLoader, SeqBiasModel> biasModels = new HashMap<ReadLoader, SeqBiasModel>();
	protected static boolean correctSeqBias=false;
	protected Map<ReadLoader, Integer> replicateIndex = new HashMap<ReadLoader, Integer>();
	//Processing
	protected double binWidth=50;
	protected double binStep = 25;
	protected double logCovThres=-9; 
	protected double genomeLen=0;
	protected double mappableGenome=0.85;
	protected double readLength=50;
	protected double connectionThreshold = 2;
	protected double junctionCoverageThreshold = 3;
	protected Map<String, Double> coverageThresholds = new HashMap<String, Double>();
	protected Map<String, Double> replicateCounts = new HashMap<String, Double>();
	protected int tUnitID=0;
	protected int tJuncID=0;
	protected int totalTUnits=0, totalTUnitJunctions=0, totalTUnitGenes=0;
	//Output files
	protected HashMap<String, FileWriter> fileHandles = new HashMap<String, FileWriter>();
	
	
	/**
	 * Constructor
	 * @param gl GenomeLoader
	 * @param expts ExptCollection
	 * @param genes Collection of known genes
	 */
	public RNAWindowScanner(String name, GenomeLoader gl, List<ReadLoader> readLoaders, Collection<AGene> genes, Map<String, List<Region>> regionsOfInterest){
		outName = name;
		loaders = readLoaders;
		scanRegions = regionsOfInterest;
		gLoad = gl;
		genomeLen = gLoad.getGenome().getGenomeLength();
		this.genes = genes;
		for(AGene g : genes){
			g.initializeCoverageCounter(loaders.size());
			for(AIsoform i : g.getIsoforms())
				i.initializeCoverageCounter(loaders.size());
		}
		
		int id=0;
		for(ReadLoader reads : loaders){
			replicateIndex.put(reads, id);
			SeqBiasModel seqBias=null;
			if(correctSeqBias){
				seqBias = new ObsExpSeqBiasModel(gLoad, 6, regionsOfInterest);
				if(!seqBias.attemptLoad(reads)){
					if(gLoad.seqIsAvailable()){
						System.err.println("ExptReplicate: Cannot find sequencing bias model for "+reads.getSourceName()+", training a new model.");
						seqBias.execute(reads, genes);
					}else{
						System.err.println("ExptReplicate: Cannot find sequencing bias model for "+reads.getSourceName()+" and genome sequence was not provided for training the model.");
						seqBias=null;
					}
				}
			}
			biasModels.put(reads, seqBias);
			id++;
		}
		
		//Initialize file handles
		openOutputFiles();
		
		initializeCoverageThresholds();
	}
	
	
	/**
	 * Execute the scanner. 
	 * For each region of interest:
	 *   For each experimental replicate:
	 *     Load alignment hits
	 *     Initialize a read coverage array & record junctions
	 *     Use a sliding window across read coverage, looking for regions with significant enrichment
	 *     Combine & trim windows
	 *   Combine all regions/junctions found with each other and gene annotated borders
	 *   While you have the hits hanging around, count coverage
	 * Print to files
	 * Done
	 */
	public void execute(int maxThreads){
		Thread[] threads = new Thread[maxThreads];
		Map[] threadRegions = new HashMap[maxThreads]; 
		int i = 0;
        for (i = 0 ; i < threads.length; i++) 
            threadRegions[i] = new HashMap<String,List<Region>>();

        for(String chr : scanRegions.keySet())
       		threadRegions[(i++) % maxThreads].put(chr, scanRegions.get(chr));
        
        for(i = 0 ; i < threads.length; i++) {
            Thread t = new Thread(new WindowScannerThread(threadRegions[i]));
            t.start();
            threads[i] = t;
        }
        boolean anyrunning = true;
        while (anyrunning) {
            anyrunning = false;
            try {
                Thread.sleep(5000);
            } catch (InterruptedException e) { }
            for (i = 0; i < threads.length; i++) {
                if (threads[i].isAlive()) {
                    anyrunning = true;
                    break;
                }
            }
        }
        //Close the output files
        closeOutputFiles();
        
        System.err.println("\nFinal results: "+totalTUnits+" TUnits, "+totalTUnitJunctions+" TUnitJunctions, "+totalTUnitGenes+" TUnitGenes, "+genes.size()+" AnnotatedGenes.");
	}
	
	/**
	 * WindowScannerThread: Thread to process reads in a given region
	 * In this version, hits for a given region are loaded via random access one experiment at a time. 
	 * We can run out of memory when hits for all  experiments are loaded simultaneously (as in the previous version).
	 * Unfortunately, we need to load the reads three distinct times for each experiment in this version. 
	 * @author Shaun Mahony
	 * @version	%I%, %G%
	 */
	class WindowScannerThread implements Runnable{
		Map<String, List<Region>> regions;
		//Results
		List<TUnit> units;
		List<TUnitJunction> junctions;
		List<TUnitGene> unitGenes;
		List<AGene> annoGenes;
		
		public WindowScannerThread(Map<String, List<Region>> r){regions = r;}
		
		public void run(){
			for(String currChr : regions.keySet()){
				units = new ArrayList<TUnit>();
				junctions = new ArrayList<TUnitJunction>();
				unitGenes = new ArrayList<TUnitGene>();
				annoGenes = new ArrayList<AGene>();
				
				for(Region currRegion : regions.get(currChr)){ //Iterate through regions
					int allHitsInRegion=0; //Counter for progress reporting only
					List<TUnit> currUnitsAllExpts = new ArrayList<TUnit>();
					List<TUnitJunction> currJunctionsAllExpts = new ArrayList<TUnitJunction>();
					
					//1. Find TUnits & TUnitJunctions in each experiment in the current region
					for(ReadLoader rl : loaders){
						String eName = rl.getSourceName();
						List<AlignHit> hits = rl.getOverlappingHits(currRegion, biasModels.get(rl));
						allHitsInRegion+= hits.size();
						if(hits.size()>0){//This probably only comes into play in test sets
							
						ArrayList<TUnit> currUnits = new ArrayList<TUnit>();
						double coverageThreshold = coverageThresholds.get(eName);
						    
						//Get a hit landscape
						double [] stackedHitCounts = makeHitLandscape(hits, currRegion, '.');
					
					    //Scan regions for enriched regions
					    int currBin=0;
					    TUnit lastUnit = null;
					    for(int i=currRegion.getStart(); i<currRegion.getEnd()-(int)binWidth; i+=(int)binStep){
							double winHits=stackedHitCounts[currBin];
							if(winHits>=coverageThreshold){
							    Region currWin = new Region(gLoad.getGenome(), currRegion.getChrom(), i, (int)(i+binWidth-1));
							    lastUnit = addUnit(currUnits, lastUnit, currWin, '.');
							}
							currBin++;
					    }
					    
					    //Trim the discovered TUnits
					    trimUnits(currUnits, hits, '.');
					    //Add to the collection
					    currUnitsAllExpts.addAll(currUnits);
					    
					    //Add all junction-mapped locations
					    List<TUnitJunction> currJunctions = findJunctions(hits, '.');
					    for(TUnitJunction tuj : currJunctions)
						if(tuj.getHitWeight()>junctionCoverageThreshold)
						    currJunctionsAllExpts.add(tuj);
					    
					    hits=null; currUnits=null; currJunctions=null; System.gc();
						}
					}
					if(allHitsInRegion>0){//This probably only comes into play in test sets
						Collections.sort(currUnitsAllExpts);
						Collections.sort(currJunctionsAllExpts);
						currJunctionsAllExpts = removeRedundantJunctions(currJunctionsAllExpts);
						
						//Get genes that overlap this region
						List<AGene> currRegGenes = AnnotationFilter.getOverlappingGenes(currRegion, genes);
						for(AGene g : currRegGenes)
							if(!annoGenes.contains(g))
								annoGenes.add(g);

						//2. Combine all TUnits with TUnitJunctions and gene annotations such that you are left with a set of non-redundant non-overlapping TUnits
						currUnitsAllExpts = combineTUnits(currUnitsAllExpts, currJunctionsAllExpts, currRegGenes);
						
						//3. Go back through the filtered TUnits to quantify coverage & startCounts in each experiment.
						for(ReadLoader rl : loaders){
							List<AlignHit> hits = rl.getOverlappingHits(currRegion, biasModels.get(rl));
							
							unitCoverageAndStarts(rl, currUnitsAllExpts, hits, '.');
							junctionCoverage(rl, currJunctionsAllExpts, hits, '.');
							hits=null; System.gc();
						}
						currUnitsAllExpts = filterLowCoverageUnits(currUnitsAllExpts, coverageThresholds);
						currJunctionsAllExpts = filterLowCoverageJunctions(currJunctionsAllExpts, coverageThresholds);
						
						//4. Quantify annotated gene hit counts and add annotations to the TUnits
						//We don't want to double count reads that span multiple exons here.
						//Do this by finding TUnits that overlap a gene's exons and adding the weights of reads that start over such units
						addCountsToAnnotatedExons(currRegGenes, currUnitsAllExpts, currJunctionsAllExpts, '.');
		
						//5. Define connected TUnits and merge into TUnitGenes
						double[][] connectionMatrix = new double[currUnitsAllExpts.size()][currUnitsAllExpts.size()];
						for(int x=0; x<currUnitsAllExpts.size(); x++){for(int y=0; y<currUnitsAllExpts.size(); y++){connectionMatrix[x][y]=0.0;}}
						for(ReadLoader rl : loaders){
							List<AlignHit> hits = rl.getOverlappingHits(currRegion, biasModels.get(rl));
							connectionMatrix = connectTUnitsWithHits(connectionMatrix, currUnitsAllExpts, hits, '.');
							//TODO: Connect TUnits with correlation across experiments?
							hits=null; System.gc();
						} 
						
						//6. Finally, add the filtered sets to the stored versions, filtering out any low coverage units that crept in during combination
						// and merge TUnits into genes (don't forget to add together starts).
						unitGenes.addAll(makeTUnitGenes(currUnitsAllExpts, connectionMatrix));
						units.addAll(currUnitsAllExpts);
						junctions.addAll(currJunctionsAllExpts);
					}//Progress
					System.err.println("WindowScannerThread: processing "+currRegion+"\t"+allHitsInRegion+" hits");
					currUnitsAllExpts=null; currJunctionsAllExpts=null; System.gc(); 
				}
				//Sort lists
				Collections.sort(units);
				Collections.sort(junctions);
				Collections.sort(unitGenes);
				
				//Process overhanging reads
				unitGenes = findBoundaryConnections(unitGenes);
				
				//Final Statistics
				System.err.println(currChr+" results: "+units.size()+" TUnits, "+junctions.size()+" TUnitJunctions, "+unitGenes.size()+" TUnitGenes, "+annoGenes.size()+" AnnotatedGenes.");
				
				//Report to files
				synchronized(fileHandles){
					printOutputFiles(units, junctions, unitGenes, annoGenes);
					totalTUnits+=units.size();
					totalTUnitJunctions+=junctions.size();
					totalTUnitGenes += unitGenes.size();
				}
				units=null; junctions=null; unitGenes=null; annoGenes=null;  System.gc(); 
			}//End of Chr loop
		}//End of run()		
	}
	
	/**
	 * Open output files
	 */
	public void openOutputFiles(){
		try {
			String outFileName = outName+".tunits.coverage";
			FileWriter outFW = new FileWriter(outFileName);
			outFW.write("#Coords\tWidth");
			for(ReadLoader rl : loaders){outFW.write("\t"+rl.getConditionName()+":"+rl.getSourceName());}
			outFW.write("\tAnnotations\n");
			fileHandles.put("TUnits", outFW);
			
			outFileName = outName+".tunitjunctions.coverage";
			outFW = new FileWriter(outFileName);
			outFW.write("#Coords\tWidth");
			for(ReadLoader rl : loaders){outFW.write("\t"+rl.getConditionName()+":"+rl.getSourceName());}
			outFW.write("\tAnnotations\n");
			fileHandles.put("TUnitJunctions", outFW);
			
			outFileName = outName+".tunitgenes.counts";
			outFW = new FileWriter(outFileName);
			outFW.write("#Coords\tWidth");
			for(ReadLoader rl : loaders){outFW.write("\t"+rl.getConditionName()+":"+rl.getSourceName());}
			outFW.write("\tAnnotations\n");
			fileHandles.put("TUnitGenes", outFW);
			
			outFileName = outName+".tunitgenes.gff";
			outFW = new FileWriter(outFileName);
			fileHandles.put("TUnitGenesGFF", outFW);
			
			outFileName = outName+".annogenes.counts";
			outFW = new FileWriter(outFileName);
			outFW.write("#Coords\tWidth");
			for(ReadLoader rl : loaders){outFW.write("\t"+rl.getConditionName()+":"+rl.getSourceName());}
			outFW.write("\tAnnotation\n");
			fileHandles.put("AnnoGenes", outFW);
			
			outFileName = outName+".annotrans.counts";
			outFW = new FileWriter(outFileName);
			outFW.write("#Coords\tWidth");
			for(ReadLoader rl : loaders){outFW.write("\t"+rl.getConditionName()+":"+rl.getSourceName());}
			outFW.write("\tAnnotation\n");
			fileHandles.put("AnnoTrans", outFW);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	/**
	 * Print output files
	 */
	public void printOutputFiles(List<TUnit> units, List<TUnitJunction> junctions, List<TUnitGene> unitGenes, List<AGene> currGenes){
		//TUnits
		try {
			FileWriter outFW = fileHandles.get("TUnits");
			for(TUnit u : units){
				outFW.write(u.getCoords().getLocationString()+"\t"+u.getCoords().getWidth());
				for(ReadLoader rl : loaders){
					outFW.write("\t"+String.format("%.2f", u.getCoverage(replicateIndex.get(rl))));
				}
				outFW.write("\t"+u.getAnnotationString()+"\n");
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//TUnitJunctions
		try {
			FileWriter outFW = fileHandles.get("TUnitJunctions");
			for(TUnitJunction u : junctions){
				outFW.write(u.getCoords().getLocationString()+"\t"+u.getCoords().getWidth());
				for(ReadLoader rl : loaders){
					outFW.write("\t"+String.format("%.2f", u.getCoverage(replicateIndex.get(rl))));
				}
				outFW.write("\t"+u.getAnnotationString()+"\n");
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//TUnitGenes
		try {
			FileWriter outFW = fileHandles.get("TUnitGenes");
			for(TUnitGene ug : unitGenes){
				outFW.write(ug.getCoords().getLocationString()+"\t"+ug.getCoords().getWidth());
				for(ReadLoader rl : loaders){
					outFW.write("\t"+String.format("%.2f", ug.getHitStarts(replicateIndex.get(rl))));
				}
				outFW.write("\t"+ug.getAnnotationString()+"\n");
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//TUnitGenes GFF
		try {
			FileWriter outFW = fileHandles.get("TUnitGenesGFF");
			for(TUnitGene ug : unitGenes){
				outFW.write(ug.getCoords().getChrom()+"\twindow_scan\tgene\t"+ug.getCoords().getStart()+"\t"+ug.getCoords().getEnd()+"\t.\t.\t.\tID:"+ug.getAnnotationString()+"\n");
				
				for(TUnit u : ug.getComponents()){
					outFW.write(u.getCoords().getChrom()+"\twindow_scan\ttunit\t"+u.getCoords().getStart()+"\t"+u.getCoords().getEnd()+"\t.\t.\t.\tID:"+u.getAnnotationString()+"\n");
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//Annotated genes
		try {
			FileWriter outFW = fileHandles.get("AnnoGenes");
			for(AGene g : currGenes){
				outFW.write(g.getCoords()+"\t"+g.getCoords().getWidth());
				for(ReadLoader rl : loaders){
					outFW.write("\t"+String.format("%.2f", g.getCoverage(replicateIndex.get(rl))));
				}outFW.write("\t"+g.getID()+":"+g.getName()+"\n");
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		//Annotated transcripts
		try {
			FileWriter outFW = fileHandles.get("AnnoTrans");
			for(AGene g : currGenes){
				for(AIsoform i : g.getIsoforms()){
					outFW.write(i.getCoords()+"\t"+i.getCoords().getWidth());
					for(ReadLoader rl : loaders){
						outFW.write("\t"+String.format("%.2f", i.getCoverage(replicateIndex.get(rl))));
					}outFW.write("\t"+i.getID()+":"+i.getName()+"\n");
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Close the output files
	 */
	protected void closeOutputFiles(){
		try {
			for(FileWriter o : fileHandles.values())
				o.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * Count all (weighted) reads in each replicate and use those counts to initialize coverage thresholds
	 * @param expts
	 */
	protected void initializeCoverageThresholds(){
		for(String chr : scanRegions.keySet()){
			for(Region reg : scanRegions.get(chr)){
				int numHits=0;
				for(ReadLoader rl : loaders){
					String eName = rl.getSourceName();
					List<AlignHit> hits = rl.getOverlappingHits(reg, biasModels.get(rl));
					if(!replicateCounts.containsKey(eName))
						replicateCounts.put(eName, 0.0);
						
					double currW = replicateCounts.get(eName);
					for(AlignHit h : hits)
						currW += h.getWeight();	
					replicateCounts.put(eName, currW);
					numHits+=hits.size();
				}
				System.err.println("Coverage threshold calc: "+reg.getLocationString()+"\t"+numHits+" hits");
			}
			System.gc();
		}
		
		for(ReadLoader rl : loaders){
			coverageThresholds.put(rl.getSourceName(), calcCoverageThreshold(replicateCounts.get(rl.getSourceName())));
			System.err.println("Coverage threshold for "+rl.getSourceName()+"\t"+coverageThresholds.get(rl.getSourceName()));
		}
	}
	
	/**
	 * Makes integer arrays corresponding to the read landscape over the current region
	 * @param hits List of AlignHits
	 * @param currReg	Current region under examination
	 * @param strand char; '+' to only examine positive reads, '-' for negative strand, '.' for both
	 * @return
	 */
    protected double[] makeHitLandscape(List<AlignHit> hits, Region currReg, char strand){
        int numBins = (int)(currReg.getWidth()/binStep);
        double[] landscape = new double[numBins+1];
        for(int i=0; i<=numBins; i++){landscape[i]=0; }
        for(AlignHit h : hits){
            if(strand=='.' || h.getStrand()==strand){
            	for(AlignBlock a : h.getAlignmentBlocks()){
            		if(a.overlaps(currReg)){
	            		int offset=inBounds(a.getReferenceStart()-currReg.getStart(),0,currReg.getWidth());
	            		int binstart = inBounds((int)((double)offset/binStep), 0, numBins);
	            		int binend = inBounds((int)((double)((a.getReferenceEnd())-currReg.getStart())/binStep), 0, numBins);
	            		for(int i=binstart; i<=binend; i++){
	            			landscape[i]+=h.getWeight();
	            		
	            		}
	            	}
            	}
            }
        }
        return(landscape);
    }
    
    /**
	 * Merge connected TUnits into TUnitGenes
	 * @param currUnits
	 * @param connections
	 * @return
	 */
	protected List<TUnitGene> makeTUnitGenes(List<TUnit> currUnits, double[][] connections) {
		List<TUnitGene> currTGenes = new ArrayList<TUnitGene>();
		
		HashMap<TUnit, TUnitGene> geneMap = new HashMap<TUnit, TUnitGene>();
		int geneID=0;
		for(int x=0; x<currUnits.size(); x++){
			if(!geneMap.containsKey(currUnits.get(x))){
				TUnitGene currGene = new TUnitGene(currUnits.get(x), geneID);
				geneMap.put(currUnits.get(x), currGene);
				recursivelyAddTUnitsToGenes(currGene, currUnits, geneMap, x, connections);
				currGene.sortComponents();
				currTGenes.add(currGene);
				geneID++;
			}
		}		
		return currTGenes;
	}
	
	/**
	 * Recursively add TUnits to the current gene
	 * @param currGene
	 * @param currUnits
	 * @param geneMap
	 * @param currX
	 * @param connections
	 */
	private void recursivelyAddTUnitsToGenes(TUnitGene currGene, List<TUnit> currUnits, HashMap<TUnit, TUnitGene> geneMap, int currX, double[][] connections){
		for(int y=0; y<currUnits.size(); y++){
			if(currX!=y && connections[currX][y]>=connectionThreshold){
				if(!geneMap.containsKey(currUnits.get(y))){
					currGene.addComponent(currUnits.get(y));
					geneMap.put(currUnits.get(y), currGene);
					recursivelyAddTUnitsToGenes(currGene, currUnits, geneMap, y, connections);
				}
			}
		}
	}
	
    /**
     * Find junction-mapping reads in a set of AlignHits
     * @param hits
     * @param strand
     * @return
     */
    protected List<TUnitJunction> findJunctions(List<AlignHit> hits, char strand){
    	List<TUnitJunction> juncs = new ArrayList<TUnitJunction>();
    	for(AlignHit h : hits){
    		if(h.getAlignmentBlocks().length>1 && (strand=='.' || h.getStrand()=='.')){
    			AlignBlock[] blocks = h.getAlignmentBlocks();
    			for(int i=0; i<blocks.length-1; i++){
    				TUnitJunction currJunc = new TUnitJunction(h.getGenome(), h.getChrom(), blocks[i].getReferenceEnd(), blocks[i+1].getReferenceStart(), tJuncID, strand);
    				currJunc.addHit(h.getWeight());
    				tJuncID++;
    				boolean jFound = false;
    				for(TUnitJunction tuj : juncs){
    					if(tuj.overlaps(currJunc)){
    						tuj.addHit(h.getWeight());
    						jFound=true;
    					}
    				}
    				if(!jFound)
    					juncs.add(currJunc);
    			}
    		}
    	}
    	Collections.sort(juncs);
    	return(juncs);
    }
    
    /**
     * Remove equivalent junctions
     * @param juncs
     * @return
     */
    protected List<TUnitJunction> removeRedundantJunctions(List<TUnitJunction> juncs){
    	List<TUnitJunction> filtered = new ArrayList<TUnitJunction>();
    	TUnitJunction lastAdded=null;
    	for(TUnitJunction j : juncs){
    		if(lastAdded==null || !j.overlaps(lastAdded)){
    			filtered.add(j);
    			j.initializeCounters(loaders.size());
    			lastAdded = j;
    		}
    	}
    	return(filtered);
    }
    
    /**
     * Remove TUnits that have coverage below a threshold
     * @param currUnits
     * @param expts
     * @param coverageThres
     * @return
     */
    protected List<TUnit> filterLowCoverageUnits(List<TUnit> currUnits, Map<String,Double> coverageThres){
    	List<TUnit> filtered = new ArrayList<TUnit>();
    	for(TUnit u : currUnits){
			boolean suffCov = false;
			for(ReadLoader rl : loaders){
				suffCov = suffCov|(u.getCoverage(replicateIndex.get(rl))>coverageThres.get(rl.getSourceName()));
			}
			if(suffCov)
				filtered.add(u);
		}
    	return filtered;
    }
    
    /**
     * Remove TUnitJunctions that have coverage below a threshold
     * @param currUnits
     * @param expts
     * @param coverageThres
     * @return
     */
    protected List<TUnitJunction> filterLowCoverageJunctions(List<TUnitJunction> currJuncs, Map<String,Double> coverageThres){
    	List<TUnitJunction> filtered = new ArrayList<TUnitJunction>();
    	for(TUnitJunction u : currJuncs){
			boolean suffCov = false;
			for(ReadLoader rl : loaders){ 
				suffCov = suffCov|(u.getCoverage(replicateIndex.get(rl))>coverageThres.get(rl.getSourceName()));
			}
			if(suffCov)
				filtered.add(u);
		}
    	return filtered;
    }
    
    
    /**
     * Combine all TUnits with TUnitJunctions and gene annotations such that you are left with a set of non-redundant non-overlapping TUnits
     * The easiest (but probably dumbest) way to do this is by painting the edges and regions of coverage onto an array representing the region
     * Assumes all lists are sorted
     * @param tunits
     * @param juncs
     * @param currGenes
     * @return
     */
    protected List<TUnit> combineTUnits(List<TUnit> tunits, List<TUnitJunction> juncs, List<AGene> currGenes){
    	List<TUnit> filtered = new ArrayList<TUnit>();
    	
    	for(int i=0; i<tunits.size(); i++){
    		//Get a region and set of overlapping TUnits
    		List<TUnit> currUnits = new ArrayList<TUnit>();
    		currUnits.add(tunits.get(i));
    		Region currReg = tunits.get(i).getCoords();
    		while(i<tunits.size()-1 && tunits.get(i+1).getCoords().overlaps(currReg)){
    			currUnits.add(tunits.get(i+1));
    			currReg = currReg.expand(0, tunits.get(i+1).getCoords().getEnd()-currReg.getEnd());
    			i++;
    		}
    		
    		//Get Junctions and Genes in this region
    		List<TUnitJunction> olapJuncs = overlappingJunctions(currReg, juncs, '.');
    		List<AGene> olapGenes = AnnotationFilter.getOverlappingGenes(currReg, currGenes);
    		
    		//Paint an array for these units
    		int[] covered = new int[currReg.getWidth()+1];
    		int[] edges = new int[currReg.getWidth()+1];
    		for(int c=0; c<covered.length; c++){
    			covered[c]=0;
    			edges[c]=0;
    		}
    		for(TUnit t : currUnits){
    			for(int c=inBounds(t.getCoords().getStart()-currReg.getStart(),0,covered.length-1); c<=inBounds(t.getCoords().getEnd()-currReg.getStart(),0,covered.length-1); c++)
    				covered[c]=1;
    			edges[inBounds(t.getCoords().getStart()-currReg.getStart(),0,covered.length-1)]=1;
    			edges[inBounds(t.getCoords().getEnd()-currReg.getStart(),0,covered.length-1)]=1;
    		}
    		for(TUnit j : olapJuncs){
    			edges[inBounds(j.getCoords().getStart()-currReg.getStart(),0,covered.length-1)]=1;
    			edges[inBounds(j.getCoords().getEnd()-currReg.getStart(),0,covered.length-1)]=1;
    		}
    		for(AGene a : olapGenes){
    			for(AIsoform iso : a.getIsoforms()){
    				for(ARegion e : iso.getExons()){
    					edges[inBounds(e.getCoords().getStart()-currReg.getStart(),0,covered.length-1)]=1;
    	    			edges[inBounds(e.getCoords().getEnd()-currReg.getStart(),0,covered.length-1)]=1;
    				}
    			}
    		}
    		
    		//Regions painted, now scan for non-overlapping TUnits
    		for(int x=0; x<currReg.getWidth(); x++){
    			if(edges[x]==1 && covered[x]==1){
    				//Start recording a possible TUnit
    				int y=x+1;
    				while(y<currReg.getWidth() && edges[y]==0 && covered[y]==1){y++;}
    				if(edges[y]==1 && covered[y]==1){
    					//New TUnit
    					TUnit unit = new TUnit(currReg.getGenome(), currReg.getChrom(), currReg.getStart()+x, currReg.getStart()+y-1, tUnitID, '.');
    					filtered.add(unit);
    					unit.initializeCounters(loaders.size());
    					tUnitID++;
    				}
    			}
    		}    		
    	}
    	Collections.sort(filtered);
    	return filtered;
    }
    
    /**
     * Add a peak to the pile
     * @param currres
     * @param last
     * @param currWin
     * @param strand
     * @return
     */
	protected TUnit addUnit(Collection<TUnit> currres, TUnit last, Region currWin, char strand){
		TUnit res =null;
		//Is this hit close to the previously added one?
		if(last!=null && (currWin.getStart() - last.getCoords().getEnd())<=binWidth){
			last.setCoords(new Region(gLoad.getGenome(), last.getCoords().getChrom(), last.getCoords().getStart(), currWin.getEnd()-1));
			res = last;
		}else{
			res = new TUnit(gLoad.getGenome(), currWin.getChrom(), currWin.getStart(), currWin.getEnd(), tUnitID, strand);
			tUnitID++;
			currres.add(res);
		}return(res);
	}
	
	/**
	 *  Trim the units back to the coordinates of the first & last align blocks 
	 *  hits must be sorted
	 *  Also counts the overlapping hits
	 */
	protected void trimUnits(List<TUnit> tunits, List<AlignHit> hits, char str){
		for(TUnit u : tunits){
			int min=u.getCoords().getStart();
			int max=u.getCoords().getEnd();
			
			//Similar to overlappingHits here, but there is no need to store in a sublist
			int i = binarySearchAlignHitBlocks(u.getCoords(), hits); 
			while(i<hits.size() && hits.get(i).getStart()<=u.getCoords().getEnd()){
				if(hits.get(i).blockOverlaps(u.getCoords()) && (str=='.' || hits.get(i).getStrand()==str)){
					for(AlignBlock b : hits.get(i).getAlignmentBlocks()){
						if(b.overlaps(u.getCoords())){
							if(b.getReferenceStart()<min)
								min = b.getReferenceStart();
							if(b.getReferenceEnd()>max)
								max = b.getReferenceEnd();
						}
					}
					//Quantify
					u.addHit(hits.get(i).getWeight());
				}
				i++;
			}
			int startOff = u.getCoords().getStart()-min;
			int endOff = max-u.getCoords().getEnd();

			u.setCoords(u.getCoords().expand(startOff, endOff));
		}
	}
	
	/**
	 * Add weighted coverage value for a defined experimental replicate to each TUnit using a list of align hits .
	 * @param expts ExptCollection (only used to find index of replicate)
	 * @param repName String replciate name (only used to find coverage array index)
	 * @param tunits List of TUnits
	 * @param hits List of AlignHits
	 * @param str char strand
	 */
	public void unitCoverageAndStarts(ReadLoader rLoader, List<TUnit> tunits, List<AlignHit> hits, char str){
		int repID = replicateIndex.get(rLoader);
		for(TUnit u : tunits){
			List<AlignHit> currBoundaryHits = new ArrayList<AlignHit>();
			Region currReg = u.getCoords();
			//Similar to overlappingHits here, but there is no need to store in a sublist
			int i = binarySearchAlignHitBlocks(currReg, hits); 
			while(i<hits.size() && hits.get(i).getStart()<=currReg.getEnd()){
				if(hits.get(i).blockOverlaps(currReg) && (str=='.' || hits.get(i).getStrand()==str)){
					//Add coverage contribution
					u.addCoverage(repID, hits.get(i).getWeight());
					//Add a start
					if(hits.get(i).getFivePrime()>=currReg.getStart() && hits.get(i).getFivePrime()<=currReg.getEnd())
						u.addHitStart(repID, hits.get(i).getWeight());
				}
				//Deal with boundary-overlapping reads
				if((!currReg.contains(hits.get(i)) && hits.get(i).getEnd()>currReg.getEnd()) || 
						(hits.get(i).getPairedRead()!=null && hits.get(i).getPairedRead().getChrom().equals(currReg.getChrom()) && hits.get(i).getPairedRead().getLocation()>currReg.getEnd()))
					currBoundaryHits.add(hits.get(i));
				i++;
			}
			if(currBoundaryHits.size()>=connectionThreshold)
				for(AlignHit h : currBoundaryHits)
					u.addBoundaryHit(h);
		}
	}
	
	/**
	 * Connect TUnits that share align hits
	 * @param expts
	 * @param repName
	 * @param tunits
	 * @param hits
	 * @param str
	 */
	public double[][] connectTUnitsWithHits(double[][] connections, List<TUnit> tunits, List<AlignHit> hits, char str){
		double[][] connect = connections;
		for(int x=0; x<tunits.size(); x++){
			TUnit u = tunits.get(x);
			Region regU = u.getCoords();
			//Similar to overlappingHits here, but there is no need to store in a sublist
			int i = binarySearchAlignHitBlocks(regU, hits); 
			while(i<hits.size() && hits.get(i).getStart()<=regU.getEnd()){
				//If the hit is overlapping the block of interest
				if(hits.get(i).blockOverlaps(regU) && (str=='.' || hits.get(i).getStrand()==str)){
					
					//Does this same hit overlap other TUnits to the right?
					int maxCoord = hits.get(i).getPairedRead()==null ?
							hits.get(i).getEnd() : Math.max(hits.get(i).getEnd(), hits.get(i).getPairedRead().getLocation());
					for(int y=x+1; y<tunits.size(); y++){
						TUnit v = tunits.get(y);
						Region regV = v.getCoords();
						  
						if(regV.getStart()>maxCoord)
							y=tunits.size();
						else{
							if(hits.get(i).blockOverlaps(regV) || 
									(hits.get(i).getPairedRead()!=null && regV.contains(hits.get(i).getPairedRead()))){
								connect[x][y]+=hits.get(i).getWeight();
								connect[y][x]+=hits.get(i).getWeight();
							}
						}
					}
				}
				i++;
			}
		}
		return(connect);
	}
	
	/**
	 * Add weighted coverage value for a defined experimental replicate to each TUnitJunction using a list of align hits  
	 * @param expts ExptCollection (only used to find index of replicate)
	 * @param repName String replciate name (only used to find coverage array index)
	 * @param juncs List of TUnitsJunctions
	 * @param hits List of AlignHits
	 * @param str char strand
	 */
	public void junctionCoverage(ReadLoader rLoader, List<TUnitJunction> juncs, List<AlignHit> hits, char str){
		int repID = replicateIndex.get(rLoader);
		for(TUnitJunction t : juncs){
			Region currReg = t.getCoords();
			//Similar to overlappingHits here, but there is no need to store in a sublist
			int i = binarySearchAlignHitBlocks(currReg, hits); 
			while(i<hits.size() && hits.get(i).getStart()<=currReg.getEnd()){
				if(hits.get(i).blockOverlaps(currReg) && (str=='.' || hits.get(i).getStrand()==str)){
					//Check for the exact overlap
					AlignBlock[] bs = hits.get(i).getAlignmentBlocks();
					for(int b=0; b<bs.length-1; b++){
						if(t.getCoords().getStart()==bs[b].getReferenceEnd() && t.getCoords().getEnd()==bs[b+1].getReferenceStart())
							t.addCoverage(repID, hits.get(i).getWeight());
					}
				}
				i++;
			}
		}
	}
	
	/**
	 * Add weighted start count for a defined experimental replicate to each ARegion using start counts stored in the TUnits.
	 * For protein_coding genes, only CDSs are examined. For other genes, all annotated "exons" are examined.
	 * This method also adds gene name annotations to TUnits & TUnitJunctions.
	 * @param expts ExptCollection (only used to find index of replicate)
	 * @param repName String replicate name (only used to find coverage array index)
	 * @param geneList List of AGenes to query
	 * @param tunits List of TUnits
	 * @param tjuncs List of TUnitJunctions
	 * @param hits List of AlignHits
	 * @param str char strand
	 * @param addAnnotation boolean add gene annotations to the TUnits 
	 */
	public void addCountsToAnnotatedExons(List<AGene> geneList, List<TUnit> tunits, List<TUnitJunction> tjuncs, char str){
		for(AGene g : geneList){
			HashMap<ARegion, AIsoform> exonIsoMap = new HashMap<ARegion, AIsoform>();  //Exon to Isoform map
			//Get list of exons
			List<ARegion> exons = new ArrayList<ARegion>();
			for(AIsoform iso :g.getIsoforms()){
				if(g.getTransType().equals("protein_coding"))
					for(ARegion e : iso.getCDSs()){
						exons.add(e);
						exonIsoMap.put(e,iso);
					}
				else
					for(ARegion e : iso.getExons()){
						exons.add(e);
						exonIsoMap.put(e,iso);
					}
			}
			
			//Iterate through the unit list for exon-overlapping TUnits
			//Add read starts that overlap an exon.
			//Because of the way that we have defined TUnits as non-overlapping units that are consistent with 
			//annotated exonic boundaries, each TUnit either wholly overlaps an exon or not -- we cannot get 
			//TUnits that span an exon and an intron, for example.
			int i=binarySearchTUnitList(g.getCoords(), tunits);
			for(int x=i; x<tunits.size() && g.getCoords().getEnd()>=tunits.get(x).getCoords().getEnd(); x++){
				boolean oLap=false;
				for(ARegion e : exons){
					if(tunits.get(x).getCoords().overlaps(e.getCoords())){
						oLap=true;
						//Add coverage to isoform
						for(ReadLoader rl : loaders){
							int repID = replicateIndex.get(rl);
							exonIsoMap.get(e).addCoverage(repID, tunits.get(x).getHitStarts(repID)); 
						}
					}
				}
				if(oLap){
					//Add coverage to gene
					for(ReadLoader rl : loaders){
						int repID = replicateIndex.get(rl);
						g.addCoverage(repID, tunits.get(x).getHitStarts(repID)); 
					}
					//Add an annotation to the TUnit	
					tunits.get(x).addAnnotation(g.getID(), g.getName(), "EXONIC");
				}else{
					//This TUnit must overlap the gene, but not an exon within it. 
					tunits.get(x).addAnnotation(g.getID(), g.getName(), "INTRONIC");
				}
			}
			
			//Iterate through the junction list for gene-contained TUnitJunctions (Annotation only)
			int j=binarySearchTUnitJunctionList(g.getCoords(), tjuncs);
			for(int x=j; x<tjuncs.size() && g.getCoords().getEnd()>=tjuncs.get(x).getCoords().getEnd(); x++){
				boolean oLap=false;
				for(AIsoform iso :g.getIsoforms()){
					ArrayList<ARegion> isoExons = iso.getExons(); 
					for(int ex=0; ex<isoExons.size()-1; ex++){
						if(tjuncs.get(x).getCoords().getStart()==isoExons.get(ex).getCoords().getEnd() && 
								tjuncs.get(x).getCoords().getEnd()==isoExons.get(ex+1).getCoords().getStart()){
							oLap=true;
						}
					}
				}
				if(oLap){
					//Add an annotation to the TUnitJunction
					tjuncs.get(x).addAnnotation(g.getID(), g.getName(), "EXONIC");
				}else{
					//This TUnitJunction must overlap the gene, but not an exon within it. 
					tjuncs.get(x).addAnnotation(g.getID(), g.getName(), "INTRONIC");
				}
			}
		}
	}
	
	/**
	 * Examine the boundary spanning reads in the TUnits in each gene and connect genes with evidence of connections
	 * @param ug
	 * @return
	 */
	protected List<TUnitGene> findBoundaryConnections(List<TUnitGene> ug){
		List<TUnitGene> processed = new ArrayList<TUnitGene>();
		int totalConnects=0;
		boolean [] added = new boolean [ug.size()]; 
		for(int i=0; i<ug.size(); i++){added[i]=false;}
		
		for(int i=0; i<ug.size(); i++){
			if(!added[i]){
				TUnitGene curr = ug.get(i);
				if(curr.getUnitBoundaryHits()!= null && curr.getUnitBoundaryHits().size()>0){
					List<AlignHit> boundaryHits = curr.getUnitBoundaryHits();
					//Assume here that only nearby right-wards genes can be connected (ug should be sorted here).
					for(int j=i+1; j<=(i+10) && j<ug.size(); j++){
						if(!added[j]){
							TUnitGene next = ug.get(j);
							for(TUnit nu : next.getComponents()){
								int connectCount=0;
								for(AlignHit h : boundaryHits){
									if(h.blockOverlaps(nu.getCoords()) || h.getPairedRead()!=null && nu.getCoords().contains(h.getPairedRead()))
										connectCount++;
								}
								if(connectCount>=connectionThreshold){
									curr = curr.joinGene(next);
									added[j]=true;
									totalConnects+=connectCount;
								}
							}
						}
					}
				}
				processed.add(curr);
				added[i]=true;
			}
		}
		System.err.println(totalConnects+" boundary-spanning connection reads processed.");
		return processed;
	}
	
	/**
	 * Calculate a coverage threshold using a Poisson assumption and a given p-value 
	 * @param totalReads
	 * @return
	 */
	protected double calcCoverageThreshold(double totalReads){
		int countThres=0;
		DRand re = new DRand();
		Poisson P = new Poisson(0, re);
		double lambda = (totalReads*(readLength/binWidth))/((genomeLen*mappableGenome)/binWidth); 
		P.setMean(lambda);
		double l=1;
		double covThres = Math.pow(10,logCovThres);
		for(int b=1; l>covThres; b++){
			l=1-P.cdf(b);
			countThres=b;
		}
		return((double)Math.max(1,countThres));		
	}
	/**
	 * Get hits from the input list that overlap a region
	 * @param currReg the Region of interest
	 * @param hits AlignHits to filter
	 * @return List of AlignHits
	 */
	private List<AlignHit> overlappingHits(Region currReg, List<AlignHit> hits, char str){
		List<AlignHit> currHits = new ArrayList<AlignHit>();
		//We cannot use the native binary search method on the list of AlignHits, because of complications arising from align blocks
		int i = binarySearchAlignHitBlocks(currReg, hits); 
		while(i<hits.size() && hits.get(i).getStart()<=currReg.getEnd()){
			if(hits.get(i).blockOverlaps(currReg) && (str=='.' || hits.get(i).getStrand()==str))
				currHits.add(hits.get(i));
			i++;
		}
		return(currHits);
	}
	
	/**
	 * Get junctions that overlap a region at one or both ends
	 * @param currReg the Region of interest
	 * @param juncs TUnitJunction List
	 * @return List of AlignHits
	 */
	private List<TUnitJunction> overlappingJunctions(Region currReg, List<TUnitJunction> juncs, char str){
		List<TUnitJunction> currJuncs = new ArrayList<TUnitJunction>();
		int startIndex = binarySearchTUnitJunctionList(currReg, juncs);
		for(int x=startIndex; x<juncs.size() && currReg.getEnd()>=juncs.get(x).getCoords().getStart(); x++){
			if(juncs.get(x).getCoords().overlaps(currReg))
				currJuncs.add(juncs.get(x));
		}
		return(currJuncs);
	}
	
	/**
	 * Search a sorted list of AlignHits for the first align hit where ANY align blocks overlap the region.
	 * Not all hits from i to the last overlapping hit are guaranteed to overlap.  
	 * @param window Region to search for
	 * @param hits List of AlignHits
	 * @return int index on list
	 */
	private int binarySearchAlignHitBlocks(Region window, List<AlignHit> hits){
		int l = 0;
        int r = hits.size();
        int windowStart = window.getStart();
        //int windowEnd = window.getEnd();
        while (r - l > 10) {
            int c = (l + r) / 2;
            
            if(hits.get(c).blockOverlaps(window)){
            	r=c;
            }else{
            	if (windowStart > hits.get(c).getStart()) {
	                l = c;
	            } else {
	                r = c;
	            }
            }
        }
        while(l<hits.size() && !hits.get(l).blockOverlaps(window)){
        	l++;
        }
        return(l);
	}
	
	/**
	 * Search a sorted list of TUnits, and return the index that this region would be inserted on the list 
	 * @param window Region to search for
	 * @param tunits List of TUnits (presumed sorted and all on the same chromosome)
	 * @return int index on list
	 */
	private int binarySearchTUnitList(Region window, List<TUnit> tunits){
		int l = 0;
        int r = tunits.size();
        int windowStart = window.getStart();
        //int windowEnd = window.getEnd();
        while (r - l > 10) {
            int c = (l + r) / 2;
            if (windowStart > tunits.get(c).getCoords().getStart()) {
                l = c;
            } else {
                r = c;
            }
        }
        while(l<tunits.size() && tunits.get(l).getCoords().getStart() < windowStart){
        	l++;
        }
        return(l);
	}
	/**
	 * Search a sorted list of TUnits, and return the index that this region would be inserted on the list 
	 * @param window Region to search for
	 * @param tunits List of TUnits (presumed sorted and all on the same chromosome)
	 * @return int index on list
	 */
	private int binarySearchTUnitJunctionList(Region window, List<TUnitJunction> tunits){
		int l = 0;
        int r = tunits.size();
        int windowStart = window.getStart();
        //int windowEnd = window.getEnd();
        while (r - l > 10) {
            int c = (l + r) / 2;
            if (windowStart > tunits.get(c).getCoords().getStart()) {
                l = c;
            } else {
                r = c;
            }
        }
        while(l<tunits.size() && tunits.get(l).getCoords().getStart() < windowStart){
        	l++;
        }
        return(l);
	}
		
    //keep the number in bounds
	protected final double inBounds(double x, double min, double max){
		if(x<min){return min;}
		if(x>max){return max;}
		return x;
	}
	protected final int inBounds(int x, int min, int max){
		if(x<min){return min;}
		if(x>max){return max;}
		return x;
	}
	
	/**
	 * Main method for program execution
	 * @param args Command-line arguments
	 */
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
		if((!ap.hasKey("fa") &&!ap.hasKey("fai")&&!ap.hasKey("species"))  || !ap.hasKey("gtf")) { 
            System.err.println("Usage:\n" +
            		           "  --species <org;genome> OR " +
                               "  --fa <FASTA file> OR --fai <FASTA index>\n" +
                               "  --gtf <GTF file>\n"+
                               "  --exptNAME <SAM file of reads>\n" +
                               "  --out <output file prefix>\n" +
                               "  --regions <regions to analyze>\n" +
                               "  --minsubsize <sub-chr size>\n" +
                               "  --minsubspace <sub-chr empty spacing>\n" +
                               "  --threads <num threads>\n" +
                               "  --genesegment [segment with genes instead of coverage]\n" +
                               "");
            return;
        }
        try {
        	//Process command-line options
        	GenomeLoader gLoad=null;
        	if(ap.hasKey("species")){
				Pair<Species, Genome> pair = Args.parseGenome(args);
				Genome currgen = pair.cdr();
				gLoad = new GenomeLoader(currgen);
        	}else if(ap.hasKey("fa")){
        		String faFile = ap.getKeyValue("fa");
        		gLoad = new GenomeLoader(new File(faFile), true);
        	}else if(ap.hasKey("fai")){
        		String faiFile = ap.getKeyValue("fai");
        		gLoad = new GenomeLoader(new File(faiFile), false);
        	}
			Args.parseFileHandles(args, "sam");
			String gtfFile = ap.getKeyValue("gtf");
			String outName = Args.parseString(args, "out", "out");
			String regionFile = Args.parseString(args, "regions", null);
			Integer minSubChrSize =Args.parseInteger(args, "minsubsize", 5000000);
			Integer minSubChrSpace = Args.parseInteger(args, "minsubspace", 500000);
			Integer threads = Args.parseInteger(args,"threads", 1);
			boolean segmentWithGenes = Args.parseArgs(args).contains("genesegment");
			
			//Load the experiments
			HashMap<String,List<File>> samFiles = new HashMap<String,List<File>>();
			ArrayList<String> conditionNames = new ArrayList<String>();
			Vector<String> exptTags=new Vector<String>();
			for(String s : args)
	        	if(s.contains("--expt"))
	        		if(!exptTags.contains(s)){
	        			exptTags.add(s);
	        			String name = s.replaceFirst("--expt", ""); 
		        		conditionNames.add(name);
	        		}
			if(exptTags.size()==0){
			    System.err.println("Error: No experiments provided.\nUse the --expt option.");
			    System.exit(1);
			}
			for(String name : conditionNames){
				if(!samFiles.containsKey(name)){
        			samFiles.put(name, new ArrayList<File>());
        		}
        		samFiles.get(name).addAll(Args.parseFileHandles(args, "expt"+name));
			}	        	
			
			//Load genes
			GTFAnnotationLoader reader = new GTFAnnotationLoader(new File(gtfFile), gLoad);
			Collection<AGene> geneSet = reader.loadGenes();
			
			//Initialize read loaders
			ArrayList<ReadLoader> loaders = new ArrayList<ReadLoader>();
			for(String name : conditionNames){
				for(File sf : samFiles.get(name)){
					//ReadLoader
					ReadLoader rl = new SAMReadLoader(sf, name, gLoad.getGenome(), gLoad.getNameTranslator());
					loaders.add(rl);
				}
			}

			//Estimate the region sets
			Map<String, List<Region>> regionsOfInterest;
			if(regionFile!=null)
				regionsOfInterest = ReadUtils.loadRegionsFromFile(gLoad.getGenome(), regionFile);
			else{
				GenomeSegmenter segmenter = new GenomeSegmenter(gLoad, minSubChrSize, minSubChrSpace);
				if(segmentWithGenes)
					regionsOfInterest = segmenter.segmentWithGenes(geneSet);
				else
					regionsOfInterest = segmenter.segmentWithReads(loaders);
			}
			ReadUtils.writeRegionsToFile(regionsOfInterest, outName+".segments.coords");
			
			System.err.println("Starting scanner");
			//Initialize the window scanner
			RNAWindowScanner scanner = new RNAWindowScanner(outName, gLoad, loaders, geneSet, regionsOfInterest);
			//Execute
			scanner.execute(threads);
				
	    } catch (NotFoundException e) {
			e.printStackTrace();
		}
	}
}
