package edu.psu.compbio.seqcode.projects.shaun;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;

import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;
import be.ac.ulg.montefiore.run.jahmm.*;
import be.ac.ulg.montefiore.run.jahmm.io.FileFormatException;
import be.ac.ulg.montefiore.run.jahmm.io.HmmReader;
import be.ac.ulg.montefiore.run.jahmm.io.HmmWriter;
import be.ac.ulg.montefiore.run.jahmm.io.OpdfGaussianReader;
import be.ac.ulg.montefiore.run.jahmm.io.OpdfGaussianWriter;
import be.ac.ulg.montefiore.run.jahmm.io.OpdfReader;
import be.ac.ulg.montefiore.run.jahmm.io.OpdfWriter;
import be.ac.ulg.montefiore.run.jahmm.learn.*;
import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Organism;
import edu.psu.compbio.seqcode.genome.location.NamedRegion;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqExptHandler;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.*;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.PointParser;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.RegionParser;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
/* Finds enriched domains in Solexa experiments using a HMM. 
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
 *  ChipSeqPeakFinder --species "Mus musculus" --genome mm8 --expt PPG_Solexa_RAR_8hr --back PPG_Solexa_WCE_ES+2d --back PPG_Solexa_WCE_2+1 --outpeak rar_8hr.peaks --outseq rar_8hr_peaks.seq --seqwin 200   
 * 
 */
public class ChipSeqHMMDomainFinder {

	private final double MIN_OBS_VAL = 0.1;
	private int chunkSize = 50000000;
	private int allowedZeroRuns=10000;
	
	//HMM variables
	private Hmm<ObservationReal> hmm;
	private ArrayList<Opdf<ObservationReal>> pdfs;
	private ArrayList<DataObs> data;
	private int numStates=2;
	private final int UNBOUNDA = 0, UNBOUNDB=1, DOMAIN = numStates-1;
	private double [] initProbs;
	private double [][] initTrans;
	private double [] initMeans;
	private double [] initCovars;
	private int trainIter=100;
	
	//General variables
	private double poissThres=-9;
	private double winWidth=500;
	private double winStep = 250;
	private int readLen=26;
	private double readLength=26;
	private int readExtension = 174;
	private int readShift = 50;
	private boolean shiftTags = true;
	private boolean needlefiltering=true;
	private Organism org;
	private Genome gen;
	private double genomeLen=0;
	private ArrayList<SeqExptHandler> IPhandles;
	protected int ipBasePoissonThres=0;
	private double iphittot;
	private ArrayList<Region> trainSet = null;
	private String trainFile =null;
	private int peakFlanks = 10000; //flanks to add onto peaks in a training set
	private String outRegName="out.domains", outGFFName="out.gff", outHMMName="out.hmm";
	public String [] annots = new String[]{
		"refGene"//,  "knownGene"//, "mgcGenes", "ensGene"
	};
	private String ucscColor = "51,102,0";
	private boolean sigNorm=false, sigSet=false; //sigmoid distribution variables calculated
	private final double sigMin=2;
	private double sigMedian=0;
	private double sigMean =0;
	private double sigStddev=0;
	private boolean printTraining=false;
	
	
	
	/* command-line driver */
	public static void main(String[] args) {
		
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("species") || !ap.hasKey("genome")||!ap.hasKey("expt")){// || !ap.hasKey("back")) { 
            System.err.println("Usage:\n " +
                               "ChipSeqHMMDomainFinder " +
                               "--species <organism name> " +
                               "--genome <genome version> "+
                               "--expt <solexa expt> " +
                               "--trainpeaks <coordinates of regions to extend & train HMM> " +
                               "--train <coordinates of regions to train HMM> " +
                               "--hmmstates <number of HMM states (2 or 3 only)> " +
                               "--out <output root> "+
                               "--binwidth <width of bins> "+
                               "--binstep <offset of bins> " +
                               "--readlen <read length> " +
                               "--readshift <read shift> " +
                               "--hmmfile <hmm file name (optional)> " +
                               "--signorm [sigmoid normalization ala ChromaSig] " +
                               "--printtraining [print the training set] ");
            return;
        }
        String species = ap.getKeyValue("species");
        String genome = ap.getKeyValue("genome");
        ArrayList<String> exptNames = new ArrayList<String>();
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--expt")) { 
                exptNames.add(args[++i]);
            }
        }
        String trainFileName=null; int ns=2; boolean extendPeaks=false;
        if(ap.hasKey("train")){trainFileName = ap.getKeyValue("train");}
        if(ap.hasKey("trainpeaks")){trainFileName = ap.getKeyValue("trainpeaks"); extendPeaks=true;}
        if(ap.hasKey("hmmstates")){ ns = new Integer(ap.getKeyValue("hmmstates")).intValue();}
        int rL=26, rE=0, rS=50;
        if(ap.hasKey("readlen")){ rL = new Integer(ap.getKeyValue("readlen")).intValue();}
        if(ap.hasKey("readext")){ rE = new Integer(ap.getKeyValue("readext")).intValue();}
        if(ap.hasKey("readshift")){ rS = new Integer(ap.getKeyValue("readshift")).intValue();}
        String hmmfilename=null;
        if(ap.hasKey("hmmfile")){hmmfilename = ap.getKeyValue("hmmfile");}
        
        //Initialize the peak finder
        ChipSeqHMMDomainFinder finder = new ChipSeqHMMDomainFinder(species, genome, exptNames, trainFileName, rL, rE, rS, ns, extendPeaks, hmmfilename);
        
        //Set the options
        if(ap.hasKey("out")){finder.setFileName(ap.getKeyValue("out"));}
        finder.setSigNorm(ap.hasKey("signorm"));
        finder.setPrintTraining(ap.hasKey("printtraining"));
        if(ap.hasKey("binwidth")){ finder.setBinWidth(new Double(ap.getKeyValue("binwidth")).doubleValue());}
        if(ap.hasKey("binstep")){ finder.setBinStep(new Double(ap.getKeyValue("binstep")).doubleValue());}
        
        //Run the peak finder
        ArrayList<Region> domains = finder.execute();
        System.out.println("Printing");
        finder.printDomainLocationsGFF(domains);
       	finder.printDomainLocations(domains);
       	finder.printHMM();
       	System.out.println("Finished!");
	}
	
	/* Constructor requires a species, genome version and lists of ip and background solexa experiments */
	public ChipSeqHMMDomainFinder(String species, String genome, ArrayList<String> ips, String trainFileName, int rL, int rE, int rS, int numHMMStates, boolean extendPeaks, String hmmFile){
        this.setReadLen(rL);
        this.setReadExtension(rE);
        this.setReadShift(rS);
		System.out.println("Initializing the ChIP-Seq HMM Domain Finder");
        iphittot=0; 
        IPhandles = new ArrayList<SeqExptHandler>();
        if(numHMMStates==2 || numHMMStates==3){numStates = numHMMStates;}        
        trainFile = trainFileName;
        try {
			//Load the reads
			org = Organism.getOrganism(species);
			gen = org.getGenome(genome);
			genomeLen = gen.getGenomeLength(); 
			
			for(String exptName : ips){
				System.out.print(String.format("%s\t", exptName));
				SeqLocator curr=null;
				if (exptName.indexOf(';') != -1) {
                    String[] pieces = exptName.split(";");
                    if(pieces.length==2)
                    	curr = new SeqLocator(pieces[0], pieces[1]);
                    else if(pieces.length==3)
                    	curr = new SeqLocator(pieces[0], pieces[1], pieces[2]);
                }else{
                	System.err.println("Experiment name incorrectly assigned");
                }
				SeqExptHandler h = new SeqExptHandler(gen, curr);
	            h.setReadLength(readLength);
                if(shiftTags){
                	h.setReadExtension(readShift*2);
                	h.setReadShift(readShift);
                }else
                	h.setReadExtension(readExtension);
                IPhandles.add(h);
				iphittot += h.getHitCount();
			
			}System.out.print(String.format("%.0f reads loaded\n", iphittot));
			
			//Load the list of training regions
			if(trainFile !=null){
				if(extendPeaks)
					trainSet = readRegionsFromFile(gen,trainFile, peakFlanks);
				else
					trainSet = readRegionsFromFile(gen,trainFile);
			}else{
				trainSet = new ArrayList<Region>();
				Iterator<NamedRegion> chroms = new ChromRegionIterator(gen);
				while (chroms.hasNext()) {
					trainSet.add(chroms.next());
				}
			}
			//count the training set 
			int totalbp=0;
			for(Region r : trainSet){
				totalbp += r.getWidth();
			}System.out.println(trainSet.size()+" training sequences loaded totaling "+totalbp+"bp");
			
			//Initialize the HMM
			if(hmmFile==null){
				if(numStates==3){
					initTrans = new double[][]{{0.98,0.01, 0.01}, {0.01,0.98, 0.01},{0.01,0.01, 0.98}};
					initMeans = new double[]{1, 5, 20.0};
					initCovars = new double[]{2.0, 5.0, 100.0};
					initProbs = new double[]{0.45, 0.45, 0.1};
				}else if(numStates==2){
					initTrans = new double[][]{{0.95,0.05}, {0.02,0.98}};
					initMeans = new double[]{2, 20.0};
					initCovars = new double[]{4.0, 100.0};
					initProbs = new double[]{0.8, 0.2};
				}
				pdfs = new ArrayList<Opdf<ObservationReal>>();
				for(int i=0; i<numStates; i++)
					pdfs.add(new OpdfGaussian(initMeans[i], initCovars[i]));
				hmm = new Hmm<ObservationReal>(initProbs, initTrans, pdfs);
			}else{
				try {
					OpdfReader opdfr = new OpdfGaussianReader();
					HmmReader reader = new HmmReader();
					FileReader fr = new FileReader(hmmFile);
					hmm = reader.read(fr, opdfr);
				} catch (FileNotFoundException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (FileFormatException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}			
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/* Run the HMM (assumes all initialization took place successfully) */
	public ArrayList<Region> execute(){
		ipBasePoissonThres = getPoissonThreshold(Math.pow(10, poissThres), iphittot, 1, genomeLen, 0.8, 1, 1);
		readLength = IPhandles.get(0).getReadLength();
		
		//Training set generation
		System.out.println("Initializing Training Set Scores");
		data = generateScores(trainSet, true);
		data = filterData(data);
		//Strip the observations out of the DataObs
		ArrayList<ArrayList<ObservationReal>> observations = new ArrayList<ArrayList<ObservationReal>>();
		int totalObsCount = 0;
		for(DataObs d : data){
			observations.add(d.obs);
			totalObsCount += d.obs.size();
		}
		
		//Train the HMM
		System.out.println("Training the "+numStates+"-state HMM with "+totalObsCount+" observations in "+observations.size()+" sequences.");
		trainHMM(observations);
		if(printTraining)
			printData();
		
		//Predict domains in the entire genome
		System.out.println("Calling bound domains in chromosome: ");
		ArrayList<Region> domains = new ArrayList<Region>(); 
		double totalDomainLen=0;
		Iterator<NamedRegion> chroms = new ChromRegionIterator(gen);
		while (chroms.hasNext()) {
			Region c = chroms.next();
			System.out.println(c.getChrom());
			ArrayList<Region> onechrom = new ArrayList<Region>();
			onechrom.add(c);
			data = generateScores(onechrom, false);
			ArrayList<Region> currDomains = callDomains(data);
			for(Region d : currDomains)
				totalDomainLen+=(double)d.getWidth();
			domains.addAll(currDomains);	
		}
		
		System.out.println(domains.size()+" Domains Called!");
		double perc = 100*totalDomainLen/genomeLen;
		System.out.println(totalDomainLen+" bases covered by domains ("+perc+"% of the genome)");		
		return(domains);
	}

	//Some option setters
	public void setNumStates(int n){numStates=n;}
	public void setSigNorm(boolean s){sigNorm=s;}
	public void setPrintTraining(boolean p){printTraining=p;}
	public void setBinWidth(double b){winWidth=b;}
	public void setBinStep(double b){winStep=b;}
	public void setTrainFileName(String f){trainFile=f;}
	public void setFileName(String f){outRegName=f+".domains"; outGFFName = f+".gff"; outHMMName = f+".hmm";}
	public void setReadLen(int x){readLength=x;}
	public void setReadExtension(int x){readExtension=x; shiftTags=false; readShift = 0;}
	public void setReadShift(int x){readShift=x; readExtension = x*2; shiftTags=true;}
	
	/*Returns the binned read-count landscape of the genome*/
	public ArrayList<DataObs> generateScores(ArrayList<Region> regions, boolean printProgress){
		ArrayList<DataObs> allData = new ArrayList<DataObs>();
		if(printProgress){System.out.println(String.format("Genome Length: %.0f", genomeLen));}
		double basesDone=0, printStep=10000000,  numPrint=0;
		
		ArrayList<Double> sigData = new ArrayList<Double>();
		double sum=0,num=0;//for simoid normalization
		for(Region currTrainReg : regions){
			//Split the job up into chunks
			for(int x=currTrainReg.getStart(); x<=currTrainReg.getEnd(); x+=chunkSize){
				int y = x+chunkSize; 
				if(y>currTrainReg.getEnd() || y>gen.getChromLength(currTrainReg.getChrom())){y=Math.min(gen.getChromLength(currTrainReg.getChrom()), currTrainReg.getEnd());}
				if(x<1){x=1;}
				Region currRegion = new Region(gen, currTrainReg.getChrom(), x, y);
				DataObs currData = new DataObs(currRegion, (int)winStep);
				
				LinkedList<StrandedRegion> ipHits = new LinkedList<StrandedRegion>();
				for(SeqExptHandler IP: IPhandles){
					if(shiftTags)
						ipHits.addAll(IP.loadShiftedExtendedHits(currRegion));
					else
						ipHits.addAll(IP.loadExtendedHits(currRegion));
				}
			
				int [] ipHitCounts = makeHitLandscape(ipHits, currRegion, ipBasePoissonThres);
				
				for(int i=currRegion.getStart(); i<currRegion.getEnd()-(int)winWidth; i+=(int)winStep){
					Region currWin = new Region(gen, currTrainReg.getChrom(), i, (int)(i+winWidth-1));
					
					int binid = (int)Math.max(0, ((double)(currWin.getStart()-currRegion.getStart())/winStep));
					double ipWinHits=(double)ipHitCounts[binid];
					
					Double count = new Double(ipWinHits);
					if(count <MIN_OBS_VAL || count.isNaN() || count.isInfinite()){
						currData.obs.add(new ObservationReal(MIN_OBS_VAL));
						sum+=MIN_OBS_VAL;
					}else{
						currData.obs.add(new ObservationReal(count));
						sum+=count;
					}num++;
					
					if(sigNorm && count>=sigMin && !sigSet)
						sigData.add(count);

					//Print out progress
					basesDone+=winStep;
					if(printProgress){
						if(basesDone > numPrint*printStep){
							if(numPrint%10==0){System.out.print(String.format("(%.0f)", (numPrint*printStep)));}
							else{System.out.print(".");}
							if(numPrint%50==0 && numPrint!=0){System.out.print("\n");}
							numPrint++;
						}
					}
				}
				allData.add(currData);
			}
		}if(printProgress){System.out.print("("+basesDone+")\n");}
		
		//Sigmoid normalization
		if(sigNorm){
			System.out.println("Sigmoid Normalizing");
			if(!sigSet){
				sigMean = sum/num;
				//find median & std dev:
				Collections.sort(sigData);
				sigMedian = sigData.get(sigData.size()/2);
				double tmpsum=0;
				for(Double d : sigData){
					tmpsum += (sigMean-d)*(sigMean-d);
				}tmpsum =tmpsum/(double)sigData.size();
				sigStddev = Math.sqrt(tmpsum);
				System.out.println("Mean: "+sigMean+" Median: "+sigMedian+" StdDev: "+sigStddev);
				sigSet=true;
			}
			//Normalize
			ArrayList<DataObs> normData = new ArrayList<DataObs>();
			for(DataObs dob : allData){
				DataObs currData = new DataObs(dob.coords, (int)winStep);
				for(ObservationReal or : dob.obs){
					if(or.value>0)
						currData.obs.add(new ObservationReal(1/(1+Math.exp(-1*(or.value-sigMedian)/sigStddev))));
					else
						currData.obs.add(new ObservationReal(MIN_OBS_VAL));
				}
				normData.add(currData);
			}allData = normData;
		}
		return(allData);
	}
	//print out the training data
	private void printData(){
		for(DataObs d : data){
			d.print();
		}
	}
	//Filters out long runs of zeros (presumably denoting runs of 'N' or other un-mappable regions)
	private ArrayList<DataObs> filterData(ArrayList<DataObs> data){
		ArrayList<DataObs> filtered = new ArrayList<DataObs>();
		
		for(DataObs d : data){
			ArrayList<Integer> zeroRunStarts = new ArrayList<Integer>();
			ArrayList<Integer> zeroRunStops = new ArrayList<Integer>();
			int zStart=0, zStop=0;
			boolean rec=false;
			for(int i=0; i<d.obs.size(); i++){
				if(d.obs.get(i).value==MIN_OBS_VAL){
					if(rec){
						zStop=i;
					}else{
						//start a new run
						rec=true;
						zStart=i; zStop=i;
					}
				}else{
					if(rec){
						//previous run long enough? 
						if((zStop-zStart)*d.stepSize >= allowedZeroRuns){
							zeroRunStarts.add(new Integer(zStart));
							zeroRunStops.add(new Integer(zStop));
						}						
					}rec=false;
				}
			}if(rec){
				if((zStop-zStart)*d.stepSize >= allowedZeroRuns){
					zeroRunStarts.add(new Integer(zStart));	zeroRunStops.add(new Integer(zStop));
				}						
			}
			
			if(zeroRunStarts.size()>0){
				for(int i=0; i<=zeroRunStarts.size(); i++){
					int rstart, rstop, ostart=0, ostop=0;
					if(i==0){//first
						rstart = d.coords.getStart();
						rstop=rstart+(zeroRunStarts.get(0)*d.stepSize);
						ostart=0; ostop=zeroRunStarts.get(0);												
					}else if(i==zeroRunStarts.size()){//last
						rstart = d.coords.getStart()+(zeroRunStops.get(i-1)*d.stepSize);
						rstop=d.coords.getEnd();
						ostart=zeroRunStops.get(i-1); ostop=d.obs.size();						
					}else{
						rstart = d.coords.getStart()+(zeroRunStops.get(i-1)*d.stepSize);
						rstop=rstart+(zeroRunStarts.get(i)*d.stepSize);
						ostart=zeroRunStops.get(i-1); ostop=zeroRunStarts.get(i);
					}
					Region currReg = new Region(d.coords.getGenome(), d.coords.getChrom(), rstart, rstop);
					DataObs obs = new DataObs(currReg, d.stepSize);
					for(int o=ostart; o<ostop; o++){obs.obs.add(d.obs.get(o));}
					if(currReg.getWidth()>=1000){
						filtered.add(obs);
					}
				}
			}else{
				filtered.add(d);
			}
		}
		
		return(filtered);
	}
	
	//Trains the HMM
	private void trainHMM(ArrayList<ArrayList<ObservationReal>> observations){
		BaumWelchScaledLearner bwsl = new BaumWelchScaledLearner();
		bwsl.setNbIterations(trainIter);
		hmm = bwsl.learn(hmm, observations);
		System.out.println(hmm.toString());
	}
	
	//Print the HMM to a file
	public void printHMM(){
		try{
			OpdfWriter opdfw = new OpdfGaussianWriter();
			FileWriter fout = new FileWriter(outHMMName);
			HmmWriter.write(fout, opdfw, hmm); 
			fout.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
	}
	
	//Uses a trained HMM to call Domains as Regions
	private ArrayList<Region> callDomains(ArrayList<DataObs> data){
		ArrayList<Region> domainCalls = new ArrayList<Region>();
		
		for(DataObs d : data) {
        	ViterbiCalculator vc = new ViterbiCalculator(d.obs, hmm);
        	int[] states = vc.stateSequence();
        	int state = states[0], start=0, i=0;
        	int genomicStart, genomicEnd;
        	for (i=0; i<states.length; i++) {
        		if (states[i]==DOMAIN && state != DOMAIN) {
        			start = i;        			
        		}
        		if (states[i]!=DOMAIN && state==DOMAIN) {
        			genomicStart = d.coords.getStart() + (start*(int)winStep); 
        			genomicEnd = d.coords.getStart() + ((i-1)*(int)winStep) +(int)winWidth;
        			domainCalls.add(new Region(d.coords.getGenome(), d.coords.getChrom(), genomicStart, genomicEnd));
        		}
        		state = states[i];
        	}
        	if (state==DOMAIN) {
        		genomicStart = d.coords.getStart() + (start*(int)winStep); 
    			genomicEnd = d.coords.getStart() + ((i-1)*(int)winStep) +(int)winWidth;
    			domainCalls.add(new Region(d.coords.getGenome(), d.coords.getChrom(), genomicStart, genomicEnd));
        	}
		}
		//Merge overlaps
        for(int i=0; i<domainCalls.size(); i++){
        	Region currR = domainCalls.get(i);
        	Region newR = currR.clone();
        	for(int j=i; j<domainCalls.size(); j++){
            	Region nextR = domainCalls.get(j);
        		if(newR.overlaps(nextR)){
        			newR = newR.combine(nextR);
        			domainCalls.remove(nextR);            			
        		}
        	}
        	domainCalls.add(newR);
        	domainCalls.remove(currR);            	
        }Collections.sort(domainCalls);
		return(domainCalls);
	}
	
	/* print out the peaks */
	public void printDomainLocations(ArrayList<Region> domains){
		try {
			FileWriter fout = new FileWriter(outRegName);
			
			for(Region r : domains){
				fout.write(r.getLocationString()+"\n");
			}
			fout.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/* print out the peaks in UCSC format*/
	public void printDomainLocationsGFF(ArrayList<Region> domains){
		try {
			FileWriter fout = new FileWriter(outGFFName);
			
			int count=1;
			fout.write("track name=\"HMM Domains\" color="+ucscColor+" itemRgb=On\n");
			for(Region r : domains){
				fout.write("chr"+r.getChrom()+"\tChIP-Seq_HMM\tDomain\t"+r.getStart()+"\t"+r.getEnd()+"\t.\t.\t.\t"+count+"\n");
				count++;
			}
			fout.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/* Binomial test for differences between two population proportions */
	//Probably not used in this class
	private double binomialSampleEquality(double X1, double X2, double n1, double n2){
		double P1 = X1/n1;
		double P2 = X2/n2;
		double P = (X1+X2)/(n1+n2);
		double Z0 = (P1-P2)/(Math.sqrt(P*(1-P)*((1/n1)+(1/n2))));
		return(Z0);
	}
	
	private ArrayList<StrandedRegion> overlappingHits(LinkedList<StrandedRegion> hits, Region window){
		ArrayList<StrandedRegion> sub = new ArrayList<StrandedRegion>();
		for(StrandedRegion r : hits){
			if(window.overlaps(r)){sub.add(r);}			
		}
		return(sub);
	}
	protected int [] makeHitLandscape(LinkedList<StrandedRegion> hits, Region currReg, int perBaseMax){
		int numBins = (int)(currReg.getWidth()/winStep);
		int [] counts = new int[currReg.getWidth()+1];
		int [] land = new int[numBins+1];
		for(int i=0; i<=numBins; i++){land[i]=0;}
		for(int i=0; i<=currReg.getWidth(); i++){counts[i]=0;}
		int count=0, sum=0;
		for(StrandedRegion r : hits){
			int offset=r.getStart()-currReg.getStart()<0 ? 0 : r.getStart()-currReg.getStart();
			counts[offset]++;
			if(!needlefiltering || (counts[offset] <= perBaseMax)){
				count++;
				int binstart = (int)Math.max(0, ((double)(offset/winStep)));
				int binend = (int)(Math.min((double)(r.getEnd()-currReg.getStart())/winStep, numBins-1));
				for(int i=binstart; i<=binend; i++){
					land[i]++;
					sum++;
				}
			}
		}
		return(land);
	}	
	
	private ArrayList<Region> readRegionsFromFile(Genome g, String filename) {
        ArrayList<Region> regions = new ArrayList<Region>();
        try {
            BufferedReader r = new BufferedReader(new FileReader(filename));
            String s;
            while ((s = r.readLine()) != null) {
                Region region = Region.fromString(g,s);
                if (region != null) {
                    regions.add(region);
                } else {
                    System.err.println("Couldn't parse " + s);
                }

            }
            r.close();
        } catch (IOException ex) {
            throw new RuntimeException("Can't read " + filename,ex);
        }
        return regions;
    }
	private ArrayList<Region> readRegionsFromFile(Genome g, String filename, int extend) {
        ArrayList<Region> regions = new ArrayList<Region>();
        ArrayList<Region> tmpregions = new ArrayList<Region>();
        try {
            BufferedReader reader = new BufferedReader(new FileReader(filename));
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                String[] words = line.split("\\s+");
                
                if(words.length>=3){
	                RegionParser rparser = new RegionParser(gen);
	            	Region p = rparser.execute(words[0]);
                	int rstart = p.getStart()-extend<1 ? 1:p.getStart()-extend;
                	int rend = p.getEnd()+extend>=gen.getChromLength(p.getChrom()) ? gen.getChromLength(p.getChrom()):p.getEnd()+extend;
                	Region r = new Region(p.getGenome(), p.getChrom(), rstart, rend);
                	tmpregions.add(r);
	            }
            }
            reader.close();
            
            //Merge overlaps
            for(int i=0; i<tmpregions.size(); i++){
            	Region currR = tmpregions.get(i);
            	Region newR = currR.clone();
            	for(int j=0; j<tmpregions.size(); j++){
                	Region nextR = tmpregions.get(j);
            		if(newR.overlaps(nextR)){
            			newR = newR.combine(nextR);
            			tmpregions.remove(nextR);            			
            		}
            	}
            	regions.add(newR);
            	tmpregions.remove(currR);            	
            }            
        } catch (IOException ex) {
            throw new RuntimeException("Can't read " + filename,ex);
        }
        return regions;
    }private ArrayList<Region> readPeaksFromFile(Genome g, String filename, int extend) {
        ArrayList<Region> regions = new ArrayList<Region>();
        ArrayList<Region> tmpregions = new ArrayList<Region>();
        try {
            BufferedReader reader = new BufferedReader(new FileReader(filename));
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                String[] words = line.split("\\s+");
                
                if(words.length>=3){
	                PointParser pparser = new PointParser(gen);
	            	Point p = pparser.execute(words[2]);
                	int rstart = p.getLocation()-extend<1 ? 1:p.getLocation()-extend;
                	int rend = p.getLocation()+extend>=gen.getChromLength(p.getChrom()) ? gen.getChromLength(p.getChrom()):p.getLocation()+extend;
                	Region r = new Region(p.getGenome(), p.getChrom(), rstart, rend);
                	tmpregions.add(r);
	            }
            }
            reader.close();
            
            //Merge overlaps
            for(int i=0; i<tmpregions.size(); i++){
            	Region currR = tmpregions.get(i);
            	Region newR = currR.clone();
            	for(int j=0; j<tmpregions.size(); j++){
                	Region nextR = tmpregions.get(j);
            		if(newR.overlaps(nextR)){
            			newR = newR.combine(nextR);
            			tmpregions.remove(nextR);            			
            		}
            	}
            	regions.add(newR);
            	tmpregions.remove(currR);            	
            }
            for(Region r : regions)
            	System.out.println(r.getLocationString());
        } catch (IOException ex) {
            throw new RuntimeException("Can't read " + filename,ex);
        }
        return regions;
    }
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
	public class DataObs{
		public ArrayList<ObservationReal> obs;
		public int stepSize;
		public Region coords;
		public DataObs(Region r, int step){
			obs=new ArrayList<ObservationReal>();
			coords=r;
			stepSize=step;
		}
		public void print(){
			System.out.println("Region:\t"+coords.getLocationString()+"\tStepSize:\t"+stepSize);
			for(ObservationReal o : obs){
				System.out.println(o.toString());
			}
		}
	}
}



