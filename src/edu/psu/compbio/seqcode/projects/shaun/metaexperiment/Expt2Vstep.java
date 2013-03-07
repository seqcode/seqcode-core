package edu.psu.compbio.seqcode.projects.shaun.metaexperiment;

import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;

import edu.psu.compbio.seqcode.gse.datasets.general.NamedRegion;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqExptHandler;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;

public class Expt2Vstep {
	private Organism org;
	private Genome gen;
	protected ArrayList<SeqExptHandler> IPhandles;
	protected ArrayList<SeqExptHandler> backhandles;
	protected int [] ipStackedHitCounts;
	protected int [] backStackedHitCounts;
	private int poissThres=-9;
	//Settings
	private static int winSize=200, winStep=100;
	private double readLength=26, readExtension=200, readShift=0;
	private boolean shiftTags=false;
	private double iphittot=0, backhittot;
	private double genomeLen=0;
	private boolean noBackground=true;
	private double ipPoissonThres=0, backPoissonThres=0;
	protected double ipBasePoissonThres=0;
	protected double backBasePoissonThres=0;
	
	
	public static void main(String[] args) throws SQLException, NotFoundException {
		Pair<Organism,Genome> pair = Args.parseGenome(args);
		if(pair==null || !Args.parseArgs(args).contains("expt")){printError();return;}
		List<SeqLocator> expts = Args.parseChipSeq(args,"expt");
		List<SeqLocator> backs = Args.parseChipSeq(args,"back");
		double rLen = Args.parseDouble(args,"readlen",SeqExptHandler.defaultReadLength);
        double rExt = Args.parseDouble(args,"readextend",SeqExptHandler.defaultReadExtension);
        double rShift =0;
        if(Args.parseArgs(args).contains("shifttags")){
        	rShift = Args.parseInteger(args,"shifttags", 0);
        }
        String outName = Args.parseString(args, "out", "out"); 
    
		//Initialize
		Expt2Vstep converter = new Expt2Vstep(pair.cdr(), expts,backs, rLen, rExt, rShift, Args.parseInteger(args,"binwidth",winSize),Args.parseInteger(args,"binstep",winStep));
        
        converter.execute(outName);
	}
	
	
	public Expt2Vstep(Genome gen, List<SeqLocator> ips,List<SeqLocator> backs, double rLen, double rExt, double rShift, int binW, int binStep) {
		System.out.println("Initializing the Converter");
		this.gen = gen;
        readLength = rLen;
        readExtension = rExt;
        readShift = rShift;
        if(readShift > 0)
        	shiftTags=true;
        iphittot=0; backhittot=0;
        winSize = binW; winStep = binStep;
        
        //Load experiments
        loadExperiments(ips, backs);
        
        //Initialize Poisson Models
        ipPoissonThres = getPoissonThreshold(Math.pow(10, poissThres), iphittot, readLength, genomeLen, 0.8, winSize, winStep);
        System.out.println("Finished");
	}
	
	public void execute(String outName){
		try {
			FileWriter fw = new FileWriter(outName+".vstep");
		
			double basesDone=0, printStep=10000000,  numPrint=0;
			
			double pt6=getPoissonThreshold(Math.pow(10, -6), iphittot, readLength, genomeLen, 0.8, winSize, winStep);
			double pt7=getPoissonThreshold(Math.pow(10, -7), iphittot, readLength, genomeLen, 0.8, winSize, winStep);
			double pt8=getPoissonThreshold(Math.pow(10, -8), iphittot, readLength, genomeLen, 0.8, winSize, winStep);
			double pt9=getPoissonThreshold(Math.pow(10, -9), iphittot, readLength, genomeLen, 0.8, winSize, winStep);
			//Print the header
			fw.write("track type=wiggle_0 name=\""+outName+"\" description=\""+outName+" summary\" totalReads="+iphittot+" poissonThres-6="+pt6+" poissonThres-7="+pt7+" poissonThres-8="+pt8+" poissonThres-9="+pt9+"\n");
			
			double ipTotHits = iphittot;
	        // set to 1 if no background so the binomial test doesn't return NaN
			double backTotHits = noBackground ? 1 : backhittot;        
			
			ChromRegionIterator chroms = new ChromRegionIterator(gen);
			while(chroms.hasNext()){
				NamedRegion currentRegion = chroms.next();
				fw.write("variableStep chrom=chr"+currentRegion.getChrom()+" span="+winSize+" step="+winStep+"\n");
				
				//Split the job up into chunks of 100Mbp
				for(int x=currentRegion.getStart(); x<=currentRegion.getEnd(); x+=100000000){
					int y = x+100000000; 
					if(y>currentRegion.getEnd()){y=currentRegion.getEnd();}
					Region currSubRegion = new Region(gen, currentRegion.getChrom(), x, y);
					
					LinkedList<StrandedRegion> ipHits = new LinkedList<StrandedRegion>();
					LinkedList<StrandedRegion> backHits = new LinkedList<StrandedRegion>();
					
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
					for(int i=currSubRegion.getStart(); i<currSubRegion.getEnd()-(int)winSize; i+=(int)winStep){
						Region currWin = new Region(gen, currentRegion.getChrom(), i, (int)(i+winSize-1));
						
						int binid = (int)Math.max(0, ((double)(currWin.getStart()-currSubRegion.getStart())/winStep));
						double ipWinHits=(double)ipStackedHitCounts[binid];
						double backWinHits= noBackground ? 0 : ((double)backStackedHitCounts[binid]);
						
						if(ipWinHits>0){
							fw.write(String.format("%d\t%.0f\n", i, ipWinHits));
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
				}
			}
			fw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
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
			//if(!needlefiltering || (counts[offset] <= perBaseMax)){
				count++;
				int binstart = (int)Math.max(0, ((double)(offset/winStep)));
				int binend = (int)(Math.min((double)(r.getEnd()-currReg.getStart())/winStep, numBins-1));
				for(int i=binstart; i<=binend; i++){
					land[i]++;
					sum++;
				}
			//}
		}
		return(land);
	}
	
	/* Binomial test for differences between two population proportions */
	protected double binomialSampleEquality(double X1, double X2, double n1, double n2){
		double P1 = X1/n1;
		double P2 = X2/n2;
		double P = (X1+X2)/(n1+n2);
		double Z0 = (P1-P2)/(Math.sqrt(P*(1-P)*((1/n1)+(1/n2))));
		return(Z0);
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
                }else
                	curr.setReadExtension(readExtension);
                IPhandles.add(curr);
				iphittot += curr.getHitCount();
			}
            System.out.print(String.format("%.0f reads loaded\n", iphittot));
			for(SeqLocator back : backs){
				System.out.print(String.format("%s\t", back.getExptName()));
				SeqExptHandler curr = new SeqExptHandler(gen, back);
                curr.setReadLength(readLength);
                if(shiftTags){
                	curr.setReadExtension(readShift*2);
                	curr.setReadShift(readShift);
                }else
                	curr.setReadExtension(readExtension);
                backhandles.add(curr);
				backhittot += curr.getHitCount();
			}
            System.out.print(String.format("%.0f reads loaded\n", backhittot));
            noBackground = backhandles.isEmpty();

		} catch (NotFoundException e) {
			e.printStackTrace();
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
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
	
	//Accessors
	public void setBinWidth(int w){winSize=w;}
	public void setBinStep(int w){winStep=w;}
////////////////////////////////////////////////////////////////////////////
	public static void printError(){
		System.err.println("Usage:\n " +
                "Expt2Vstep \n " +
                "Required:\n  " +
                "--species <organism name;genome version> \n  " +
                "--expt <experiment>\n  " +
                "--readlen --readextend --shifttags " +
                "--binwidth --binstep\n" +
                "--out <output file name>\n  "
                );
	}
}
