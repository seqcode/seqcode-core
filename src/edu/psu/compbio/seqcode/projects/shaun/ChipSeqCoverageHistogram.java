package edu.psu.compbio.seqcode.projects.shaun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.chipseq.ChipSeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.general.NamedRegion;
import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.deepseq.BackgroundCollection;
import edu.psu.compbio.seqcode.gse.deepseq.DeepSeqExpt;
import edu.psu.compbio.seqcode.gse.deepseq.PoissonBackgroundModel;
import edu.psu.compbio.seqcode.gse.deepseq.ReadHit;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.PointParser;
import edu.psu.compbio.seqcode.gse.ewok.verbs.RegionParser;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;

public class ChipSeqCoverageHistogram {

	public Genome gen;
	public DeepSeqExpt signal=null;
	public int readLen, readExt;
	public int binWidth, binStep;
	public double confLevel=-1;
	public boolean printHisto = false;
	public static int MAXSECTION=100000000;
	private double[] landscape=null;
	private double[] startcounts=null;
	protected BackgroundCollection signalPerBaseBack=new BackgroundCollection();
	protected double perBaseLogConf=-7;
	protected double mappableGenome = 0.8;
	protected boolean needlefiltering=true;
	public int maxRC = 1000;
	public Integer[] readCounts = new Integer[maxRC+1];
	private ArrayList<Region> towers = new ArrayList<Region>();	
	
	public static void main(String[] args) {
		try{
			ArgParser ap = new ArgParser(args);
			Pair<Organism, Genome> pair = Args.parseGenome(args);
			
	        List<ChipSeqLocator> expts = Args.parseChipSeq(args,"expt");
	        int rL = 36, rE=164, ws=50, wo=25;
	        double cl=-1;
	        if(ap.hasKey("readlen")){ rL = new Integer(ap.getKeyValue("readlen")).intValue();}
	        if(ap.hasKey("readext")){ rE = new Integer(ap.getKeyValue("readext")).intValue();}
	        if(ap.hasKey("win")){ ws = new Integer(ap.getKeyValue("win")).intValue();}
	        if(ap.hasKey("off")){ wo = new Integer(ap.getKeyValue("off")).intValue();}
	        if(ap.hasKey("conflevel")){ cl = new Double(ap.getKeyValue("conflevel")).doubleValue();}
	        String towerFile = ap.hasKey("towers") ? ap.getKeyValue("towers") : null;
	    	boolean hasTowers = towerFile==null ? false : true;
	        if (expts.size() == 0) {
	            System.err.println("Usage:\n " +
	                               "ChipSeqCoverageHistogram \n" +
	                               "--species <organism name;genome version> \n"+
	                               "--expt <solexa expt> \n" +
	                               "--win <winsize> \n" +
	                               "--off <winoffset>\n" +
	                               "--readlen <readlength>\n" +
	                               "--readext <readext>\n" +
	                               "--printhisto [flag]\n" +
	                               "--conflevel <fraction of bins>\n" +
	                               "--towers <file containing towers in Shaun's format>\n");
	            return;
	        }else{
	        	ChipSeqCoverageHistogram csch = new ChipSeqCoverageHistogram(pair.cdr(), expts, ws, wo, rL, rE);
	        	if(ap.hasKey("printhisto"))
	        		csch.setPrintHisto(true);
	        	if(ap.hasKey("conflevel"))
	        		csch.setConfLevel(cl);
	        	if(hasTowers)
					csch.loadTowers(towerFile);
	        	
	        	csch.execute();
	        	
	        	csch.close();
	        }
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public ChipSeqCoverageHistogram(Genome g, List<ChipSeqLocator> e, int ws, int wo, int rL, int rE){
		gen = g;
		signal = new DeepSeqExpt(gen, e, "readdb", readLen);
		signal.setThreePrimeExt(rE);
		readLen = rL;
		readExt = rE;
		binWidth = ws;
		binStep = wo;
		signalPerBaseBack=new BackgroundCollection();
		signalPerBaseBack.addBackgroundModel(new PoissonBackgroundModel(-1, perBaseLogConf, signal.getWeightTotal(), 1, gen.getGenomeLength(), mappableGenome, 1, 1, '.', 1, true));
	}
	public void setPrintHisto(boolean ph){printHisto=ph;}
	public void setConfLevel(double cl){confLevel=cl;}
	
	public void execute(){
		//initialize count array
		for(int i=0; i<=maxRC; i++){
			readCounts[i]=0;
		}
		double totalBins=0;
		
		//Scan chromosomes
		Iterator<NamedRegion> chroms = new ChromRegionIterator(gen);
		while (chroms.hasNext()) {
			Region currentRegion = chroms.next();
			for(int x=currentRegion.getStart(); x<=currentRegion.getEnd(); x+=MAXSECTION){
				int y = x+MAXSECTION; 
				if(y>currentRegion.getEnd()){y=currentRegion.getEnd();}
				Region currSubRegion = new Region(gen, currentRegion.getChrom(), x, y);
				ArrayList<ReadHit> ipHits = new ArrayList<ReadHit>();
				ipHits.addAll(signal.loadExtHits(currSubRegion));
				
				ArrayList<Region> currTowers = new ArrayList<Region>();
				for(Region t : towers)
					if(currSubRegion.overlaps(t))
						currTowers.add(t);
				
				makeHitLandscape(ipHits, currSubRegion, signalPerBaseBack.getMaxThreshold('.'), '.');
				double ipStackedHitCounts[] = landscape.clone();
				//Scan regions
	            int currBin=0;
				for(int i=currSubRegion.getStart(); i<currSubRegion.getEnd()-(int)binWidth; i+=(int)binStep){
					boolean isTower=false;
					Region currWin = new Region(currSubRegion.getGenome(), currSubRegion.getChrom(), i, i+binWidth-1);
					for(Region t : currTowers)
						if(currWin.overlaps(t))
							isTower=true;
					if(!isTower){
						Integer ipWinHits=(int)ipStackedHitCounts[currBin];
						if(ipWinHits>maxRC)
							ipWinHits=maxRC;
						readCounts[ipWinHits]++;
						totalBins++;
					}
					currBin++;					
				}
			}
		}
		if(printHisto){
			for(int i=0; i<=maxRC; i++){
				System.out.println(i+"\t"+readCounts[i]);
			}
		}
		if(confLevel>0 && confLevel<=1){
			double sumBins=0;
			for(int i=0; i<=maxRC; i++){
				sumBins+=readCounts[i];
				if(sumBins/totalBins >confLevel){
					System.out.println("ConfLevel_"+confLevel+"\t"+i);
					i=maxRC+1;
				}
			}
		}
	}
	
	public void close(){
		if(signal!=null){
			signal.closeLoaders();
		}
	}
	//Makes integer arrays corresponding to the read landscape over the current region
	protected void makeHitLandscape(ArrayList<ReadHit> hits, Region currReg, int perBaseMax, char strand){
		int numBins = (int)(currReg.getWidth()/binStep);
		int [] counts = new int[currReg.getWidth()+1];
		//double [] land = new double[numBins+1];
		landscape = new double[numBins+1];
		startcounts = new double[numBins+1];
		for(int i=0; i<=numBins; i++){landscape[i]=0; startcounts[i]=0;}
		for(int i=0; i<=currReg.getWidth(); i++){counts[i]=0;}
		for(ReadHit r : hits){
			if(strand=='.' || r.getStrand()==strand){
				int offset=inBounds(r.getStart()-currReg.getStart(),0,currReg.getWidth());
				counts[offset]++;//small issue here... counts will not be corrected for scalingFactor in control
				if(!needlefiltering || (counts[offset] <= perBaseMax)){
					int binstart = inBounds((int)((double)offset/binStep), 0, numBins);
					int binend = inBounds((int)((double)(r.getEnd()-currReg.getStart())/binStep), 0, numBins);
					for(int i=binstart; i<=binend; i++){
						landscape[i]+=r.getWeight();
					}
					if(r.getStrand()=='+')
						startcounts[binstart]+=r.getWeight();
					else
						startcounts[binend]+=r.getWeight();
				}
			}
		}
		//return(land);
	}
	protected final int inBounds(int x, int min, int max){
		if(x<min){return min;}
		if(x>max){return max;}
		return x;
	}
	public void loadTowers(String fname){
		towers = loadRegionsFromPeakFile(fname, -1);
	}
	//Load a set of regions from a peak file
	public ArrayList<Region> loadRegionsFromPeakFile(String filename, int win){
		ArrayList<Region> regs = new ArrayList<Region>();
		try{
			File pFile = new File(filename);
			if(!pFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(pFile));
	        String line = reader.readLine(); //Ignore first line
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            
	            if(words.length>=3 && win!=-1){
	                PointParser pparser = new PointParser(gen);
	            	Point p = pparser.execute(words[2]);
	            	int rstart = p.getLocation()-(win/2)<1 ? 1:p.getLocation()-(win/2);
                	int rend = p.getLocation()+(win/2)>gen.getChromLength(p.getChrom()) ? gen.getChromLength(p.getChrom()):p.getLocation()+(win/2)-1;
                	Region r = new Region(p.getGenome(), p.getChrom(), rstart, rend);
                	regs.add(r);
                }else if(words.length>=1){
	            	RegionParser parser = new RegionParser(gen);
	            	Region q = parser.execute(words[0]);
	            	if(win!=-1){
	            		int rstart = q.getMidpoint().getLocation()-(win/2)<1 ? 1:q.getMidpoint().getLocation()-(win/2);
	                	int rend = q.getMidpoint().getLocation()+(win/2)>gen.getChromLength(q.getChrom()) ? gen.getChromLength(q.getChrom()):q.getMidpoint().getLocation()+(win/2)-1;
	                	Region r = new Region(q.getGenome(), q.getChrom(), rstart, rend);
	                	if(r!=null){regs.add(r);}
	            	}else{
	            		if(q!=null){regs.add(q);}
	            	}
	            }
	        }reader.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return(regs);
	}
}
