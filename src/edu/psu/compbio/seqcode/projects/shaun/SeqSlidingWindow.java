package edu.psu.compbio.seqcode.projects.shaun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.chipseq.ChipSeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.deepseq.DeepSeqExpt;
import edu.psu.compbio.seqcode.gse.deepseq.ReadHit;
import edu.psu.compbio.seqcode.gse.ewok.verbs.RegionParser;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;

public class SeqSlidingWindow {

	protected DeepSeqExpt signal;
	protected DeepSeqExpt control;
	protected int regionExtension = 1000;
	protected Genome gen=null;
	protected boolean noControl=true;
	protected boolean dbconnected=true;
	protected int binWidth=500;
	protected int binStep=250;
	protected String regFile=null;
	protected int pbMax = 2;
	protected double[] landscape =null;
	protected double[] startcounts =null;
	protected boolean needleFiltering = true;
	
	public SeqSlidingWindow(String[] args){
		
		if(args.length==0){
			System.out.println("SeqSlidingWindow:\n" +
					"--species <org;species>\n" +
					"--binwidth <width>\n" +
					"--binstep <step>\n" +
					"--regfile <region file>\n" +
					"--rdbexpt <ReadDB identifier>\n" +
					"--rdbctrl <ReadDB identifier>\n");
		}else{
			ArgParser ap = new ArgParser(args);
			try {
				Pair<Organism, Genome> pair = Args.parseGenome(args);
				gen = pair.cdr();
			} catch (NotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			binWidth = Args.parseInteger(args,"binwidth",binWidth);
	        binStep=Args.parseInteger(args,"binstep",binStep);
	        regFile = Args.parseString(args, "regfile", "NONE");
			List<ChipSeqLocator> dbexpts = Args.parseChipSeq(args,"dbexpt");
	        List<ChipSeqLocator> dbctrls = Args.parseChipSeq(args,"dbctrl");
	        List<ChipSeqLocator> rdbexpts = Args.parseChipSeq(args,"rdbexpt");
	        List<ChipSeqLocator> rdbctrls = Args.parseChipSeq(args,"rdbctrl");
	        List<File> expts = Args.parseFileHandles(args, "expt");
	        List<File> ctrls = Args.parseFileHandles(args, "ctrl");
	        boolean nonUnique = ap.hasKey("nonunique") ? true : false;
	        String fileFormat = Args.parseString(args, "format", "ELAND");
	    	int readLength = Args.parseInteger(args,"readlen",32);
	        if(expts.size()>0 && dbexpts.size() == 0 && rdbexpts.size()==0){
	        	signal = new DeepSeqExpt(gen, expts, nonUnique, fileFormat, readLength);
	        }else if(dbexpts.size()>0 && expts.size() == 0){
	        	signal = new DeepSeqExpt(gen, dbexpts, "db", readLength);
	        	dbconnected=true;
	        }else if(rdbexpts.size()>0 && expts.size() == 0){
	            	signal = new DeepSeqExpt(gen, rdbexpts, "readdb", readLength);
	            	dbconnected=true;
	        }else{System.err.println("Must provide either an aligner output file or Gifford lab DB experiment name for the signal experiment (but not both)");System.exit(1);}
	        
	        if(ctrls.size()>0 && dbctrls.size() == 0 && dbctrls.size()==0){
	        	control = new DeepSeqExpt(gen, ctrls, nonUnique, fileFormat, readLength); noControl=false;
	        }else if(dbctrls.size()>0 && ctrls.size() == 0){
	        	control = new DeepSeqExpt(gen, dbctrls, "db", readLength); noControl=false;
	        	dbconnected=true;
	        }else if(rdbctrls.size()>0 && ctrls.size() == 0){
	        	control = new DeepSeqExpt(gen, rdbctrls, "readdb", readLength); noControl=false;
	        	dbconnected=true;
	        }else{
	        	if(dbctrls.size()>0 && ctrls.size()>0){
	        		System.err.println("Cannot mix files and db loading yet...");;System.exit(1);
	        	}else{
	        		noControl=true; control=null;
	        	}
	        }
	       
		
			//Print some info
			System.err.println("Signal hit count: "+(int)signal.getHitCount()+", weight: "+(int)signal.getWeightTotal());
			if(!noControl)
				System.err.println("Control hit count: "+(int)control.getHitCount()+", weight: "+(int)control.getWeightTotal());
		}
   		
	}
	
	public void execute(){
		//Naive scaling factor for now 
		if(!noControl)
			control.setScalingFactor(scalingRatioByMedian(10000));
		
		Iterator<Region> testRegions=loadRegionsFromFile(regFile).iterator();
		
		while (testRegions.hasNext()) {
			Region currentRegion = testRegions.next();
			System.out.println("#"+currentRegion.getLocationString());
			
			List<ReadHit> ipHits = signal.loadHits(currentRegion);
			makeHitLandscape(ipHits, currentRegion, pbMax, '.');
			double ipStackedHitCounts[] = landscape.clone();
			double backStackedHitCounts[] = null;
			if (!noControl) {
				List<ReadHit> backHits = control.loadHits(currentRegion);
            	makeHitLandscape(backHits, currentRegion, pbMax, '.');
            	backStackedHitCounts = landscape.clone();
            }
            
			int currBin=0;
			for(int i=currentRegion.getStart(); i<currentRegion.getEnd()-(int)binWidth; i+=(int)binStep){
				Region currReg = new Region(gen, currentRegion.getChrom(), i, i+binWidth);
				
				double sigweight=ipStackedHitCounts[currBin];
				double ctrlweight= noControl ? 0 : backStackedHitCounts[currBin]*control.getScalingFactor();
				
				double sigavg = sigweight>0 ? sigweight/(double)currReg.getWidth() : 0;
				double ctrlavg = ctrlweight>0 ? ctrlweight/(double)currReg.getWidth() : 0;
				double fold = ctrlweight>0 ? sigweight/ctrlweight : sigweight;
				double logfold = fold>0 ? Math.log(fold) : 0;
				
				//Print here
				if(!noControl){
					System.out.println(currReg.getLocationString()+"\t"+sigweight+"\t"+ctrlweight+"\t"+fold+"\t"+logfold);
				}else{
					System.out.println(currReg.getLocationString()+"\t"+sigweight);
				}
				currBin++;
			}	
		}
	}
	
	//Load a set of regions from a file
	public ArrayList<Region> loadRegionsFromFile(String filename){
		ArrayList<Region> regs = new ArrayList<Region>();
		try{
			File pFile = new File(filename);
			if(!pFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
			BufferedReader reader = new BufferedReader(new FileReader(pFile));
			String line;
			while ((line = reader.readLine()) != null) {
			    line = line.trim();
			    String[] words = line.split("\\s+");
			    if(words.length>=1 && words[0].contains(":")){
				RegionParser parser = new RegionParser(gen);
				Region r = parser.execute(words[0]);
				if(r!=null){regs.add(r);}
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
	
	//Makes integer arrays corresponding to the read landscape over the current region
	protected void makeHitLandscape(List<ReadHit> ipHits, Region currReg, int perBaseMax, char strand){
		int numBins = (int)(currReg.getWidth()/binStep);
		int [] counts = new int[currReg.getWidth()+1];
		//double [] land = new double[numBins+1];
		landscape = new double[numBins+1];
		startcounts = new double[numBins+1];
		for(int i=0; i<=numBins; i++){landscape[i]=0; startcounts[i]=0;}
		for(int i=0; i<=currReg.getWidth(); i++){counts[i]=0;}
		for(ReadHit r : ipHits){
			if(strand=='.' || r.getStrand()==strand){
				int offset=inBounds(r.getStart()-currReg.getStart(),0,currReg.getWidth());
				counts[offset]++;//small issue here... counts will not be corrected for scalingFactor in control
				if(!needleFiltering ||  (counts[offset] <= perBaseMax)){
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
	
	public static void main(String[] args){
		SeqSlidingWindow ssw = new SeqSlidingWindow(args);
		ssw.execute();
	}
	
	public double scalingRatioByMedian(int windowSize){
		double scalingRatio=1;
		
	    ArrayList<Float> ratios = new ArrayList<Float>();
		for(String chrom:gen.getChromList()) {
            int chrlen = gen.getChromLength(chrom);
            for (int start = 0; start  < chrlen - windowSize; start += windowSize) {
                Region r = new Region(gen, chrom, start, start + windowSize);
                double countA = signal.sumWeights(r);
                double countB = control.sumWeights(r);
                ratios.add((float)(countA / countB));
            }
        }
        Collections.sort(ratios);
		scalingRatio = ratios.get(ratios.size() / 2);
        System.err.println(String.format("Scaling ratio estimated as %.3f based on %d regions of size %d",
        		scalingRatio, ratios.size(), windowSize));
		
		return(scalingRatio);
	}
	
}
