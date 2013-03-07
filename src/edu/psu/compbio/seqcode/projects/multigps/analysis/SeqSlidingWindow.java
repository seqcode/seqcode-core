package edu.psu.compbio.seqcode.projects.multigps.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.deepseq.DeepSeqExpt;
import edu.psu.compbio.seqcode.gse.deepseq.ReadHit;
import edu.psu.compbio.seqcode.gse.ewok.verbs.ChromosomeGenerator;
import edu.psu.compbio.seqcode.gse.ewok.verbs.RegionParser;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ExperimentSet;
import edu.psu.compbio.seqcode.projects.multigps.framework.Config;
import edu.psu.compbio.seqcode.projects.multigps.framework.StrandedBaseCount;

public class SeqSlidingWindow {

	protected ExperimentManager exptMan;
	protected ExperimentSet eset;
	protected Config config;
	protected int regionExtension = 1000;
	protected Genome gen=null;
	protected int binWidth=500;
	protected int binStep=250;
	protected String regFile=null;
	protected int pbMax = 2;
	protected double[] landscape =null;
	protected double[] startcounts =null;
	protected boolean needleFiltering = true;
	
	public SeqSlidingWindow(String[] args){
		
		if(args.length==0){
			System.err.println("SeqSlidingWindow:\n" +
					"--species <org;species>\n" +
					"\tOR\n" +
					"\t--geninfo <genome info file>\n" +
					"--binwidth <width>\n" +
					"--binstep <step>\n" +
					"--regfile <region file>\n" +
					"--rdbexpt <ReadDB identifier>\n" +
					"--rdbctrl <ReadDB identifier>\n");
		}else{
			config = new Config(args);
			gen = config.getGenome();
			exptMan = new ExperimentManager(config);
			eset = exptMan.getExperimentSet();
			
			if(eset.getConditions().size()==1 && eset.getConditions().get(0).getReplicates().size()==1){
				
				System.err.println("Conditions:\t"+eset.getConditions().size());
				for(ExperimentCondition c : eset.getConditions()){
					System.err.println("Condition "+c.getName()+":\t#Replicates:\t"+c.getReplicates().size());
				}
				for(ExperimentCondition c : eset.getConditions()){
					for(ControlledExperiment r : c.getReplicates()){
						System.err.println("Condition "+c.getName()+":\tRep "+r.getName());
						if(r.getControl()==null)
							System.err.println("\tSignal:\t"+r.getSignal().getHitCount());
						else
							System.err.println("\tSignal:\t"+r.getSignal().getHitCount()+"\tControl:\t"+r.getControl().getHitCount());
					}
				}
	
				
				binWidth = Args.parseInteger(args,"binwidth",binWidth);
		        binStep=Args.parseInteger(args,"binstep",binStep);
		        regFile = Args.parseString(args, "regfile", "NONE");
			}else{
				System.err.println("Only one condition & one replicate handled in this version");
				System.exit(1);
			}
		}
   		
	}
	
	public void execute(){
		ControlledExperiment expt = eset.getConditions().get(0).getReplicates().get(0);
		Iterator<Region> testRegions;
		
		if(!regFile.equals("NONE")){
			testRegions=loadRegionsFromFile(regFile).iterator();
		}else{
			testRegions = new ChromosomeGenerator().execute(gen);
		}

		while (testRegions.hasNext()) {
			Region currentRegion = testRegions.next();
			System.out.println("#"+currentRegion.getLocationString());
			//Split the job up into large chunks
            for(int x=currentRegion.getStart(); x<=currentRegion.getEnd(); x+=config.MAXSECTION){
                int y = x+config.MAXSECTION; 
                if(y>currentRegion.getEnd()){y=currentRegion.getEnd();}
                Region currSubRegion = new Region(gen, currentRegion.getChrom(), x, y);
			
				List<StrandedBaseCount> ipHits = expt.getSignal().getUnstrandedBases(currSubRegion);
				makeHitStartLandscape(ipHits, currSubRegion, pbMax, '.');
				double ipStackedHitCounts[] = landscape.clone();
				double backStackedHitCounts[] = null;
				if (expt.hasControl()) {
					List<StrandedBaseCount> backHits = expt.getControl().getUnstrandedBases(currSubRegion);
	            	makeHitStartLandscape(backHits, currSubRegion, pbMax, '.');
	            	backStackedHitCounts = landscape.clone();
	            }
	            
				int currBin=0;
				for(int i=currSubRegion.getStart(); i<currSubRegion.getEnd()-(int)binWidth; i+=(int)binStep){
					Region currReg = new Region(gen, currSubRegion.getChrom(), i, i+binWidth);
					
					double sigweight=ipStackedHitCounts[currBin];
					double ctrlweight= expt.hasControl() ? backStackedHitCounts[currBin]*expt.getControlScaling() : 0 ;
					
					double sigavg = sigweight>0 ? sigweight/(double)currReg.getWidth() : 0;
					double ctrlavg = ctrlweight>0 ? ctrlweight/(double)currReg.getWidth() : 0;
					double fold = ctrlweight>0 ? sigweight/ctrlweight : sigweight;
					double logfold = fold>0 ? Math.log(fold) : 0;
					
					//Print here
					if(expt.hasControl()){
						System.out.println(currReg.getLocationString()+"\t"+sigweight+"\t"+ctrlweight+"\t"+fold+"\t"+logfold);
					}else{
						System.out.println(currReg.getLocationString()+"\t"+sigweight);
					}
					currBin++;
				}
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
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return(regs);
	}
	
	//Makes integer arrays corresponding to the read landscape over the current region
	protected void makeHitStartLandscape(List<StrandedBaseCount> ipHits, Region currReg, int perBaseMax, char strand){
		int numBins = (int)(currReg.getWidth()/binStep);
		int [] counts = new int[currReg.getWidth()+1];
		landscape = new double[numBins+1];
		for(int i=0; i<=numBins; i++){landscape[i]=0; }
		for(int i=0; i<=currReg.getWidth(); i++){counts[i]=0;}
		for(StrandedBaseCount r : ipHits){
			if(strand=='.' || r.getStrand()==strand){
				int offset=inBounds(r.getCoordinate()-currReg.getStart(),0,currReg.getWidth());
				int binstart = inBounds((int)((double)(offset)/binStep), 0, numBins);
				
				float weight=0;
				if(needleFiltering){
					weight = counts[offset]+r.getCount()<=perBaseMax ? r.getCount() : perBaseMax - counts[offset];
					counts[offset]+=r.getCount();
					if(counts[offset]>perBaseMax)
						counts[offset]=perBaseMax;
				}else{
					weight = r.getCount();
				}
				
				landscape[binstart]+=weight;
			}
		}
	}
	protected final int inBounds(int x, int min, int max){
		if(x<min){return min;}
		if(x>max){return max;}
		return x;
	}
	
	public void close(){
		exptMan.close();
	}
	
	public static void main(String[] args){
		SeqSlidingWindow ssw = new SeqSlidingWindow(args);
		ssw.execute();
		ssw.close();
	}
		
}
