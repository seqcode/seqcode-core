package edu.psu.compbio.seqcode.projects.sequtils;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import edu.psu.compbio.seqcode.deepseq.StrandedBaseCount;
import edu.psu.compbio.seqcode.deepseq.experiments.ControlledExperiment;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.NamedRegion;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.projects.multigps.utilities.Utils;

/**
 * This class aims to count tags overlapping a set of regions. 
 *   
 * @author mahony
 *
 */
public class RegionTagCounts {

	private GenomeConfig gConfig;
	private ExptConfig eConfig;
	private List<Region> testRegs;
	private ExperimentManager manager=null;
	private String outName = "out";
	private boolean totalTagNorm=false; //normalize to total tags
	private boolean sigPropNorm=false; //normalize to signal proportion
	
	public RegionTagCounts(GenomeConfig gcon, ExptConfig econ, List<Region> regs, String out){
		gConfig = gcon;
		eConfig = econ;
		testRegs = regs;
		manager = new ExperimentManager(eConfig);
		outName = out;
	}
	
	public void setTotalTagNorm(boolean n){totalTagNorm=n;}
	public void setSigPropNorm(boolean n){sigPropNorm=n;}
	
	
	public void execute(){
		for(ExperimentCondition c : manager.getConditions()){
			for(ControlledExperiment rep : c.getReplicates()){
				System.err.println("Condition "+c.getName()+":\tRep "+rep.getName());
				double scaling = rep.getControlScaling();
				double sigStrength = 1-(scaling/(rep.getSignal().getHitCount()/rep.getControl().getHitCount()));
				double sigCount = sigStrength * rep.getSignal().getHitCount();
				
				try {
					FileWriter fw = new FileWriter(outName+".region-counts.txt");
													       			
					ChromRegionIterator chroms = new ChromRegionIterator(gConfig.getGenome());
					while(chroms.hasNext()){
						NamedRegion currentRegion = chroms.next();
						
						//Split the job up into chunks of 100Mbp
						for(int x=currentRegion.getStart(); x<=currentRegion.getEnd(); x+=100000000){
							int y = x+100000000; 
							if(y>currentRegion.getEnd()){y=currentRegion.getEnd();}
							Region currSubRegion = new Region(gConfig.getGenome(), currentRegion.getChrom(), x, y);
							
							List<StrandedBaseCount> hits = rep.getSignal().getBases(currSubRegion);
		                    double stackedTagStarts[] = makeTagStartLandscape(hits, currSubRegion);
		                    
		                    //Get coverage of points that lie within the current region
		                    for(Region r : testRegs){
		                    	if(currSubRegion.contains(r)){
			                    
									int offsetStart = inBounds(r.getStart()-currSubRegion.getStart(), 0, currSubRegion.getWidth()-1);
									int offsetEnd =inBounds(r.getEnd()-currSubRegion.getStart(), 0, currSubRegion.getWidth()-1);
									double sum=0;
									for(int o=offsetStart; o<=offsetEnd; o++){
										sum+=stackedTagStarts[o];
									}
									
									if(totalTagNorm)
										fw.write(r.getLocationString()+"\t"+String.format("%e", sum/rep.getSignal().getHitCount()) +"\n");
									else if(sigPropNorm)
										fw.write(r.getLocationString()+"\t"+String.format("%e", sum/sigCount) +"\n");
									else
										fw.write(r.getLocationString()+"\t"+String.format("%.0f", sum) +"\n");
									
		                    	}
							}
						}
					}
					System.out.print("\n");
					fw.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	protected double[] makeTagStartLandscape(List<StrandedBaseCount> hits, Region currReg){
		double[] counts = new double[(int)currReg.getWidth()+1];
        for(int i=0; i<=currReg.getWidth(); i++){counts[i]=0;}
        for(StrandedBaseCount r : hits){
        	if(r.getCoordinate()>=currReg.getStart() && r.getCoordinate()<=currReg.getEnd()){
	        	int offset=inBounds(r.getCoordinate()-currReg.getStart(),0,currReg.getWidth());
	        	counts[offset]+=r.getCount();
        	}
        }
        return(counts);
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
	
	public void close(){
		if(manager !=null)
			manager.close();
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
		if(args.length==0 || ap.hasKey("h")){
			System.err.println("RegionTagCounts:");
			System.err.println("Genome:" +
					"\t--species <Species;Genome>\n" +
					"\tOR\n" +
					"\t--geninfo <genome info file> AND --seq <fasta seq directory>\n" +
					"Coverage Testing:\n" +
					"\t--reg <region coords>\n" +
					"\t--out <output file root>\n" +
					"\t--signormcounts [flag to normalize counts by inferred signal proportion]\n" +
					"\t--totalnormcounts [flag to normalize counts by total tags]\n"
					);
		}else{
			
			GenomeConfig gcon = new GenomeConfig(args);
			ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		
			List<Region> testSites = new ArrayList<Region>();
			Collection<String> regFiles = Args.parseStrings(args, "reg");
			for(String rf : regFiles)
				testSites.addAll(Utils.loadRegionsFromFile(rf, gcon.getGenome(), -1));
			String outName = Args.parseString(args, "out", "out");
			boolean signorm = Args.parseFlags(args).contains("signormcounts");
			boolean totnorm = Args.parseFlags(args).contains("totalnormcounts");
			
			RegionTagCounts rct = new RegionTagCounts(gcon, econ, testSites, outName);
			rct.setSigPropNorm(signorm); rct.setTotalTagNorm(totnorm);
			rct.execute();
			rct.close();
		}	
		
	}

}
