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
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.ChromRegionIterator;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.projects.multigps.utilities.Utils;

/**
 * This class aims to test the sequencing read coverage at a set of points (e.g. coverage overlapping a set of SNPs). 
 *   
 * @author mahony
 *
 */
public class ReadCoverageTester {

	private GenomeConfig gConfig;
	private ExptConfig eConfig;
	private List<Point> testSites;
	private ExperimentManager manager=null;
	private String outName = "out";
	private Integer readLen = 36;
	
	public ReadCoverageTester(GenomeConfig gcon, ExptConfig econ, List<Point> sites, int readLength){
		gConfig = gcon;
		eConfig = econ;
		testSites = sites;
		manager = new ExperimentManager(eConfig);
		readLen = readLength;
		
	}
	
	
	public void execute(){
		for(ExperimentCondition c : manager.getConditions()){
			for(ControlledExperiment rep : c.getReplicates()){
				System.err.println("Condition "+c.getName()+":\tRep "+rep.getName());
				try {
					FileWriter fw = new FileWriter(outName+"."+c.getName()+"."+rep.getName()+".point-coverage.txt");
													       			
					ChromRegionIterator chroms = new ChromRegionIterator(gConfig.getGenome());
					while(chroms.hasNext()){
						NamedRegion currentRegion = chroms.next();
						
						//Split the job up into chunks of 100Mbp
						for(int x=currentRegion.getStart(); x<=currentRegion.getEnd(); x+=100000000){
							int y = x+100000000; 
							if(y>currentRegion.getEnd()){y=currentRegion.getEnd();}
							Region currSubRegion = new Region(gConfig.getGenome(), currentRegion.getChrom(), x, y);
							
							List<StrandedBaseCount> hits = rep.getSignal().getBases(currSubRegion);
		                    double stackedReadCoverage[] = makeReadCoverageLandscape(hits, currSubRegion);
		                    
		                    //Get coverage of points that lie within the current region
		                    for(Point pt : testSites){
		                    	if(currSubRegion.contains(pt)){
			                    
									int offset = pt.getLocation()-currSubRegion.getStart();
									double sum=stackedReadCoverage[offset];
									
									fw.write("chr"+pt.getChrom()+":"+pt.getLocation()+"\t"+String.format("%.0f", sum) +"\n");
									
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
	
	protected double[] makeReadCoverageLandscape(List<StrandedBaseCount> hits, Region currReg){
		double[] counts = new double[(int)currReg.getWidth()+1];
        for(int i=0; i<=currReg.getWidth(); i++){counts[i]=0;}
        for(StrandedBaseCount r : hits){
        	if(r.getCoordinate()>=currReg.getStart() && r.getCoordinate()<=currReg.getEnd()){
	        	if(r.getStrand()=='+'){
	        		int offsetStart=inBounds(r.getCoordinate()-currReg.getStart(),0,currReg.getWidth());
	        		int offsetEnd=inBounds(r.getCoordinate()+readLen-currReg.getStart(),0,currReg.getWidth());
	        		for(int o=offsetStart; o<offsetEnd; o++)
	        			counts[o]+=r.getCount();
	            }else{
	            	int offsetStart=inBounds(r.getCoordinate()-readLen+1-currReg.getStart(),0,currReg.getWidth());
	            	int offsetEnd=inBounds(r.getCoordinate()-currReg.getStart()+1,0,currReg.getWidth());
	            	for(int o=offsetStart; o<offsetEnd; o++)
	            		counts[o]+=r.getCount();
	            }
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
			System.err.println("ReadCoverageTester:");
			System.err.println("Genome:" +
					"\t--species <Species;Genome>\n" +
					"\tOR\n" +
					"\t--geninfo <genome info file> AND --seq <fasta seq directory>\n" +
					"Coverage Testing:\n" +
					"\t--sites <test site coords>\n" +
					"\t--readlen <read length>\n"
					);
		}else{
			
			GenomeConfig gcon = new GenomeConfig(args);
			ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		
			List<Point> testSites = new ArrayList<Point>();
			Collection<String> siteFiles = Args.parseStrings(args, "sites");
			for(String sf : siteFiles)
				testSites.addAll(Utils.loadPointsFromFile(sf, gcon.getGenome()));
			
			Integer readL = Args.parseInteger(args, "readlen", 36);
			
			ReadCoverageTester rct = new ReadCoverageTester(gcon, econ, testSites, readL);
			rct.execute();
			rct.close();
		}	
		
	}

}
