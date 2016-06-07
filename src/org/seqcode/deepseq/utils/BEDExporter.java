package org.seqcode.deepseq.utils;

import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.NamedRegion;
import org.seqcode.genome.location.Region;
import org.seqcode.gsebricks.verbs.location.ChromRegionIterator;
import org.seqcode.gseutils.Args;
import org.seqcode.gseutils.NotFoundException;


public class BEDExporter {
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected ExperimentManager manager;
	protected int [] stackedHitCountsPos;
	protected int [] stackedHitCountsNeg;
	final int MAXSECTION=50000000;
	private int readLength=1;
	private String outName="out";
	
	public static void main(String[] args) throws SQLException, NotFoundException {
		GenomeConfig gconfig = new GenomeConfig(args);
		ExptConfig econfig = new ExptConfig(gconfig.getGenome(), args);
		if(args.length==0 || gconfig.helpWanted()){
			System.err.println("BEDExporter usage:\n" +
					GenomeConfig.getArgsList()+"\n"+
					ExptConfig.getArgsList()+"n"+
					"BEDExporter:\n"+
					"\t--readlen <read length>\n" +
					"\t--out <output file name>");
			System.exit(1);
		}else{
			String outName = Args.parseString(args,"out","out");
			int readLen = Args.parseInteger(args,"readlen",40);
			
			BEDExporter exporter = new BEDExporter(gconfig, econfig, outName, readLen);
			exporter.execute();
			exporter.close();
		}
	}
	
	
	public BEDExporter(GenomeConfig gcon, ExptConfig econ, String out, int rL) {
		gconfig = gcon;
		econfig = econ;
		manager = new ExperimentManager(econfig);
		
		readLength = rL;
	}
	
	public void execute(){
		for(ExperimentCondition c : manager.getConditions()){
			for(ControlledExperiment rep : c.getReplicates()){
				System.err.println("Condition "+c.getName()+":\tRep "+rep.getName());
				try {
					FileWriter fw = new FileWriter(outName+"."+c.getName()+"."+rep.getName()+".idx");
					fw.write("chrom\tindex\tforward\treverse\tvalue\n");
					
					double basesDone=0, printStep=10000000,  numPrint=0;
								       			
					ChromRegionIterator chroms = new ChromRegionIterator(gconfig.getGenome());
					while(chroms.hasNext()){
						NamedRegion currentRegion = chroms.next();
						
						//Split the job up into chunks of 100Mbp
						for(int x=currentRegion.getStart(); x<=currentRegion.getEnd(); x+=100000000){
							int y = x+100000000; 
							if(y>currentRegion.getEnd()){y=currentRegion.getEnd();}
							Region currSubRegion = new Region(gconfig.getGenome(), currentRegion.getChrom(), x, y);
							
							List<StrandedBaseCount> hits = rep.getSignal().getBases(currSubRegion);
		                    double stackedHitCountsPos[] = make5PrimeLandscape(hits, currSubRegion, '+');
		                    double stackedHitCountsNeg[] = make5PrimeLandscape(hits, currSubRegion, '-');
		                    
		                    //Scan regions
							for(int i=currSubRegion.getStart(); i<currSubRegion.getEnd(); i++){
								int offset = i-currSubRegion.getStart();
								double posHits=stackedHitCountsPos[offset];
								double negHits=stackedHitCountsNeg[offset];
								
								if(posHits>0 || negHits>0){
									for(int hitc=0; hitc<(int)posHits; hitc++){
										fw.write("chr"+currSubRegion.getChrom()+"\t"+i+"\t"+(i+readLength)+"\t"+"+"+"\n");
									}
									for(int hitc=0; hitc<(int)negHits; hitc++){
										fw.write("chr"+currSubRegion.getChrom()+"\t"+(i+1-readLength)+"\t"+(i+1)+"\t"+"-"+"\n");
									}
								}
								//Print out progress
								basesDone++;
								if(basesDone > numPrint*printStep){
									if(numPrint%10==0){System.out.print(String.format("(%.0f)", (numPrint*printStep)));}
									else{System.out.print(".");}
									if(numPrint%50==0 && numPrint!=0){System.out.print("\n");}
									numPrint++;
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
	
	public void close(){
		manager.close();
	}
	
	protected double[] make5PrimeLandscape(List<StrandedBaseCount> hits, Region currReg, char strand){
		double[] startcounts = new double[(int)currReg.getWidth()+1];
        for(int i=0; i<=currReg.getWidth(); i++){startcounts[i]=0;}
        for(StrandedBaseCount r : hits){
            if(strand=='.' || r.getStrand()==strand){
            	if(r.getCoordinate()>=currReg.getStart() && r.getCoordinate()<=currReg.getEnd()){
	                int offset5=inBounds(r.getCoordinate()-currReg.getStart(),0,currReg.getWidth());
	                startcounts[offset5]+=r.getCount();
            	}
            }
        }
        return(startcounts);
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

}
