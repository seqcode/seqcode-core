package org.seqcode.deepseq.utils;

import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import org.seqcode.deepseq.StrandedBaseCount;
import org.seqcode.deepseq.StrandedBaseCountFilterByBase;
import org.seqcode.deepseq.experiments.ControlledExperiment;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.NamedRegion;
import org.seqcode.genome.location.Region;
import org.seqcode.gsebricks.verbs.location.ChromRegionIterator;
import org.seqcode.utils.Args;
import org.seqcode.utils.NotFoundException;


/**
 * Outputs a GeneTrack index format file for a deep-seq experiment.
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class IDXExporter {
	protected ExperimentManager manager;
	protected GenomeConfig gcon=null;
	protected ExptConfig econ=null;
	protected int [] stackedHitCountsPos;
	protected int [] stackedHitCountsNeg;
	protected String outName="out";
	protected char baseLimit='.'; //Only use tags with this character at baseLimitRelPosition relative to 5' end (. = use all tags)
	protected int baseLimitRelPosition=0;
	protected StrandedBaseCountFilterByBase sbcFilter;
	protected boolean filterByBase=false;
	
	
	public static void main(String[] args) throws SQLException, NotFoundException {
		
		IDXExporter exporter = new IDXExporter(args);
		exporter.execute();
		exporter.close();
			
	}
	
	
	public IDXExporter(String [] args) {
		if(args.length==0){
			System.err.println("IDXExporter usage:\n" +
					GenomeConfig.getArgsList()+"\n"+
					ExptConfig.getArgsList()+"n"+
					"IDXExporter:\n"+
					"\t--out <output file name>\n"+
					"\t--baselimit <./A/C/G/T: only use tags with this base at below position>\n" +
					"\t--baselimitposition <only use tags with above base at this position>\n");
			System.exit(1);
		}else{
			gcon = new GenomeConfig(args);
			econ = new ExptConfig(gcon.getGenome(), args);
			manager = new ExperimentManager(econ);
			
			outName = Args.parseString(args,"out",outName);
			
			baseLimit = Args.parseString(args, "baselimit", ".").charAt(0);
			baseLimitRelPosition = Args.parseInteger(args, "baselimitposition", 0);
			if(baseLimit!='.'){
				sbcFilter = new StrandedBaseCountFilterByBase(gcon, baseLimit, baseLimitRelPosition);
				filterByBase=true;
			}
		}
	}
	
	public void execute(){
		for(ExperimentCondition c : manager.getConditions()){
			for(ControlledExperiment rep : c.getReplicates()){
				System.err.println("Condition "+c.getName()+":\tRep "+rep.getName());
				try {
					FileWriter fw = new FileWriter(outName+"."+c.getName()+"."+rep.getName()+".idx");
					fw.write("chrom\tindex\tforward\treverse\tvalue\n");
					
					double basesDone=0, printStep=10000000,  numPrint=0;
								       			
					ChromRegionIterator chroms = new ChromRegionIterator(gcon.getGenome());
					while(chroms.hasNext()){
						NamedRegion currentRegion = chroms.next();
						
						//Split the job up into chunks of 100Mbp
						for(int x=currentRegion.getStart(); x<=currentRegion.getEnd(); x+=100000000){
							int y = x+100000000; 
							if(y>currentRegion.getEnd()){y=currentRegion.getEnd();}
							Region currSubRegion = new Region(gcon.getGenome(), currentRegion.getChrom(), x, y);
							
							List<StrandedBaseCount> hits = rep.getSignal().getBases(currSubRegion);
							if(filterByBase)
								hits = sbcFilter.execute(currSubRegion, hits);
		                    double stackedHitCountsPos[] = make5PrimeLandscape(hits, currSubRegion, '+');
		                    double stackedHitCountsNeg[] = make5PrimeLandscape(hits, currSubRegion, '-');
		                    
		                    //Scan regions
							for(int i=currSubRegion.getStart(); i<currSubRegion.getEnd(); i++){
								int offset = i-currSubRegion.getStart();
								double posHits=stackedHitCountsPos[offset];
								double negHits=stackedHitCountsNeg[offset];
								double sum = posHits+negHits;
								
								if(posHits>0 || negHits>0){
									fw.write("chr"+currSubRegion.getChrom()+"\t"+i+"\t"+String.format("%.0f\t%.0f\t%.0f", posHits, negHits, sum) +"\n");
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
		if(manager!=null)
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
