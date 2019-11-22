package org.seqcode.projects.galaxyexo;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.deepseq.experiments.Sample;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.gseutils.ArgParser;
import org.seqcode.gseutils.Args;

/**
 * Utility to sort locations by tag counts.  
 * If multiple samples are provided, experiments are normalized by total tags and locations are sorted based on average tag counts. 
 * 
 * input : Signal experiment
 * 		 : Points
 *       : Window size
 * output: Regions sorted by tag counts
 * 
 * @author naomi yamada
 */

public class RegionCountSorter {
	protected GenomeConfig gconf;
	protected Genome genome;
	protected ExperimentManager manager;
	protected ExptConfig econf;
	protected List<Point> points;
	protected List<StrandedPoint> spoints;
	protected List<Region> regions; 
	protected List<StrandedRegion> strandedReg;
	protected int window;
	
	public RegionCountSorter(GenomeConfig gcon, ExperimentManager man, ExptConfig econ, int win){
		gconf = gcon;
		genome = gconf.getGenome();
		manager= man;
		econf = econ;
		window = win;
	}
	
	public void setPoint(List<Point> p ){ points = p;}
	public void setStrandedPoint(List<StrandedPoint> p){spoints=p;}
	public void setStrandedRegion(List<StrandedRegion> sreg){ strandedReg=sreg; }
	public void setRegion(List<Region> reg){ regions=reg;}
	
	public void execute(){
		
		// get the regions to count reads
		List<Region> countRegions = new ArrayList<Region>();
		if (points != null){
			for (Point p : points)
				countRegions.add(p.expand(window/2));
		}else if (spoints !=null){
			for (StrandedPoint p : spoints)
				countRegions.add(p.expand(window/2));
		}else if (strandedReg != null){
			if(window >0){
				for (StrandedRegion sreg : strandedReg)					
					countRegions.add(sreg.resize(window));
			}else{
				countRegions.addAll(strandedReg);
			}
		}else if(regions != null){
			if (window>0){
				for (Region reg : regions)
					countRegions.add(reg.resize(window));
			}else{
				countRegions = regions;
			}
		}
		
		double sumReads = 0;
		// summing total tag counts for all experiments
		for (Sample sample : manager.getSamples())
			sumReads += sample.getHitCount();
		
		// get normalization constant for the total tag normalization
		double normCONST = sumReads/manager.getSamples().size();
		
		RegionCounts[] regionCounts = new RegionCounts[countRegions.size()];
		for (int i = 0 ; i < countRegions.size(); i++){
			float counts = 0;
			for (Sample sample : manager.getSamples())
				counts += (float) (normCONST*sample.countHits(countRegions.get(i))/sample.getHitCount());
			// add counts
			regionCounts[i] = new RegionCounts(counts);
			
			//add peaks or regions
			if (points !=null)
				regionCounts[i].setCoord(points.get(i));
			else if (spoints !=null)
				regionCounts[i].setStrandedCoord(spoints.get(i));
			else if (strandedReg != null)
				regionCounts[i].setStrandedRegion(strandedReg.get(i));
			else if (regions != null)
				regionCounts[i].setRegion(regions.get(i));		
		}	
		Arrays.sort(regionCounts);
		
		// outputting the list of regions in descending order of counts
		for (int i = 0; i < countRegions.size(); i++){
			if (regionCounts[i].getCoord()!=null)
				System.out.println(regionCounts[i].getCoord());
			else if (regionCounts[i].getStrandedCoord()!=null)
				System.out.println(regionCounts[i].getStrandedCoord());
			else if (regionCounts[i].getRegion()!=null)
				System.out.println(regionCounts[i].getRegion());
			else
				System.out.println(regionCounts[i].getStrandedRegion());
		}
		manager.close();
	}
	
	private class RegionCounts implements Comparable<RegionCounts>{	
		
		protected Point coord=null;  //Event coordinate
		protected StrandedPoint spoints=null;
		protected Region reg=null; //
		protected StrandedRegion sreg=null;
		protected float count;
		
		public RegionCounts(float count){
			this.count=count;
		}
		
		public void setCoord(Point coord){ this.coord=coord;}
		public void setStrandedCoord(StrandedPoint spoints){this.spoints=spoints;}
		public void setRegion(Region reg){this.reg=reg;}
		public void setStrandedRegion(StrandedRegion sreg){this.sreg=sreg;}
		
		public Point getCoord(){ return coord;}
		public StrandedPoint getStrandedCoord(){return spoints;}
		public Region getRegion(){return reg;}
		public StrandedRegion getStrandedRegion(){return sreg;}
		
		public int compareTo(RegionCounts regcounts){
			float compareVal = regcounts.count;
			if (compareVal > this.count){
				return 1;
			}else{
				return -1;
			}
		}
	}

	
	public static void main(String[] args){
		ArgParser ap = new ArgParser(args);		
        if((!ap.hasKey("peak") && !ap.hasKey("speak") && !ap.hasKey("region") && !ap.hasKey("sregion") )) { 
        	System.err.println("Please provide peak, region, sregion file !");
            System.err.println("Usage:\n" +
                               "RegionCountSorter\n" +
                               "--species <organism;genome> OR\n" +
                               "--geninfo <genome info file> AND --seq <path to seqs>\n" +
                               "--peak <list of peaks> OR --speak <list of stranded peaks> OR --region <list of regions> OR --sregion <list of stranded regions> \n" +
                               "--expt <experiments> \n" +
                               "\nOPTIONS:\n" +
                               "--win <window size around peaks> \n" +
                               "");
            System.exit(0);
        }
		
		GenomeConfig gconf = new GenomeConfig(args);
		ExptConfig econf = new ExptConfig(gconf.getGenome(), args);		
		econf.setPerBaseReadFiltering(false);		
		ExperimentManager manager = new ExperimentManager(econf);
		
		// parsing command line arguments	
		int win = Args.parseInteger(args, "win", -1);
		RegionCountSorter sorter = new RegionCountSorter(gconf,manager,econf,win);		
		List<Region> reg = new ArrayList<Region>() ;
		List<StrandedRegion> sreg ;
		if (ap.hasKey("peak")){
			List<Point> points = RegionFileUtilities.loadPointsFromFile(ap.getKeyValue("peak"), gconf.getGenome());
			sorter.setPoint(points);
		}else if (ap.hasKey("speak")){
			List<StrandedPoint> points = RegionFileUtilities.loadStrandedPointsFromFile(gconf.getGenome(), ap.getKeyValue("speak"));
			sorter.setStrandedPoint(points);
		}else if (ap.hasKey("sregion")){
			sreg = RegionFileUtilities.loadStrandedRegionsFromMotifFile(gconf.getGenome(), ap.getKeyValue("sregion"), -1);
			sorter.setStrandedRegion(sreg);
		}else if (ap.hasKey("region")){
			reg = RegionFileUtilities.loadRegionsFromFile(ap.getKeyValue("region"), gconf.getGenome(), -1);
			sorter.setRegion(reg);
		}			
		sorter.execute();
	}
}
