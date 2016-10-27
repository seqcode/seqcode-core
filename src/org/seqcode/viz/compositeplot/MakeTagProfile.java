package org.seqcode.viz.compositeplot;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.List;

import org.seqcode.data.io.RegionFileUtilities;
import org.seqcode.deepseq.composite.CompositeTagDistribution;
import org.seqcode.deepseq.experiments.ExperimentCondition;
import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.gseutils.Args;

public class MakeTagProfile {

	
	protected ExperimentManager manager;
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected int winSize;
	protected List<StrandedPoint> points;
	protected CompositeTagDistribution tagDist;
	protected String outName;
	
	public MakeTagProfile(GenomeConfig gcon, ExptConfig econ, List<StrandedPoint> pts, int win, String out){
		gconfig = gcon;
		econfig = econ;
		manager = new ExperimentManager(econ);
		winSize = win;
		points = pts;
		outName=out;
		
		tagDist = new CompositeTagDistribution(points, manager, winSize, true);
		
	}
	
	public void execute(){
		try{	
			double Ymax=0;
			//full composites
			for(ExperimentCondition cond : manager.getConditions()){
				TagProfile full = new TagProfile(tagDist.getCompositeWatson(cond), 
						tagDist.getCompositeCrick(cond), 
						tagDist.getCenterOffset());
				TagProfilePaintable fullPainter = new TagProfilePaintable(full);
				fullPainter.autoYmax(true);
				fullPainter.setFilledColumns(true);
				fullPainter.setWatsonColor(new Color(122,88,143));
				fullPainter.setCrickColor(new Color(154,114,179));
				Ymax = Math.max(Ymax,  fullPainter.getYmax());
				String fullImageFileName = outName+"_compositeFull."+cond.getName()+".png";
				fullPainter.saveImage(new File(fullImageFileName), 1060, 500, true);
				fullPainter.setProfileLeftLimit(-50);
				fullPainter.setProfileRightLimit(+49);
				String fullZoomImageFileName = outName+"_compositeFullZoom."+cond.getName()+".png";
				fullPainter.saveImage(new File(fullZoomImageFileName), 1060, 500, true);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void close(){
		manager.close();
	}
	
	//Main method to make new composite plots
	public static void main(String[] args){
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		if(args.length==0){
			System.err.println("MakeTagProfile:\n"+
					"\t--cpoints <stranded point file>\n"+
					"\t--cwin <window around points>\n"+
					"\t--out <outfile base name>\n"+
					"Genome:\n" +
					"\t--species <Species;Genome>\n" +
					"\tOR\n" +
					"\t--geninfo <genome info file> AND --seq <fasta seq directory>\n" +
					"Experiment Design File:\n" +
					"\t--design <file name>\n");			
		}else{
			int w = Args.parseInteger(args, "cwin", 400);
			String pFile = Args.parseString(args, "cpoints", null);
			String out = Args.parseString(args, "out", "out");
			List<StrandedPoint> pts = RegionFileUtilities.loadStrandedPointsFromFile(gcon.getGenome(), pFile);
			
			MakeTagProfile maker = new MakeTagProfile(gcon, econ, pts, w, out);
			
			maker.execute();
			
			maker.close();
		}
	}
}
