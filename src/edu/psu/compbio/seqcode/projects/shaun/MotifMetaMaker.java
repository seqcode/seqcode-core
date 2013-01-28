package edu.psu.compbio.seqcode.projects.shaun;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Vector;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.motifs.MarkovBackgroundModel;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.io.motifs.BackgroundModelIO;
import edu.psu.compbio.seqcode.gse.viz.metagenes.BinningParameters;
import edu.psu.compbio.seqcode.gse.viz.metagenes.MetaNonFrame;
import edu.psu.compbio.seqcode.gse.viz.metagenes.MetaProfileHandler;
import edu.psu.compbio.seqcode.gse.viz.metagenes.MotifProfiler;
import edu.psu.compbio.seqcode.gse.viz.metagenes.PointProfiler;
import edu.psu.compbio.seqcode.gse.viz.metagenes.swing.MetaFrame;

public class MotifMetaMaker {
	private static boolean batchRun = false;
	private static boolean cluster = false;
	private static boolean usingColorQuanta=false;
	private static double[] colorQuanta=null;
	
	public static void main(String[] args) throws IOException, ParseException {
		try {
			if(args.length < 2){ printError();}
			
			Pair<Organism, Genome> pair = Args.parseGenome(args);
			Genome gen = pair.cdr();
			int winLen = Args.parseInteger(args,"win", 10000);
			int bins = Args.parseInteger(args,"bins", 100);
			String profilerType = Args.parseString(args, "profiler", "motif");	
			String motifName = Args.parseString(args,"motif", null);
			String backName = Args.parseString(args,"mback", null);
			double minthres = Args.parseDouble(args, "mthres", 0);
			String peakFile = Args.parseString(args, "peaks", null);
			String outName = Args.parseString(args, "out", "meta");
			if(Args.parseFlags(args).contains("batch")){batchRun=true;}
			if(Args.parseFlags(args).contains("cluster")){cluster=true;}
			if(Args.parseArgs(args).contains("quanta")){
				usingColorQuanta=true; 
				for(int a=0; a<args.length; a++){
					if(args[a].equals("--quanta")){
						int numQ=new Integer(args[a+1]);
						colorQuanta=new double[numQ];
						for(int q=0; q<numQ && (a+2+q)<args.length; q++){
							colorQuanta[q]=new Double(args[a+2+q]);
						}
					}
				}
				minthres = colorQuanta[0];
			}
			Color c = Color.blue;
			String newCol = Args.parseString(args, "color", "blue");
			if(newCol.equals("red"))
				c=Color.red;
			if(newCol.equals("green"))
				c=new Color(0,153,0);
			if(newCol.equals("black"))
				c=Color.black;
		
			
			if(gen==null || motifName==null){printError();}
	
			BinningParameters params = new BinningParameters(winLen, bins);
			System.out.println("Binding Parameters:\tWindow size: "+params.getWindowSize()+"\tBins: "+params.getNumBins());
		
			PointProfiler profiler=null;
			boolean normalizeProfile=false;
			if(profilerType.equals("motif")){
				ArrayList<WeightMatrix> motifs = new ArrayList<WeightMatrix>();
				//Load the background
				MarkovBackgroundModel back = BackgroundModelIO.parseMarkovBackgroundModel(backName, gen);
		    	//Load the motifs
		    	FreqMatrixImport motifImport = new FreqMatrixImport();
				motifImport.setBackground(back);
				for(WeightMatrix wm : motifImport.readTransfacMatrices(motifName)){
					motifs.add(wm);
				}
				System.out.println("Loading data...");
				profiler = new MotifProfiler(params, gen, motifs.get(0), minthres);
			}
			
			if(batchRun){
				System.out.println("Batch running...");
				MetaNonFrame nonframe = new MetaNonFrame(gen, params, profiler, normalizeProfile);
				if(usingColorQuanta)
					nonframe.setLinePanelColorQuanta(colorQuanta);
				nonframe.setColor(c);
				MetaProfileHandler handler = nonframe.getHandler();
				if(peakFile != null){
					Vector<Point> points = nonframe.getUtils().loadPoints(new File(peakFile));
					handler.addPoints(points);
				}else{
					Iterator<Point> points = nonframe.getUtils().loadTSSs();
					handler.addPoints(points);
				}
				while(handler.addingPoints()){}
				if(cluster)
					nonframe.clusterLinePanel();
				//Set the panel sizes here...
				nonframe.setLineMin(0);
				nonframe.saveImages(outName);
				nonframe.savePointsToFile(outName);
				System.out.println("Finished");
			}else{
				System.out.println("Initializing Meta-point frame...");
				MetaFrame frame = new MetaFrame(gen, params, profiler, normalizeProfile);
				if(usingColorQuanta)
					frame.setLinePanelColorQuanta(colorQuanta);
				frame.setColor(c);
				frame.startup();
				MetaProfileHandler handler = frame.getHandler();
				if(peakFile != null){
					Vector<Point> points = frame.getUtils().loadPoints(new File(peakFile));
					handler.addPoints(points);
				}
			}
		} catch (NotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private static void printError(){
		System.err.println("Usage: MotifMetaMaker --species <organism;genome> \n" +
				"--win <profile width> --bins <num bins> \n" +
				"--profiler <motif> \n" +
				"--motif <motif names> --mback <background model name> --mthres <threshold>\n" +
				"--peaks <peaks file name> --out <output root name> \n" +
				"--color <red/green/blue> \n" +
				"--cluster [flag to cluster in batch mode] \n" +
				"--batch [a flag to run without displaying the window]");
		System.exit(1);
	}
}
