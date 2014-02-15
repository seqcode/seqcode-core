package edu.psu.compbio.seqcode.projects.shaun;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedPoint;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.viz.metagenes.BinningParameters;
import edu.psu.compbio.seqcode.gse.viz.metagenes.EventProfiler;
import edu.psu.compbio.seqcode.gse.viz.metagenes.MetaNonFrame;
import edu.psu.compbio.seqcode.gse.viz.metagenes.MetaProfileHandler;
import edu.psu.compbio.seqcode.gse.viz.metagenes.MotifProfiler;
import edu.psu.compbio.seqcode.gse.viz.metagenes.PointProfiler;
import edu.psu.compbio.seqcode.gse.viz.metagenes.swing.MetaFrame;

public class EventMetaMaker {
	private static boolean batchRun = false;
	private static boolean cluster = false;
	private static Genome gen;
	
	public static void main(String[] args) {
		try {
			if(args.length < 2){ printError();}
			
			Pair<Organism, Genome> pair = Args.parseGenome(args);
			gen = pair.cdr();
			int winLen = Args.parseInteger(args,"win", 10000);
			int bins = Args.parseInteger(args,"bins", 100);
			String profilerType = Args.parseString(args, "profiler", "events");	
			String eventFile = Args.parseString(args,"events", null);
			String peakFile = Args.parseString(args, "peaks", null);
			String outName = Args.parseString(args, "out", "meta");
			double lineMax = Args.parseDouble(args,"linemax", 100);
			if(Args.parseFlags(args).contains("batch")){batchRun=true;}
			if(Args.parseFlags(args).contains("cluster")){cluster=true;}
			Color c = Color.blue;
			String newCol = Args.parseString(args, "color", "blue");
			if(newCol.equals("red"))
				c=Color.red;
			if(newCol.equals("orange"))
				c=Color.orange;
			if(newCol.equals("green"))
				c=new Color(0,153,0);
			if(newCol.equals("black"))
				c=Color.black;
		
			
			if(gen==null || eventFile==null){printError();}
	
			BinningParameters params = new BinningParameters(winLen, bins);
			System.out.println("Binding Parameters:\tWindow size: "+params.getWindowSize()+"\tBins: "+params.getNumBins());
		
			PointProfiler profiler=null;
			boolean normalizeProfile=false;
			if(profilerType.equals("events")){
				ArrayList<Point> events = loadPoints(new File(eventFile), gen);
				System.out.println("Loading data...");
				profiler = new EventProfiler(params, gen, events);
			}
			
			if(batchRun){
				System.out.println("Batch running...");
				MetaNonFrame nonframe = new MetaNonFrame(gen, params, profiler, normalizeProfile, false);
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
				nonframe.setLineMax(lineMax);
				nonframe.saveImages(outName);
				nonframe.savePointsToFile(outName);
				System.out.println("Finished");
			}else{
				System.out.println("Initializing Meta-point frame...");
				MetaFrame frame = new MetaFrame(gen, params, profiler, normalizeProfile);
				frame.setColor(c);
				frame.startup();
				MetaProfileHandler handler = frame.getHandler();
				if(peakFile != null){
					Vector<Point> points = frame.getUtils().loadPoints(new File(peakFile));
					handler.addPoints(points);
				}frame.setLineMax(lineMax);
			}
		} catch (NotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private static void printError(){
		System.err.println("Usage: EventMetaMaker --species <organism;genome> \n" +
				"--win <profile width> --bins <num bins> \n" +
				"--profiler <events> \n" +
				"--events <point file>\n" +
				"--linemax <max>\n" +
				"--peaks <peaks file name> --out <output root name> \n" +
				"--color <red/green/blue/black> \n" +
				"--cluster [flag to cluster in batch mode] \n" +
				"--batch [a flag to run without displaying the window]");
		System.exit(1);
	}
	public static ArrayList<Point> loadPoints(File f, Genome g) throws IOException {
		System.out.println("Loading points");
		ArrayList<Point> pts = new ArrayList<Point>();
		BufferedReader br = new BufferedReader(new FileReader(f));
		Pattern ptpatt = Pattern.compile("([^:\\s]+):(\\d+)");
		Pattern strptpatt = Pattern.compile("([^:\\s]+):(\\d+):([^:\\s]+)");
		String line;
		while((line = br.readLine()) != null) {
			Matcher m = ptpatt.matcher(line);
			if(m.find()) { 
				String chrom = m.group(1);
				int location = Integer.parseInt(m.group(2));
				char strand = '?';
				Matcher sm = strptpatt.matcher(line);
				if(sm.find()){
					String strandstr = sm.group(3);
					if(strandstr.length() > 0) { strand = strandstr.charAt(0); }
				}
				Point pt = null;
				if(strand == '+') { 
					pt = new StrandedPoint(g, chrom, location, strand);
				} else if (strand == '-') { 
					pt = new StrandedPoint(g, chrom, location, strand);					
				} else { 
					pt = new Point(g, chrom, location);
				}
				pts.add(pt);
			} else { 
				System.err.println(String.format("Couldn't find point in line \"%s\"", line));
			}
		}
		br.close();
		System.err.println(pts.size()+" points loaded");
		return pts;
	}
}
