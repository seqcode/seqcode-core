package edu.psu.compbio.seqcode.gse.viz.metagenes;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.FlowLayout;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.locators.ChipChipLocator;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.deepseq.DeepSeqExpt;
import edu.psu.compbio.seqcode.gse.ewok.verbs.chipseq.SeqExpander;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.viz.metagenes.swing.MetaFrame;
import edu.psu.compbio.seqcode.projects.multigps.experiments.Sample;

public class MetaMaker {
	private static boolean batchRun = false;
	private static boolean cluster = false;
	
	public static void main(String[] args) {
		try {
			if(args.length < 2){ printError();}
			
			Pair<Organism, Genome> pair = Args.parseGenome(args);
			Genome gen = pair.cdr();
			int winLen = Args.parseInteger(args,"win", 10000);
			int bins = Args.parseInteger(args,"bins", 100);
			int readExt = Args.parseInteger(args,"readext", 0);
			double lineMin = Args.parseDouble(args,"linemin", 0);
			double lineMax = Args.parseDouble(args,"linemax", 100);
			int lineThick = Args.parseInteger(args,"linethick", 1);
			double pbMax = Args.parseDouble(args,"pbmax", 100);
			char strand = Args.parseString(args, "strand", "/").charAt(0);
			boolean drawColorBar = !Args.parseFlags(args).contains("nocolorbar");
			boolean saveSVG = Args.parseFlags(args).contains("svg");
			boolean transparent = Args.parseFlags(args).contains("transparent");
			boolean printMatrix = Args.parseFlags(args).contains("printMatrix");
			String profilerType = Args.parseString(args, "profiler", "simplechipseq");	
			String profileStyle = Args.parseString(args, "style", "Line");	
			List<String> expts = (List<String>) Args.parseStrings(args,"expt");
			Collection<String> exptFilenames = Args.parseStrings(args, "exptfile");
			String format = Args.parseString(args, "format","SAM");
			List<String> backs = (List<String>) Args.parseStrings(args,"back");
			List<String> peakFiles = (List<String>)Args.parseStrings(args, "peaks");
			String outName = Args.parseString(args, "out", "meta");
			if(Args.parseFlags(args).contains("batch")){batchRun=true;}
			if(Args.parseFlags(args).contains("cluster")){cluster=true;}
			Color c = Color.blue;
			String newCol = Args.parseString(args, "color", "blue");
			File file_mat = new File(outName+"_matrix.peaks");
			if(!file_mat.exists()){
				file_mat.createNewFile();
			}
			FileWriter fw_mat = new FileWriter(file_mat.getAbsoluteFile());
			BufferedWriter br_mat = new BufferedWriter(fw_mat);
			
			if(newCol.equals("red"))
				c=Color.red;
			if(newCol.equals("green"))
				c=Color.green;
			for(int s=0; s<args.length; s++){
				if(args[s].equals("--color4")){
					Integer R = new Integer(args[s+1]);
					Integer G = new Integer(args[s+2]);
					Integer B = new Integer(args[s+3]);
					Integer A = new Integer(args[s+4]);
					c= new Color(R,G,B,A);
				}
			}
		
			
			if(gen==null || (expts.size()==0 && exptFilenames.size()==0)){printError();}
	
			BinningParameters params = new BinningParameters(winLen, bins);
			System.out.println("Binding Parameters:\tWindow size: "+params.getWindowSize()+"\tBins: "+params.getNumBins());
		
			PointProfiler profiler=null;
			boolean normalizeProfile=false;
			if(profilerType.equals("simplechipseq")){
				List<SeqLocator> exptlocs = Args.parseSeqExpt(args,"expt");
				if(exptlocs.size()>0){
					ArrayList<SeqExpander> exptexps = new ArrayList<SeqExpander>();
					for(SeqLocator loc : exptlocs){
						System.out.println(loc.getExptName()+"\t"+loc.getAlignName());
						exptexps.add(new SeqExpander(loc));
					}
					System.out.println("Loading data...");
					profiler = new ChipSeqProfiler(params, exptexps, readExt, pbMax,strand);
				}else{
					List<File> exptFiles = new ArrayList<File>();
					for(String s : exptFilenames)
						exptFiles.add(new File(s));
					DeepSeqExpt dse = new DeepSeqExpt(gen, exptFiles, false, format, 40);
					profiler = new ChipSeqProfiler(params, dse, readExt, pbMax,strand);
				}
			}else if(profilerType.equals("fiveprime")){
				List<SeqLocator> exptlocs = Args.parseSeqExpt(args,"expt");
				if(exptlocs.size()>0){
					ArrayList<SeqExpander> exptexps = new ArrayList<SeqExpander>();
					for(SeqLocator loc : exptlocs)
						exptexps.add(new SeqExpander(loc));
					System.out.println("Loading data...");
					profiler = new Stranded5PrimeProfiler(params, exptexps, strand, pbMax);
				}else{
					List<File> exptFiles = new ArrayList<File>();
					for(String s : exptFilenames)
						exptFiles.add(new File(s));
					DeepSeqExpt dse = new DeepSeqExpt(gen, exptFiles, false, format, 40);
					profiler = new Stranded5PrimeProfiler(params, dse, strand, pbMax);
				}
			}
			
			if(batchRun){
				System.out.println("Batch running...");
				
				if(peakFiles.size()==1 || peakFiles.size()==0){
					MetaNonFrame nonframe = new MetaNonFrame(gen, params, profiler, normalizeProfile, saveSVG);
					nonframe.setColor(c);
					nonframe.setDrawColorBar(drawColorBar);
					nonframe.setTransparent(transparent);
					MetaProfileHandler handler = nonframe.getHandler();
					if(peakFiles.size()==1){
						System.out.println("Single set mode...");
						String peakFile = peakFiles.get(0);
						Vector<Point> points = nonframe.getUtils().loadPoints(new File(peakFile));
						if(printMatrix){
							double[][] mat_out = null;
							for(int k=0; k<points.size(); k++){
								if(k==0){
									PointProfile temp = (PointProfile) profiler.execute(points.get(k));
									mat_out = new double[points.size()][temp.length()];
									for(int j=0; j< temp.length(); j++){
										mat_out[k][j] = temp.value(j);
									}
								}
								else{
									PointProfile temp = (PointProfile) profiler.execute(points.get(k));
									for(int j=0; j< temp.length(); j++){
										mat_out[k][j] = temp.value(j);
									}
								}
							}
							for(int k =0; k< mat_out.length; k++ ){
								br_mat.write(points.get(k).getLocationString()+"\t");
								for (int j=0; j< mat_out[k].length; j++){
									br_mat.write(mat_out[k][j]+"\t");
								}
								br_mat.write("\n");
							}
						}
						handler.addPoints(points);
					}else{
						System.out.println("All TSS mode...");
						Iterator<Point> points = nonframe.getUtils().loadTSSs();
						handler.addPoints(points);
					}
					while(handler.addingPoints()){}
					if(cluster)
						nonframe.clusterLinePanel();
					//Set the panel sizes here...
					nonframe.setStyle(profileStyle);
					nonframe.setLineMin(lineMin);
					nonframe.setLineMax(lineMax);
					nonframe.setLineThick(lineThick);
					nonframe.saveImages(outName);
					nonframe.savePointsToFile(outName);
				}else if(peakFiles.size()>1){
					System.out.println("Multiple set mode...");
					MetaNonFrameMultiSet multinonframe = new MetaNonFrameMultiSet(peakFiles, gen, params, profiler, true);
					for(int x=0; x<peakFiles.size(); x++){
						String pf = peakFiles.get(x);
						Vector<Point> points = multinonframe.getUtils().loadPoints(new File(pf));
						List<MetaProfileHandler> handlers = multinonframe.getHandlers();
						handlers.get(x).addPoints(points);
						while(handlers.get(x).addingPoints()){}
					}
					multinonframe.saveImage(outName);
					multinonframe.savePointsToFile(outName);
				}
				System.out.println("Finished");
				if(profiler!=null)
					profiler.cleanup();
			}else{
				System.out.println("Initializing Meta-point frame...");
				MetaFrame frame = new MetaFrame(gen, params, profiler, normalizeProfile);
				frame.setColor(c);
				frame.setLineMax(lineMax);
				frame.setLineMin(lineMin);
				frame.setLineThick(lineThick);
				frame.startup();
				if(peakFiles.size() > 0){
					MetaProfileHandler handler = frame.getHandler();
					for(String pf : peakFiles){
						Vector<Point> points = frame.getUtils().loadPoints(new File(pf));
						handler.addPoints(points);
					}
				}
				frame.setLineMax(lineMax);
				frame.setLineMin(lineMin);
			}
		} catch (SQLException e) {
			e.printStackTrace();
		} catch (NotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	private static void printError(){
		System.err.println("Usage: MetaMaker --species <organism;genome> \n" +
				"--win <profile width> --bins <num bins> \n" +
				"--readext <read extension> \n" +
				"--linemin <min>  --linemax <max> \n" +
				"--pbmax <per base max>\n" +
				"--profiler <chipseq/fiveprime> \n" +
				"--expt <experiment names> OR --exptfile <file names> AND --format <SAM> \n" +
				"--back <control experiment names (only applies to chipseq)> \n" +
				"--peaks <peaks file name> --out <output root name> \n" +
				"--color <red/green/blue> or --color4 <R G B A>\n" +
				"--strand <+-/>\n" +
				"--printMatrix [flag to print the matrix of tags] \n"+
				"--cluster [flag to cluster in batch mode] \n" +
				"--batch [a flag to run without displaying the window]\n" +
				"--nocolorbar [flag to turn off colorbar in batch mode]\n" +
				"--transparent [flag for transparent background]\n" +
				"--style <Line/Histo>\n" +
				"--svg [flag to save an SVG image]\n");
		System.exit(1);
	}
}
