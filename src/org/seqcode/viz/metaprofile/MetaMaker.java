package org.seqcode.viz.metaprofile;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;

import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.viz.metaprofile.swing.MetaFrame;
import org.seqcode.viz.metaprofile.swing.MetaNonFrame;
import org.seqcode.viz.metaprofile.swing.MetaNonFrameMultiSet;


public class MetaMaker {
	
	private GenomeConfig gconfig;
	private MetaConfig mconfig;
	private ExptConfig econfig;
	private ExperimentManager manager;
	
	public MetaMaker(GenomeConfig g, MetaConfig m, ExptConfig e){
		gconfig = g;
		mconfig = m;
		econfig = e;
		if(mconfig.profilerType.equals("nucleosome"))
			econfig.setLoadPairs(true);
		manager = new ExperimentManager(econfig, true);
	}
	
	public void run(){
		try {
			if(mconfig.printHelp){
				System.err.println("MetaMaker:\n" +
						gconfig.getArgsList()+"\n"+
						mconfig.getArgsList()+"\n"+
						econfig.getArgsList()+"\n");
			}else{
				Genome gen = gconfig.getGenome();
				
				File file_mat = new File(mconfig.outName+"_matrix.peaks");
				if(!file_mat.exists()){
					file_mat.createNewFile();
				}
				FileWriter fw_mat = new FileWriter(file_mat.getAbsoluteFile());
				BufferedWriter br_mat = new BufferedWriter(fw_mat);
				
				BinningParameters params = new BinningParameters(mconfig.winLen, mconfig.bins);
				System.out.println("Binding Parameters:\tWindow size: "+params.getWindowSize()+"\tBins: "+params.getNumBins());
				System.out.println("LineMin: "+mconfig.lineMin+"\tLineMax"+mconfig.lineMax);
			
				PointProfiler profiler=null;
				boolean normalizeProfile=false;
				if(mconfig.profilerType.equals("simplechipseq")){
					profiler = new ChipSeqProfiler(gen, params, manager, mconfig.readExt, mconfig.strand);
				}else if(mconfig.profilerType.equals("fiveprime")){
					profiler = new Stranded5PrimeProfiler(gconfig, params, manager, mconfig.strand, mconfig.fivePrimeShift, mconfig.baseLimit, mconfig.baseLimitRelPosition, true, mconfig.cpmNorm);
				}else if(mconfig.profilerType.equals("nucleosome")){
					profiler = new PairedSeqMidpointProfiler(gen, params, manager);
				}
				
				if(mconfig.batchRun){
					System.out.println("Batch running...");
					System.setProperty("java.awt.headless", "true");
					if(mconfig.peakFiles.size()==1 || mconfig.peakFiles.size()==0){
						MetaNonFrame nonframe = new MetaNonFrame(gen, params, profiler, normalizeProfile, mconfig.saveSVG);
						nonframe.setColor(mconfig.color);
						nonframe.setDrawColorBar(mconfig.drawColorBar);
						nonframe.setTransparent(mconfig.transparent);
						nonframe.setDrawBorder(mconfig.drawBorder);
						MetaProfileHandler handler = nonframe.getHandler();
						if(mconfig.peakFiles.size()==1){
							System.out.println("Single set mode...");
							String peakFile = mconfig.peakFiles.get(0);
							Vector<Point> points = nonframe.getUtils().loadPoints(new File(peakFile));
							if(mconfig.printMatrix){
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
							Iterator<Point> points = nonframe.getUtils().loadTSSs("refGene");
							handler.addPoints(points);
						}
						while(handler.addingPoints()){}
						if(mconfig.cluster)
							nonframe.clusterLinePanel();
						//Set the panel sizes here...
						nonframe.setStyle(mconfig.profileStyle);
						nonframe.setLineMax(mconfig.lineMax);
						nonframe.setLineMin(mconfig.lineMin);
						nonframe.setLineThick(mconfig.lineThick);
						nonframe.saveImages(mconfig.outName);
						nonframe.savePointsToFile(mconfig.outName);
					}else if(mconfig.peakFiles.size()>1){
						System.out.println("Multiple set mode...");
						MetaNonFrameMultiSet multinonframe = new MetaNonFrameMultiSet(mconfig.peakFiles, gen, params, profiler, true);
						for(int x=0; x<mconfig.peakFiles.size(); x++){
							String pf = mconfig.peakFiles.get(x);
							Vector<Point> points = multinonframe.getUtils().loadPoints(new File(pf));
							List<MetaProfileHandler> handlers = multinonframe.getHandlers();
							handlers.get(x).addPoints(points);
							while(handlers.get(x).addingPoints()){}
						}
						multinonframe.saveImage(mconfig.outName);
						multinonframe.savePointsToFile(mconfig.outName);
					}
					System.out.println("Finished");
					if(profiler!=null)
						profiler.cleanup();
				}else{
					System.out.println("Initializing Meta-point frame...");
					MetaFrame frame = new MetaFrame(gen, params, profiler, normalizeProfile);
					frame.setColor(mconfig.color);
					frame.setLineMax(mconfig.lineMax);
					frame.setLineMin(mconfig.lineMin);
					frame.setLineThick(mconfig.lineThick);
					frame.startup();
					if(mconfig.peakFiles.size() > 0){
						MetaProfileHandler handler = frame.getHandler();
						for(String pf : mconfig.peakFiles){
							Vector<Point> points = frame.getUtils().loadPoints(new File(pf));
							handler.addPoints(points);
						}
					}
					frame.setLineMax(mconfig.lineMax);
					frame.setLineMin(mconfig.lineMin);
				}
				
				br_mat.close();
				fw_mat.close();
				manager.close();
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		GenomeConfig gconfig = new GenomeConfig(args);
		ExptConfig econfig = new ExptConfig(gconfig.getGenome(), args);
		MetaConfig mconfig = new MetaConfig(args);
		
		MetaMaker maker = new MetaMaker(gconfig, mconfig, econfig);
		
		maker.run();
	}
	
}
