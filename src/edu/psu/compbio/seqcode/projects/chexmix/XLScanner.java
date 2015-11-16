package edu.psu.compbio.seqcode.projects.chexmix;

import java.io.File;
import java.util.List;

import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Point;

public class XLScanner {

	protected ExperimentManager manager;
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected ChExMixConfig cconfig;
	protected ProteinDNAInteractionModel model;
	protected WindowedTagDistributions tagsDistribs;
	protected List<Point> scanPoints;
	protected int scanWindow=50;
	
	
	public XLScanner(GenomeConfig gcon, ExptConfig econ, ChExMixConfig ccon){
		gconfig = gcon;
		econfig = econ;
		cconfig = ccon;
		cconfig.makeChExMixOutputDirs(true);
		manager = new ExperimentManager(econfig);
	}
	
	public void execute(){
		//Load appropriate data
		model = ProteinDNAInteractionModel.loadFromFile(cconfig, new File(cconfig.getModelFilename()));
		int winSize = model.getWidth();
		scanPoints = cconfig.getScanPoints();
		
		tagsDistribs = new WindowedTagDistributions(scanPoints, manager, winSize+scanWindow, true);
		
		
	}
	
	
	
	public static void main(String[] args){
		System.setProperty("java.awt.headless", "true");
		System.err.println("ChExMix version: "+ChExMixConfig.version);
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		//econ.setLoadRead2(false);//Enforce for chip-exo
		ChExMixConfig ccon = new ChExMixConfig(gcon, args);
		if(ccon.helpWanted()){
			System.err.println(ccon.getArgsList());
		}else if(ccon.getModelFilename()==null){
			System.err.println("Error: no ChExMix model provided. Use --model");
		}else{
			XLScanner xlScan = new XLScanner(gcon, econ, ccon);
			xlScan.execute();
		}
		
	}
}
