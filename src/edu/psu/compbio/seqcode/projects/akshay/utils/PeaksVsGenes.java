package edu.psu.compbio.seqcode.projects.akshay.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.StrandedPoint;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.PointParser;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.RegionParser;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;

public class PeaksVsGenes {
	
	private Genome gen;
	private List<Point> peaks;
	private List<StrandedPoint> genes;
	private HashMap <String, List<String>> startToGenename;
	
	private ExperimentManager manager = null;
	private boolean fromDB = false;
	private ExptConfig econfig = null;
	
	
	private int radius = 50000; // Default is 50kb
	
	private HashMap<String,Double> peak_attributes;
	private HashMap<String, Double> gene_attributes;
	
	
	
	
	
	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Loaders and Settors
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	public void laodpeaks(String peaksfile){
		peaks = RegionFileUtilities.loadPeaksFromPeakFile(gen, peaksfile, (Integer) null);
	}
	
	public void laodgenes(String genefile, boolean cuffdiff){
		if(!cuffdiff){
			genes = RegionFileUtilities.loadStrandedPointFromRefTssFile(gen, genefile);
		}else{
			try{
				File gFile = new File(genefile);
				if(!gFile.isFile()){System.err.println("Invalid cuffdiff filename");System.exit(1);}
				BufferedReader reader = new BufferedReader(new FileReader(gFile));
				String line;
				while ((line = reader.readLine()) != null) {
					line = line.trim();
					String[] words = line.split("\t");
					if(!words[0].contains("test_id")){
						if(words[3].contains("-")){
							RegionParser rparser = new RegionParser(gen);
			            	Region q = rparser.execute(words[3]);
			            	genes.add(new StrandedPoint(q.getMidpoint(),'+'));
			            	if(startToGenename.containsKey(genes.get(genes.size()-1).getLocationString())){ 
			            		startToGenename.get(genes.get(genes.size()-1).getLocationString()).add(words[2]);
			            	}else{
			            		List<String> names = new ArrayList<String>();
			            		names.add(words[2]);
			            		startToGenename.put(genes.get(genes.size()-1).getLocationString(), names);
			            	}
						}else{
							System.err.println("Invalid cuffdiff file format");System.exit(1);
						}
						
					}
				}
				reader.close();
			}catch(Exception e){
				e.printStackTrace();
			}
		}
	}
	
	
	/**
	 * 
	 * @param peaksfile
	 * @param readDB
	 */
	public void loadpeakattributes(String peaksfile, boolean readDB){ // Add readdb loading utility later
		try{
			File pFile = new File(peaksfile);
			if(!pFile.isFile()){System.err.println("Invalid peaks filename");System.exit(1);}
			BufferedReader reader = new BufferedReader(new FileReader(pFile));
			String line;
			while ((line = reader.readLine()) != null) {
				line = line.trim();
				String[] words = line.split("\t");
				if(words.length >=1 && !words[0].contains("#") &&  !words[0].equals("Region") && !words[0].equals("Position")){
					if(words[0].contains("-")){
	                	RegionParser rparser = new RegionParser(gen);
		            	Region q = rparser.execute(words[0]);
		            	peak_attributes.put(q.getMidpoint().getLocationString(), Double.parseDouble(words[1]));		            	
	            	}else{
	            		PointParser pparser = new PointParser(gen);
	            		Point p = pparser.execute(words[0]);
	            		peak_attributes.put(p.getLocationString(), Double.parseDouble(words[1]));
	            	}
				}
					
			}
			reader.close();
		}catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public void loadgeneattributes(){
		
	}
	
	

}
