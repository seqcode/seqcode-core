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
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.StrandedPoint;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.PointParser;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.RegionParser;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;
import edu.psu.compbio.seqcode.gse.utils.strings.StringUtils;

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
	
	
	
	
	public PeaksVsGenes(Genome g){
		gen =g;
	}
	
	public static void main(String[] args){
		GenomeConfig gconfig = new GenomeConfig(args);
		ArgParser ap = new ArgParser(args);
		int radius = Args.parseInteger(args, "radius", 50000);
		String gene_attribute = Args.parseString(args, "gAttribute", "tstat");
		boolean cuffdiff = ap.hasKey("cuffdiff");
		
		boolean printNearestGene = ap.hasKey("printNearestGene");
		boolean printCloseHighAttributeGenes  = ap.hasKey("printCloseHighAttributeGenes");
		double threshold = 0.0;
		if(printCloseHighAttributeGenes){
			threshold = Args.parseDouble(args, "threshold", 2);
		}
		boolean printAllCloseGenes = ap.hasKey("printAllCloseGenes");
		
		String peaksfile = Args.parseString(args, "peaks", null);
		if(peaksfile == null){
			System.err.println("Provide ChIP-Seq peaks file!!");
			return;
		}
		String genefile = Args.parseString(args, "genes", null);
		if(genefile == null){
			System.err.println("Provide genes file!!");
			return;
		}
		
		// For now formDB is always false .. with add that utility later
		boolean fromDB=false;
		PeaksVsGenes analyzer = new PeaksVsGenes(gconfig.getGenome());
		
		analyzer.laodpeaks(peaksfile, fromDB);
		analyzer.laodgenes(genefile, cuffdiff, fromDB, gene_attribute);
		analyzer.setRadius(radius);
		
		if(printNearestGene){
			analyzer.printNearestGene();
		}else if(printCloseHighAttributeGenes){
			analyzer.printCloseHighAttributeGenes(threshold);
		}else if(printAllCloseGenes){
			analyzer.printAllCloseGenes();
		}
		
		
	}
	
	
	
	
	public void printNearestGene(){
		HashMap<String, List<StrandedPoint>> genesbyChrs = hashbychrom(genes);
		for(Point p: peaks){
			int MINDISTANCE = Integer.MAX_VALUE;
			StrandedPoint nearestGene= null;
			boolean hasgene = false;
			if(genesbyChrs.containsKey(p.getChrom())){
				for(StrandedPoint gene : genesbyChrs.get(p.getChrom())){
					int distance = gene.distance(p);
					if(distance < MINDISTANCE && distance < radius){
						MINDISTANCE = distance;
						nearestGene = gene;
						hasgene = true;
					}
				}
			}
			
			if(hasgene){
				String gene_names = StringUtils.join(startToGenename.get(nearestGene.getLocationString()), "\t"); 
				System.out.println(p.getLocationString()+"\t"+gene_names+"\t"+Integer.toString(MINDISTANCE));
			}else{
				System.out.println(p.getLocationString()+"\tNULL\tNULL");
			}
			
		}
		
	}
	
	public void printCloseHighAttributeGenes(double threshold){
		HashMap<String, List<StrandedPoint>> genesbyChrs = hashbychrom(genes);
		for(Point p: peaks){
			List<StrandedPoint> highattgenes = new ArrayList<StrandedPoint>();
			List<Integer> distances = new ArrayList<Integer>();
			boolean hasgene = false;
			if(genesbyChrs.containsKey(p.getChrom())){
				for(StrandedPoint sp : genesbyChrs.get(p.getChrom())){
					int distance = sp.distance(p);
					if(distance < radius  && Math.abs(gene_attributes.get(sp.getLocationString())) > threshold){
						highattgenes.add(sp);
						hasgene = true;
						distances.add(distance);
					}
				}
			}
			if(hasgene){
				List<String> gene_names = new ArrayList<String>();
				for(StrandedPoint sp : highattgenes){
					gene_names.addAll(startToGenename.get(sp.getLocationString()));
				}
				String gene_names_string = StringUtils.join(gene_names, "\t");
				String distance_string = StringUtils.join(distances, "\t");
				System.out.println(p.getLocationString()+"\t"+gene_names_string+"\t"+distance_string);
			}else{
				System.out.println(p.getLocationString()+"NULL\tNULL");
			}
		}
	}
	
	
	public void printAllCloseGenes(){
		this.printCloseHighAttributeGenes(0);
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////
	///Slave methods
	//////////////////////////////////////////////////////
	
	public HashMap<String, List<StrandedPoint>> hashbychrom(List<StrandedPoint> pts){
		HashMap<String, List<StrandedPoint>> byChr = new HashMap<String, List<StrandedPoint>>();
		for(StrandedPoint p : pts){
			if(!byChr.containsKey(p.getChrom()))
				byChr.put(p.getChrom(), new ArrayList<StrandedPoint>());
			byChr.get(p.getChrom()).add(p);
		}
		return byChr;
	}
	
	
	
	
	
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Loaders and Settors
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	public void setRadius(int r){ radius = r;}
	
	public void laodpeaks(String peaksfile, boolean fromDB){
		peaks = RegionFileUtilities.loadPeaksFromPeakFile(gen, peaksfile, (Integer) null);
		if(!fromDB){
			try{
				File pFile = new File(peaksfile);
				if(!pFile.isFile()){System.err.println("Invalid peaks filename");System.exit(1);}
				peak_attributes = new HashMap<String, Double>();
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
		
	}
	
	public void laodgenes(String genefile, boolean cuffdiff, boolean fromDB, String attribute){
		if(!cuffdiff){
			genes = RegionFileUtilities.loadStrandedPointFromRefTssFile(gen, genefile);
			if(!fromDB){
				try{
					gene_attributes = new HashMap<String, Double>();
					File gFile = new File(genefile);
					if(!gFile.isFile()){System.err.println("Invalid genes list filename");System.exit(1);}
					BufferedReader reader = new BufferedReader(new FileReader(gFile));
					String line;
					while ((line = reader.readLine()) != null) {
						line = line.trim();
						String[] words = line.split("\t");
						if(words.length >=1 && !words[0].contains("#") &&  !words[0].equals("Region") && !words[0].equals("Position")){
							String[] subwords = words[0].split(":");
							PointParser pparser = new PointParser(gen);
							Point p = pparser.execute(subwords[0]+":"+subwords[1]);
							StrandedPoint sp = new StrandedPoint(p,subwords[2].charAt(0));
							gene_attributes.put(sp.getLocationString(), Double.parseDouble(words[1]));
						}
							
					}
					reader.close();
				}catch(Exception e){
					e.printStackTrace();
				}
			}
		}else{
			try{
				genes = new ArrayList<StrandedPoint>();
				gene_attributes = new HashMap<String, Double>();
				startToGenename = new HashMap<String, List<String>>();
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
			            	if(!fromDB){
			            		if(attribute == null){
			            			throw new Exception("Provide which gene attribute to load from cuffdiff file!!");
			            		}
			            		if(attribute.toLowerCase().contains("foldchange")){
			            			gene_attributes.put(genes.get(genes.size()-1).getLocationString(), Double.parseDouble(words[9]));
			            		}else if(attribute.toLowerCase().contains("tstat")){
			            			gene_attributes.put(genes.get(genes.size()-1).getLocationString(), Double.parseDouble(words[10]));
			            		}else if(attribute.toLowerCase().contains("pvalue")){
			            			gene_attributes.put(genes.get(genes.size()-1).getLocationString(), Double.parseDouble(words[11]));
			            		}else if(attribute.toLowerCase().contains("qvalue")){
			            			gene_attributes.put(genes.get(genes.size()-1).getLocationString(), Double.parseDouble(words[12]));
			            		}
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
	

}
