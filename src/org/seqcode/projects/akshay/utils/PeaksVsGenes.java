package org.seqcode.projects.akshay.utils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.seqcode.deepseq.experiments.ExperimentManager;
import org.seqcode.deepseq.experiments.ExptConfig;
import org.seqcode.genome.Genome;
import org.seqcode.genome.GenomeConfig;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.gse.gsebricks.verbs.location.PointParser;
import org.seqcode.gse.gsebricks.verbs.location.RegionParser;
import org.seqcode.gse.tools.utils.Args;
import org.seqcode.gse.utils.ArgParser;
import org.seqcode.gse.utils.io.RegionFileUtilities;
import org.seqcode.gse.utils.strings.StringUtils;
import org.seqcode.projects.shaun.rnaseq.GTFReader;
import org.seqcode.projects.shaun.rnaseq.genemodels.GeneTUnit;
import org.seqcode.projects.shaun.rnaseq.genemodels.SplicedTUnit;


public class PeaksVsGenes {
	
	private Genome gen;
	private List<Point> peaks;
	private List<StrandedPoint> genes;
	private HashMap <String, List<String>> startToGenename;
	
	private ExperimentManager manager = null;
	private boolean fromDB = false;
	private ExptConfig econfig = null;
	
	private boolean gtf = false;
	
	private int radius = 50000; // Default is 50kb
	
	private HashMap<String,Double> peak_attributes;
	private HashMap<String, Double> gene_attributes;
	
	
	
	
	public PeaksVsGenes(Genome g, boolean gtf_loading, boolean DB_loading){
		gen =g;
		gtf = gtf_loading;
		fromDB = DB_loading;
		
	}
	
	public static void main(String[] args){
		GenomeConfig gconfig = new GenomeConfig(args);
		ArgParser ap = new ArgParser(args);
		int radius = Args.parseInteger(args, "radius", 50000);
		
		
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
		
		boolean cuffdiff = ap.hasKey("cuffdiff");
		boolean gtf = !cuffdiff;
		String gene_attribute = "";
		String genefile = "";
		if(cuffdiff){
			gene_attribute = Args.parseString(args, "gAttribute", "tstat");
			genefile = Args.parseString(args, "genes", null);
			if(genefile == null){
				System.err.println("Provide genes file!!");
				return;
			}
		}else{
			genefile = Args.parseString(args, "gtf", null);
			if(genefile == null){
				System.err.println("Provide gtf file!!");
				return;
			}
			// When loading from a gtf file you can provide gene attributes separately (StrandedPoint	attribute-value)
			if(ap.hasKey("gAttributeFile")){
				gene_attribute = ap.getKeyValue("gAttributeFile");
			}
		}

		// For now formDB is always false .. will add that utility later
		boolean fromDB=false;
		PeaksVsGenes analyzer = new PeaksVsGenes(gconfig.getGenome(), gtf,fromDB);
		
		analyzer.laodpeaks(peaksfile);
		analyzer.laodgenes(genefile, gene_attribute);
		analyzer.setRadius(radius);
		
		if(printNearestGene){
			analyzer.printNearestGene();
		}else if(printCloseHighAttributeGenes){
			if(gtf && !fromDB && !ap.hasKey("gAttributeFile")){
				System.err.println("This option only works when using cuffdiff files or when a gene attribute file is provided");
				return;
			}
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
				String gene_names = StringUtils.join(startToGenename.get(nearestGene.getLocationString()), ","); 
				System.out.println(p.getLocationString()+"\t"+gene_names+"\t"+Integer.toString(MINDISTANCE));
			}else{
				System.out.println(p.getLocationString()+"\tNULL\tNULL");
			}
			
		}
		
	}
	
	public void printCloseHighAttributeGenes(double threshold){
		try{
			if(gtf && !fromDB && gene_attributes.size() == 0){
				throw new Exception("This option needs geneattrbutes (Fold-change, diff-pvalue ...) via cuffdiff or any readDB loading or a gene attribute file!!");
			}
		}catch(Exception e){
			e.printStackTrace();
		}
		HashMap<String, List<StrandedPoint>> genesbyChrs = hashbychrom(genes);
		for(Point p: peaks){
			List<StrandedPoint> highattgenes = new ArrayList<StrandedPoint>();
			List<Integer> distances = new ArrayList<Integer>();
			boolean hasgene = false;
			if(genesbyChrs.containsKey(p.getChrom())){
				for(StrandedPoint sp : genesbyChrs.get(p.getChrom())){
					int distance = sp.distance(p);
					if(distance < radius  && gene_attributes.containsKey(sp.getLocationString()) && gene_attributes.get(sp.getLocationString()) > threshold){
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
				String gene_names_string = StringUtils.join(gene_names, ",");
				String distance_string = StringUtils.join(distances, ",");
				System.out.println(p.getLocationString()+"\t"+gene_names_string+"\t"+distance_string);
			}else{
				System.out.println(p.getLocationString()+"NULL\tNULL");
			}
		}
	}
	
	
	public void printAllCloseGenes(){
		if(!gtf || fromDB){
			this.printCloseHighAttributeGenes(0);
		}else{
			HashMap<String, List<StrandedPoint>> genesbyChrs = hashbychrom(genes);
			
			for(Point p: peaks){
				List<Integer> distances = new ArrayList<Integer>();
				List<StrandedPoint> closeGenes= new ArrayList<StrandedPoint>();
				boolean hasclosegenes = false;
				if(genesbyChrs.containsKey(p.getChrom())){
					for(StrandedPoint gene : genesbyChrs.get(p.getChrom())){
						int distance = gene.distance(p);
						if( distance < radius){
							closeGenes.add(gene);
							distances.add(distance);
							hasclosegenes = true;
						}
					}
				}
				
				if(hasclosegenes){
					List<String> gene_names = new ArrayList<String>();
					for(StrandedPoint sp : closeGenes){
						gene_names.addAll(startToGenename.get(sp.getLocationString()));
					}
					String gene_names_string = StringUtils.join(gene_names, ",");
					String distance_string = StringUtils.join(distances, ",");
					System.out.println(p.getLocationString()+"\t"+gene_names_string+"\t"+distance_string);
				}else{
					System.out.println(p.getLocationString()+"\tNULL\tNULL");
				}
				
			}
			
		}
		
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
	
	public void laodpeaks(String peaksfile){
		peaks = RegionFileUtilities.loadPeaksFromPeakFile(gen, peaksfile, -1);
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
					if(words.length >=2){
						if(words.length >=1 && !words[0].contains("#") &&  !words[0].equals("Region") && !words[0].equals("Position")){
							if(words[0].contains("-")){
								RegionParser rparser = new RegionParser(gen);
								Region q = rparser.execute(words[0]);
								peak_attributes.put(q.getMidpoint().getLocationString(), Double.parseDouble(words[1]));		            	
							}else{
								PointParser pparser = new PointParser(gen);
								Point p = pparser.execute(words[0]);
								if(words[1].contains("-Infinity") || words[1].contains("-inf")){
									words[1] = Double.toHexString(-1*Double.MAX_VALUE);
								}else if(words[1].contains("Infinity") || words[1].contains("inf")){
									words[1] = Double.toHexString(Double.MAX_VALUE);
								}
								peak_attributes.put(p.getLocationString(), Double.parseDouble(words[1]));
							}
						}
						
					}
				}
				reader.close();
			}catch(Exception e){
				e.printStackTrace();
			}
		}
		
	}
	
	// If gtf is false, it assumes loading from a cuffdiff file
	// will add options to load gene attributes from other readDB expts later
	/**
	 * 
	 * @param genefile
	 * @param attribute Is attribute name when loading from cuffdiff, and file name when loading from gtf file
	 */
	public void laodgenes(String genefile, String attribute){
		if(gtf){
			//genes = RegionFileUtilities.loadStrandedPointFromRefTssFile(gen, genefile);
			GTFReader gffreader = new GTFReader(new File(genefile), gen);
			List<GeneTUnit> geneObjs = gffreader.loadGenes();
			genes = new ArrayList<StrandedPoint>();
			startToGenename = new HashMap<String, List<String>>();
			for (GeneTUnit gu: geneObjs){
				genes.add(new StrandedPoint(gu.getTSS(),gu.getStrand()));
				if(startToGenename.containsKey(genes.get(genes.size()-1).getLocationString())){
					startToGenename.get(genes.get(genes.size()-1).getLocationString()).add(gu.getName());
				}else{
					List<String> g_names = new ArrayList<String>();
					g_names.add(gu.getName());
					startToGenename.put(genes.get(genes.size()-1).getLocationString(), g_names);
				}
			}
			
			// Loading gene attributes form the attribute file
			try {
				File attFile = new File(attribute);
				if(!attFile.isFile()){System.err.println("Invalid gene attribute filename");System.exit(1);}
				BufferedReader attReader = new BufferedReader(new FileReader(attFile));
				String line;
				gene_attributes = new HashMap<String, Double>();
				while ((line = attReader.readLine()) != null){
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
				attReader.close();
				
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
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
			            			if(words[6].contains("OK") && words[9].contains("-inf")){words[9] = Double.toHexString(-1*Double.MAX_VALUE);}
			            			if(words[6].contains("OK") && words[9].contains("inf")){words[9] = Double.toHexString(Double.MAX_VALUE);}
			            			if(words[6].contains("NOTEST") || words[6].contains("FAIL") || words[6].contains("LOWDATA") || words[6].contains("HIDATA")){words[9] = Double.toHexString(0);}
			            			gene_attributes.put(genes.get(genes.size()-1).getLocationString(), Double.parseDouble(words[9]));
			            			
			            		}else if(attribute.toLowerCase().contains("tstat")){
			            			if(words[6].contains("OK") && words[9].contains("-inf") && words[10].contains("nan")){words[10] = Double.toHexString(-1*Double.MAX_VALUE);}
			            			if(words[6].contains("OK") && words[9].contains("inf") && words[10].contains("nan")){words[10] = Double.toHexString(Double.MAX_VALUE);}
			            			if(words[6].contains("NOTEST") || words[6].contains("FAIL") || words[6].contains("LOWDATA") || words[6].contains("HIDATA")){words[10] = Double.toHexString(0);}
			            			gene_attributes.put(genes.get(genes.size()-1).getLocationString(), Double.parseDouble(words[10]));
			            		}else if(attribute.toLowerCase().contains("pvalue")){
			            			gene_attributes.put(genes.get(genes.size()-1).getLocationString(), -1*Math.log10(Double.parseDouble(words[11])));
			            		}else if(attribute.toLowerCase().contains("qvalue")){
			            			gene_attributes.put(genes.get(genes.size()-1).getLocationString(), -1*Math.log10(Double.parseDouble(words[12])));
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
