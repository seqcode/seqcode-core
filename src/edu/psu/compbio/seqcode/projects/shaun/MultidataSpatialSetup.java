package edu.psu.compbio.seqcode.projects.shaun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Pattern;

import javax.swing.JFrame;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Species;
import edu.psu.compbio.seqcode.genome.location.ExonicGene;
import edu.psu.compbio.seqcode.genome.location.Gene;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.PointParser;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.RegionParser;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.StrandedRegionParser;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.viz.paintable.PaintableFrame;

public class MultidataSpatialSetup {

	private MultidataSpatialPaintable painter;
	private PaintableFrame plotter;
	private int screenSizeX=1000, screenSizeY=900;
	private String inputFile;
	private ArrayList<String> times=new ArrayList<String>();
	private ArrayList<ExonicGene> genes = new ArrayList<ExonicGene>();
	private HashMap<String, ArrayList<Double>> expression = new HashMap<String, ArrayList<Double>>();
	private HashMap<String, ArrayList<Region>> sites = new HashMap<String, ArrayList<Region>>();
	private ArrayList<Pair<String, StrandedRegion>> motifs = new ArrayList<Pair<String, StrandedRegion>>();
	private ArrayList<Region> lits = new ArrayList<Region>();
	private HashMap<String, SeqLocator> locs = new HashMap<String, SeqLocator>();
	private int rstart=0;
	private int rstop = 0;
	private String chr;
	private Genome gen;
	private Species org;
	

	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("data")) { 
            System.err.println("Usage:\n " +
                               "MultidataSpatialSetup " +
                               "--data <file name> ");
            return;
        }
        String dfile = ap.getKeyValue("data");
        MultidataSpatialSetup setup = new MultidataSpatialSetup(dfile);		
	}
	
	public MultidataSpatialSetup(String df){
		inputFile = df;
		try {
			org = Species.getSpecies("Mus musculus");
			gen = Genome.findGenome("mm8");
			
			//Load the file contents
			loadFile(inputFile);
			
			Region gRegion = new Region(gen, chr, rstart, rstop);
			
			//Paint the picture
			MultidataSpatialPaintable painter = new MultidataSpatialPaintable(times, gRegion, genes, expression, sites, motifs, lits, locs);
			plotter = new PaintableFrame("Genomic Data", painter);
			plotter.setSize(screenSizeX, screenSizeY);
			plotter.setVisible(true);
			plotter.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
	}
	
	private void loadFile(String inF){
		try{
			File aFile = new File(inF);
			if(aFile.isFile()){
				BufferedReader reader;
				reader = new BufferedReader(new FileReader(aFile));
				String line;
				while((line= reader.readLine())!=null){
					String [] tokens = line.split("\\t");
					if(tokens[0].equals("timepoints")){
						for(int i=1; i<tokens.length; i++){
							times.add(tokens[i]);
						}
					}else if(tokens.length==2){//Setup annotations
						if(tokens[0].equals("regionstart")){
							rstart = new Integer(tokens[1]).intValue();
						}
						if(tokens[0].equals("regionstop")){
							rstop = new Integer(tokens[1]).intValue();
						}
						if(tokens[0].equals("regionchr")){
							chr = new String(tokens[1]);
						}						
					}else if(tokens.length>=15){//A gene
						int start = new Integer(tokens[4]).intValue();
						int stop = new Integer(tokens[5]).intValue();
						String name = tokens[12];
						String id = tokens[1];
						char str = tokens[3].charAt(0);
						String [] estarts = tokens[9].split(",");
						String [] estops = tokens[10].split(",");
						ExonicGene e = new ExonicGene(gen, chr, start, stop, name, id, str, "Manual");
						for(int i=0; i<estarts.length; i++){
							Region r = new Region(gen, chr, new Integer(estarts[i]).intValue(), new Integer(estops[i]).intValue());
							e.addExon(r);
						}
						genes.add(e);System.out.println(name);
					}else if(tokens[0].equals("Expr")){//An expression set
						String name = tokens[1];
						ArrayList<Double> vals = new ArrayList<Double>();
						for(int i=2; i<tokens.length; i++){
							vals.add(new Double(tokens[i]));
						}
						expression.put(name, vals);
					}else if(tokens[0].equals("Site")){//A binding site (timepoint, position, zscore)
						String time = tokens[1];
						//PointParser pparser = new PointParser(gen);
		            	//Point p = pparser.execute(tokens[2]);
						RegionParser rparser = new RegionParser(gen);
						Region r = rparser.execute(tokens[2]);
		            	//Double zscore = new Double(tokens[3]);
		            	if(!sites.containsKey(time)){
		            		sites.put(time, new ArrayList<Region>());
		            	}sites.get(time).add(r);
					}else if(tokens[0].equals("Motif")){//Motif hits
						StrandedRegionParser parser = new StrandedRegionParser(gen);
						StrandedRegion s = parser.execute(tokens[2]);
		            	String type = tokens[1];
		            	if(s!=null)
		            		motifs.add(new Pair<String, StrandedRegion>(type, s));
					}else if(tokens[0].equals("Lit")){//Known sites
						Region site = new Region(gen, tokens[3], new Integer(tokens[4]).intValue(), new Integer(tokens[5]).intValue());
						lits.add(site);
					}else if(tokens[0].equals("Expt")){//Experiments
						Set<String> reps = new TreeSet<String>();
						SeqLocator l = new SeqLocator(tokens[2], reps, tokens[3]);
						locs.put(tokens[1], l);
					}
					
				}
				reader.close();
			}
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
