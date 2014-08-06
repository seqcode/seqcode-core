package edu.psu.compbio.seqcode.projects.shaun.viz;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import edu.psu.compbio.seqcode.genome.location.ExonicGene;
import edu.psu.compbio.seqcode.genome.location.Gene;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.RefGeneGenerator;
import edu.psu.compbio.seqcode.gse.utils.Pair;

public class RNASeqFigurePaintable extends FigurePaintable{

	private int topBorder=50, bottomBorder=50;
	private int leftBorder=25, rightBorder=25;
	private ArrayList<String> exptNames;
	protected RefGeneGenerator<Region> geneGen=null; 
	protected ArrayList<ArrayList<Gene>> geneSets=new ArrayList<ArrayList<Gene>>();
	protected ArrayList<String> geneSetNames = new ArrayList<String>();
	protected HashMap<String, ArrayList<Pair<Gene,Double>>> junctions = new HashMap<String, ArrayList<Pair<Gene,Double>>>();
	
	public RNASeqFigurePaintable(FigureOptions opts){
		options = opts;
		reverseIt=options.reverseOrder;
		chr = options.gRegion.getChrom();
		rstart = options.gRegion.getStart();
		rstop = options.gRegion.getEnd();
		exptNames = options.exptNames;
		
		//Initialize genes
		if(options.useDBGenes){
			geneSetNames.add("Reference");
			ArrayList<Gene> geneSet = new ArrayList<Gene>();
			geneGen = new RefGeneGenerator<Region>(options.genome, "refGene");
			geneGen.retrieveExons(true);
			geneGen.setWantAlias(true);
			Iterator<Gene> gi = geneGen.execute(options.gRegion);
			while(gi.hasNext()) {
				geneSet.add(gi.next());
			}geneSets.add(geneSet);
		}if(options.transcriptGTF != null){//Load GTF
			geneSetNames.add(options.transcriptGTF.getName());
			ArrayList<Gene> geneSet = new ArrayList<Gene>();
			geneSet.addAll(loadGenes(options.transcriptGTF, options.gRegion));
			geneSets.add(geneSet);
		} 
		
		//Initialize junction
		for(String t : exptNames){
			if(options.experiments.get(t).junctionsFile!=null){
				ArrayList<Pair<Gene,Double>> scoredJunctions = loadScoredBED(options.experiments.get(t).junctionsFile, options.gRegion);
				junctions.put(t, scoredJunctions);
			}
		}
	    
		//Initialize experiment painters
			
	
	}
	
	public void paintItem (Graphics g, int x1, int y1, int x2, int y2){
		Graphics2D g2d = (Graphics2D)g;
		FontMetrics metrics = g2d.getFontMetrics();
		int screenSizeX = x2-x1;
		int screenSizeY = y2-y1;
		//g2d.setColor(Color.white);
		//g2d.fillRect(0, 0, screenSizeX, screenSizeY);
		topBound = topBorder;
		bottomBound = screenSizeY-bottomBorder;
		leftBound = leftBorder+30;
		rightBound = screenSizeX-rightBorder;
		baseLine = topBound;
			

		//Experiment tracks
		int offset=0;
		for(String t : exptNames){
			System.err.println("Processing: "+t);
    		//Experiment label
   			if(options.drawExptLabels){
   				g2d.setColor(Color.gray);
   	    		g2d.setFont(new Font("Ariel", Font.BOLD, options.labelFontSize));
   	    		metrics = g2d.getFontMetrics();
   				AffineTransform oldtrans = g2d.getTransform();
   				AffineTransform newtrans = new AffineTransform();
   		        newtrans.translate(leftBound-(metrics.getHeight()), baseLine+offset+metrics.getHeight()+(options.experiments.get(t).exptTrackHeight/2)); 
   		        newtrans.rotate(Math.toRadians(-90));
   		        g2d.setTransform(newtrans);
   		        g2d.drawString(t,(-1*metrics.stringWidth(t))/2,0);
   		        g2d.setTransform(oldtrans);
   			}    		
    		//RNA-seq paired reads for this experiment (if any)
    		
    		//RNA-seq read depth for this experiment (if any)
    		
    		//Junctions for this experiment (if any)
    		if(junctions.containsKey(t)){
    			offset += drawJunctions(g2d, x1, baseLine+offset, x2, junctions.get(t));
    		}
		}
		
		//Draw some coordinates & axis
		drawCoordinates(g2d, x1, baseLine+offset, x2);
		offset+=20;
		
		//Gene tracks
		for(int i=0; i<geneSetNames.size(); i++){
			offset += drawGenes(g2d, x1, baseLine+offset, x2, geneSets.get(i), options.drawGeneLabels, geneSetNames.get(i));
		}
 	}
	
	
	public ArrayList<Pair<Gene,Double>> loadScoredBED(File jFile, Region reg){
		ArrayList<Pair<Gene,Double>> sJuncts = new ArrayList<Pair<Gene,Double>>();
		try {
			if(!jFile.isFile()){System.err.println("Invalid junctions filename");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(jFile));
	        String line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\t");
	            if(words.length==12){
	            	String chr = words[0].replaceAll("chr", "");
		            if(chr.equals("MT"))
		            	chr = "M";
		            Integer start = new Integer(words[1]);
		            Integer stop = new Integer(words[2]);
		            String name = words[3];
		            Double score = new Double(words[4]);
		            char strand = words[5].charAt(0);
		            String [] blockSizes = words[10].split(",");
		            String [] blockStarts = words[11].split(",");
		            
		            ExonicGene currGene = new ExonicGene(options.genome, chr, start, stop, name, name, strand, "user");
		    		if(currGene.overlaps(options.gRegion)){
		    			for(int b=0; b<blockSizes.length; b++){
		    				Integer bSize = new Integer(blockSizes[b]);
		    				Integer bStartOff = new Integer(blockStarts[b]);
		    				Region block = new Region(options.genome, chr, start+bStartOff, start+bStartOff+bSize);
		    				currGene.addExon(block);
		    			}sJuncts.add(new Pair(currGene, score));
		    		}
	            }
	        }
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return(sJuncts);
	}
}
