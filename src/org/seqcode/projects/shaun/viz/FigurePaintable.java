package org.seqcode.projects.shaun.viz;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import javax.swing.Scrollable;

import org.seqcode.genome.location.ExonicGene;
import org.seqcode.genome.location.Gene;
import org.seqcode.genome.location.Region;
import org.seqcode.gse.datasets.seqdata.SeqLocator;
import org.seqcode.gse.projects.gps.DeepSeqExpt;
import org.seqcode.gse.projects.gps.ExtReadHit;
import org.seqcode.gse.projects.gps.ReadHit;
import org.seqcode.gse.seqview.paintable.NonOverlappingLayout;
import org.seqcode.gse.utils.Pair;
import org.seqcode.gse.viz.paintable.AbstractPaintable;


public abstract class FigurePaintable extends AbstractPaintable{

	protected FigureOptions options;
	protected NonOverlappingLayout<Gene> layout;
	protected boolean reverseIt;
	protected int topBound, bottomBound, leftBound, rightBound, baseLine;
	protected int rstart, rstop;
	protected String chr;
	
	protected List<Region> getReads(SeqLocator loc, Region reg, int perBaseMax, boolean isPaired){
		List<Region> reads = new ArrayList<Region>();
		List<SeqLocator> loclist = new ArrayList<SeqLocator>();
		loclist.add(loc);
		DeepSeqExpt dse = new DeepSeqExpt(options.genome, loclist, "readdb", -1);
		dse.setPairedEnd(isPaired);
		dse.setThreePrimeExt(options.readExt);
		
		List<ExtReadHit> hits = dse.loadExtHits(reg);
		HashMap<ExtReadHit, Double> readCounter = new HashMap<ExtReadHit, Double>();
		for(ExtReadHit h : hits){
			if(reg.contains(h)){ //(some paired end hits may not be in region...}
				if(readCounter.containsKey(h))
					readCounter.put(h, readCounter.get(h)+h.getWeight());
				else
					readCounter.put(h, h.getWeight());
				
				if(readCounter.get(h)<= perBaseMax)
					reads.add(h);
			}
		}
		return reads;
	}
	
	protected List<Region> getPairs(SeqLocator loc, Region reg, int perBaseMax, boolean isPaired){
		List<Region> pairs = new ArrayList<Region>();
		List<SeqLocator> loclist = new ArrayList<SeqLocator>();
		loclist.add(loc);
		DeepSeqExpt dse = new DeepSeqExpt(options.genome, loclist, "readdb", -1);
		dse.setPairedEnd(isPaired);
		dse.setThreePrimeExt(options.readExt);
		
		List<ReadHit> hits = dse.loadPairsAsSingle(reg);
		HashMap<ReadHit, Double> readCounter = new HashMap<ReadHit, Double>();
		for(ReadHit h : hits){
			if(reg.contains(h)){ //(some paired end hits may not be in region...}
				if(readCounter.containsKey(h))
					readCounter.put(h, readCounter.get(h)+h.getWeight());
				else
					readCounter.put(h, h.getWeight());
				
				if(readCounter.get(h)<= perBaseMax)
					pairs.add(h);
			}
		}
		return pairs;
	}
	
	protected void drawCoordinates(Graphics2D g2d, int x1, int y, int x2){
		g2d.setColor(Color.black);
		g2d.drawLine(leftBound, y, rightBound, y);
		g2d.setFont(new Font("Ariel", Font.PLAIN, options.fontSize));
		FontMetrics metrics = g2d.getFontMetrics();
		AffineTransform oldtrans = g2d.getTransform();
        AffineTransform newtrans = new AffineTransform();
        String text = reverseIt ? new String("chr"+chr+":"+rstop) : new String("chr"+chr+":"+rstart);
        newtrans.translate(leftBound, y+metrics.stringWidth(text)+metrics.getHeight()+10);
        newtrans.rotate(Math.toRadians(-90));
        g2d.setTransform(newtrans);
        g2d.drawString(text,0,0);
        g2d.setTransform(oldtrans);
        newtrans = new AffineTransform();
        text = reverseIt ? new String("chr"+chr+":"+rstart) : new String("chr"+chr+":"+rstop);
        newtrans.translate(rightBound+(metrics.getHeight()), y+metrics.stringWidth(text)+metrics.getHeight()+10); 
        newtrans.rotate(Math.toRadians(-90));
        g2d.setTransform(newtrans);
        g2d.drawString(text,0,0);
        g2d.setTransform(oldtrans);
	}
	protected int drawGenes(Graphics2D g2d, int x1, int y, int x2, List<Gene> genes, boolean drawNames){return(drawGenes(g2d, x1, y, x2, genes, drawNames, ""));}
	protected int drawGenes(Graphics2D g2d, int x1, int y, int x2, List<Gene> genes, boolean drawNames, String trackName){
		layout = new NonOverlappingLayout<Gene>();
		layout.setRegions(genes);
		int[] a = new int[7];
        int[] b = new int[7];
        FontMetrics metrics = g2d.getFontMetrics();
        for(Gene g : genes){
			if(g instanceof ExonicGene) { 
				int track = 0;
	            if(!layout.hasTrack(g)) { 
	                System.err.println("No track assigned to gene: " + g.getName());
	            } else { 
	                track = layout.getTrack(g);
	            }
	            
                ExonicGene gene = (ExonicGene)g;
	            int gy1 = y+ (2 * options.geneHeight * track);
				int halfGeneHeight = options.geneHeight / 2;
	            int gmy = gy1 + halfGeneHeight;
	            
	            int geneStart = gene.getStart(), geneEnd = gene.getEnd();
	            boolean strand = reverseIt ? (gene.getStrand() == '-'):(gene.getStrand() == '+');
	            int gx1 = reverseIt ? xcoord(geneEnd):xcoord(geneStart);
	            int gx2 = reverseIt? xcoord(geneStart):xcoord(geneEnd);
	            gx1=Math.max(leftBound, gx1);
	            gx2=Math.min(rightBound, gx2);
	            
	            g2d.setColor(Color.black);
	            
	            
            	g2d.drawLine(gx1, gmy, gx2, gmy);
	            if((g.getStrand()!='.') && (g.getFivePrime()>=options.gRegion.getStart() && g.getFivePrime()<=options.gRegion.getEnd())){
	            	arrangeArrow(a, b, strand, halfGeneHeight, gx1, gx2, gmy);
	            	g2d.drawPolyline(a, b, 7);
	            }
	            
	            Iterator<Region> exons = gene.getExons();
	            while(exons.hasNext()) { 
	                Region exon = exons.next();
	                if(exon.overlaps(options.gRegion)){
		                int ex1 = reverseIt ? xcoord(exon.getEnd()) : xcoord(exon.getStart());
		                int ex2 = reverseIt?xcoord(exon.getStart()): xcoord(exon.getEnd());
		                int eleft = Math.max(leftBound, ex1);
		                int eright = Math.min(rightBound, ex2);
		
		                int rectwidth = eright - eleft + 1;
		                
		                if(g.getName().equals(options.specialGeneName)){
		                	g2d.setColor(options.specialGeneColor);
		                }else{
		                	g2d.setColor(options.geneColor);
		                }
		                g2d.fillRect(eleft, gy1, rectwidth, options.geneHeight);
		                g2d.setColor(Color.black);
		                g2d.drawRect(eleft, gy1, rectwidth, options.geneHeight);
	                }
	            }
	            //Gene name
	            if(drawNames){
	            	g2d.setColor(Color.gray);
	            	g2d.setFont(new Font("Ariel", Font.BOLD, options.geneFontSize));
	   	    		metrics = g2d.getFontMetrics();
	   				AffineTransform oldtrans = g2d.getTransform();
	   				AffineTransform newtrans = new AffineTransform();
	   				int nx = strand ? gx1 : gx2;
	   				int ny = gmy+(options.geneHeight)+(metrics.getHeight()*2);
	   				int str = strand ? 1 : 0;
	   				newtrans.translate(nx, ny);
	   		        newtrans.rotate(Math.toRadians(-90));
	   		        g2d.setTransform(newtrans);
	   		        g2d.drawString(gene.getName(),-1*metrics.stringWidth(gene.getName()),str*(metrics.getHeight()));
	   		        g2d.setTransform(oldtrans);
	            }
			}
		}
        
        int finalOff = layout.getNumTracks()*2*options.geneHeight+options.geneFontSize;
        //Track name
        if(options.drawTrackNames){
			g2d.setColor(Color.lightGray);
    		g2d.setFont(new Font("Ariel", Font.BOLD, options.geneFontSize));
    		metrics = g2d.getFontMetrics();
			g2d.drawString(trackName,(leftBound+rightBound)/2-metrics.stringWidth(trackName)/2,y+finalOff);
	        g2d.setStroke(new BasicStroke(1.0f));
    		//g2d.drawLine(leftBound, y+finalOff, rightBound, y+finalOff);
		}
        return(finalOff+10);
	}
	
	protected int drawJunctions(Graphics2D g2d, int x1, int y, int x2, List<Pair<Gene, Double>> scoredJunctions){
		ArrayList<Gene> juncts = new ArrayList<Gene>();
		double maxScore = 0;
		for(Pair<Gene,Double> p : scoredJunctions){ 
			juncts.add(p.car());
			if(p.cdr()>maxScore)
				maxScore = p.cdr();
		}
		layout = new NonOverlappingLayout<Gene>();
		layout.setRegions(juncts);
		g2d.setColor(Color.black);
		g2d.setFont(new Font("Ariel", Font.PLAIN, options.geneFontSize));
        FontMetrics metrics = g2d.getFontMetrics();
        int maxY = y+ (2 * options.geneHeight * layout.getNumTracks());
        for(Pair<Gene,Double> p : scoredJunctions){
			Gene g = p.car();
			Double s = p.cdr();
        	if(g instanceof ExonicGene) { 
				int track = 0;
	            if(!layout.hasTrack(g)) { 
	                System.err.println("No track assigned to gene: " + g.getName());
	            } else { 
	                track = layout.getTrack(g);
	            }
	            
                ExonicGene gene = (ExonicGene)g;
	            int gy1 = maxY - (2 * options.geneHeight * track) - options.geneHeight;
				int halfGeneHeight = options.geneHeight / 2;
	            int gmy = gy1 + halfGeneHeight;
	            
	            int geneStart = gene.getStart(), geneEnd = gene.getEnd();
	            boolean strand = reverseIt ? (gene.getStrand() == '-'):(gene.getStrand() == '+');
	            int gx1 = reverseIt ? xcoord(geneEnd):xcoord(geneStart);
	            int gx2 = reverseIt? xcoord(geneStart):xcoord(geneEnd);
	            gx1=Math.max(leftBound, gx1);
	            gx2=Math.min(rightBound, gx2);
	            
	            g2d.setColor(Color.black);
            	g2d.drawLine(gx1, gmy, gx2, gmy);
	            
            	//Score
            	String scoreStr = String.format("%.0f", s);
            	g2d.drawString(scoreStr, (gx1+gx2)/2-metrics.stringWidth(scoreStr)/2, gy1);
            	
	            Iterator<Region> exons = gene.getExons();
	            while(exons.hasNext()) { 
	                Region exon = exons.next();
	                if(exon.overlaps(options.gRegion)){
		                int ex1 = reverseIt ? xcoord(exon.getEnd()) : xcoord(exon.getStart());
		                int ex2 = reverseIt?xcoord(exon.getStart()): xcoord(exon.getEnd());
		                int eleft = Math.max(x1, ex1);
		                int eright = Math.min(x2, ex2);
		
		                int rectwidth = eright - eleft + 1;
		                
		                double percScore = s/maxScore;
		                Color cx = new Color((int)(255-(255*percScore)), (int)(255-(255*percScore)), (int)(255-(255*percScore)));
		                g2d.setColor(cx);
			           
		                g2d.fillRect(eleft, gy1, rectwidth, options.geneHeight);
		                g2d.setColor(Color.black);
		                g2d.drawRect(eleft, gy1, rectwidth, options.geneHeight);
	                }
	            }
			}
		}
        int finalOff = layout.getNumTracks()*2*options.geneHeight;
        return(finalOff+10);
	}
	
	protected int xcoord(int coord) { 
		double frac = ((double)(coord-rstart)/(double)(rstop-rstart));
		int pos;
		if(reverseIt)
			pos = leftBound+(int)((double)(rightBound-leftBound)*(1-frac));
		else
			pos= leftBound+(int)((double)(rightBound-leftBound)*frac);
		return(pos);
	}
	
	protected void arrangeArrow(int[] a, int[] b, boolean strand, int geneHeight, int gx1, int gx2, int my) { 
        double arrowHt = 0.15 * geneHeight; 
        double arrowWd = 1;
        int a1, a2, a3;
        
                
        // forward arrow
        if(strand) {
            int startX = gx1;
            a1 = startX; 
            a2 = (int) Math.round(startX + (arrowWd * 6)); 
            a3 = (int) Math.round(startX + (arrowWd * 10)); 
            
        } else { 
            // backward arrow
            int startX = gx2+1;
            a1 = startX; 
            a2 = (int) Math.round(startX - (arrowWd * 6)); 
            a3 = (int) Math.round(startX - (arrowWd * 10)); 
        }
        
        a[0] = a1;
        a[1] = a1;
        a[2] = a2;
        a[3] = a2;
        a[4] = a3;
        a[5] = a2;
        a[6] = a2;

        int b1 = (int) Math.round(my);
        int b2 = (int) Math.round(my + (arrowHt * 13));
        int b3 = (int) Math.round(my + (arrowHt * 10));
        int b4 = (int) Math.round(my + (arrowHt * 16));
        
        b[0] = b1;
        b[1] = b2;
        b[2] = b2;
        b[3] = b3;
        b[4] = b2;
        b[5] = b4;
        b[6] = b2;
    }
	public ArrayList<Gene> loadGenes(File gtf, Region reg){
		ArrayList<Gene> theGenes = new ArrayList<Gene>();
		
		try {
			if(!gtf.isFile()){System.err.println("Invalid GTF filename");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(gtf));
	        String currGeneName=null;
	        String currTransName=null;
	        String currTrans=""; 
	        ArrayList<Region> exons=new ArrayList<Region>();
	        Integer geneStart=-1, geneStop=-1; String geneChr=""; char geneStrand='+';
	        String line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\t");
	            if(words[2].equals("exon")){
		            String chr = words[0].replaceAll("chr", "");
		            if(chr.equals("MT"))
		            	chr = "M";
		            Integer start = new Integer(words[3]);
		            Integer stop = new Integer(words[4]);
		            char strand = words[6].charAt(0);
		            String transID = ""; String transName=null; String geneName=null;
		            String[] dwords = line.split(";");
		            for(String d : dwords){
		            	String [] dxwords = d.split("\\s+");
		            	for(int s=0; s<dxwords.length-1; s++){
		            		if(dxwords[s].equals("transcript_id"))
		            			transID = dxwords[s+1].replaceAll("\"", "");
		            		if(dxwords[s].equals("transcript_name"))
		            			transName = dxwords[s+1].replaceAll("\"", "");
		            		if(dxwords[s].equals("gene_name"))
		            			geneName = dxwords[s+1].replaceAll("\"", "");
		            	}
		            }
		            Region exon = new Region(options.genome, chr, start, stop);
		            if(transID.equals(currTrans)){  //Add a new exon 
		            	exons.add(exon);
		            	if(stop>geneStop)
		            		geneStop = stop;
		            	if(start<geneStart)
		            		geneStart=start;
		            }else{ //Check if recorded gene should be saved. Record a new gene.
		            	if(geneStart != -1){
		            		ExonicGene currGene = new ExonicGene(options.genome, geneChr, geneStart, geneStop, currGeneName==null ? currTrans : currGeneName, currTrans, geneStrand, "user");
		            		if(currGene.overlaps(options.gRegion)){
		            			for(Region e : exons)
		            				currGene.addExon(e);
		            			theGenes.add(currGene);	
		            			System.out.println(currGene.getName()+"\t"+currGene.getLocationString());
		            		}
		            	}	
		            	//New gene
		            	geneStart = start; geneStop = stop; geneChr = chr; geneStrand = strand;
		            	exons=new ArrayList<Region>();
		            	exons.add(exon);
		            	currTrans = transID;
		            	currGeneName = geneName;
		            	currTransName = transName;
		            	
		            }
	            }
	        }
	        if(exons.size()>0){
		        ExonicGene currGene = new ExonicGene(options.genome, geneChr, geneStart, geneStop, currGeneName==null ? currTrans : currGeneName, currTrans, geneStrand, "user");
	    		if(currGene.overlaps(options.gRegion)){
	    			for(Region e : exons)
	    				currGene.addExon(e);
	    			theGenes.add(currGene);	
	    			System.out.println(currGene.getName()+"\t"+currGene.getLocationString());
	    		}
	        }
	        
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return theGenes;
	}
}
