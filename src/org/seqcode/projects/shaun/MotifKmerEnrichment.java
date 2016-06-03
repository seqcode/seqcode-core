package org.seqcode.projects.shaun;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.geom.AffineTransform;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;

import javax.imageio.ImageIO;

import org.seqcode.gse.tools.utils.Args;


public class MotifKmerEnrichment {
	protected static int sortingExpt=0;
	protected ArrayList<HashMap<String, Double>> kmerCounts = new ArrayList<HashMap<String, Double>>();
	protected ArrayList<Kmer> instances = new ArrayList<Kmer>();
	protected ArrayList<String> exptNames=new ArrayList<String>();
	protected ArrayList<Double> totalCounts = new ArrayList<Double>();
	protected int exptCount=0;
	protected String outFilename="";
	
	protected int motifWidth = 200;
	protected int motifSpacing= 20;
	protected int imageHeight = 600;
	protected int imageBorder = 20;
	
	protected Color overMaxColor = Color.yellow;
	protected Color overMidColor = Color.black;
	protected Color overMinColor = Color.blue;
	protected double maxOver=2.0, midOver=0, minOver=-2;
	
	//Constructor
	public MotifKmerEnrichment(Collection<String> filenames, String out){		
		outFilename=out;
		
		//Load the kmer counts from the files
		for(String f : filenames){
			exptNames.add(f);
			kmerCounts.add(loadKmersFromFile(f));
			exptCount++;
		}
		
	}
	
	//execute the analysis
	public void execute(){
		HashMap<String, Integer> addedIndex = new HashMap<String, Integer>();
		int addedCount = 0;
		//Put the contents of the HashMaps into the Kmer ArrayList
		int eCount=0;
		for(HashMap<String, Double> count : kmerCounts){
			totalCounts.add(eCount, 0.0);
			for(String k : count.keySet()){
				int index;
				if(addedIndex.containsKey(k))
					index = addedIndex.get(k);
				else{
					index = addedCount;
					addedIndex.put(k, index);
					Kmer mer = new Kmer(k, exptCount);
					instances.add(mer);
					addedCount++;					
				}
				Kmer mer =instances.get(index);
				mer.counts[eCount] = count.get(k);
				totalCounts.set(eCount, totalCounts.get(eCount)+count.get(k));
			}
			eCount++;
		}
		
		sortingExpt=0;
		//Test sort
		Collections.sort(instances);
		
		//Test freq conversion
		System.out.print("Kmer");
		for(String e : exptNames)
			System.out.print("\t"+e);
		System.out.println("");
		for(Kmer mer : instances){
			System.out.print(mer.kmer);
			for(int x=0; x<mer.counts.length; x++){
				double freq = mer.counts[x]/totalCounts.get(x);
				System.out.print("\t"+freq);
			}
			System.out.println("");
		}
		
		//Draw the motif k-mers
		drawKmers(outFilename);
		
	}
	
	//Draw the k-mers
	public void drawKmers(String imageOut){
		String filename = imageOut + ".png";
		File f = new File(filename);
		int w = (imageBorder*2) + (exptCount*motifWidth) + (motifSpacing*(exptCount-1));
		int h = imageHeight;
	    BufferedImage im = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
	    Graphics g = im.getGraphics();
	    Graphics2D g2 = (Graphics2D)g;
	    g2.setRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
	    //g2.setColor(Color.lightGray);
	    g2.setColor(Color.white);
	    g2.fillRect(0, 0, w, h);
	    
	    int currH = imageBorder+20;
	    int currW = imageBorder;
	    
	    //Draw names
	    for(int e=0; e <exptCount; e++){
	    	g2.setColor(Color.black);
	    	g2.drawString(exptNames.get(e), currW, currH);
	    	currW += motifWidth+motifSpacing;
	    }
	    currW = imageBorder;
	    currH = imageBorder+20+10;
	    //Draw kmers
	    double availHeight = (double)(imageHeight-currH-imageBorder);
	    Font baseFont = new Font("Lucida Sans",Font.BOLD,(int)(availHeight));
	    g2.setFont(baseFont);
	    FontMetrics fontmetrics = g2.getFontMetrics();
	    FontMetrics newFontMet;
	    for(int e=0; e <exptCount; e++){
	    	currH = imageBorder+20+10;
	    	sortingExpt=e;
			Collections.sort(instances);
			
			for(Kmer mer : instances){
				double freq = (mer.counts[e]/totalCounts.get(e));
				double currKPixels = freq*availHeight;
				g2.setColor(Color.black);
				
				double stringWidth = fontmetrics.stringWidth(mer.kmer);
				
				AffineTransform transform = new AffineTransform((1/(stringWidth/motifWidth)),0,0,freq,0,0);
                Font thisFont = baseFont.deriveFont(transform);
                g2.setFont(thisFont);
                newFontMet = g2.getFontMetrics();
                
                //Calculate over-rep
                double sum = 0;
                for(int z=0; z <exptCount; z++)
                	if(z!=e)
                		sum+=mer.counts[z]/totalCounts.get(z);
                double mean = sum/(double)(exptCount-1);
                double over = freq/mean;

                Color c = overColor(Math.log(over)/Math.log(2));
                g2.setColor(c);
                
                g2.drawString(mer.kmer, currW, currH+newFontMet.getHeight()-(newFontMet.getLeading()));
                currH += currKPixels;
			}
			currW+= motifWidth+motifSpacing;
	    }
	    
	    try{
	    	ImageIO.write(im, "png", f);
	    }
	    catch(IOException e){
	    	System.err.println("Error in printing file "+filename);
	    }
	}
	
	//convert color
	private Color overColor(double v){
		Color c;
		if(v>midOver){
			Color maxColor = overMaxColor;
			Color minColor = overMidColor;
			double sVal = v>maxOver ? 1 : (v-midOver)/(maxOver-midOver);
			
			int red = (int)(maxColor.getRed() * sVal + minColor.getRed() * (1 - sVal));
		    int green = (int)(maxColor.getGreen() * sVal + minColor.getGreen() * (1 - sVal));
		    int blue = (int)(maxColor.getBlue() *sVal + minColor.getBlue() * (1 - sVal));
		    c = new Color(red, green, blue);
		}else{
			Color maxColor = overMidColor;
			Color minColor = overMinColor;
			double sVal = v<minOver ? 0 : ((midOver-minOver)-(midOver-v))/(midOver-minOver);
			int red = (int)(maxColor.getRed() * sVal + minColor.getRed() * (1 - sVal));
	        int green = (int)(maxColor.getGreen() * sVal + minColor.getGreen() * (1 - sVal));
	        int blue = (int)(maxColor.getBlue() *sVal + minColor.getBlue() * (1 - sVal));
	        c = new Color(red, green, blue);				
		}
		return(c);
	}
	
	//load k-mer instances
	public HashMap<String, Double> loadKmersFromFile(String fname){
		HashMap<String, Double> kmers = new HashMap<String, Double>();
		try{
			File pFile = new File(fname);
			if(!pFile.isFile()){System.err.println("Invalid file name: "+fname);System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(pFile));
	        String line = reader.readLine();
	        while ((line = reader.readLine()) != null) {
	        	line = line.trim();
	        	String[] words = line.split("\\s+");
	        	
	        	if(words.length>=3){
		        	String kmer = words[2].toUpperCase();
		        	if(kmers.containsKey(kmer)){
		        		kmers.put(kmer, kmers.get(kmer)+1.0);
		        	}else{
		        		kmers.put(kmer, 1.0);
		        	}
	        	}
	        }reader.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return(kmers);
	}
	
	//main
	public static void main(String[] args){
		if(args.length==0){
			System.err.println("MotifKmerEnrichment:\n\t--kmers <file containing motif hits>");
			System.exit(1);
		}else{
			Collection<String> files = Args.parseStrings(args, "kmers");
			String outName = Args.parseString(args, "out", "out");
			MotifKmerEnrichment mke = new MotifKmerEnrichment(files, outName);
			mke.execute();
		}		
	}
	
	public class Kmer implements Comparable<Kmer>{
		public String kmer;
		public Double[] counts;
		
		public Kmer(String k, int numExpt){
			kmer = k;
			counts = new Double[numExpt];
			for(int i=0; i<numExpt; i++){ counts[i]=0.0;}
		}
		public int compareTo(Kmer k){
			if(this.counts[sortingExpt] > k.counts[sortingExpt])
				return -1;
			else if (this.counts[sortingExpt] < k.counts[sortingExpt])
				return 1;
			return 0;
		}
	}
}
