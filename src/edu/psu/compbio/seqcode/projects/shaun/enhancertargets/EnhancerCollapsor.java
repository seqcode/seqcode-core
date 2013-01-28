package edu.psu.compbio.seqcode.projects.shaun.enhancertargets;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;

import edu.psu.compbio.seqcode.gse.datasets.general.NamedPoint;
import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.ewok.verbs.PointParser;
import edu.psu.compbio.seqcode.gse.ewok.verbs.RegionParser;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;

public class EnhancerCollapsor {
	private Genome gen =null;
	private HashMap<String,ArrayList<Region>> peakRegions = new HashMap<String,ArrayList<Region>>();
	private HashMap<String,ArrayList<Point>> peakPoints = new HashMap<String,ArrayList<Point>>();
	private HashMap<String, ArrayList<BindingRegion>> bindingRegions = new HashMap<String,ArrayList<BindingRegion>>();
	private int window=200;
	private boolean limitRegionsToWin=false;
	private int pointCount=0;
	
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("species") || !ap.hasKey("peaks")) { 
            System.err.println("Usage:\n " +
                               "EnhancerCollapsor \n" +
                               " Required: \n" +
                               "  --species <organism;genome> " +
                               "  --peaks <files containing coordinates of primary peaks> \n" +
                               " More Info: \n"+
                               "  --win <window of sequence around positive/negative points> \n" +
                               "  --out <output filename> \n"+
                               "");
            return;
        }
        try {
    		Pair<Organism, Genome> pair = Args.parseGenome(args);
    		Collection<String> peakSet = ap.hasKey("peaks") ? Args.parseStrings(args, "peaks") : null;
    		int win = ap.hasKey("win") ? new Integer(ap.getKeyValue("win")).intValue():200;
    		String out = ap.hasKey("out") ? new String(ap.getKeyValue("out")):"out.txt";
        
        	//initialize
			EnhancerCollapsor collapsor = new EnhancerCollapsor(pair.cdr());
			collapsor.setWin(win);
			
			//load positive & negative sets
			collapsor.loadPeaks(peakSet);

			//collapse 
			collapsor.collapse();
			
			//print
			collapsor.printBindingRegions(out);
        } catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public EnhancerCollapsor(Genome g){
		gen=g;
	}
	public void setWin(int w){window=w;}

	//print
	public void printBindingRegions(){
		for(String chr : bindingRegions.keySet()){
			for(BindingRegion b : bindingRegions.get(chr)){
				System.out.println(b);
			}
		}
	}
	//print to file
	public void printBindingRegions(String outName){
		try{
			FileWriter fw = new FileWriter(outName);
			for(String chr : bindingRegions.keySet()){
				for(BindingRegion b : bindingRegions.get(chr)){
					fw.write(b.toString()+"\n");
				}
			}
			fw.close();
			System.out.println("Results written to "+outName);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	//Load peaks
	public void loadPeaks(Collection<String> peakF){
		if(peakF!=null){
			for(String f : peakF){
				String n = f.replaceAll(".peaks", "").replaceAll(".txt", "").replaceAll("_signal", "");
				String[] x = n.split("/");
				peakRegions.put(x[x.length-1], loadRegionsFromPeakFile(f, window));
				peakPoints.put(x[x.length-1], loadPeaksFromPeakFile(f, window));
			}
		}
		
		//Put this peak into the binding regions list
		for(String name : peakPoints.keySet()){
			for(Point p : peakPoints.get(name)){
				if(!bindingRegions.containsKey(p.getChrom())){
					bindingRegions.put(p.getChrom(), new ArrayList<BindingRegion>());
				}
				bindingRegions.get(p.getChrom()).add(new BindingRegion(new NamedPoint(p,name)));
				pointCount++;
			}
		}
		
		System.out.println("Initial peak count: "+pointCount);
		
		//Sort the binding region collections
		for(String chr : bindingRegions.keySet()){
			Collections.sort(bindingRegions.get(chr));
		}
		System.out.println("BindingRegions sorted");
	}
	
	//Load a set of regions from a peak file
	public ArrayList<Region> loadRegionsFromPeakFile(String filename, int win){
		ArrayList<Region> regs = new ArrayList<Region>();
		try{
			File pFile = new File(filename);
			if(!pFile.isFile()){System.err.println("Invalid positive file name: "+filename);System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(pFile));
	        String line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            if(!words[0].equals("Region")){
		            if(words.length>=3 && win!=-1 && limitRegionsToWin){
		                PointParser pparser = new PointParser(gen);
		            	Point p = pparser.execute(words[2]);
		            	int rstart = p.getLocation()-(win/2)<1 ? 1:p.getLocation()-(win/2);
	                	int rend = p.getLocation()+(win/2)>gen.getChromLength(p.getChrom()) ? gen.getChromLength(p.getChrom()):p.getLocation()+(win/2)-1;
	                	Region r = new Region(p.getGenome(), p.getChrom(), rstart, rend);
	                	regs.add(r);
	                }else if(words.length>=1){
		            	RegionParser parser = new RegionParser(gen);
		            	Region r = parser.execute(words[0]);
		            	if(r!=null){regs.add(r);}
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
		return(regs);
	}
	//Load a set of peaks from a peak file
	public ArrayList<Point> loadPeaksFromPeakFile(String filename, int win){
		ArrayList<Point> points = new ArrayList<Point>();
		try{
			File pFile = new File(filename);
			if(!pFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(pFile));
	        String line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            if(!words[0].equals("Region")){
		            if(words.length>=3 && win!=-1){
		                PointParser pparser = new PointParser(gen);
		            	Point p = pparser.execute(words[2]);
		            	points.add(p);
	                }else if(words.length>=1){
		            	RegionParser parser = new RegionParser(gen);
		            	Region r = parser.execute(words[0]);
		            	if(r!=null){points.add(new Point(r.getGenome(), r.getChrom(), (r.getStart()+r.getEnd())/2));}
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
		return(points);
	}

	//Collapse
	public void collapse(){
		int currCount = pointCount, lastCount = pointCount;
		do{
			lastCount=currCount;
			for(String chr : bindingRegions.keySet()){
				if(bindingRegions.get(chr).size()>0){
					for(int i=0; i<bindingRegions.get(chr).size()-1; i++){
						BindingRegion curr = bindingRegions.get(chr).get(i);
						BindingRegion next = bindingRegions.get(chr).get(i+1);
						if(curr.reg.distance(next.reg)<window){
							curr.addPoint(next);
							bindingRegions.get(chr).remove(i+1);
							i--;
							pointCount--;
						}
					}
				}
			}currCount = pointCount; 
			System.out.println("After collapsing: "+pointCount+" points");
		}while(currCount<lastCount);
	}
}
