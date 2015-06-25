package edu.psu.compbio.seqcode.projects.shaun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.genome.Species;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.PointParser;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.location.RegionParser;
import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.RealValuedHistogram;

public class ChipSeqMetaPeakMaker {

	private Species org;
	private Genome gen;
	private ArrayList<Region> regions;
	private ArrayList<Point> peaks;
	private RealValuedHistogram histo;
	private int windowSize;
	
	//main
	public static void main(String[] args) throws SQLException, NotFoundException {
		int readLen=26, readExt=0;
		Pair<Species,Genome> pair = Args.parseGenome(args);
		if(pair==null){return;}
		List<SeqLocator> expts = Args.parseSeqExpt(args,"expt");
		String peaks = Args.parseString(args,"peaks","error");
		int histwin = Args.parseInteger(args,"histwin", 1000);
		int histbins = Args.parseInteger(args,"bins", 2000);
		
		ChipSeqMetaPeakMaker metamaker = new ChipSeqMetaPeakMaker(pair.car(), pair.cdr(), histwin, histbins);
		metamaker.loadPeaks(peaks);
		metamaker.execute(expts, readLen, readExt);
		metamaker.printMeta();
	}
	
	//Constructor
	public ChipSeqMetaPeakMaker(Species o, Genome g, int histwindow, int numBins){
		org=o;
		gen =g;
		regions = new ArrayList<Region>();
		peaks = new ArrayList<Point>();
		histo = new RealValuedHistogram(-1*histwindow, histwindow, numBins);
		windowSize = histwindow;
	}
	
	//Load the positive hits
	public void loadPeaks(String inFile){
	    File pFile = new File(inFile);
	    try {
			if(!pFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
			BufferedReader reader = new BufferedReader(new FileReader(pFile));
			
		    String line;
		    while ((line = reader.readLine()) != null) {
		        line = line.trim();
		        String[] words = line.split("\\s+");
		        PointParser pparser = new PointParser(gen);
		    	Point p = pparser.execute(words[2]);
		    	peaks.add(p);
		        RegionParser parser = new RegionParser(gen);
		       	Region r = parser.execute(words[0]);
		       	regions.add(r);
		    }
	    } catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public RealValuedHistogram execute(String [][] exp){
		return(execute(exp, 26.0, 174.0));
	}
	public RealValuedHistogram execute(String [][] exp, double readLength, double readExtension){
		ArrayList<SeqLocator> locs = new ArrayList<SeqLocator>();
		Set<String> reps = new TreeSet<String>();
		for(String[] e : exp){
			locs.add(new SeqLocator(e[0], reps, e[1]));					
		}
		return(execute(locs, readLength, readExtension));
	}
	public RealValuedHistogram execute(List<SeqLocator> locs, double readLength, double readExtension){
		double hittot = 0;
		ArrayList<SeqData> handles = new ArrayList<SeqData>();
		try {
			//Load experiments
            for(SeqLocator l : locs){
				System.err.print(String.format("%s\t", l.getExptName()));
				SeqData curr = new SeqData(gen, l);
			    curr.setReadLength(readLength);
                curr.setReadExtension(readExtension);
				handles.add(curr);
				hittot += curr.getHitCount();	
			}System.err.print(String.format("%.0f reads loaded\n", hittot));
	            
            //Iterate through peaks
            for(Point p : peaks) { 
            	//System.out.println(p.toString());
            	String c = p.getChrom();
				int st = p.getLocation();
				Region currRegion = new Region(gen, c, st-windowSize, st+windowSize);
				LinkedList<StrandedRegion> ipHits = new LinkedList<StrandedRegion>();
				for(SeqData IP: handles){
					ipHits.addAll(IP.loadExtendedHits(currRegion));									
				}
				for(StrandedRegion r : ipHits){
					histo.addValueRange(r.getStart()-p.getLocation(), r.getEnd()-p.getLocation(), 1/hittot);
				}
				
            }	
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return(histo);
		
	}
	public void printMeta(){
		histo.printContents();
	}
}
