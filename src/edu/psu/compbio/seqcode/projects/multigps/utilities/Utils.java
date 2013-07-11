package edu.psu.compbio.seqcode.projects.multigps.utilities;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedPoint;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.verbs.PointParser;
import edu.psu.compbio.seqcode.gse.ewok.verbs.RegionParser;
import edu.psu.compbio.seqcode.gse.ewok.verbs.StrandedPointParser;
import edu.psu.compbio.seqcode.gse.ewok.verbs.StrandedRegionParser;
import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.gse.utils.io.parsing.GFFEntry;
import edu.psu.compbio.seqcode.gse.utils.io.parsing.ParseGFF;

/**
 * Utils: a collection of general methods
 * 
 * @author mahony
 *
 */
public class Utils {


	/**
	 * Loads a set of points from the third or first column of a file
	 * (Suitable for GPS & StatisticalPeakFinder files)
	 * @param filename
	 * @return
	 */
	public static List<Point> loadPointsFromFile(String filename, Genome gen){
		List<Point> points = new ArrayList<Point>();

		try{
			File pFile = new File(filename);
			if(!pFile.isFile()){System.err.println("Invalid file name: "+filename);System.exit(1);}
			BufferedReader reader = new BufferedReader(new FileReader(pFile));
			String line;
			while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            
	            if(words.length>0 && !words[0].contains("#") && !words[0].equals("Region") && !words[0].equals("Position")){
	                if(words.length>=3 && words[2].contains(":")){
		                PointParser pparser = new PointParser(gen);
		            	Point p = pparser.execute(words[2]);
		            	if(p!=null)
		            		points.add(p);		                
	                }else if(words.length>=1 && words[0].contains(":")){
		            	if(words[0].contains("-")){
		                	RegionParser rparser = new RegionParser(gen);
			            	Region q = rparser.execute(words[0]);
			            	if(q!=null)
			            		points.add(q.getMidpoint());			            	
		            	}else{
		            		PointParser pparser = new PointParser(gen);
			            	Point p = pparser.execute(words[0]);
			            	if(p!=null)
			            		points.add(p);
		            	}
		            }
                }
	        }reader.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		System.err.println("Loaded "+points.size()+" points from "+filename);

		return(points);
	}

	/**
	 * Loads a set of points from GFF file
	 * @param filename
	 * @return
	 */
	public static List<Point> loadPointsFromGFFFile(String filename, Genome gen){
		List<Point> points = new ArrayList<Point>();
		try {
			ParseGFF parser = new ParseGFF(new File(filename));
			while(parser.hasNext()){
				GFFEntry site = parser.next();
				Point currPt = new Point(gen, site.getChr(), site.getMidPoint());
				points.add(currPt);
			}
		} catch (IOException e) {
			//Silent exceptions
		}
		
		System.err.println("Loaded "+points.size()+" points from "+filename);

		return(points);
	}

	
	public static List<Pair<Point,Point>> loadIntersFromFile(String filename, Genome gen) {
		List<Pair<Point,Point>> inters = new ArrayList<Pair<Point,Point>>();
		BufferedReader r;
		String s;
		String[] split = {""};
		try {
			r = new BufferedReader(new FileReader(filename));

			while ((s = r.readLine()) != null) {
				split = s.split("\t");
				Point tmp1 = Point.fromString(gen, split[0]);
				Point tmp2 = Point.fromString(gen, split[1]);
				if (tmp1==null || tmp2==null) {
					System.err.println(s);
				} else {
					inters.add(new Pair<Point,Point>(tmp1,tmp2));
				}
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return inters;
	}

	/**
	 * Loads a set of regions from the third or first column of a file
	 * (Suitable for GPS & StatisticalPeakFinder files
	 * @param filename String
	 * @param win integer width of region to impose (-1 leaves region width alone)
	 * @return
	 */
	public static List<Region> loadRegionsFromFile(String filename, Genome gen, int win){
		List<Region> regs = new ArrayList<Region>();

		try{
			File pFile = new File(filename);
			if(!pFile.isFile()){System.err.println("Invalid file name: "+filename);System.exit(1);}
			BufferedReader reader = new BufferedReader(new FileReader(pFile));
			String line;
			while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            
	            if(words.length>0 && !words[0].contains("#") && !words[0].equals("Region") && !words[0].equals("Position")){
	                if(words.length>=3 && words[2].contains(":")){
		                PointParser pparser = new PointParser(gen);
		            	Point p = pparser.execute(words[2]);
		                if(win==-1 && words[0].contains(":") && words[0].contains("-")){
		                	RegionParser rparser = new RegionParser(gen);
			            	Region q = rparser.execute(words[0]);
			            	regs.add(q);
		                }else{
		                	regs.add(p.expand(win/2));
		                }
	                }else if(words.length>=1 && words[0].contains(":")){
		            	if(words[0].contains("-")){
		                	RegionParser rparser = new RegionParser(gen);
			            	Region q = rparser.execute(words[0]);
			            	if(win==-1){
			                	if(q!=null){regs.add(q);}
			                }else
			                	regs.add(q.getMidpoint().expand(win/2));
		            	}else{
		            		PointParser pparser = new PointParser(gen);
			            	Point p = pparser.execute(words[0]);
			            	regs.add(p.expand(win/2));
		            	}
		            }
                }
	        }reader.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return(regs);
	}

	//Load a set of stranded points from a file (stranded point in first column)
	public static List<StrandedPoint> loadStrandedPointsFromFile(Genome gen, String filename){
		List<StrandedPoint> points = new ArrayList<StrandedPoint>();
		try{
			File pFile = new File(filename);
			if(!pFile.isFile()){System.err.println("Invalid file name: "+filename);System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(pFile));
	        String line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            char strand = '+';
	            if(words.length>0 && !words[0].contains("#") && !words[0].equals("Region") && !words[0].equals("Position")){
	                if(words.length>=1 && words[0].contains(":")){
		            	if(words[0].split(":").length>2){
		            		StrandedPointParser pparser = new StrandedPointParser(gen);
		            		StrandedPoint sq = pparser.execute(words[0]);
			            	if(sq!=null){points.add(sq);}
			                
		            	}else{
		            		PointParser pparser = new PointParser(gen);
			            	Point p = pparser.execute(words[0]);
			            	StrandedPoint sp = new StrandedPoint(p, strand);
			            	points.add(sp);
		            	}
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
}
