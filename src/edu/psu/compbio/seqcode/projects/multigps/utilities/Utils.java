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
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.verbs.PointParser;
import edu.psu.compbio.seqcode.gse.ewok.verbs.RegionParser;
import edu.psu.compbio.seqcode.gse.utils.Pair;

/**
 * Utils: a collection of general methods
 * 
 * @author mahony
 *
 */
public class Utils {


	/**
	 * Loads a set of points from the third or first column of a file
	 * (Suitable for GPS & StatisticalPeakFinder files
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

				if(words.length>=3 && words[2].contains(":")){
					PointParser pparser = new PointParser(gen);
					Point p = pparser.execute(words[2]);
					points.add(p);
				}else if(words.length>=1 && words[0].contains(":")){
					if(words[0].contains("-")){
						RegionParser parser = new RegionParser(gen);
						Region q = parser.execute(words[0]);
						points.add(q.getMidpoint());
					}else{
						PointParser pparser = new PointParser(gen);
						Point p = pparser.execute(words[0]);
						points.add(p);
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

				if(words.length>=3 && win!=-1 && words[2].contains(":")){
					PointParser pparser = new PointParser(gen);
					Point p = pparser.execute(words[2]);
					int rstart = p.getLocation()-(win/2)<1 ? 1:p.getLocation()-(win/2);
					int rend = p.getLocation()+(win/2)>gen.getChromLength(p.getChrom()) ? gen.getChromLength(p.getChrom()):p.getLocation()+(win/2)-1;
					Region r = new Region(p.getGenome(), p.getChrom(), rstart, rend);
					regs.add(r);
				}else if(words.length>=1 && words[0].contains(":")){
					if(words[0].contains("-")){
						RegionParser parser = new RegionParser(gen);
						Region q = parser.execute(words[0]);
						if(win!=-1){
							int rstart = q.getMidpoint().getLocation()-(win/2)<1 ? 1:q.getMidpoint().getLocation()-(win/2);
							int rend = q.getMidpoint().getLocation()+(win/2)>gen.getChromLength(q.getChrom()) ? gen.getChromLength(q.getChrom()):q.getMidpoint().getLocation()+(win/2)-1;
							Region r = new Region(q.getGenome(), q.getChrom(), rstart, rend);
							if(r!=null){regs.add(r);}
						}else{
							if(q!=null){regs.add(q);}
						}
					}else{
						if(win==-1)
							win=200; 
						PointParser pparser = new PointParser(gen);
						Point p = pparser.execute(words[0]);
						int rstart = p.getLocation()-(win/2)<1 ? 1:p.getLocation()-(win/2);
						int rend = p.getLocation()+(win/2)>gen.getChromLength(p.getChrom()) ? gen.getChromLength(p.getChrom()):p.getLocation()+(win/2)-1;
						Region r = new Region(p.getGenome(), p.getChrom(), rstart, rend);
						regs.add(r);
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

}
