package org.seqcode.data.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import org.seqcode.data.io.parsing.BEDLine;
import org.seqcode.data.io.parsing.BEDParser;
import org.seqcode.data.io.parsing.GFFEntry;
import org.seqcode.data.io.parsing.ParseGFF;
import org.seqcode.genome.Genome;
import org.seqcode.genome.location.NamedRegion;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedPoint;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.genome.sequence.SequenceGenerator;
import org.seqcode.genome.sequence.SequenceUtils;
import org.seqcode.gsebricks.verbs.location.ChromRegionIterator;
import org.seqcode.gsebricks.verbs.location.PointParser;
import org.seqcode.gsebricks.verbs.location.RegionParser;
import org.seqcode.gsebricks.verbs.location.StrandedPointParser;
import org.seqcode.gsebricks.verbs.location.StrandedRegionParser;
import org.seqcode.utils.Pair;


public class RegionFileUtilities {

	public RegionFileUtilities(){}
	
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
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return(points);
	}
	/**
	 * Load a set of regions from a peak file
	 * @param gen
	 * @param filename
	 * @param win
	 * @return
	 */
	public static List<Region> loadRegionsFromPeakFile(Genome gen, String filename, int win){
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
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return(regs);
	}
	
	/**
	 * Load Stranded points from stranded peak file, Usually used to load TSS
	 * @param gen
	 * @param filename
	 * @return
	 */
	public static List<StrandedPoint> loadStrandedPointFromRefTssFile(Genome gen, String filename){
		List<StrandedPoint> pts = new ArrayList<StrandedPoint>();
		try{
			File pFile = new File(filename);
			if(!pFile.isFile()){System.err.println("Invalid file name: "+filename);System.exit(1);}
			BufferedReader reader = new BufferedReader(new FileReader(pFile));
			String line;
			while((line = reader.readLine()) != null){
				line=line.trim();
				String[] words = line.split("\\s+");
				if(words.length >=1 && !words[0].contains("#") &&  !words[0].equals("Region") && !words[0].equals("Position")){
					String[] subwords = words[0].split(":");
					PointParser pparser = new PointParser(gen);
					Point p = pparser.execute(subwords[0]+":"+subwords[1]);
					StrandedPoint sp = new StrandedPoint(p,subwords[2].charAt(0));
					pts.add(sp);
				}
			}
			reader.close();
		} catch(FileNotFoundException e){
			e.printStackTrace();
		} catch(IOException e){
			e.printStackTrace();
		}
		
		return pts;
	}
	
	/**
	 * Load a set of stranded regions from a file
	 * @param gen
	 * @param filename
	 * @param win
	 * @return
	 */
	public static List<StrandedRegion> loadStrandedRegionsFromMotifFile(Genome gen, String filename, int win){
		List<StrandedRegion> regs = new ArrayList<StrandedRegion>();
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
	                if(words.length>=3 && words[2].contains(":")){
		                PointParser pparser = new PointParser(gen);
		            	Point p = pparser.execute(words[2]);
		                if(win==-1 && words[0].contains(":") && words[0].contains("-")){
		                	StrandedRegionParser rparser = new StrandedRegionParser(gen);
			            	StrandedRegion sq = rparser.execute(words[0]);
			            	regs.add(sq);
		                }else{
		                	StrandedRegion sp = new StrandedRegion(p.expand(win/2), strand);
		                	regs.add(sp);
		                }
	                }else if(words.length>=1 && words[0].contains(":")){
	                	String[] subwords = words[0].split(":");
		            	if(subwords[1].contains("-")){
		            		StrandedRegionParser rparser = new StrandedRegionParser(gen);
		            		StrandedRegion sq = rparser.execute(words[0]);
			            	if(win==-1){
			                	if(sq!=null){regs.add(sq);}
			                }else
			                	regs.add(new StrandedRegion(sq.getMidpoint().expand(win/2), sq.getStrand()));
		            	}else{
		            		PointParser pparser = new PointParser(gen);
			            	Point p = pparser.execute(subwords[0]+":"+subwords[1]);
			            	StrandedRegion sp=null;
			            	if(subwords.length>=3)
			            		sp = new StrandedRegion(p.expand(win/2), subwords[2].charAt(0));
			            	else
			            		sp = new StrandedRegion(p.expand(win/2), strand);
			            	regs.add(sp);
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
	
	/**
	 * Load a set of stranded regions from a BED file
	 * @param gen
	 * @param filename
	 * @param win
	 * @return
	 */
	public static List<StrandedRegion> loadStrandedRegionsFromBEDFile(Genome gen, String filename, int win){
		List<StrandedRegion> regs = new ArrayList<StrandedRegion>();
		try{
			BEDParser parser = new BEDParser(new File(filename));
			BEDLine line;
	        while (parser.hasNext()) {
	            line = parser.next();
	            StrandedRegion sq = new StrandedRegion(
	            		gen,
	            		line.getChrom(),
	            		line.getChromStart()+1,
	            		line.getChromEnd(),
	            		line.getStrand());
	            if(win==-1)
	            	regs.add(sq);
	            else
	            	regs.add(new StrandedRegion(sq.getMidpoint().expand(win), sq.getStrand()));
	        }
	        parser.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return(regs);
	}
	
	/**
	 * Load a set of stranded points from a file
	 * @param gen
	 * @param filename
	 * @param win
	 * @return
	 */
	public static List<StrandedPoint> loadStrandedPointsFromMotifFile(Genome gen, String filename, int win){
		List<StrandedRegion> regs = new ArrayList<StrandedRegion>();
		List<StrandedPoint> peaks = new ArrayList<StrandedPoint>();
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
	                if(words.length>=3 && words[2].contains(":")){
		                PointParser pparser = new PointParser(gen);
		            	Point p = pparser.execute(words[2]);
		                if(win==-1 && words[0].contains(":") && words[0].contains("-")){
		                	StrandedRegionParser rparser = new StrandedRegionParser(gen);
			            	StrandedRegion sq = rparser.execute(words[0]);
			            	regs.add(sq);
		                }else{
		                	StrandedRegion sp = new StrandedRegion(p.expand(win/2), strand);
		                	regs.add(sp);
		                }
	                }else if(words.length>=1 && words[0].contains(":")){
	                	String[] subwords = words[0].split(":");
		            	if(subwords[1].contains("-")){
		            		StrandedRegionParser rparser = new StrandedRegionParser(gen);
		            		StrandedRegion sq = rparser.execute(words[0]);
			            	if(win==-1){
			                	if(sq!=null){regs.add(sq);}
			                }else
			                	regs.add(sq.expand(win/2, win/2));
		            	}else{
		            		PointParser pparser = new PointParser(gen);
			            	Point p = pparser.execute(subwords[0]+":"+subwords[1]);
			            	StrandedRegion sp=null;
			            	if(subwords.length>=3)
			            		sp = new StrandedRegion(p.expand(win/2), subwords[2].charAt(0));
			            	else
			            		sp = new StrandedRegion(p.expand(win/2), strand);
			            	regs.add(sp);
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
		for(StrandedRegion r : regs)
			peaks.add(new StrandedPoint(r.getGenome(), r.getChrom(), r.getMidpoint().getLocation(), r.getStrand()));
		return(peaks);
	}
	
	/**
	 * Load a set of regions from a peak file
	 * @param gen
	 * @param filename
	 * @param win
	 * @return
	 */
	public static List<Point> loadPeaksFromPeakFile(Genome gen, String filename, int win){
		List<Point> peaks = new ArrayList<Point>();
		try{
			File pFile = new File(filename);
			if(!pFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(pFile));
	        String line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            
	            if(words.length>0 && !words[0].contains("#") && !words[0].equals("Region") && !words[0].equals("Position")){
	                if(words.length>=3 && words[2].contains(":")){
		                PointParser pparser = new PointParser(gen);
		            	Point p = pparser.execute(words[2]);
		            	peaks.add(p);		                
	                }else if(words.length>=1 && words[0].contains(":")){
		            	if(words[0].contains("-")){
		                	RegionParser rparser = new RegionParser(gen);
			            	Region q = rparser.execute(words[0]);
			            	peaks.add(q.getMidpoint());			            	
		            	}else{
		            		PointParser pparser = new PointParser(gen);
			            	Point p = pparser.execute(words[0]);
			            	peaks.add(p);
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
		return(peaks);
	}
	
	public static List<String> loadLinesFromFile(String filename){
		List<String> lines = new ArrayList<String>();
		try{
			File inFile = new File(filename);
			if(!inFile.isFile()){System.err.println("Invalid file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(inFile));
	        String line;//Ignore first line;
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            if(words.length>0 && !words[0].contains("#") && !words[0].equals("Region") && !words[0].equals("Position")){
	            	lines.add(line);
	            }
	        }reader.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        return(lines);
	}
	
	
	/**
	 * Get sequences for a set of regions
	 * @param regions
	 * @param seqgen
	 * @return
	 */
	public static List<String> getSequencesForRegions(List<Region> regions, SequenceGenerator seqgen){
		ArrayList<String> seqs = new ArrayList<String>(); 
		if(seqgen==null)
			seqgen = new SequenceGenerator();
		for(Region r : regions){
			seqs.add(seqgen.execute(r).toUpperCase());
		}return(seqs);
	}
	
	//Get sequences for a set of regions
	public static List<String> getSequencesForStrandedRegions(List<StrandedRegion> regions, SequenceGenerator seqgen){
		ArrayList<String> seqs = new ArrayList<String>(); 
		if(seqgen==null)
			seqgen = new SequenceGenerator();
		for(StrandedRegion r : regions){
			String seq = seqgen.execute(r).toUpperCase();
			if(r.getStrand()=='-')
				seq = SequenceUtils.reverseComplement(seq);
			seqs.add(seq);
		}return(seqs);
	}
	
	//Randomly pick a set of Regions
	public static List<Region> randomRegionPick(Genome gen, List<Region> blackList, int numSamples, int sampleSize){
		List<Region> regs = new ArrayList<Region>();
		Random rand = new Random();
		int validSamples=0;
		
		//First see how big the genome is:
		int numChroms=0;
		long genomeSize=0;
		long [] chromoSize = new long[gen.getChromList().size()];
		String [] chromoNames = new String[gen.getChromList().size()];
		Iterator<NamedRegion> chroms = new ChromRegionIterator(gen);
		while (chroms.hasNext()) {
			NamedRegion currentChrom = chroms.next();
			genomeSize += (double)currentChrom.getWidth();
			chromoSize[numChroms]=currentChrom.getWidth();
			chromoNames[numChroms]=currentChrom.getChrom();
			numChroms++;				
		}

		//Now, iteratively generate random positions and check if they are valid and not overlapping repeats. 
		while(validSamples<numSamples){
			Region potential;				
			long randPos = (long)(1+(rand.nextDouble()*genomeSize));
			//find the chr
			boolean found=false;
			long total=0;
			for(int c=0; c<numChroms && !found; c++){
				if(randPos<total+chromoSize[c]){
					found=true;
					if(randPos+sampleSize<total+chromoSize[c]){
						potential = new Region(gen, chromoNames[c], (int)(randPos-total), (int)(randPos+sampleSize-total));
						
						//is this region in the blacklist? 
						boolean valid=true;
						if(blackList!=null){
							for(Region r : blackList){
								if(potential.overlaps(r)){valid=false;}
							}
						}
						if(valid){
							validSamples++;
							regs.add(potential);
						}
					}
				}total+=chromoSize[c];
			}
		}
		return(regs);
	}
	/**
	 * Regions to midpoints
	 * @param regs
	 * @return
	 */
	public static List<Point> regions2midpoints(List<Region> regs){
		List<Point> p = new ArrayList<Point>();
		for(Region r : regs)
			p.add(r.getMidpoint());
		return p;
	}
	
	
	/**
	   * Convert a base to an int value
	   * 
	   * @param base
	   * @return
	   */
	  public static int base2int(char base) {
	    int intVal = -1;
	    switch (base) {
	      case 'A':
	        intVal = 0;
	        break;
	      case 'C':
	        intVal = 1;
	        break;
	      case 'G':
	        intVal = 2;
	        break;
	      case 'T':
	        intVal = 3;
	        break;
	      default:
	        throw new IllegalArgumentException("Invalid character: " + base);
	    }
	    return intVal;
	  }


	  /**
	   * Return a base for the specified integer
	   * 
	   * @param x
	   * @return
	   */
	  public static char int2base(int x) {
	    char base;
	    switch (x) {
	      case 0:
	        base = 'A';
	        break;
	      case 1:
	        base = 'C';
	        break;
	      case 2:
	        base = 'G';
	        break;
	      case 3:
	        base = 'T';
	        break;
	      default:
	        throw new IllegalArgumentException("Invalid int: " + x);
	    }
	    return (base);
	  }


	  /**
	   * Convert a nucleotide sequence to an integer value
	   * 
	   * @param seq
	   * @return
	   */
	  public static int seq2int(String seq) {
	    int intVal = 0;
	    int len = seq.length();

	    for (int i = 0; i < len; i++) {
	      long currInt = base2int(seq.charAt(i));
	      if (currInt == -1) {
	        return -1;
	      }
	      intVal = intVal << 2;
	      intVal += currInt;
	    }
	    return intVal;
	  }


	  /**
	   * 
	   * @param x
	   * @return
	   */
	  public static String int2seq(long x, int kmerLen) {
	    /**
	     * check that the x is valid for the specified maxKmerLen. Note: 4 << (2 *
	     * (kmerLen - 1)) = 4^kmerLen
	     */
	    if (x > ((4 << (2 * (kmerLen - 1))) - 1)) {
	      throw new IllegalArgumentException("Invalid int value, " + x + ", for kmerLen " + kmerLen);
	    }
	    StringBuffer seq = new StringBuffer(kmerLen);
	    for (int i = 0; i < kmerLen; i++) {
	      int baseVal = (int) (x % 4);
	      seq.append(int2base(baseVal));
	      x = x >> 2;
	    }
	    return seq.reverse().toString();
	  }

	  /**
	   * Get all possible k-mers of a given length k
	   * @param kmerLen
	   * @return
	   */
	  public static List<String> getAllKmers(int kmerLen){
		  List<String> kmers = new ArrayList<String>();
		  int numK = (int) Math.pow(4, kmerLen);
		  for(int k=0; k<numK; k++)
			  kmers.add(int2seq(k, kmerLen));
		  return kmers;
	  }
}
