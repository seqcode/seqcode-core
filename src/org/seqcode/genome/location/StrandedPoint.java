package org.seqcode.genome.location;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.regex.Matcher;

import org.seqcode.genome.Genome;


public class StrandedPoint extends Point implements Stranded {
	
	private char strand;
	
	public int compareToStranded(StrandedPoint p) {
		int compare = this.compareTo(p);
		if (compare==0) {
			return this.strand - p.strand;
		} else {
			return compare;
		}
	}

	public StrandedPoint(Genome g, String chrom, int start, char str) {
		super(g, chrom, start);
		strand = str;
	}
	
	public StrandedPoint(Point p, char str) {
		this(p.getGenome(), p.getChrom(), p.getLocation(), str);
	}

	public char getStrand() {
		return strand;
	}

    public int hashCode() { 
        int code = super.hashCode();
        code += strand; code *= 37;
        return code;
    }

    public boolean equals(Object o) { 
        if(!(o instanceof StrandedPoint)) { return false; }
        StrandedPoint sp = (StrandedPoint)o;
        if(strand != sp.strand) { return false; }
        return super.equals(sp);
    }
    
    public String getLocationString() {
    	if (strand==' ') {
    		return super.getLocationString();
    	} else {
    		return super.getLocationString()+":"+strand;
    	}
    }
    
    public StrandedRegion expand(int upstream, int downstream) {
        if (strand == '+') {
            int ns = getLocation() - upstream;
            int ne = getLocation() + downstream;
            if (ns < 1) {ns = 1;}
            return new StrandedRegion(getGenome(),getChrom(),ns,ne,strand);
        } else if (strand == '-') {
            int ns = getLocation() - downstream;
            int ne = getLocation() + upstream;
            if (ns < 1) {ns = 1;}
            return new StrandedRegion(getGenome(),getChrom(),ns,ne,strand);                
        } else {
            throw new IllegalArgumentException("Strand isn't + or - so I don't know what to do");
        }

    }
    
    /**
     * parses the input String into a StrandedPoint. Understands abbreviates in the
     * coordinates such as k and m. <br>
     * The method accepts the form: <blockquote>
     * 
     * <pre>
     * chromosome:start:strand
     * </pre>
     * 
     * </blockquote>
     */
    public static StrandedPoint fromString(Genome genome, String input) throws NumberFormatException {
      String pieces[] = input.split(":");
      char strand = ' ';
      if (pieces.length == 3 && pieces[2].length() == 1) {
        strand = pieces[2].charAt(0);
        input = pieces[0] + ":" + pieces[1];
      }

      String trimmed = input.trim();
      
      Matcher pointmatcher = POINT_PATTERN.matcher(trimmed);

      StrandedPoint output = null;
      
      if (pointmatcher.find()) {
        if (pointmatcher.groupCount() != 2) {
          return null;
        }
        String chromStr = pointmatcher.group(1);
        String locStr = pointmatcher.group(2);
        if (chromStr.startsWith("chr")) {
          chromStr = chromStr.substring(3, chromStr.length());
        }
        locStr = locStr.replaceAll(",", "");
        output = new StrandedPoint(genome, chromStr, Region.stringToNum(locStr), strand);
      }

      return output;
    }
    /**
     * load from file a list of StrandedPoints. 
     * Understands abbreviates in the coordinates such as k and m. <br>
     * The method accepts the form: <blockquote>
     * 
     * <pre>
     * chromosome:start:strand
     * </pre>
     * 
     * </blockquote>
     */
    public static ArrayList<StrandedPoint> fromFile(Genome genome, String fileName)  {
    	ArrayList<String> strs = new ArrayList<String>();
		try {	
			BufferedReader bin = new BufferedReader(new InputStreamReader(new FileInputStream(new File(fileName))));
	        String line;
	        while((line = bin.readLine()) != null) { 
	            line = line.trim();
	            if (line.length()!=0)
	            	strs.add(line);
	        }			
	        if (bin != null) {
	            bin.close();
	        }
        } catch (IOException e) {
        	System.err.println("Error when processing "+fileName);
            e.printStackTrace(System.err);
        }   
		
      ArrayList<StrandedPoint> points = new ArrayList<StrandedPoint>();
      for (String s:strs){
    	  StrandedPoint p=null;
    	  try{
    		  p = StrandedPoint.fromString(genome, s);
    	  }
    	  catch(NumberFormatException e){
    		  e.printStackTrace();
    		  continue;
    	  }
    	  if (p==null)
    		  continue;
    	  points.add(p);
      }
      return points;
    }
}
