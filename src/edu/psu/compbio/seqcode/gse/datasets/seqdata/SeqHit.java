package edu.psu.compbio.seqcode.gse.datasets.seqdata;

import edu.psu.compbio.seqcode.gse.datasets.general.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;

public class SeqHit extends StrandedRegion {

  private double weight = 1.0;

  public SeqHit(Genome g, String chrom, int start, int end, char strand, double weight) {
    super(g,chrom,start,end,strand);
    this.weight = weight;
  }

    public double getWeight() {
    return weight;
  }


  public void setWeight(double weight) {
    this.weight = weight;
  }

  public String toString() { 
    return String.format("%s:%c", getLocationString(), getStrand());
  }

  public SeqHit extendHit(int ext) { 
    if(getStrand() == '+') { 
      return new SeqHit(getGenome(), getChrom(), getStart(), getEnd() + ext, getStrand(), weight);
    } else { 
      return new SeqHit(getGenome(), getChrom(), getStart() - ext, getEnd(), getStrand(), weight);
    }
  }
  public SeqHit fivePrime() { 
    if(getStrand() == '+') { 
	      return new SeqHit(getGenome(), getChrom(), getStart(), getStart(), getStrand(), weight);
	    } else { 
	      return new SeqHit(getGenome(), getChrom(), getEnd(), getEnd(), getStrand(), weight);
	    }
  }
  
  public SeqHit shiftExtendHit(int ext, int shift) { 
    if(getStrand() == '+') { 
      return new SeqHit(getGenome(), getChrom(), getStart()+shift-(ext/2), getEnd()+shift+(ext/2), getStrand(), weight);
    } else { 
      return new SeqHit(getGenome(), getChrom(), getStart()-shift-(ext/2), getEnd()-shift+(ext/2), getStrand(), weight);
    }
  }

  public int hashCode() { 
    int code = 17;
    code += super.hashCode(); code *= 37;
    return code; 
  }

  public boolean equals(Object o) { 
      return this == o;
  }  
}
