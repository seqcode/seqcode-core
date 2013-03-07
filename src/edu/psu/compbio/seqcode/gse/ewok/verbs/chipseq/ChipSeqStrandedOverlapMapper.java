package edu.psu.compbio.seqcode.gse.ewok.verbs.chipseq;

import java.sql.*;
import java.util.*;
import java.io.*;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.general.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.*;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Mapper;
import edu.psu.compbio.seqcode.gse.utils.Closeable;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.utils.database.*;

public class ChipSeqStrandedOverlapMapper implements Closeable, Mapper<Region, WeightedRunningOverlapSum[]> {

  public static final int TOTAL_SUM_INDEX = 0;
  public static final int POS_SUM_INDEX = 1;
  public static final int NEG_SUM_INDEX = 2;
  
  private SeqDataLoader loader;
  private SeqLocator locator;
  private Genome lastGenome;
  private LinkedList<SeqAlignment> alignments;

  private int extension;
  private int shift = 0;
  private boolean shifting = false;


  public ChipSeqStrandedOverlapMapper(SeqLocator loc, int extension)  throws SQLException, IOException {
    this.extension = extension;
    loader = new SeqDataLoader();
    locator = loc;
    alignments = null;
    lastGenome = null;
  }
  private void getAligns(Genome genome) throws SQLException {
      if (alignments != null && genome.equals(lastGenome)) {
          return;
      }
      alignments = new LinkedList<SeqAlignment>();
      lastGenome = genome;
      try {
          alignments.addAll(locator.loadAlignments(loader, genome));
      } catch (SQLException e) {
            e.printStackTrace(System.err);
      } catch (NotFoundException e) {
          e.printStackTrace();
      }
  }

  public WeightedRunningOverlapSum[] execute(Region a) {
    try {
      Genome g = a.getGenome();
      try {
          getAligns(g);
      } catch (SQLException e) {
          throw new DatabaseException(e.toString(),e);
      }
      String chrom = a.getChrom();
      WeightedRunningOverlapSum[] sums = { new WeightedRunningOverlapSum(g, chrom), 
          new WeightedRunningOverlapSum(g, chrom), 
          new WeightedRunningOverlapSum(g, chrom) };
      
      Collection<SeqHit> hits = loader.loadByRegion(alignments, a);
      for (SeqHit hit : hits) {
        if (hit.getStrand() == '+') {
          if (shifting) {
            int intervalStart = hit.getStart() + shift - (extension / 2);
            int intervalEnd = hit.getEnd() + shift + (extension / 2); 
            sums[POS_SUM_INDEX].addWeightedInterval(intervalStart, intervalEnd, hit.getWeight());
            sums[TOTAL_SUM_INDEX].addWeightedInterval(intervalStart, intervalEnd, hit.getWeight());            
          }
          else {
            sums[POS_SUM_INDEX].addWeightedInterval(hit.getStart(), hit.getEnd() + extension, hit.getWeight());
            sums[TOTAL_SUM_INDEX].addWeightedInterval(hit.getStart(), hit.getEnd() + extension, hit.getWeight());
          }
        }
        else {
          if (shifting) {
            int intervalStart = hit.getStart() - shift - (extension / 2);
            int intervalEnd = hit.getEnd() - shift + (extension / 2);
            sums[NEG_SUM_INDEX].addWeightedInterval(intervalStart, intervalEnd, hit.getWeight());
            sums[TOTAL_SUM_INDEX].addWeightedInterval(intervalStart, intervalEnd, hit.getWeight());
          }
          else {
            sums[NEG_SUM_INDEX].addWeightedInterval(hit.getStart() - extension, hit.getEnd(), hit.getWeight());
            sums[TOTAL_SUM_INDEX].addWeightedInterval(hit.getStart() - extension, hit.getEnd(), hit.getWeight());
          }
        }
      }
      return sums;
    }
    catch (Exception sqlex) {
      throw new DatabaseException(sqlex.toString(), sqlex);
    }
  }


  public void setExtension(int e) {
    extension = e;
  }


  public void setShift(int s) {
    shift = s;
    if (s > 0)
      shifting = true;
    else shifting = false;
  }


  public void close() {
    if (loader != null) {
      loader.close();
      loader = null;
      alignments.clear();
    }
  }


  public boolean isClosed() {
    return loader == null;
  }
}