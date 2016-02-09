package edu.psu.compbio.seqcode.projects.akshay.SeqUnwinder.loadData;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.utils.io.RegionFileUtilities;
import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;

public class KmerCounter {
	
	protected GenomeConfig gcon;
	/** List of ChIP-Seq peaks */
	protected List<Point> peaks;
	/** ChIP-Seq peaks expanded into regions */
	protected List<Region> regs;
	/** Binding site labels (multiple labels are seperated by a ";")*/
	protected List<String> labels;
	/** File name to write the output */
	protected String outfilename = "tmpCounts.mat";
	
	SequenceGenerator<Region> seqgen = null;
	
	@SuppressWarnings("unchecked")
	public KmerCounter(GenomeConfig gc) {
		gcon = gc;
		seqgen = gcon.getSequenceGenerator();
	}
	
	public void setPeaks(List<Point> ps){peaks  = ps;}
	public void setRegions(List<Region> rs){regs = rs;}
	
	
	// Print the the kmer counts (for a given range of value k) in the sequences
	// for each peak
	public void printPeakSeqKmerRange(int kmin, int kmax) throws IOException {
		
		int numK = 0;
		for (int k = kmin; k <= kmax; k++) {
			numK = numK + (int) Math.pow(4, k);
		}
		int[] kmerCounts = new int[numK];
		
		FileWriter of = new FileWriter(outfilename);
		BufferedWriter bw = new BufferedWriter(of);

		// Printing the header line
		//System.out.print("Region");
		for (int k = kmin; k <= kmax; k++) {
			int N = (int) Math.pow(4, k);
			for (int i = 0; i < N; i++)
				bw.write("\t" + RegionFileUtilities.int2seq(i, k));
		}
		bw.write("Label"+"\n");
		

		//for (Region r : regs) {
		for(int r=0; r<regs.size(); r++){
			for (int i = 0; i < numK; i++)
				kmerCounts[i] = 0;

			String seq = seqgen.execute(regs.get(r)).toUpperCase();
			// Check if the sequence (seq) contains any N's if present ignore
			// them
			if (seq.contains("N"))
				continue;

			int ind = 0;
			for (int k = kmin; k <= kmax; k++) {
				for (int i = 0; i < (seq.length() - k + 1); i++) {
					String currK = seq.substring(i, i + k);
					String revCurrK = SequenceUtils.reverseComplement(currK);
					int currKInt = RegionFileUtilities.seq2int(currK);
					int revCurrKInt = RegionFileUtilities.seq2int(revCurrK);
					int kmer = currKInt < revCurrKInt ? currKInt : revCurrKInt;
					kmerCounts[ind + kmer]++;
				}
				ind = ind + (int) Math.pow(4, k);
			}
			//System.out.print(r.getLocationString());
			for (int i = 0; i < numK; i++)
				bw.write("\t" + kmerCounts[i]);
			bw.write(labels.get(r)+"\n");
		}
		bw.close();
	}
	

}
