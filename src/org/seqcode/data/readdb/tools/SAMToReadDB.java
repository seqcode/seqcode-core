package org.seqcode.data.readdb.tools;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.cli.*;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;

/**
 * Reads SAM or BAM data on stdin. Produces a file on stdout in the format
 * expected by ImportHits. The weight for a hit is 1/(# of hits for that read)
 * Ignores secondary & supplementary (i.e. chimeric) alignments.
 * 
 * Options: --uniquehits (flag to only print 1:1 read to hit mappings)
 * --pairedend (flag to print pairs) --junctions (flag to print junction mapping
 * reads as pairs) --read1 (flag to print only read 1 hits) --read2 (flag to
 * print only read 2 hits)
 */

public class SAMToReadDB {

	public static boolean uniqueOnly;
	public static boolean inclPairedEnd;
	public static boolean inclJunction;
	public static boolean read1, read2;

	public static void main(String args[]) throws IOException, ParseException {
		Options options = new Options();
		options.addOption("u", "uniquehits", false, "only output hits with a single mapping");
		options.addOption("p", "pairedend", false, "output paired-end hits");
		options.addOption("j", "junctions", false, "output junction mapping reads (reads with a single gap)");
		options.addOption("1", "read1", false, "output only read 1 hits");
		options.addOption("2", "read2", false, "output only read 2 hits");
		CommandLineParser parser = new GnuParser();
		CommandLine cl = parser.parse(options, args, false);
		uniqueOnly = cl.hasOption("uniquehits");
		inclPairedEnd = cl.hasOption("pairedend");
		inclJunction = cl.hasOption("junctions");
		read1 = cl.hasOption("read1");
		read2 = cl.hasOption("read2");
		SamReaderFactory factory = SamReaderFactory.makeDefault()
				.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS,
						SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
				.validationStringency(ValidationStringency.SILENT);
		SamReader reader = factory.open(SamInputResource.of(System.in));

		CloseableIterator<SAMRecord> iter = reader.iterator();
		List<SAMRecord> records = new ArrayList<SAMRecord>();
		String lastName = "";
		while (iter.hasNext()) { // Group neighboring reads by name
			SAMRecord record = iter.next();
			if (record.getReadUnmappedFlag()) {
				continue;
			}
			if (record.getSupplementaryAlignmentFlag()) {
				continue;
			}
			if (!record.getReadName().equals(lastName)) {
				if (records.size() > 0) {
					processRecord(records);
					records.clear();
				}
			}
			records.add(record);
			lastName = record.getReadName();
		}
		if (records.size() > 0)
			processRecord(records);

		iter.close();
		reader.close();
	}

	public static void processRecord(List<SAMRecord> records) {

		boolean currUnique = false;
		int lcount = 0, rcount = 0; // Have to figure out something for BWA when
									// reporting multiple alignments
		for (SAMRecord record : records) // get weights for L & R reads
											// separately
			if (!record.getReadPairedFlag() || record.getFirstOfPairFlag())
				lcount++;
			else
				rcount++;

		for (SAMRecord record : records) {
			int count = lcount;
			if (record.getReadPairedFlag() && record.getSecondOfPairFlag())
				count = rcount;

			if (record.getIntegerAttribute("NH") != null) // This tag overrides
															// the record count
				count = record.getIntegerAttribute("NH");
			if (count == 1 && record.getMappingQuality() != 0) // Second clause
																// for BWA
				currUnique = true;
			float weight = 1 / (float) count; // Fix this if using to produce
												// multiple mappings for each
												// read

			if (inclPairedEnd || inclJunction) {
				/*
				 * Only accept proper, congruent pairs. It also assumes that the
				 * left and right mates have the same length, and that there are
				 * no gaps in the second mate alignment (SAM doesn't store the
				 * paired read's end) Note: if you change this, you may have to
				 * change the SAMStats output also
				 */
				if (inclPairedEnd) {
					// May need to revisit this section if laoding multiple
					// mapping pairs
					if (record.getReadPairedFlag() && record.getFirstOfPairFlag() && record.getProperPairFlag()) {
						if (!uniqueOnly || currUnique) {
							// Print
							boolean neg = record.getReadNegativeStrandFlag();
							boolean mateneg = record.getMateNegativeStrandFlag();
							String len = record.getReadLength() + "\t";
							System.out.println(record.getReferenceName() + "\t"
									+ (neg ? record.getAlignmentEnd() : record.getAlignmentStart()) + "\t"
									+ (neg ? "-\t" : "+\t") + len + record.getMateReferenceName() + "\t"
									+ (mateneg ? record.getMateAlignmentStart() + record.getReadLength() - 1
											: record.getMateAlignmentStart())
									+ "\t" + (mateneg ? "-\t" : "+\t") + len + weight + "\t" + 1);
						}
					}
				}

				/*
				 * Outputs as paired alignments those reads that are aligned in
				 * >2 blocks Note: if you change this, you may have to change
				 * the SAMStats output also
				 */
				if (inclJunction) {
					if (!uniqueOnly || currUnique) {
						List<AlignmentBlock> blocks = record.getAlignmentBlocks();
						if (blocks.size() >= 2) {
							for (int ab = 0; ab < blocks.size() - 1; ab++) {
								AlignmentBlock lBlock = blocks.get(ab);
								int lStart = lBlock.getReferenceStart();
								int lEnd = lStart + lBlock.getLength() - 1;
								int lLen = lBlock.getLength();
								AlignmentBlock rBlock = blocks.get(ab + 1);
								int rStart = rBlock.getReferenceStart();
								int rEnd = rStart + rBlock.getLength() - 1;
								int rLen = rBlock.getLength();
								boolean neg = record.getReadNegativeStrandFlag();
								String refname = record.getReferenceName() + "\t";
								System.out.println(refname + (neg ? lEnd : lStart) + "\t" + (neg ? "-\t" : "+\t") + lLen
										+ "\t" + refname + (neg ? rEnd : rStart) + "\t" + (neg ? "-\t" : "+\t") + rLen
										+ "\t" + weight + "\t" + 0);
							}
						}
					}
				}
			} else { // Just output reads (ignore alignment blocks for now)

				if (uniqueOnly && !currUnique) {
					return;
				}
				if ((!read1 && !read2) || (!record.getReadPairedFlag()) || (read1 && record.getFirstOfPairFlag())
						|| (read2 && record.getSecondOfPairFlag())) {
					System.out.println(String.format("%s\t%d\t%s\t%d\t%f", record.getReferenceName(),
							record.getReadNegativeStrandFlag() ? record.getAlignmentEnd() : record.getAlignmentStart(),
							record.getReadNegativeStrandFlag() ? "-" : "+", record.getReadLength(), weight));
				}
			}
		}
	}
}