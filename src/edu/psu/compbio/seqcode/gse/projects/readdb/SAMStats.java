package edu.psu.compbio.seqcode.gse.projects.readdb;

import java.io.IOException;
import java.util.List;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloseableIterator;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import edu.psu.compbio.seqcode.gse.utils.RealValuedHistogram;

public class SAMStats {

	private Double totalReads = 0.0, totalHits=0.0, LHits=0.0, RHits=0.0; 
	private Double singleEnd=0.0, properPair=0.0, pairMapped=0.0, notPrimary=0.0;
	private Double junctions=0.0;
	private Double uniquelyMapped=0.0;
	private Double weight=0.0;
	private Integer readLen =-1;
	private Double pairedEndSameChr=0.0, pairedEndDiffChr=0.0;
	private Boolean bowtie1=false, bowtie2=false;
	private RealValuedHistogram histo;
	
	public static void main(String args[]) throws IOException, ParseException {
        Options options = new Options();
        options.addOption("i","insdist",false,"print insert size distribution");
        options.addOption("s","stats",false,"print mapping stats");
        options.addOption("l","readlen",false,"read length");
        options.addOption("c","readcount",false,"read count");
        options.addOption("bt2",false,"input is from bowtie2");
        options.addOption("bt1",false,"input is from bowtie1 (tested for --best --strata -m 1)");
        CommandLineParser parser = new GnuParser();
        CommandLine cl = parser.parse( options, args, false );       
        
        SAMStats s = new SAMStats(cl.hasOption("bt1"), cl.hasOption("bt2"));
        if(cl.hasOption("insdist"))
        	s.printInsertDistrib();
        if(cl.hasOption("stats"))
        	s.printStats();
        if(cl.hasOption("readlen"))
        	s.printReadLength();
        if(cl.hasOption("readcount"))
        	s.printReadCount();
    }
	
	public SAMStats(boolean bowtie1, boolean bowtie2){
		this.bowtie1 = bowtie1;
		this.bowtie2 = bowtie2;
		histo = new RealValuedHistogram(0, 1000, 100);
		SAMFileReader reader = new SAMFileReader(System.in);
		reader.setValidationStringency(ValidationStringency.SILENT);
        CloseableIterator<SAMRecord> iter = reader.iterator();
        while (iter.hasNext()) {
            SAMRecord record = iter.next();
            if(readLen==-1)
            	readLen = record.getReadLength();
            if(bowtie1)
            	processBT1SAMRecord(record);
            else if(bowtie2)
            	processBT2SAMRecord(record);
            else
            	processSAMRecord(record);
        }
        iter.close();
        reader.close();        
    }
	
	public void processSAMRecord(SAMRecord r){
		totalReads++;
		if(!r.getReadUnmappedFlag()){
			totalHits++;
			int count = 1;  //Have to figure out something for BWA when reporting multiple alignments
			if(r.getIntegerAttribute("NH")!=null)
				count = r.getIntegerAttribute("NH");
			if(count==1 && r.getMappingQuality()!=0) //Second clause for BWA
				uniquelyMapped++;
			
			weight += 1/(float)count;

			if(r.getReadPairedFlag()){
				if(r.getMateUnmappedFlag()){
					singleEnd++;
				}else{
					pairMapped++;
					if(r.getMateReferenceName().equals(r.getReferenceName())){
						pairedEndSameChr++;
					}else{
						pairedEndDiffChr++;
					}
				}
			}else{
				singleEnd++;
			}
			
			List<AlignmentBlock> blocks = r.getAlignmentBlocks();
    		if(blocks.size()>=2){
    			junctions+=blocks.size()-1;
    		}
			
			if(!r.getNotPrimaryAlignmentFlag()){
				if(!r.getReadPairedFlag() || r.getFirstOfPairFlag()){
					LHits++;
					if(r.getReadPairedFlag() && r.getProperPairFlag()){
						properPair++;
						if(!r.getReadNegativeStrandFlag() && r.getMateNegativeStrandFlag()){
							double dist = (r.getMateAlignmentStart()+r.getReadLength())-r.getAlignmentStart();
							histo.addValue(dist);
						}
					}
				}else if(r.getSecondOfPairFlag()){
					RHits++;
				}
			}else{
				notPrimary++;
			}
		}
	}
	
	public void processBT1SAMRecord(SAMRecord r){
		totalReads++;
		if(!r.getReadUnmappedFlag())
			if(r.getIntegerAttribute("XM")!=null){
				int xm = r.getIntegerAttribute("XM");
				if(xm!=0)
					weight++;
			}
		else{
			totalHits++;
			int count =1; //TODO: Fix this if using bowtie for multi-mapping reads
			boolean currUnique = true;
	    	
			if(count==1 && currUnique){
				uniquelyMapped++;
			}
			
			weight += 1/(float)count;
			
			if(r.getReadPairedFlag()){
				if(r.getMateUnmappedFlag()){
					singleEnd++;
				}else{
					pairMapped++;
					if(r.getMateReferenceName().equals(r.getReferenceName())){
						pairedEndSameChr++;
					}else{
						pairedEndDiffChr++;
					}
				}
			}else{
				singleEnd++;
			}
			
			List<AlignmentBlock> blocks = r.getAlignmentBlocks();
    		if(blocks.size()>=2){
    			junctions+=blocks.size()-1;
    		}
    		
			if(!r.getNotPrimaryAlignmentFlag()){
				if(!r.getReadPairedFlag() || r.getFirstOfPairFlag()){
					LHits++;
					if(r.getReadPairedFlag() && r.getProperPairFlag()){
						properPair++;
						if(!r.getReadNegativeStrandFlag() && r.getMateNegativeStrandFlag()){
							double dist = (r.getMateAlignmentStart()+r.getReadLength())-r.getAlignmentStart();
							histo.addValue(dist);
						}
					}
				}else if(r.getSecondOfPairFlag()){
					RHits++;
				}
			}else{
				notPrimary++;
			}
		}
	}
	
	public void processBT2SAMRecord(SAMRecord r){
		totalReads++;
		if(!r.getReadUnmappedFlag()){
			totalHits++;
			int count =1; //TODO: Fix this if using bowtie2 for multi-mapping reads
			
			int primAScore = r.getIntegerAttribute("AS");
	    	int secAScore=-1000000;
	    	if(r.getIntegerAttribute("XS")!=null)
	    		secAScore = r.getIntegerAttribute("XS");
	    	boolean currUnique = primAScore > secAScore ? true : false;
	    	
			if(count==1 && currUnique){
				uniquelyMapped++;
			}
			
			weight += 1/(float)count;
			
			if(r.getReadPairedFlag()){
				if(r.getMateUnmappedFlag()){
					singleEnd++;
				}else{
					pairMapped++;
					if(r.getMateReferenceName().equals(r.getReferenceName())){
						pairedEndSameChr++;
					}else{
						pairedEndDiffChr++;
					}
				}
			}else{
				singleEnd++;
			}
			
			List<AlignmentBlock> blocks = r.getAlignmentBlocks();
    		if(blocks.size()>=2){
    			junctions+=blocks.size()-1;
    		}
    		
			if(!r.getNotPrimaryAlignmentFlag()){
				if(!r.getReadPairedFlag() || r.getFirstOfPairFlag()){
					LHits++;
					if(r.getReadPairedFlag() && r.getProperPairFlag()){
						properPair++;
						if(!r.getReadNegativeStrandFlag() && r.getMateNegativeStrandFlag()){
							double dist = (r.getMateAlignmentStart()+r.getReadLength())-r.getAlignmentStart();
							histo.addValue(dist);
						}
					}
				}else if(r.getSecondOfPairFlag()){
					RHits++;
				}
			}else{
				notPrimary++;
			}
		}
	}
	
	public void printReadLength(){
		System.out.println(readLen);
	}
	
	public void printReadCount(){
		System.out.println(String.format("%.0f",+totalReads));
	}
	
	public void printInsertDistrib(){
		histo.printContents();
	}

	public void printStats(){
		System.out.println("\nTotalHits:\t"+String.format("%.0f",+totalHits));
		System.out.println("MappedSeq:\t"+String.format("%.0f",weight));
		System.out.println("UniquelyMapped:\t"+String.format("%.0f",uniquelyMapped));
		System.out.println("SingleEndMapped:\t"+String.format("%.0f",singleEnd));
		if(pairMapped>0){
			System.out.println("LeftHits:\t"+String.format("%.0f",LHits));
			System.out.println("RightHits:\t"+String.format("%.0f",RHits));
			System.out.println("PairedEndMapped:\t"+String.format("%.0f",pairMapped));
			System.out.println("ProperPairs:\t"+String.format("%.0f",properPair));
			System.out.println("Junctions:\t"+String.format("%.0f",junctions));
			System.out.println("PairedEndMapped_SameChr:\t"+String.format("%.0f",pairedEndSameChr));
			System.out.println("PairedEndMapped_DiffChr:\t"+String.format("%.0f",pairedEndDiffChr));
			System.out.println("NotPrimary:\t"+String.format("%.0f",notPrimary));				
		}
	}
}
