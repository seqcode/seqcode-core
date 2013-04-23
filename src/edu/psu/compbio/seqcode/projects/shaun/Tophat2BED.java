package edu.psu.compbio.seqcode.projects.shaun;
import java.io.File;
import java.util.List;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloseableIterator;

import edu.psu.compbio.seqcode.gse.utils.ArgParser;
import edu.psu.compbio.seqcode.projects.multigps.framework.*;


public class Tophat2BED {

	protected File samFile;
	
    public Tophat2BED(File f) {
    	samFile = f;
    }

    /**
	 * Convert SAM/BAM alignment blocks to BED
	 */
	public void convert() {
		SAMFileReader reader = new SAMFileReader(samFile);
		reader.setValidationStringency(ValidationStringency.SILENT);
		CloseableIterator<SAMRecord> iter = reader.iterator();
		while (iter.hasNext()) {
		    SAMRecord record = iter.next();
		    
		    if (record.getReadUnmappedFlag()) {continue; }
		    float weight = 1/(float)record.getIntegerAttribute("NH");
		    
		    List<AlignmentBlock> blocks = record.getAlignmentBlocks();
		    for(int a=0; a<blocks.size(); a++){ //Iterate over alignment blocks
		    	AlignmentBlock currBlock = blocks.get(a);
		    	int aStart = currBlock.getReferenceStart();
		    	int aEnd = aStart + currBlock.getLength()-1;
		    	int aLen = currBlock.getLength();
		    	
		    	/*boolean nearbyBlocks=true;
		    	while(nearbyBlocks && a<blocks.size()-1){
		    		if(blocks.get(a+1).getReferenceStart() - currBlock.getReferenceStart() < record.getReadLength()){
		    			aEnd = blocks.get(a+1).getReferenceStart() + blocks.get(a+1).getLength()-1;
		    			aLen += blocks.get(a+1).getLength();
		    			a++;
		    		}else{
		    			nearbyBlocks=false;
		    		}
		    	}*/
		    	
		    	System.out.println(record.getReferenceName()+"\t"+aStart+"\t"+aEnd+"\t"+record.getReadName()+"\t"+weight+"\t"+(record.getReadNegativeStrandFlag() ? '-' : '+'));
			}
		}
		iter.close();
		reader.close();
    }
	
	/**
	 * 
	 * @param args
	 */
	public static void main(String[] args){
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("bam")) { 
            System.err.println("Tophat2BED:\n" +
                               "\tUsage:\n" +
                               "\t--bam <BAM/SAM file>\n");
        }else{
        	String bamFile = ap.getKeyValue("bam");
        	
        	Tophat2BED converter = new Tophat2BED(new File (bamFile));
        	converter.convert();
        	
        }
	}
}
