package edu.psu.compbio.seqcode.gse.tools.seqdata;

import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;

import edu.psu.compbio.seqcode.gse.datasets.general.MetadataLoader;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqAlignment;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqDataLoader;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqExpt;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;

/**
 * Export existing SeqExpt & SeqAlignment descriptions to a file. 
 * This exports a file like deepseq.list for backing up the seqdata mysql database. 
 *  
 * @author mahony
 *
 * Usage: ExportSeqAlignmentList --list "filename"
 * 
 * The assumed output file is in the deepseq.list format, with the following fields:
 * 
 *0) ReadDBID
 *1) ExptType
 *2) Lab
 *3) ExptCondition
 *4) ExptTarget
 *5) CellLine
 *6) Replicate
 *7) Aligner
 *8) Genome
 *9) Permissions
 *10) PubSource
 *11) PublicDBID
 *12) CollabExptID
 *13) CollabAlignID
 *14) ReadType
 *15) AlignType
 *16) ReadLength
 *17) TotalReads
 *18) AlignedHits
 *19) UniquelyAlignedHits
 *20) DBLoadedHits
 *21) DBLoadedWeight
 *22) DBLoadedPairs
 *23) DBLoadedPairWeight
 *24) ReadsFile
 *25) AlignDir
 *26) AlignFile
 *27) IDXFile
 *28) AlignParamFile
 *29) ExptNote
 *30) LoadDate
 *31) ExptName
 * 
 */
public class ExportSeqAlignmentList {
	public static void main(String args[]) throws SQLException, IOException{
		String outFile="out";
		ArgParser ap = new ArgParser(args);
		if(ap.hasKey("h")){
			System.err.println("ExportSeqAlignmentList:\n" +
					"\t--out <output filename>\n"
					);
		}
		if(ap.hasKey("out"))
			outFile = ap.getKeyValue("out");
		
		SeqDataLoader loader = new SeqDataLoader();
        MetadataLoader core = new MetadataLoader();
        
        try{
	        FileWriter fw = new FileWriter(outFile);
	        
	        for(SeqExpt expt : loader.loadAllExperiments()){
	        	for(SeqAlignment align : loader.loadAllAlignments(expt)){
	        		String permissions ="";
	        		for(String p : align.getPermissions())
	        			permissions = permissions+p+";";
	        		
	        		fw.write(align.getDBID()+"\t"+
	        				expt.getExptType().getName()+"\t"+
	        				expt.getLab().getName()+"\t"+
	        				expt.getExptCondition().getName()+"\t"+
	        				expt.getExptTarget().getName()+"\t"+
	        				expt.getCellLine().getName()+"\t"+
	        				expt.getReplicate()+"\t"+
	        				align.getName()+"\t"+
	        				align.getGenome().getName()+"\t"+
	        				permissions+"\t"+
	        				expt.getPublicSource()+"\t"+
	        				expt.getPublicDBID()+"\t"+
	        				expt.getCollabID()+"\t"+
	        				align.getCollabAlignID()+"\t"+
	        				expt.getReadType().getName()+"\t"+
	        				align.getAlignType().getName()+"\t"+
	        				expt.getReadLength()+"\t"+
	        				expt.getNumRead()+"\t"+
	        				//alignedHits missing
	        				//uniquelyAlignedHits missing
	        				align.getNumHits()+"\t"+
	        				align.getTotalWeight()+"\t"+
	        				align.getNumPairs()+"\t"+
	        				align.getTotalPairWeight()+"\t"+
	        				expt.getFQFile()+"\t"+
	        				align.getAlignDir()+"\t"+
	        				align.getAlignFile()+"\t"+
	        				align.getIDXFile()+"\t"
	        				);
	        		
	        		
	        	}
	        }
	        fw.close();
        } catch (IOException e) {
			e.printStackTrace();
		}

	}
}
