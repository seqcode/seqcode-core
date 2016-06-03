package org.seqcode.gse.tools.pughdb;

import java.sql.SQLException;
import java.util.Collection;

import org.seqcode.gse.datasets.pughdb.PughDBLoader;
import org.seqcode.gse.datasets.pughdb.PughLabSampleInfo;


public class ExportPughDBSampleList {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		PughDBLoader loader = new PughDBLoader();
		Collection<PughLabSampleInfo> info;
		try {
			info = loader.loadAllSamples();
		
			for(PughLabSampleInfo i : info){
				System.out.println(i.buildName()+"\t"+i.getUniqID()+"\t"+i.convertAssayCode()+"\t"+i.getSeqDataConditionName()+"\t"+i.getSeqDataTargetName()+"\t"+i.getSeqDataCelllineName()+"\t"+i.getSeqDataReplicateName());
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
	}

}
