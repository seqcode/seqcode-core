package org.seqcode.data.pegr.tools;

import java.sql.SQLException;
import java.util.Collection;

import org.seqcode.data.pegr.PughDBLoader;
import org.seqcode.data.pegr.PughLabSampleInfo;

public class ExportPughDBSampleList {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		PughDBLoader loader = new PughDBLoader();
		Collection<PughLabSampleInfo> info;
		try {
			info = loader.loadAllSamples();

			for (PughLabSampleInfo i : info) {
				System.out.println(i.buildName() + "\t" + i.getUniqID() + "\t" + i.convertAssayCode() + "\t"
						+ i.getSeqDataConditionName() + "\t" + i.getSeqDataTargetName() + "\t"
						+ i.getSeqDataCelllineName() + "\t" + i.getSeqDataReplicateName());
			}
		} catch (SQLException e) {
			e.printStackTrace();
		}
	}

}
