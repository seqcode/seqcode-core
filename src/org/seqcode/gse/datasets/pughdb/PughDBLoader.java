package org.seqcode.gse.datasets.pughdb;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import org.seqcode.gse.utils.database.DatabaseException;
import org.seqcode.gse.utils.database.DatabaseFactory;


/**
 * @author mahony
 * 
 * PughDBLoader serves as a clearinghouse for query interactions with the Pugh lab database
 */
public class PughDBLoader {
	public static String role = "pughdb";
	private java.sql.Connection cxn=null;

	public static void main(String[] args) throws Exception{
		PughDBLoader loader = new PughDBLoader();
		Collection<PughLabSampleInfo> info = loader.loadAllSamples();
		for(PughLabSampleInfo i : info){
			System.out.println(i.toString());
		}
	}

	public PughDBLoader(){
		
	}
	
	public java.sql.Connection getConnection() {
        if (cxn == null) {
            try {
                cxn = DatabaseFactory.getConnection(role);
            } catch (SQLException e) {
                throw new DatabaseException(e.toString(),e);
            }
        }
        return cxn;
    }
	
	/**
	 * Load all samples in Pugh DB
	 * @return
	 * @throws SQLException
	 */
	public Collection<PughLabSampleInfo> loadAllSamples() throws SQLException {
		PreparedStatement ps = PughLabSampleInfo.createLoadAll(getConnection());
		ps.setFetchSize(1000);
		List<PughLabSampleInfo> samples = new ArrayList<PughLabSampleInfo>();
		ResultSet rs = ps.executeQuery();
		while (rs.next()) {
			samples.add(new PughLabSampleInfo(rs));
		}
		Collections.sort(samples);
		rs.close();
		ps.close();
		return samples;
	}
	
}
