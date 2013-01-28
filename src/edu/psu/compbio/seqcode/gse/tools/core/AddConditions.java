package edu.psu.compbio.seqcode.gse.tools.core;

import java.sql.*;

import edu.psu.compbio.seqcode.gse.utils.*;
import edu.psu.compbio.seqcode.gse.utils.database.DatabaseException;
import edu.psu.compbio.seqcode.gse.utils.database.DatabaseFactory;
import edu.psu.compbio.seqcode.gse.utils.database.Sequence;
import edu.psu.compbio.seqcode.gse.utils.database.UnknownRoleException;

/**
 * adds the listed conditions to the DB:
 * AddConditions YPD culture tissue
 */

public class AddConditions {

    public static void main(String args[]) throws Exception {
        java.sql.Connection cxn = DatabaseFactory.getConnection("core");        
        PreparedStatement insert = cxn.prepareStatement("insert into conditions(id,name) values(" +
                                                        Sequence.getInsertSQL(cxn,"condition_id") + ",?)");
        for (int i = 0; i < args.length; i++) {
            insert.setString(1,args[i]);
            try {
                insert.execute();
            } catch (SQLException e) {
                System.err.println(e.toString());
            }

        }

    }

}