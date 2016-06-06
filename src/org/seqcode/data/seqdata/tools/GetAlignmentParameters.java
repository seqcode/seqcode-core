package org.seqcode.data.seqdata.tools;

import java.sql.*;
import java.util.*;
import java.io.IOException;

import org.seqcode.data.seqdata.*;
import org.seqcode.utils.*;


/**
 * Dumps data from the alignmentparameters table for a single alignment 
 * 
 * GetAlignmentParameters --id ID
 */

public class GetAlignmentParameters {
    public static void main(String args[]) throws SQLException, NotFoundException, IOException {
        
        SeqDataLoader loader = new SeqDataLoader();
        
        Integer id = Args.parseInteger(args,"id", -1);
        if(id != -1){
        	SeqAlignment align = loader.loadAlignment(id);
        	Map<String,String> params = loader.getAlignmentParameters(align);
        	
        	for(String s : params.keySet()){
        		System.out.println(s+"="+params.get(s));
        	}
        }
        loader.close();
    }
}