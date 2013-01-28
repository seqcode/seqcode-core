/*
 * Author: tdanford
 * Date: Aug 25, 2008
 */
package edu.psu.compbio.seqcode.gse.datasets.alignments;

import java.sql.*;
import java.util.Collection;

import edu.psu.compbio.seqcode.gse.datasets.DBID;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.utils.database.Sequence;

/*
create sequence alignment_version_id;	 
create table alignment_version (
	id number(10) constraint cst_alignment_version_id unique not null,
	name varchar2(1000) not null
);
*/
	
public interface AlignmentVersion { 
	public String getName();
	public DBID getDBID();
	
	public Collection<AlignBlock> getAlignBlocks(Region r);
}

