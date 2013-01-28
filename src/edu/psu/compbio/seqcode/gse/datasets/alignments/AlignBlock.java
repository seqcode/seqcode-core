/*
 * Author: tdanford
 * Date: Aug 25, 2008
 */
package edu.psu.compbio.seqcode.gse.datasets.alignments;

import edu.psu.compbio.seqcode.gse.datasets.DBID;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;

/*
create sequence align_block_id;
create table align_block (
	id number(10) constraint cst_align_block_id unique not null,
	alignment number(10) constraint cst_align_block_align references alignment(id),
	chromosome number(10) not null,
	start_pos number(10) not null,
	end_pos number(10) not null,
	strand char(1) not null,
	bit_string clob,
	gapped_length number(10) not null
);
 */

public interface AlignBlock {
	public DBID getDBID();
	public Alignment getAlignment();
	public Genome getGenome();
	public String getChrom();
	public int getStartPos();
	public int getStopPos();
	public char getStrand();
	public char[] getBitString();
	public int getGappedLength();
}

