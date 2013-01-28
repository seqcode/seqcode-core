/*
 * Author: tdanford
 * Date: Aug 12, 2008
 */
package edu.psu.compbio.seqcode.gse.viz.metagenes;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Filter;

public interface PointProfiler<PointClass extends Point, ProfileClass extends Profile> extends Filter<PointClass,ProfileClass> { 
	public BinningParameters getBinningParameters();
	public void cleanup();
}
