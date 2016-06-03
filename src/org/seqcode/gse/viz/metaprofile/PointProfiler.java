/*
 * Author: tdanford
 * Date: Aug 12, 2008
 */
package org.seqcode.gse.viz.metaprofile;

import org.seqcode.genome.location.Point;
import org.seqcode.gse.gsebricks.verbs.Filter;

public interface PointProfiler<PointClass extends Point, ProfileClass extends Profile> extends Filter<PointClass,ProfileClass> { 
	public BinningParameters getBinningParameters();
	public void cleanup();
}
