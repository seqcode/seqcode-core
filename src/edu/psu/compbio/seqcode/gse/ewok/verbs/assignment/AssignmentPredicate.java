/*
 * Author: tdanford
 * Date: Jun 9, 2008
 */
package edu.psu.compbio.seqcode.gse.ewok.verbs.assignment;

import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Mapper;

public interface AssignmentPredicate<X extends Region> {

	public boolean isValidAssignment(Region item, X event);
	public Mapper<Region,Region> assignmentZoneMapper();
	
	public static class Filter<Y extends Region> 
		implements edu.psu.compbio.seqcode.gse.ewok.verbs.Filter<Y,Y> {
		
		private AssignmentPredicate<Y> predicate;
		private Region baseRegion;
		
		public Filter(AssignmentPredicate<Y> ap, Region r) {
			predicate = ap;
			baseRegion = r;
		}
		
		public Y execute(Y val) { 
			return predicate.isValidAssignment(baseRegion, val) ? val : null;
		}
	}
}


