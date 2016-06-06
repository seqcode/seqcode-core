/*
 * Author: tdanford
 * Date: Jun 9, 2008
 */
package org.seqcode.gsebricks.verbs.assignment;

import org.seqcode.genome.location.Region;
import org.seqcode.gsebricks.verbs.Mapper;

public interface AssignmentPredicate<X extends Region> {

	public boolean isValidAssignment(Region item, X event);
	public Mapper<Region,Region> assignmentZoneMapper();
	
	public static class Filter<Y extends Region> 
		implements org.seqcode.gsebricks.verbs.Filter<Y,Y> {
		
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


