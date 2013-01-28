/*
 * Author: tdanford
 * Date: Dec 4, 2008
 */
package edu.psu.compbio.seqcode.gse.utils.models.bns;

import edu.psu.compbio.seqcode.gse.utils.ArrayUtils;
import edu.psu.compbio.seqcode.gse.utils.models.*;
import edu.psu.compbio.seqcode.gse.utils.models.data.DataFrame;

public class NaiveBayes<X extends Model> extends BN<X> {

	public NaiveBayes(DataFrame<X> data, String classField, String... attrs) { 
		super(data, ArrayUtils.prepend(classField, attrs));
		for(int i = 0; i < attrs.length; i++) { 
			graph.addEdge(classField, attrs[i]);
		}
		learnCPDs();
	}
	
	
}
