/*
 * Author: tdanford
 * Date: Dec 4, 2008
 */
package org.seqcode.ml.bayesnets;

import org.seqcode.gseutils.ArrayUtils;
import org.seqcode.gseutils.models.*;
import org.seqcode.ml.regression.DataFrame;

public class NaiveBayes<X extends Model> extends BN<X> {

	public NaiveBayes(DataFrame<X> data, String classField, String... attrs) {
		super(data, ArrayUtils.prepend(classField, attrs));
		for (int i = 0; i < attrs.length; i++) {
			graph.addEdge(classField, attrs[i]);
		}
		learnCPDs();
	}

}
