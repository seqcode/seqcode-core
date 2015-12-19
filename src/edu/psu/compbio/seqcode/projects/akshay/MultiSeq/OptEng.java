package edu.psu.compbio.seqcode.projects.akshay.MultiSeq;

import edu.psu.compbio.seqcode.projects.akshay.MultiSeq.Optimizer.OptObject;
import weka.core.Optimization;
import weka.core.RevisionUtils;

public class OptEng extends Optimization {
	
	OptObject m_oO = null;
	
	public OptEng(OptObject oO){
		m_oO = oO;
	}
	 
	@Override
	public String getRevision() {
		return RevisionUtils.extract("$Revision: 11247 $");
	}

	@Override
	protected double objectiveFunction(double[] x) throws Exception {
		return m_oO.objectiveFunction(x);
	}

	@Override
	protected double[] evaluateGradient(double[] x) throws Exception {
		return m_oO.evaluateGradient(x);
	}

}

