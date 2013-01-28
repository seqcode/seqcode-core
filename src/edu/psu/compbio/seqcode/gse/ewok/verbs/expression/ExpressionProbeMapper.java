package edu.psu.compbio.seqcode.gse.ewok.verbs.expression;

import java.sql.SQLException;

import edu.psu.compbio.seqcode.gse.datasets.expression.Experiment;
import edu.psu.compbio.seqcode.gse.datasets.expression.ExprMeasurement;
import edu.psu.compbio.seqcode.gse.datasets.expression.ExpressionLoader;
import edu.psu.compbio.seqcode.gse.datasets.expression.Probe;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Filter;
import edu.psu.compbio.seqcode.gse.utils.Closeable;

public class ExpressionProbeMapper 
	implements Filter<Probe,ExprMeasurement>, Closeable {
	
	private ExpressionLoader loader;
	private boolean shouldCloseLoader;
	private Experiment expt;
	
	public ExpressionProbeMapper(ExpressionLoader el, String exptName) throws SQLException { 
		loader = el;
		shouldCloseLoader = false;
		expt = loader.loadExperiment(exptName);
	}
	
	public ExpressionProbeMapper(String exptName) throws SQLException { 
		loader = new ExpressionLoader();
		shouldCloseLoader = true;
		expt = loader.loadExperiment(exptName);
	}
	
	public ExpressionLoader getLoader() { return loader; }
	public Experiment getExperiment() { return expt; }

	public ExprMeasurement execute(Probe a) {
		try {
			return loader.loadMeasurement(expt, a);
		} catch (SQLException e) {
			e.printStackTrace();
			return null;
		}
	}

	public void close() {
		if(shouldCloseLoader && !loader.isClosed()) { 
			loader.close();
		}
		loader = null;
	}

	public boolean isClosed() {
		return loader == null || loader.isClosed();
	}
}
