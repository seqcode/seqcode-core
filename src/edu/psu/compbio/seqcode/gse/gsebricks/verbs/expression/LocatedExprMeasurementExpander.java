package edu.psu.compbio.seqcode.gse.gsebricks.verbs.expression;

import java.sql.SQLException;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedList;

import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.datasets.expression.Experiment;
import edu.psu.compbio.seqcode.gse.datasets.expression.ExprMeasurement;
import edu.psu.compbio.seqcode.gse.datasets.expression.ExpressionLoader;
import edu.psu.compbio.seqcode.gse.datasets.expression.LocatedExprMeasurement;
import edu.psu.compbio.seqcode.gse.datasets.expression.LocatedProbe;
import edu.psu.compbio.seqcode.gse.datasets.expression.ProbePlatform;
import edu.psu.compbio.seqcode.gse.gsebricks.iterators.EmptyIterator;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.Expander;
import edu.psu.compbio.seqcode.gse.utils.Closeable;

public class LocatedExprMeasurementExpander 
	implements Expander<Region,LocatedExprMeasurement>, Closeable {
	
	private ExpressionLoader loader;
	private boolean shouldCloseLoader;
    private Experiment expt;
    private ProbePlatform plat;
	
	public LocatedExprMeasurementExpander(ExpressionLoader l, 
			String exptName, String platName) throws SQLException { 
		shouldCloseLoader = false;
		loader = l;
        expt = loader.loadExperiment(exptName);
        plat = loader.loadPlatform(platName);
	}

	public Iterator<LocatedExprMeasurement> execute(Region a) {        
        try {
            Collection<LocatedExprMeasurement> exprs = loader.loadMeasurementsInRegion(a, plat, expt);
            return exprs.iterator();        
        } catch (SQLException e) {
            e.printStackTrace();
            return new EmptyIterator<LocatedExprMeasurement>();
        }
	}

	public void close() {
		if(shouldCloseLoader && !loader.isClosed()) { 
			loader.close();
		}
		loader.close();
	}

	public boolean isClosed() {
		return loader == null || loader.isClosed();
	}
}
