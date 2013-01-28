package edu.psu.compbio.seqcode.gse.ewok.verbs.expression;

import java.sql.SQLException;
import java.util.Collection;
import java.util.Iterator;

import edu.psu.compbio.seqcode.gse.datasets.expression.ExpressionLoader;
import edu.psu.compbio.seqcode.gse.datasets.expression.LocatedProbe;
import edu.psu.compbio.seqcode.gse.datasets.expression.ProbePlatform;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.ewok.verbs.Expander;
import edu.psu.compbio.seqcode.gse.utils.Closeable;
import edu.psu.compbio.seqcode.gse.utils.iterators.EmptyIterator;

public class ExpressionRegionExpander 
	implements Expander<Region,LocatedProbe>, Closeable {
	
	private boolean shouldCloseLoader;
	private ExpressionLoader loader;
	private ProbePlatform platform;

	public ExpressionRegionExpander(String pp) throws SQLException {
		shouldCloseLoader = true;
		loader = new ExpressionLoader();
		platform = loader.loadPlatform(pp);
	}
	
	public ExpressionRegionExpander(ExpressionLoader el, String pp) throws SQLException { 
		shouldCloseLoader = false;
		loader = el;
		platform = loader.loadPlatform(pp);
	}
	
	public ExpressionLoader getLoader() { return loader; }

	public Iterator<LocatedProbe> execute(Region a) {
		try {
			Collection<LocatedProbe> list = loader.loadProbesInRegion(a, platform);
			return list.iterator();
			
		} catch (SQLException e) {
			e.printStackTrace();
			return new EmptyIterator<LocatedProbe>();
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
