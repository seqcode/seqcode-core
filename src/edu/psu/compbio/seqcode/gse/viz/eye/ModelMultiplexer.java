/*
 * Author: tdanford
 * Date: Sep 18, 2008
 */
package edu.psu.compbio.seqcode.gse.viz.eye;

import java.util.*;

import edu.psu.compbio.seqcode.gse.utils.models.Model;

public class ModelMultiplexer {

	private Vector<ModelPaintable> paintables;
	
	public ModelMultiplexer(ModelPaintable... ps) { 
		paintables = new Vector<ModelPaintable>();
		for(int i = 0; i < ps.length; i++) { 
			paintables.add(ps[i]);
		}
	}
	
	public ModelMultiplexer(Collection<ModelPaintable> ps) { 
		paintables = new Vector<ModelPaintable>();
	}
	
	public void addModelPaintable(ModelPaintable mp) { 
		paintables.add(mp);
	}
	
	public void addModel(Model m) { 
		for(ModelPaintable p : paintables) { 
			p.addModel(m);
		}
	}
	
	public <X extends Model> void addModels(Iterator<X> itr) { 
		while(itr.hasNext()) { 
			addModel(itr.next());
		}
	}
}
