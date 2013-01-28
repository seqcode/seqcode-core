/*
 * Author: tdanford
 * Date: Sep 16, 2008
 */
package edu.psu.compbio.seqcode.gse.viz.eye;

import java.util.Iterator;

import edu.psu.compbio.seqcode.gse.utils.models.Model;
import edu.psu.compbio.seqcode.gse.viz.paintable.Paintable;

public interface ModelPaintable extends Paintable {

	public void addModels(Iterator<? extends Model> itr);
	public void addModel(Model m);
	public void addValue(Object v);
	public void clearModels();
	
	public ModelPaintable setProperty(ModelPaintableProperty p);
	public ModelPaintable setProperty(String key, ModelPaintableProperty p);
	public ModelPaintable synchronizeProperty(String k, ModelPaintable p);
}
