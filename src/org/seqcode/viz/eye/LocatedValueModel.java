/*
 * Author: tdanford
 * Date: Aug 7, 2009
 */
package org.seqcode.viz.eye;

import java.awt.Color;

import org.seqcode.gseutils.models.Model;


public class LocatedValueModel extends Model {
	
	public Integer location;
	public Double value; 
	public Color color;
	
	public LocatedValueModel() {}
	
	public LocatedValueModel(int loc, double val, Color c) { 
		location = loc; value = val;
		color = c;
	}
}