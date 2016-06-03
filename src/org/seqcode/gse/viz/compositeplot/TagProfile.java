package org.seqcode.gse.viz.compositeplot;

/**
 * TagProfile: per-base container for tag counts
 * @author mahony
 *
 */
public class TagProfile {

	protected double[] watsonTags; //binned watson tag counts
	protected double[] crickTags;  //binned crick tag counts
	protected int center; //The index at the center of the profile (rel coord = 0)
	protected int width;
	protected boolean stranded;
	
	
	/**
	 * TagProfile
	 * @param watson
	 * @param crick : if null, this is an unstranded tag density
	 * @param centerBin
	 */
	public TagProfile(double [] watson, double[] crick, int centerBin){
		watsonTags = watson;
		crickTags = crick;
		this.center = centerBin;
		width = watson.length;
		stranded = crickTags==null ? false : true;
	}
	
	//Accessors
	public double[] getWatsonTags(){return watsonTags;}
	public double[] getCrickTags(){return crickTags;}
	public int getWidth(){return width;}
	public int getCenterBin(){return center;}
	public int getLeftRelCoord(){return -center;}
	public int getRightRelCoord(){return (width-center-1);}
	public boolean isStranded(){return stranded;}
	
	protected int getBinAtRelCoord(int coord){
		int x=-1;
		if(coord>=getLeftRelCoord() && coord<=getRightRelCoord())
			x= center+coord;
		return x;
	}
	
	public double getWatsonTagsAtRelCoord(int coord){return watsonTags[getBinAtRelCoord(coord)];}
	public double getCrickTagsAtRelCoord(int coord){return crickTags[getBinAtRelCoord(coord)];}

}
