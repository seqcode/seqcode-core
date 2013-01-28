package edu.psu.compbio.seqcode.projects.multigps.framework;

/**
 * <tt>ReadHit</tt> represents the stranded coordinates of a mapped alignment hit
 * Read and ReadHit are only used as a convenient way to represent multiply mapping reads when reading from certain file formats.
 * <u>Note</u>: <br>
 *              1) For a single read, we can have multiple (read) hits <br>
 *              2) The <tt>start < end</tt> of the <tt>ReadHit</tt> ALWAYS   <br>
 *                 strand <tt>str</tt> is used to indicate whether the read is
 *                 forward or reverse.
 * @author mahony
 *
 */
public class ReadHit{

	protected String chrom;
	protected int start, end;
	protected char strand;
	protected float weight;
		
	public ReadHit(String c, int s, int e, char str){this(c,s,e,str, 1);}
	public ReadHit(String c, int s, int e, char str, float w){
		chrom = c;
		start = s;
		end = e;
		strand = str;
		weight = w;
	}
	
	//Accessors
	public String getChrom(){return chrom;}
	public int getStart(){return start;}
	public int getEnd(){return end;}
	public int getFivePrime(){return(strand=='+' ? start :end);}
	public char getStrand(){return strand;}
	public float getWeight(){return(weight);}
	public void setWeight(float w){weight=w;}
	public double getReadLength(){return(this.getEnd()-this.getStart()+1);}
}
