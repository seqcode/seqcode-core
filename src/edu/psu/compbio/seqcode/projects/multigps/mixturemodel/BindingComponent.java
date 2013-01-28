package edu.psu.compbio.seqcode.projects.multigps.mixturemodel;

import edu.psu.compbio.seqcode.gse.datasets.general.Point;
import edu.psu.compbio.seqcode.projects.multigps.framework.BindingModel;
import edu.psu.compbio.seqcode.projects.multigps.framework.ReadHit;
import edu.psu.compbio.seqcode.projects.multigps.framework.StrandedBaseCount;

/**
 * BindingComponents are used in mixture models to represent potential binding events.
 * 
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class BindingComponent implements Comparable<BindingComponent>{

	protected Point coord;  //Event coordinate
	protected int position; //Position without the chromosome name (for computational convenience)
	protected double pi; //Emission probability
	protected double sum_resp; //Sum of read responsibilities
	protected double[][] readProfile_plus;  //Read profiles are left uninitialized, as they are only used to store read distributions at the end of training. 
	protected double[][] readProfile_minus; //The read distributions stored in readProfiles are centered on the binding position. 
	protected int index=0;
	
	public BindingComponent(Point pos, int numReps){
		coord=pos;
		position = coord.getLocation();
		sum_resp = 0;
		pi = 1;
		readProfile_plus= new double[numReps][];
		readProfile_minus= new double[numReps][];
		for(int r=0; r<numReps; r++){
			readProfile_plus[r]=null;
			readProfile_minus[r]=null;
		}	
	}
	
	//Accessors
	public double scoreHit(ReadHit h, BindingModel model){
		int dist = h.getStrand()=='+' ? h.getFivePrime()-coord.getLocation():coord.getLocation()-h.getFivePrime();
		return model.probability(dist);
	}
	public double scoreBase(StrandedBaseCount b, BindingModel model){
		int dist = b.getStrand()=='+' ? b.getCoordinate()-coord.getLocation():coord.getLocation()-b.getCoordinate();
		return model.probability(dist);
	}
	public double score(int dist,  BindingModel model){
		return model.probability(dist);
	}
	public double getPi(){return pi;}
	public Point getCoord(){return coord;}
	public int getPosition(){return position;}
	public int getIndex(){return index;}
	
	public boolean isNonZero(){return pi>0;}

	public double getSumResponsibility(){return sum_resp;}
	
	public double[] getReadProfile_plus(int repIndex){
		return readProfile_plus[repIndex];
	}
	
	public double[] getReadProfile_minus(int repIndex){
		return readProfile_minus[repIndex];
	}
	
	
	public double getReadProfile(int repIndex, int i, char strand){
		double result=0;
		if (strand=='+')
			result=readProfile_plus[repIndex][i];
		else if (strand=='-')
			result=readProfile_minus[repIndex][i];
		return result;
	}
	
	//Mutators
	public void setPi(double p){pi=p;}
	public void setPosition(int p){position = p; updateCoordFromLocation();}
	public void setCoord(Point p){coord=p; position=p.getLocation();}
	public void updateCoordFromLocation(){Point newCoord = new Point(coord.getGenome(), coord.getChrom(), position); coord=newCoord;}
	public void setIndex(int i){index=i;}
	public void setSumResponsibility(double sum_resp) { this.sum_resp = sum_resp; }
		
	public void uniformInit(double initValue){
		pi=initValue;
	}
	
	public void setReadProfile(int repIndex, double[] profile, char strand){
		if (readProfile_plus[repIndex]==null){
			readProfile_plus[repIndex] = new double[profile.length];
			readProfile_minus[repIndex] = new double[profile.length];
		}
		if (strand=='+')
			readProfile_plus[repIndex] = profile;
		else if (strand=='-')
			readProfile_minus[repIndex] = profile;
	}

	//Comparable default method
	public int compareTo(BindingComponent m) {
		return getCoord().compareTo(m.getCoord());
	}
	
	//Compare by responsibility
	public int compareByResp(BindingComponent m) {
		return Double.compare(sum_resp, m.sum_resp);
	}
	
	public String toString(){
		return "B\t"+coord.getLocationString()+"\t"+String.format("%.3f",pi)+"\t"+String.format("%.3f", sum_resp)+"\t"+index;
	}
}//end of BindingComponent class
