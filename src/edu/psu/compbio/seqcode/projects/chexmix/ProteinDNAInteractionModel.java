package edu.psu.compbio.seqcode.projects.chexmix;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * ProteinDNAInteractionModel: defines a collection of composite model components 
 * 
 *   Coordinates are 0-based, but the midpoint of the model is defined by the centerOffset variable. 
 * @author mahony
 *
 */
public class ProteinDNAInteractionModel {

	protected int modelWidth;
	protected int centerOffset;
	protected ChExMixConfig cmConfig;
	protected List<CompositeModelComponent> allComponents;
	protected List<CompositeModelComponent> XLComponents;
	protected CompositeModelComponent CSComponent;
	protected CompositeModelComponent backgroundComponent;
	protected int numXLComponents;
	protected TagProbabilityDensity XLTagDist, CSTagDist, backTagDist; //TODO: Convert these to condition/replicate specific?
	
	public ProteinDNAInteractionModel(ChExMixConfig config, int width, TagProbabilityDensity initXLdistrib, TagProbabilityDensity initCSdistrib, 
			TagProbabilityDensity initBackDistrib, double noisePi){
		modelWidth = width;
		centerOffset = modelWidth/2;
		cmConfig = config;
		
		//Background component
		backTagDist = initBackDistrib;
		backgroundComponent = new CompositeModelComponent(backTagDist, centerOffset, 0, "Back",  false, true);

		//ChIP-seq component
		CSTagDist = initCSdistrib;
		CSComponent = new CompositeModelComponent(CSTagDist, centerOffset, 1, "CS", false, true);
				
		//XL components
		XLTagDist = initXLdistrib;
		XLComponents = new ArrayList<CompositeModelComponent>();
		numXLComponents = (initCSdistrib.getInfluenceRange()/2)/cmConfig.getXLComponentSpacing();
		int xlPos = centerOffset-(numXLComponents*cmConfig.getXLComponentSpacing())/2;
		for(int i=0; i<numXLComponents; i++){
			XLComponents.add(new CompositeModelComponent(XLTagDist, xlPos, i+2, "XL", true, true));
			xlPos += cmConfig.getXLComponentSpacing(); 
		}
		
		//All components
		allComponents = new ArrayList<CompositeModelComponent>();
		allComponents.add(backgroundComponent);
		allComponents.add(CSComponent);
		allComponents.addAll(XLComponents);
		
		//Set initial pi values
		setInitialPi(noisePi);
	}
	
	//Accessors
	public int getWidth(){return modelWidth;}
	public int getCenterOffset(){return centerOffset;}
	public TagProbabilityDensity getXLTagDistribution(){return XLTagDist;}
	public TagProbabilityDensity getCSTagDistribution(){return CSTagDist;}
	public TagProbabilityDensity getBackTagDistribution(){return backTagDist;}
	public int getNumComponents(){return allComponents.size();}
	public List<CompositeModelComponent> getAllComponents(){return allComponents;}
	public List<CompositeModelComponent> getXLComponents(){return XLComponents;}
	public CompositeModelComponent getCSComponent(){return CSComponent;}
	public CompositeModelComponent getBackgroundComponent(){return backgroundComponent;}
	
	/**
	 * Return the non-zero components
	 * @return
	 */
	public List<CompositeModelComponent> getNonZeroComponents(){
		List<CompositeModelComponent> comps = new ArrayList<CompositeModelComponent>();
		for(CompositeModelComponent c : allComponents)
			if(c.isNonZero())
				comps.add(c);
		return comps;
	}
	
	/**
	 * Set the initial pi values for all components. 
	 * Assumes model is initialized. 
	 */
	protected void setInitialPi(double noisePi){
		backgroundComponent.setPi(noisePi);
		double bindingPi = 1-noisePi;
		double initCS = Math.max(bindingPi*cmConfig.INIT_CS_TO_XL_RATIO, cmConfig.MIN_CS_PI);
		CSComponent.setPi(initCS);
		for(CompositeModelComponent xl : XLComponents){
			double xoPi = bindingPi*(1-initCS)/(double)numXLComponents;
			xl.setPi(xoPi);
		}
	}
	
	/**
	 * toString
	 */
	public String toString(){
		Collections.sort(XLComponents);
		String out = "ProteinDNAInteractionModel:\n";
		for(CompositeModelComponent xl : XLComponents){
			if(xl.isNonZero())
				out = out+"\t"+xl.getIndex()+"\tXL:\t"+xl.toString(centerOffset)+"\n";
		}
		out = out+"\t"+CSComponent.getIndex()+"\tCS:\t"+CSComponent.toString(centerOffset)+"\n";
		out = out+"\t"+backgroundComponent.getIndex()+"\tBack:\t"+backgroundComponent.toString(centerOffset)+"\n";
		return out;
	}
}
