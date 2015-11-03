package edu.psu.compbio.seqcode.projects.chexmix;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
	protected Map<CompositeModelComponent, TagProbabilityDensity> XLComponentDensities;
	protected CompositeModelComponent CSComponent;
	protected CompositeModelComponent backgroundComponent;
	protected int numXLComponents; 
	
	public ProteinDNAInteractionModel(ChExMixConfig config, int width, TagProbabilityDensity initXLdistrib, TagProbabilityDensity initCSdistrib, 
			TagProbabilityDensity initBackDistrib, double noisePi){
		modelWidth = width;
		centerOffset = modelWidth/2;
		cmConfig = config;
		
		//Background component
		backgroundComponent = new CompositeModelComponent(initBackDistrib, centerOffset, 0, "Back",  false, true);

		//ChIP-seq component
		CSComponent = new CompositeModelComponent(initCSdistrib, centerOffset, 1, "CS", false, true);
				
		//XL components
		XLComponents = new ArrayList<CompositeModelComponent>();
		XLComponentDensities = new HashMap<CompositeModelComponent, TagProbabilityDensity>();
		numXLComponents = (initCSdistrib.getInfluenceRange()/2)/cmConfig.getXLComponentSpacing();
		int xlPos = centerOffset-(numXLComponents*cmConfig.getXLComponentSpacing())/2;
		for(int i=0; i<numXLComponents; i++){
			TagProbabilityDensity XLTagDist = initXLdistrib.clone();
			CompositeModelComponent newComp = new CompositeModelComponent(XLTagDist, xlPos, i+2, "XL", true, true);
			XLComponents.add(newComp);
			XLComponentDensities.put(newComp, XLTagDist);
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
	public TagProbabilityDensity getXLTagDistribution(CompositeModelComponent comp){return XLComponentDensities.get(comp);}
	public TagProbabilityDensity getCSTagDistribution(){return CSComponent.getTagDistribution();}
	public TagProbabilityDensity getBackTagDistribution(){return backgroundComponent.getTagDistribution();}
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
