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
	protected TagDistribution XLTagDist, CSTagDist, backTagDist; //TODO: Convert these to condition/replicate specific?
	
	public ProteinDNAInteractionModel(ChExMixConfig config, int width, TagDistribution initXLdistrib, TagDistribution initCSdistrib, 
			TagDistribution initBackDistrib, double noisePi){
		modelWidth = width;
		centerOffset = modelWidth/2;
		cmConfig = config;
		
		//Background component
		backTagDist = initBackDistrib;
		backgroundComponent = new CompositeModelComponent(backTagDist, centerOffset, 0, false, true);

		//ChIP-seq component
		CSTagDist = initCSdistrib;
		CSComponent = new CompositeModelComponent(CSTagDist, centerOffset, 1, false, true);
				
		//XL components
		XLTagDist = initXLdistrib;
		XLComponents = new ArrayList<CompositeModelComponent>();
		numXLComponents = (initCSdistrib.getInfluenceRange()/4)/cmConfig.XLComponentSpacing;
		int xlPos = centerOffset-(numXLComponents*cmConfig.XLComponentSpacing)/2;
		for(int i=0; i<numXLComponents; i++){
			XLComponents.add(new CompositeModelComponent(XLTagDist, xlPos, i+2, true, true));
			xlPos += cmConfig.XLComponentSpacing; 
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
	public TagDistribution getXLTagDistribution(){return XLTagDist;}
	public TagDistribution getCSTagDistribution(){return CSTagDist;}
	public TagDistribution getBackTagDistribution(){return backTagDist;}
	public int getNumComponents(){return allComponents.size();}
	public List<CompositeModelComponent> getAllComponents(){return allComponents;}
	public List<CompositeModelComponent> getXLComponents(){return XLComponents;}
	public CompositeModelComponent getCSComponent(){return CSComponent;}
	public CompositeModelComponent getBackgroundComponent(){return backgroundComponent;}
	
	/**
	 * Set the initial pi values for all components. 
	 * Assumes model is initialized. 
	 */
	protected void setInitialPi(double noisePi){
		backgroundComponent.setPi(noisePi);
		double bindingPi = 1-noisePi;
		CSComponent.setPi(bindingPi*cmConfig.INIT_CS_TO_XO_RATIO);
		for(CompositeModelComponent xl : XLComponents){
			double xoPi = bindingPi*(1-cmConfig.INIT_CS_TO_XO_RATIO)/(double)numXLComponents;
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
