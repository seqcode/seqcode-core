package edu.psu.compbio.seqcode.projects.akshay.MultiSeq;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import weka.classifiers.AbstractClassifier;
import weka.classifiers.functions.Logistic;
import weka.core.*;
import weka.core.pmml.PMMLProducer;
import weka.filters.unsupervised.attribute.NominalToBinary;
import weka.filters.unsupervised.attribute.RemoveUseless;
import weka.filters.unsupervised.attribute.ReplaceMissingValues;

public class MultiLogistic extends AbstractClassifier implements OptionHandler, WeightedInstancesHandler, TechnicalInformationHandler, PMMLProducer, Aggregateable<MultiLogistic> {

	// Model parameters compatible with a normal Multinomial logit model
	
	/** The coefficients of the optimized Structured Multinomial logit model. However, this only stores the leaf node parameters*/
	protected double[][] m_Par; // [NumPredictors+Intercept][NumClasses]
	
	/** The training data saved as a matrix */
	protected double[][] m_Data; // [NumInstances][NumPredictors+Intercept]
	
	/** The number of class lables or the number of leaf nodes */
	protected int m_NumClasses;
	
	/** The regularization parameter, Is multiplied to the logit part of the loss function */
	protected double m_Ridge = 100;
	
	/** The index of the class attribute. Usually the last attribute */
	protected int m_ClassIndex;
	
	/** An attribute filter */
	private RemoveUseless m_AttFilter;
	
	/** The filter used to make attributes numeric. */
	private NominalToBinary m_NominalToBinary;
	
	/** The filter used to get rid of missing values. */
	private ReplaceMissingValues m_ReplaceMissingValues;
	 
	/** Log-likelihood of the searched model */
	protected double m_LL;
	
	/** The maximum number of iterations. */
	private int m_MaxIts = -1;
	 
	/** Wether to use conjugate gradient descent rather than BFGS updates. */
	private boolean m_useConjugateGradientDescent = false;
	  
	private Instances m_structure;
	
	
	// Model parameters of the structured multinomial logit not compatible with the weka logit class
	
	/** All model parameters including the internal non-leaf nodes */
	protected double[][] sm_Par; //[NumPredictors+Intercept][NumNodes]
	
	/** Total number of nodes in the structured representation of the classes. Including the root node */
	protected int sm_NumNodes;
	
	//private ClassRelation
	  
	  
	  
	
	
	
	@Override
	public void buildClassifier(Instances data) throws Exception {
		// TODO Auto-generated method stub
		
	}

	@Override
	public MultiLogistic aggregate(MultiLogistic toAggregate) throws Exception {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void finalizeAggregation() throws Exception {
		// TODO Auto-generated method stub
		
	}

	@Override
	public String toPMML(Instances train) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public TechnicalInformation getTechnicalInformation() {
		// TODO Auto-generated method stub
		return null;
	}
	
	
	private class OptObject {
		
	}
	
	
	protected class ClassRelationStructure{
		
		protected List<Node> leafs = new ArrayList<Node>();
		protected Node root;
		protected List<List<Node>> layers = new ArrayList<List<Node>>();
		protected List<Node> allNodes = new ArrayList<Node>();
		
		
		public ClassRelationStructure(File structureFile) throws IOException {
			BufferedReader br = new BufferedReader(new FileReader(structureFile));
			String line;
			List<String> fileContents = new ArrayList<String>();
			HashMap<Integer,Node> addedNodes = new HashMap<Integer,Node>();
	        while ((line = br.readLine()) != null) {
	        	if(!line.startsWith("#")){
	        		fileContents.add(line);
	        		String[] pieces = line.split("\t");
	        		if(pieces.length <4){System.err.println("Incorrect class structure format!!!"); System.exit(1);}
	        	}
	        }br.close();
			for(int n=0; n<fileContents.size(); n++){
				String[] pieces = fileContents.get(n).split("\t");
				
				// Add the current node if not already added
        		Node currNode = new Node(Integer.parseInt(pieces[0]), pieces[1],Integer.parseInt(pieces[2]) == 1? true : false);
        		if(!addedNodes.containsKey(currNode.nodeIndex)){
        			addedNodes.put(currNode.nodeIndex, currNode);
        		}
        		
        		if(!currNode.nodeName.equals("Root")){
        			String[] parentInds = pieces[3].split(",");
        			for(String s :  parentInds){
        				int pind = Integer.parseInt(s);
        				String pString = fileContents.get(pind);
        				String[] Ppieces = pString.split("\t");
        				if(addedNodes.containsKey(pind)){
        					addedNodes.get(pind).addChild(currNode);
        				}else{
        					addedNodes.put(pind, new Node(pind,Ppieces[1],false));
        					addedNodes.get(pind).addChild(currNode);
        				}
        			}
        		}
			}
			
			
			
			
			
			
			
			
        
		}
		
		
		
		protected class Node implements Comparable<Node>{
			
			protected int nodeIndex;
			protected String nodeName;
			protected boolean isLeaf;
			protected int level;
			protected boolean isRoot;
			protected List<Node> parents = new ArrayList<Node>();
			protected List<Node> children = new ArrayList<Node>();
			
			protected Node(int nInd, String nName, boolean leaf) {
				nodeIndex = nInd;
				nodeName = nName;
				isLeaf = leaf;
			}
			
			protected void addParent(Node p){
				boolean add=true;
				for(Node ps: parents){
					if(ps.nodeIndex == p.nodeIndex){
						add = false;
					}
				}
				if(add)
					parents.add(p);
			}
			protected void addChild(Node c){
				boolean add = true;
				for(Node cs : children){
					if(cs.nodeIndex == c.nodeIndex){
						add = false;
					}
				}
				if(add)
					children.add(c);
			}
			
			@Override
			public int compareTo(Node o) {
				if(o.nodeIndex == nodeIndex){return 0;}
				else if(nodeIndex > o.nodeIndex){return 1;}
				else{return -1;}
			}
			
		}
	}
	
	
	
	
	
	
	
	

}
