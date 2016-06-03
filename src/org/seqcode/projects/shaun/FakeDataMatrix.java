package org.seqcode.projects.shaun;

import java.util.Random;

import cern.jet.random.Normal;
import cern.jet.random.engine.RandomEngine;

public class FakeDataMatrix {

	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int numCells=120;
		int numTypes = 3;
		int numGenes = 1000;
		double gaussSigma=3;
		Random random = new Random();
		RandomEngine randomEngine = new cern.jet.random.engine.MersenneTwister();
		
		//Make the signatures
		double [][] typeSigs = new double[numTypes][numGenes];
		for(int i=0; i<numTypes; i++){
			for(int j=0; j<numGenes; j++){
				double rand = random.nextDouble();
				if(rand < 0.33333)
					typeSigs[i][j]=-10;
				else if(rand < 0.66666)
					typeSigs[i][j]=0;
				else
					typeSigs[i][j]=+10;
			}
		}
		
		//Make the noisy data matrix
		Normal gauss = new Normal(0, 1, randomEngine);
		double [][] data = new double[numCells][numGenes];
		for(int c=0; c<numCells; c++){
			int cellType = c%numTypes;
			
			for(int j=0; j<numGenes; j++){
				double x = typeSigs[cellType][j];
				double rand = gauss.nextDouble(x, gaussSigma);
				
				data[c][j]=rand;
			}
		}
		
		//Print header
		for(int c=0; c<numCells; c++){
			int cellType = c%numTypes;
			System.out.print("\t"+cellType+"_"+c);
		}
		System.out.println("");
		
		//Print data
		for(int j=0; j<numGenes; j++){
			System.out.print(j);
			for(int c=0; c<numCells; c++){
				System.out.print("\t"+data[c][j]);
			}System.out.println("");
		}
	}

}
