package edu.psu.compbio.seqcode.projects.kunz.isotop;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Random;
/*
 * Figure out if unclustered (random) points will look different than clustered point
 * 		-All prototypes seem to line up so even random points may end up looking clustered
 * Need to add a relative scale to grid
 * Vector quantization: not working (all vectors filled with zero maybe?)
 *       -Need a way to tell what bins (names) are covered by each prototype
 * 
 */

import edu.psu.compbio.seqcode.projects.kunz.chromeSOM.Node;

public class Isochrome 
{
	public double epochs;
	public ArrayList<Prototype> protos;
	public int neighborNumber;
	public double[][] simMat;
	public double[][] delta;
	public int protoNum;
	public Isochrome()
	{
		protoNum = 30;
		epochs = 1000;
		neighborNumber = 6;
		protos = new ArrayList<Prototype>();
		simMat = new double[0][0];
	}
	public void go()
	{
		populateProtos();
		System.out.println("Quantizing");
		vectorQuant();
		pairWise();
		neighborhood();
		System.out.println("Graphing");
		graphDistances();
		System.out.println("Training Started");
		train();
		writeFile();
		System.out.println(protos.size());
	}
	public void populateProtos()
	{
		String ff = System.getProperty("user.dir")+"/MatrixLanding.txt";
		Reader reader = new Reader(ff);
		protos = reader.matrixReader();
		for(int i = 0; i<protos.size(); i++)
		{
			Prototype p = protos.get(i);
			p.myIndex=i;
			if(p.pmag == 0)
				System.out.println(p.name);
		}
		
	}
	public void vectorQuant() //k-means
	{
		ArrayList<Prototype> ks = new ArrayList<Prototype>();
		//init
		for(int i = 0; i < protoNum; i++)
		{
			Prototype pep = protos.get((int)(Math.random()*protos.size()));
			ks.add(new Prototype(pep.P_hiDVec, pep.name));
		}
		for(double i = 100; i>0;i--)
		{
			//find closest node for each prototype
			for(int j = 0 ; j < protos.size(); j++)
			{
				Prototype pop = protos.get(j);
				Prototype closest = ks.get(0);
				for(int k = 0; k< ks.size(); k++)
				{
					if(cosineSim(ks.get(k), pop)>cosineSim(pop,closest))
						closest=ks.get(k);
				}
				closest.neighbors.add(closest);
			}
			//find average value among neighbors for each node
			for(int j = 0; j< ks.size(); j++)
			{
				Prototype pop = ks.get(j);
				double[] num = new double[pop.P_hiDVec.length];
				double denom = 0;
				for(Prototype neighb : pop.neighbors)
				{
					for(int k = 0; k<neighb.P_hiDVec.length; k++)
					{
						num[k] += (i/1000) * neighb.P_hiDVec[k];
					}
					denom+= i/1000;
				}
				for (int k = 0; k < num.length; k++)
				{
					num[k]=num[k]/denom;
				}
				pop.P_hiDVec = num;
			}
		}
		for(int k = 0; k< ks.size(); k++)
		{
			ks.get(k).neighbors.removeAll(ks.get(k).neighbors);
		}
		protos.removeAll(protos);
		protos.addAll(ks);

		
		/*
		 * 
		 * Need a way to tell what bins (names) are covered by each prototype
		 * 
		 */
	}
	public void pairWise()
	{
		simMat = new double[protos.size()][protos.size()];
		for(int i = 0; i<protos.size(); i++)
		{
			Prototype p = protos.get(i);
			
			for(int j = 0; j<protos.size(); j++)
			{
				if(i == j)
					simMat[i][j] = 0;
				else
				{
					Prototype o = protos.get(j);
					simMat[i][j] = 1-cosineSim(p,o);
				}	
			}
		}
		/*
		for(int i = 0; i<protos.size(); i++)
		{
			for(int j = 0; j<protos.size(); j++)
			{
					System.out.print(simMat[i][j]+ "-");
			}
			System.out.println();
		}*/
	}
	public void neighborhood()
	{
		for(int i = 0; i<simMat.length; i++)
		{
			int[] inds = new int[neighborNumber]; for(int pops = 0; pops < neighborNumber; pops++){inds[pops] = 0;}
			double[] vals = new double[neighborNumber]; for(int pops = 0; pops < neighborNumber; pops++){vals[pops] = 1;}
			for(int j = 0; j < simMat[i].length; j++)
			{
				if(i != j)
				{
					int loc = findMax(vals);
					if(simMat[i][j]<vals[loc])
					{
						vals[loc] = simMat[i][j];
						inds[loc] = j;
					}
				}
			}
			protos.get(i).neighborIndexes = inds;
			protos.get(i).neighborDists = new double[inds.length];
			for(int y = 0; y<inds.length; y++)
			{
				protos.get(i).neighbors.add(protos.get(inds[y]));
				protos.get(i).neighborDists[y] = simMat[protos.get(i).myIndex][protos.get(i).neighbors.get(y).myIndex];
			}
		}
		/*for(int h = 0; h < protos.size(); h++)
		{
			Prototype fit = protos.get(h);
			System.out.print(fit.name+":");
			for(int j = 0; j< fit.neighborIndexes.length; j++)
			{
				System.out.print("   "+fit.neighbors.get(j).name + " = " + fit.neighborDists[j]);
			}
			System.out.println();
		}*/
		
	}
	//Dijkstra algorithm repeated for each Prototype as a source -> All-Pairs shortest path problem
	public void graphDistances()
	{
		double inf = 99999;
		delta = new double[protos.size()][protos.size()];
		for(int i = 0; i<protos.size();i++)
		{	for(int j = 0; j<protos.size(); j++)
			{
				delta[i][j]=inf;
			}
		}
		for(int i = 0; i<protos.size(); i++)
		{
			ArrayList<Prototype> visited = new ArrayList<Prototype>();
			ArrayList<Prototype> unvisited = new ArrayList<Prototype>();
			for(Prototype cray : protos)
			{
				unvisited.add(cray);
			}
			Prototype source = protos.get(i);
			Prototype current = protos.get(i);
			delta[i][source.myIndex]=0;
			while(unvisited.size()>0)
			{
				unvisited.remove(current);
				visited.add(current);
				//check distances for all pairwise distances of current vertex
				for(int j = 0; j< current.neighbors.size(); j++)
				{
					if(delta[i][current.neighbors.get(j).myIndex] > delta[i][current.myIndex]+current.neighborDists[j])
						delta[i][current.neighbors.get(j).myIndex] = delta[i][current.myIndex]+current.neighborDists[j];
				}
				//find closest unvisited vertex
				double min = inf;
				for(int j = 0; j < delta[i].length; j++)
				{
					if(delta[i][j] <= min && unvisited.contains(protos.get(j)))
					{
						min = delta[i][j];
						current = protos.get(j);
					}
				}
				if(delta[i][current.myIndex] == inf)
				{
					min = inf;
					for(int j = 0; j < simMat[i].length; j++) //find closest unvisited node to source if all unvisited nodes are disconnected from network
					{
						if(simMat[i][j] < min && unvisited.contains(protos.get(j)))
						{
							min = simMat[i][j];
							current = protos.get(j);
						}
					}
					
					int find = current.myIndex;
					min = simMat[find][i];
					Prototype close = source;
					for(int j = 0; j < simMat[find].length; j++) //find closest visited node to previously found current vertex
					{
						if(simMat[find][j] < min && visited.contains(protos.get(j)))
						{
							min = simMat[find][j];
							close = protos.get(j);
						}
					}
					delta[i][find] = delta[i][close.myIndex] + simMat[close.myIndex][find]; //sets delta for previously unseen current
				}
			}
		}
		/*for(int i = 0; i<delta.length; i++)
		{
			for(int j = 0; j < delta[i].length; j++)
			{
				delta[i][j] *= 100; //make it bigger
			}
		}*/
	}
	public int findMax(double[] u)
	{
		int ree = 0;
		double max = u[0];
		for(int i = 1; i<u.length;i++)
		{
			if(u[i]>max)
			{
				max = u[i];
				ree=i;
			}
		}
		return ree;
	}
	public void train()
	{
		populateMSpace();
		for(int i = 0; i<epochs; i++)
		{
			epoch(i);
		}
	}
	public void populateMSpace()
	{
		Prototype pepper;
		for(int i = 0; i<protos.size(); i++)
		{
			pepper = protos.get(i);
			for(int j = 0; j<pepper.m.length;j++)
			{
				Random randy = new Random();
				pepper.m[j] = randy.nextGaussian();
			}
		}
	}
	public void epoch(int iter)
	{
		double alpha = (1-(iter/epochs));
		for(int i = 0; i<protos.size(); i++)
		{
			Prototype stimulator = protos.get((int)(Math.random()*protos.size()));
			double[] stimulus = new double[stimulator.m.length];
			
			//generates a stimulus
			for(int j = 0; j<stimulus.length;j++)
			{
				Random randy = new Random();
				stimulus[j] = randy.nextGaussian();//-stimulator.m[j]; //generate a stimulus distributed around each prototype
			}
			//System.out.println(stimulus[0] + ", " + stimulus[1]);
			//finds the best matched unit for the stimulus
			Prototype BMU = protos.get(0); 
			double bm = euclidean(BMU.m, stimulus);
			for(int k = 1; k<protos.size(); k++)
			{
				double compare = euclidean(protos.get(i).m, stimulus);
				if (compare<bm)
				{
					bm = compare;
					BMU = protos.get(k);
				}
			}
			//updates protype m's
			for(Prototype updating: protos)
			{
				double nu = nuFactor(BMU, updating);
				double weightFactor = alpha*nu;
				//System.out.println(nu+ " "+ alpha + " " + weightFactor);
				for(int j = 0; j<updating.m.length; j++)
				{
					updating.m[j] += weightFactor*(stimulus[j]-updating.m[j]);
				}
			}
			
		}
	}
	public double nuFactor(Prototype BMU, Prototype updating)
	{
		double d = delta[BMU.myIndex][updating.myIndex]*delta[BMU.myIndex][updating.myIndex];
		double lambda = 5; //neighborhood width (compare to sigma in SOM)
		double mu = 0;
		for(int i = 0; i< BMU.neighborDists.length; i++) /////**** should it be average distance from BMU to neighbors || updating to neighbors
		{
			mu += BMU.neighborDists[i];
		}
		mu = mu / neighborNumber;
		mu = mu * mu;
		lambda = lambda * lambda;
		
		return Math.exp((-.5)*(d/(lambda*mu)));
	}
	//returns euclidean distance for comparing projection dimensionality similarity
	public double euclidean(double[] a, double[] b)
	{
		if(a.length == b.length)
		{
			double dist = 0;
			for(int i = 0; i<a.length; i++)
			{
				dist += ((a[i]-b[i])*(a[i]-b[i]));
			}
			return Math.sqrt(dist);
		}
		else
			return -1;
	}
	//Cosine similarity
	public double cosineSim(Prototype nope, Prototype dope)							//Cosine Similarity = (A  B)/(||A||*||B||)
	{																			//Higher values indicate more similarity (ranges from [-1, 1])
		double dot = 0;
		double magN = nope.pmag;
		double magD = dope.pmag;
		
		for(int i = 0; i < nope.P_hiDVec.length; i++)
		{
			dot += (nope.P_hiDVec[i]) * (dope.P_hiDVec[i]);
		}
		dot = dot/(magN*magD);
		if (dot> 1&& dot<1.001)
			dot=1;
		return dot;
	}

	public void writeFile()
	{
		BufferedWriter b;
		try 
		{
			String g = "ISO ("+epochs+", "+ neighborNumber+", " + protos.size() +").txt";
			System.out.println(g);
			
			FileWriter ff = new FileWriter(g,true);
			b = new BufferedWriter(ff);											//ff set in constructor = the file name of the SOMlander
			PrintWriter printer = new PrintWriter(b);
			printer.println("ISO");
			for(int i = 0; i<protos.size(); i++)
			{
				printer.print(protos.get(i).name + " [" + protos.get(i).m[0] + "; "+protos.get(i).m[1] + "]");	
				printer.print("\n");
			}
			printer.print("\n");
			printer.close();
		}
		catch (IOException e) 
		{
			e.printStackTrace();
			System.out.print("no way");
		}
	}
}
