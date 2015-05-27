package edu.psu.compbio.seqcode.projects.kunz.chromeSOM;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;


public class BatchMap
{
	public NodeSystem n;
	public Reader reader;
	public int xNodes, yNodes, iterations, count;
	public double sgm;
	public ArrayList<DataPoint> points;
	public String lander;
	public BatchMap(int xNode, int yNode, int sigma,int it)
	{
		xNodes = xNode;
		yNodes = yNode;
		sgm = sigma;
		//System.out.println(xNodes+","+yNodes);
		points = new ArrayList<DataPoint>();
		n = new NodeSystem(xNode,yNode);
		iterations = it;
		String ff = System.getProperty("user.dir")+"/MatrixLanding.txt";
		reader = new Reader(ff);
	}
	
	//initialize -> iterate -> terminate
	public void go() 
	{
		generateDataPoints();
		initialize();
		iterater();
		finishTraining();
	}
	
	//Reader class used to generate data points from file 
	public void generateDataPoints() 
	{
		reader.matrixReader(); 									//reader object and file location set in constructor
		points = reader.points;
	}
	
	//Node vectors are each set equal to a randomly selected data point's vector - Also calls for neighbor assignment
	public void initialize() 
	{
		for(int i = 0; i<n.size(); i++)
			n.get(i).initialize(points.get((int)(Math.random()*points.size())));
		assignNeighbors();
		neighborhooder();
		assignInitialNodes();
		System.out.println("Initialization Done");
	}
	
	//Assigns data points to nodes based on cosine similarity metric (cosineSim() method)
	public void assignInitialNodes()  
	{
		for(int i = 0; i<points.size(); i++) 					//Data Points and Nodes have a magnitude field that refers to their vectors' magnitudes
			points.get(i).updateMag();							//Magnitude is found and saved because they are used so many times over in cosine similarity method
		for(int i = 0; i<n.size(); i++)
			n.get(i).updateMag();
		Node noodle; DataPoint doodle;
		for(int i = 0; i<points.size(); i++)					//Goes through each data points
		{
			doodle = points.get(i);
			noodle = n.get((int)(Math.random()*n.size())); 		//initial node is randomly selected
			doodle.myNode = noodle;
			double cosineSim = cosineSim(noodle, doodle);
			Node compare;
			for(int j = 0; j < n.size(); j++) 					//Each data point goes through each Node looking for most similar
			{
				compare = n.get(j);
				double dist = cosineSim(compare,doodle);
				if(dist>cosineSim)  							//Node with highest cosine similarity value to data point = most similar
				{
					cosineSim = dist;
					noodle = compare;
				}
			}
			doodle.myNode.dataPoints.remove(doodle); 			//keeping track of which data points and Nodes are assigned to eachother
			doodle.myNode = noodle;
			noodle.dataPoints.add(doodle);	
		}
	}
	
	//Assigns data points to nodes based on cosine similarity metric (cosineSim() method)
	public void assignNodes()
	{
		for(int i = 0; i<n.size(); i++) 						//Data point vectors do not change -> only Node update is needed aside from initial  
			n.get(i).updateMag();
		for(int i = 0; i<points.size(); i++)
		{
			DataPoint doodle = points.get(i);
			Node noodle = n.get((int)(Math.random()*n.size()));
			double cosineSim = cosineSim(noodle, doodle);
			Node compare;
			for(int j = 0; j < n.size(); j++)
			{
				compare = n.get(j);
				double sim = cosineSim(compare,doodle);
				if(sim>cosineSim)
				{
					cosineSim = sim;
					noodle = compare;
				}
			}
			doodle.myNode.dataPoints.remove(doodle);
			doodle.myNode = noodle;
			noodle.dataPoints.add(doodle);
		}
	}
	
	//Cosine similarity
	public double cosineSim(Node nope, DataPoint dope)							//Cosine Similarity = (A · B)/(||A||*||B||)
	{																			//Higher values indicate more similarity (ranges from [-1, 1])
		double dot = 0;
		double magN = nope.mag;
		double magD = dope.mag;
		
		for(int i = 0; i < nope.g.length; i++)
		{
			dot += (nope.g[i]) * (dope.g[i]);
		}
		dot = dot/(magN*magD);
		return dot;
	}
	
	//Calls iterate methods the number of times dictated by the parameter
	public void iterater()
	{
		double dd = iterations;
		for(double i = 0; i < dd; i++)
		{
			iterate(i);
		}
	}
	
	//Finds the weight factor of a given data point based on learning rate and neighborhood function
	public double nFactor(int o, double itercount) 
	{                  
		double sig = sgm-1;
		sig = (sig - (sig*(itercount/iterations)))+1; 							//sigma = width of Guassian kernel (sgm ->1)
		double r = 0;
		double lr = (1-(itercount/iterations)); 								//lr = learning rate (1 -> 0) 
		r = lr * Math.exp(-1*((o*o)/(2*sig*sig)));
		return r;
	}
	
	//Where the magic happens
	public void iterate(double itercount)
	{
		for(int i = 0; i< n.size(); i++) 										//Updates each Node; neighborhood shells -> nodes -> data points
		{
			Node nodey = n.get(i);
			double[] numVec = new double[nodey.g.length]; for(int m = 0; m<nodey.g.length; m++) {numVec[m] = 0.0;}  //initialize a matrix for sum of all datapoints*weight for weighted average
			double sumDenom = 0;																					//SumDenom = sum of weight factors used as denominator in weighted average
			
			for(int o = 0; o < nodey.neighborHoods.size(); o++)        			//Cycles through neighborhood shells, starting with most similar
			{
				for(int k = 0; k<nodey.neighborHoods.get(o).size(); k++)		//Cycles through Nodes in each shell
				{
					Node nono = nodey.neighborHoods.get(o).get(k);
					
					double nj = nono.dataPoints.size();
					double h = nFactor(o,itercount); 							//o = degree of separation between nodes = how many shells apart					
					//double weight = nj*h;
					double weight = h;
					
					for(int j = 0; j<nj; j++)									//Cycles through the data points assigned to each Node in each shell
					{
						DataPoint p = nono.dataPoints.get(j);
						for(int l = 0; l<p.g.length; l++)
						{
							double doo = ((p.g[l])*weight)+(numVec[l]);
							numVec[l] = doo;
						}
						sumDenom += weight;
					}
				}
			}
			for(int l = 0; l<numVec.length; l++) 						//Takes weighted average = updated vector for node
			{
				double doo = numVec[l]/sumDenom;
				numVec[l] = doo;
			}
			nodey.g = numVec;									//Node vectors are updated in order, but reassignment happens after all are updated -> batched SOM
		}
		count++;
		assignNodes();
	}
	
	//Assigns the first shell of neighbors - hard coded to follow hex grid rules
	public void assignNeighbors() 								
	{
		for(int i = 0; i<xNodes; i++)
		{
			for(int j = 0; j<yNodes; j++)
			{
				Node a = n.get(i,j);
				if(j%2==0)
				{
					a.neighbors.add(n.get(i-1,j-1));
					a.neighbors.add(n.get(i,j-1));
					
					a.neighbors.add(n.get(i-1,j));
					a.neighbors.add(n.get(i+1,j));
					
					a.neighbors.add(n.get(i-1,j+1));
					a.neighbors.add(n.get(i,j+1));
				}
				else
				{
					a.neighbors.add(n.get(i,j-1));
					a.neighbors.add(n.get(i+1,j-1));
					
					a.neighbors.add(n.get(i-1,j));
					a.neighbors.add(n.get(i+1,j));
					
					a.neighbors.add(n.get(i,j+1));
					a.neighbors.add(n.get(i+1,j+1));
				}
			}
		}
	}
	
	//Each node has every other node (including self) assigned to a shell based on the initial shell from assignNeighbors() method
	public void neighborhooder()						
	{
		for(int i = 0; i<n.size(); i++)
		{
			Node a = n.get(i);
			a.neighborHoods = new ArrayList<ArrayList<Node>>();
			a.neighborHoods.add(new ArrayList<Node>());
			a.neighborHoods.get(0).add(a);							//adds self to the neighborhood neighborhood[0]
			a.neighborHoods.add(a.neighbors);						//adds 6 closest neighbors as found by assignNeighbors() method to neighborhood[1]
		}
		int shell = 1;
		int maxShell = (xNodes + yNodes)/4;
		while(shell < maxShell)
		{
			for(int i = 0; i<n.size(); i++)
			{
				Node a = n.get(i);
				a.neighborHoods.add(new ArrayList<Node>());
				for(int k = 0; k<a.neighborHoods.get(shell).size();k++)
				{
					Node b = a.neighborHoods.get(shell).get(k);
					for(int j = 0; j < b.neighborHoods.get(shell).size(); j++)
					{
						Node add = b.neighborHoods.get(shell).get(j);
						if(!a.neighborHoodContains(add)) 							//checks Node a's neighborHoods to see if Node add is already present
						{
							a.neighborHoods.get(shell+1).add(add);
						}
					}
				}
			}
			shell++;
		}
		/*
		 	For Diagnosing neighborhood issues
		 */
		/*Node bb = n.get(0);
		int p = 0; 
		for(int i = 0; i< bb.neighborHoods.size(); i++)
		{
			System.out.println(i + " " +bb.neighborHoods.get(i).size());
			for(int k = 0; k<bb.neighborHoods.get(i).size(); k++)
			{
				p++;
			}
		}
		System.out.println(p);
		*/
	}
	
	//Called by go() after iterations are complete - not very useful, just calls another methods and outprints
	public void finishTraining()       
	{
		writeFile();
		System.out.println("Training done");
	}
	
	//Saves trained SOM in text file
	public void writeFile()
	{
		BufferedWriter b;
		try 
		{
			String g = "SOM ("+xNodes+"x"+yNodes+"), "+(int)sgm+", " + iterations +".txt";
			System.out.println(g);
			
			FileWriter ff = new FileWriter(g,true);
			b = new BufferedWriter(ff);											//ff set in constructor = the file name of the SOMlander
			PrintWriter printer = new PrintWriter(b);
			printer.print(xNodes+"x"+yNodes+ "\nSigma:" + sgm +"\n");
			for(int i = 0; i<yNodes; i++)
			{
				for(int j = 0; j<xNodes; j++)
				{
					Node node = n.get(j,i);
					printer.print("["+j +","+i+"]\t"+"(");
					for(int o = 0; o<node.dataPoints.size(); o++)
					{
						printer.print(node.dataPoints.get(o).name);
						if(o<node.dataPoints.size()-1)
							printer.print(", ");
					}
					printer.print(")" + "\t");
				}
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
