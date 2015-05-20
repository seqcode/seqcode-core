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
	public BatchMap(int xNode, int yNode, int sigma)
	{
		String ffs = System.getProperty("user.dir")+"/SOMlander.txt";
		lander = ffs;
		xNodes = xNode;
		yNodes = yNode;
		sgm = sigma;
		//System.out.println(xNodes+","+yNodes);
		points = new ArrayList<DataPoint>();
		n = new NodeSystem(xNode,yNode);
		iterations = 1000;
		String ff = System.getProperty("user.dir")+"/MatrixLanding.txt";
		reader = new Reader(ff);
	}
	public void go()
	{
		generateDataPoints();
		initialize();
		iterater();
		finishTraining();
	}
	public void generateDataPoints()
	{
		reader.matrixReader();
		points = reader.points;
	}
	public void initialize() //random
	{
		for(int i = 0; i<n.size(); i++)
			n.get(i).initialize(points.get((int)(Math.random()*points.size())));
		assignNeighbors();
		neighborhooder();
		assignInitialNodes();
		System.out.println("Initialization Done");
	}
	public void assignInitialNodes()  // reimplement to check nodes more efficiently 
	{
		for(int i = 0; i<points.size(); i++)
			points.get(i).updateMag();
		for(int i = 0; i<n.size(); i++)
			n.get(i).updateMag();
		Node noodle; DataPoint doodle;
		for(int i = 0; i<points.size(); i++)
		{
			doodle = points.get(i);
			noodle = n.get((int)(Math.random()*n.size()));
			doodle.myNode = noodle;
			double cosineSim = cosineSim(noodle, doodle);
			Node compare;
			for(int j = 0; j < n.size(); j++)
			{
				compare = n.get(j);
				double dist = cosineSim(compare,doodle);
				if(dist>cosineSim)
				{
					cosineSim = dist;
					noodle = compare;
				}
			}
			doodle.myNode.dataPoints.remove(doodle);
			doodle.myNode = noodle;
			noodle.dataPoints.add(doodle);	
		}
	}
	public void assignNodes()
	{
		for(int i = 0; i<n.size(); i++)
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
	
	public double cosineSim(Node nope, DataPoint dope)
	{
		//cosine similarity
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
	
	public void iterater()
	{
		double dd = iterations;
		for(double i = 0; i < dd; i++)
		{
			iterate(i);
		}
	}
	public double nFactor(int o, double itercount) 
	{                  
		double sig = sgm-1;
		sig = (sig - (sig*(itercount/iterations)))+1; //sigma = width of Guassian kernel
		double r = 0;
		//System.out.println((1-(itercount/iterations)));
		r = (1-(itercount/iterations)) * Math.exp(-1*((o*o)/(2*sig*sig)));
		//System.out.println(o + "  " + r);
		return r;
	}
	public void iterate(double itercount)
	{
		for(int i = 0; i< n.size(); i++)
		{
			Node nodey = n.get(i);
			double[] numVec = new double[nodey.g.length]; for(int m = 0; m<nodey.g.length; m++) {numVec[m] = 0.0;}
			double sumDenom = 0;
			
			for(int o = 0; o < nodey.neighborHoods.size(); o++)
			{
				for(int k = 0; k<nodey.neighborHoods.get(o).size(); k++)
				{
					Node nono = nodey.neighborHoods.get(o).get(k);
					
					double nj = nono.dataPoints.size();
					double h = nFactor(o,itercount); //o = degree of separation between nodes (0-3)
					
					//double weight = nj*h;
					double weight = h;
					
					for(int j = 0; j<nj; j++)
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
			for(int l = 0; l<numVec.length; l++)
			{
				double doo = numVec[l]/sumDenom;
				numVec[l] = doo;
			}
			nodey.g = numVec;
		}
		count++;
		assignNodes();
	}
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
	public void neighborhooder()
	{
		for(int i = 0; i<n.size(); i++)
		{
			Node a = n.get(i);
			a.neighborHoods = new ArrayList<ArrayList<Node>>();
			a.neighborHoods.add(new ArrayList<Node>());
			a.neighborHoods.get(0).add(a);//adds self to the neighborhood 
			a.neighborHoods.add(a.neighbors);//adds 6 closest neighbors as found by neighborhood method above
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
						if(!a.neighborHoodContains(add)) //checks Node a's neighborHoods to see if Node add is already present
						{
							a.neighborHoods.get(shell+1).add(add);
						}
					}
				}
			}
			shell++;
		}
		/*
		 * For Diagnosing neighborhood issues
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
	public void finishTraining()
	{
		writeFile();
		System.out.println("Training done");
	}
	public void writeFile()
	{
		BufferedWriter b;
		try 
		{
			FileWriter ff = new FileWriter(lander,true);
			b = new BufferedWriter(ff);
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
