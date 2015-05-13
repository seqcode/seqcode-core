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
	public BatchMap(int xNode, int yNode)
	{
		lander = "src/kunzSOMstuff/SOMlander.txt";
		xNodes = xNode;
		yNodes = yNode;
		System.out.println(xNodes+","+yNodes);
		points = new ArrayList<DataPoint>();
		//reader = new Reader("src/kunzSOMstuff/cities1000_dist-matrix.txt", "src/kunzSOMstuff/Look Up.txt");
		reader = new Reader("src/kunzSOMstuff/MatrixLanding.txt");
		n = new NodeSystem(xNode,yNode);
		iterations = 1000;
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
		//System.out.println(dot);
		return dot;
	}
	public double cosineSim(Node nope, Node dope)
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
		//System.out.println(dot);
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
		//System.out.println(o);
		sgm = 175;
		double sig = sgm - (sgm*(itercount/iterations)); //sigma = width of Guassian kernel... what is this supposed to start at?
		
		double r = 0;
		r = (iterations-itercount) * Math.exp(-1*((o*o)/(2*sig*sig)));
		/*if(count % 100 == 0)
		{
			System.out.println(o + " - " + r + " - "+Math.exp(-1*((o*o)/(2*sig*sig))));
		}*/
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
		neighborhood2();
	}
	public void neighborhood2()
	{
		Node b; Node a;
		for(int i = 0; i<n.size(); i++)
		{
			a = n.get(i);
			for(int k = 0; k<a.neighbors.size();k++)
			{
				b = a.neighbors.get(k);
				for(int j = 0; j < b.neighbors.size(); j++)
				{
					if(a.neighbors.indexOf(b.neighbors.get(j))==-1&&b.neighbors.get(j)!=a)
					{
						if(a.neighbors2.indexOf(b.neighbors.get(j))==-1)
							a.neighbors2.add(b.neighbors.get(j));
					}
				}
			}
		}
		neighborhood3();
	}
	public void neighborhood3()
	{
		Node b; Node a;
		for(int i = 0; i<n.size(); i++)
		{
			a = n.get(i);
			for(int k = 0; k<a.neighbors2.size();k++)
			{
				b = a.neighbors2.get(k);
				for(int j = 0; j < b.neighbors.size(); j++)
				{
					if(a.neighbors2.indexOf(b.neighbors.get(j))==-1 && a.neighbors.indexOf(b.neighbors.get(j))==-1 && b.neighbors.get(j)!=a)
					{
						if(a.neighbors3.indexOf(b.neighbors.get(j))==-1)
							a.neighbors3.add(b.neighbors.get(j));
					}
				}
			}
		}
		assignNeighborHoods();
	}
	public void assignNeighborHoods()
	{
		for(int i = 0; i<n.size(); i++)
		{
			n.get(i).neighborHoods();
		}
	}
	public void finishTraining()
	{
		System.out.println("Training done");
		writeFile();
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
	
	//To be cut and moved/redone in draw hex
	public void findCountedDPS(ArrayList<String> theWord)
	{
		for(int i = 0; i<n.size(); i++)
		{
			n.get(i).countedPoints.clear();
			for(int j = 0; j < n.get(i).dataPoints.size(); j++)
			{
				if(theWord.contains(n.get(i).dataPoints.get(j).chrome))
				{
					n.get(i).countedPoints.add(n.get(i).dataPoints.get(j));
				}
			}
		}
	}
}
