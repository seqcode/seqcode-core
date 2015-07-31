package edu.psu.compbio.seqcode.projects.kunz.chromeSOM;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Scanner;

import javax.swing.JPanel;

//fix heat mapping?
//make heat mapping option between overall point distribution & normal vs expected (heated true vs false)
public class DrawHex extends JPanel// implements Grid2D
{
	private static final long serialVersionUID = 1L;
	public int nodes,xs,ys,xNodes,yNodes, colorNum, winW, winH;
	public double maxDataPoints, minDataPoints, totalWeight, counter;
	public ArrayList<Node> coords;
	public ArrayList<DataPoint> dataPoints;
	public ArrayList<MiniNode> nodeList;
	public ArrayList<Color> colors;
	MiniSystem nodeSystem;
	public boolean weighting, heat;
	public DrawHex(String s)
	{
		weighting = false;
		String ffs = System.getProperty("user.dir")+"/"+s;
		yNodes=0;
		xNodes=0;
		reader(ffs);
		nodes = yNodes*xNodes;
		heat = true;
		
    	colors = new ArrayList<Color>();
    	colorNum=100;
	}
	public void search(String file)
	{

		totalWeight = 0;
		for(int i = 0; i < nodeSystem.size(); i++)
		{
			MiniNode mini = nodeSystem.get(i);
			mini.counting.clear();
			mini.weight = 0;
		}
		ArrayList<String> strings = inputRead(file);
		for(String whole: strings)
		{
			
			int chr = 0;
			int locus1 = 0; int locus2 = 0;
			double weight = 0;
			counter = 0;
			//System.out.println(whole);
			if(whole.contains("\t") && whole.contains(":")&&whole.contains("chr")&&whole.contains("-"))
			{
				weighting = true;
				chr  = Integer.parseInt(whole.substring(whole.indexOf("chr")+3,whole.indexOf(":")));
				locus1 = Integer.parseInt(whole.substring(whole.indexOf(":")+1,whole.indexOf("-")));
				locus2 = Integer.parseInt(whole.substring(whole.indexOf("-")+1,whole.indexOf("\t")));
				int locus = (locus1 +locus2)/2;
				weight = Double.parseDouble(whole.substring(whole.indexOf("\t")+1,whole.length()));
				for(int i = 0; i< dataPoints.size(); i++)
				{
					if(dataPoints.get(i).chrome == chr && dataPoints.get(i).minLocus <= locus && dataPoints.get(i).maxLocus>=locus)// && weight > 14)
					{
						//System.out.println(dataPoints.get(i).name + " = " + chr + ":" + locus +  "     " + weight);
						dataPoints.get(i).myMini.counting.add(dataPoints.get(i));
						dataPoints.get(i).myMini.weight += weight;
						totalWeight += weight;
					}
				}
			}
			else if (whole.contains(":")&&whole.contains("chr")&&whole.contains("-"))
			{
				weighting = false;
				chr  = Integer.parseInt(whole.substring(whole.indexOf("chr")+3,whole.indexOf(":")));
				locus1 = Integer.parseInt(whole.substring(whole.indexOf(":")+1,whole.indexOf("-")));
				locus2 = Integer.parseInt(whole.substring(whole.indexOf("-")+1,whole.length()));
				int locus = (locus1 +locus2)/2;
				for(int i = 0; i< dataPoints.size(); i++)
				{
					if(dataPoints.get(i).chrome == chr && dataPoints.get(i).minLocus <= locus && dataPoints.get(i).maxLocus>=locus)
					{
						counter ++;
						dataPoints.get(i).myMini.counting.add(dataPoints.get(i));
					}
				}
			}
			
		}
		for(int i = 0; i < nodeSystem.size(); i++)
		{
			MiniNode mini = nodeSystem.get(i);
			mini.notCounting.addAll(mini.dataPoints);
			mini.notCounting.removeAll(mini.counting);
		}
		heatMapping();
		repaint();
		
	}
	public ArrayList<String> inputRead(String file)
	{
		ArrayList<String> StringMat = new ArrayList<String>();
		try 
		{
			Scanner in = new Scanner(new FileReader(file));
			String sizer = in.next();
			//System.out.println(xo + "  x  "+ yo);
			in.next();
			while(in.hasNextLine())
			{
				StringMat.add(in.nextLine());
			}
		} catch (FileNotFoundException e) {e.printStackTrace();}
		return StringMat;

	}
	public void countingDPS(int chr)
	{
		counter = 0;
		weighting = false;
		for(int i =0; i<nodeSystem.size();i++)
		{
			MiniNode mini = nodeSystem.get(i);
			mini.counting.clear();
			for(int j =0; j<mini.dataPoints.size(); j++)
			{
				if(chr <= -1)
				{
					if(Math.random()>.9)
					{
						mini.counting.add(mini.dataPoints.get(j));
						counter++;
					}
				}
				else if (chr == 0)
				{
					mini.counting.add(mini.dataPoints.get(j));
					counter++;
				}
				else if (mini.dataPoints.get(j).chrome == chr)
				{
					mini.counting.add(mini.dataPoints.get(j));
					counter++;
				}
			}
			//System.out.println(mini.counting.size());
		}
	    maxDataPoints = 0;
	    minDataPoints = 0;
		heatMapping();
		repaint();
	}
	public void reader(String f)
	{
		ArrayList<String> StringMat = new ArrayList<String>();
		try 
		{
			Scanner in = new Scanner(new FileReader(f));
			String sizer = in.next();
			int xo = Integer.parseInt(sizer.substring(0,sizer.indexOf("x")));
			int yo = Integer.parseInt(sizer.substring(sizer.indexOf("x")+1,sizer.length()));
			//System.out.println(xo + "  x  "+ yo);
			in.next();
			while(in.hasNext())
			{
				StringMat.add(in.next());
			}
			createNoder(StringMat, xo, yo);
		} catch (FileNotFoundException e) {e.printStackTrace();}
	}
	public void createNoder(ArrayList<String> strings, int x, int y)
	{
		nodeList = new ArrayList<MiniNode>(x*y);
		dataPoints = new ArrayList<DataPoint>();
		MiniNode current = null;
		xNodes = x;
		yNodes = y;
		for(String whole: strings)
		{
			if(whole.contains("["))
			{
				current = new MiniNode(whole);
				nodeList.add(current);
			}
			else
				current.addDP(whole);
		}
		//Data point management
		for(int i = 0; i<nodeList.size(); i++)
		{
			for(int j = 0; j < nodeList.get(i).dataPoints.size(); j++)
			{
				dataPoints.add(nodeList.get(i).dataPoints.get(j));
			}
		}
		
		//System.out.println(dataPoints.size());
		nodeSystem = new MiniSystem(nodeList,xNodes,yNodes);
	}
	public void paintComponent(Graphics g)
	  {
	    super.paintComponent(g);
	    colorBar(g);
	    setBackground(Color.GRAY);
	    for(int i = 0; i<nodeList.size(); i++)
	    {
	    	g.setColor(nodeList.get(i).color);
	    	g.fillPolygon(nodeList.get(i).p);
	    	g.setColor(Color.BLACK);
	    	g.setFont(new Font("Serif", 5, 9));
	    	g.drawString(""+nodeList.get(i).counting.size(),nodeList.get(i).xLoc, nodeList.get(i).yLoc);
	    	g.drawPolygon(nodeList.get(i).p);
	    }
	}
	public void nodeBuild(int winW, int winH)
	{  
		int xMin = (int)(.1 * winW);
	    int yMin = (int)(.05 * winH);
	    int xMax = (int)(.90 * winW);
	    int yMax = (int)(.85 * winH);
		
	    double winWidth = xMax-xMin;
	    double winHeight = yMax-yMin;
	    
	    double height = winHeight/yNodes;
	    double width = winWidth/xNodes;
	    height = height/3;
	    width = width/2;
	    int x = (int) (width);
	    int y = (int) (height);
	    int w = (int) (winWidth/width);
	    int h = (int)(winHeight/height);
	    for(int i = 0; i<h; i++)
		{
	    	for(int m = 0; m<w; m++)
			{
				if((i%6==1 && m%2==0)||(i%6==4 && m%2 == 1))
				{
					MiniNode p = nodeSystem.get(m/2,i/3);
					p.xLoc = m*x+xMin;
					p.yLoc = i*y+yMin;
				}
			}
		}
	   for(int c = 0; c<nodeList.size(); c++)
	   {
		  MiniNode m = nodeList.get(c);
		  
		  int[] xPoints = new int[6];
		  int[] yPoints = new int[6];
		  
		  xPoints[0] = m.xLoc;
		  yPoints[0] = m.yLoc-(2*y);
		  
		  xPoints[1] = m.xLoc+x;
		  yPoints[1] = m.yLoc-y;
		  
		  xPoints[2] = m.xLoc+x;
		  yPoints[2] = m.yLoc+y;
			
		  xPoints[3] = m.xLoc;
		  yPoints[3] = m.yLoc+(2*y);

		  xPoints[4] = m.xLoc-x;
		  yPoints[4] = m.yLoc+y;
			
		  xPoints[5] = m.xLoc-x;
		  yPoints[5] = m.yLoc-y;
		  
		  m.polygonMaker(xPoints, yPoints, 6);
	   }
	}
	public void colors()
	{
		for(int k = 0; k<colorNum;k++)
		{
			int blue;
			int red;
			if(k<colorNum/2)
			{
				red= (int)(k*255/colorNum)*2; 
				if(red >255) 
					red = 255-(red-255);
			}
			else
				red = 255;
			if(k<colorNum/2)
				blue = 255;
			else
			{
				blue = (int)(255-(k*255/colorNum))*2;
				if(blue >255)
					blue = 255-(blue-255);
			}
			int green = (int)(255-(k*255/colorNum))*2;if(green >255) green = 255-(green-255);
			
			Color c = new Color(red, green, blue);
			colors.add(c);
		}
	}
	public void colorBar(Graphics g)
	{
		g.drawRect((int)(getWidth()*.2), (int)(getHeight()*.96), (int)(getWidth()*.55), (int)(getHeight()*.02));
		g.setColor(Color.BLACK);
		g.drawString(""+minDataPoints, (int)(getWidth()*.17),(int)(getHeight()*.95));
		g.drawString(""+maxDataPoints, (int)(getWidth()*.2+(int)(getWidth()*.55)),(int)(getHeight()*.95));
		for(int k = 0; k<colorNum;k++)
		{
			int blue;
			int red;
			if(k<colorNum/2)
			{
				red= (int)(k*255/colorNum)*2; 
				if(red >255) 
					red = 255-(red-255);
			}
			else
				red = 255;
			if(k<colorNum/2)
				blue = 255;
			else
			{
				blue = (int)(255-(k*255/colorNum))*2;
				if(blue >255)
					blue = 255-(blue-255);
			}
			int green = (int)(255-(k*255/colorNum))*2;if(green >255) green = 255-(green-255);
			
			Color c = new Color(red, green, blue);
			g.setColor(c);
			g.fillRect((int)(getWidth()*.2+((k*getWidth()*.55)/colorNum)),(int)(getHeight()*.96),(int)(getWidth()*.55/colorNum)+3,(int)(getHeight()*.02));
		}
		
	}
	public void heatMapping()
	{
		
		//Is this coloring scheme by ratio a good idea? Does this properly compare actual to expected?
		//Switch dataPoints.size() to nonCounted.size() & dataPoints.size()-counter
		if(weighting == true)
		{
			
			minDataPoints = 0;
			maxDataPoints = 0;
			for(int i = 0; i<nodeList.size(); i++)
			{
				
				double scalar = ((double)nodeList.get(i).dataPoints.size())/((double)dataPoints.size()); //compares ratios
				double w = nodeList.get(i).weight/totalWeight;
				w = w - scalar;
				if(w>maxDataPoints)
					maxDataPoints = w;
				if(w<minDataPoints)
					minDataPoints =  w;
			}
			if(maxDataPoints <= minDataPoints) maxDataPoints = minDataPoints + 1;
			for(int i = 0; i<nodeSystem.size(); i++)
			{	   
				MiniNode p = nodeSystem.get(i);
				double scalar = ((double)p.dataPoints.size())/((double)dataPoints.size()); //compares ratios
				double w = p.weight/totalWeight;
				if(scalar == 0) 
					scalar = .000001;
				w = w - scalar;
				p.color = colors.get((int) ((w-minDataPoints)*(colors.size()-1)/(maxDataPoints-minDataPoints)));
			}
		}
		else
		{
			minDataPoints = 0;
			maxDataPoints = 0;

			for(int i = 0; i<nodeList.size(); i++)
			{
				double scalar = ((double)nodeList.get(i).dataPoints.size())/((double)dataPoints.size()); //compares ratios
				double w = nodeList.get(i).counting.size()/counter;
				w = w - scalar;
				if(w>maxDataPoints)
					maxDataPoints = w;
				if(w<minDataPoints)
					minDataPoints = w;
			}
			if(maxDataPoints == minDataPoints) maxDataPoints++;
			for(int i = 0; i<nodeSystem.size(); i++)
			{	

				MiniNode p = nodeSystem.get(i);
				double scalar = ((double)p.dataPoints.size())/((double)dataPoints.size()); if(scalar == 0 ) scalar =.000001; //compares ratios
				double w = p.counting.size()/counter;
				w = w - scalar;
				p.color = colors.get((int) ((w-minDataPoints)*(colors.size()-1)/(maxDataPoints-minDataPoints)));
			}
		}
	}
}
