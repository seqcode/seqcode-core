package edu.psu.compbio.seqcode.projects.kunz.isotop;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Scanner;

import javax.swing.JPanel;

import edu.psu.compbio.seqcode.projects.kunz.chromeSOM.MiniNode;


public class DrawIso extends JPanel
{
	public ArrayList<Prototype> protos; 
	public ArrayList<Color> colors;
	public int colorNum, minDataPoints, maxDataPoints;
	public double minX, maxX, minY, maxY, scalarX, scalarY;
	public boolean coded;
	public Color look; public int coder;
	public boolean weighting;
	public DrawIso(String file)
	{
		protos = new ArrayList<Prototype>();
		protos = read(file);
		minX = 0; maxX =0; minY = 0; maxY = 0; 
		coded = false;
		colorNum = 100;
		look = Color.BLUE;
	    colorCoder();
	}
	public ArrayList<Prototype> read(String file)
	{
		String ffs = System.getProperty("user.dir")+"/"+file;
		ArrayList<String> StringMat = new ArrayList<String>();
		try 
		{
			Scanner in = new Scanner(new FileReader(ffs));
			//System.out.println(xo + "  x  "+ yo);
			in.nextLine();
			while(in.hasNextLine())
			{
				StringMat.add(in.nextLine());
			}
		} catch (FileNotFoundException e) {e.printStackTrace();}
		for(String sr : StringMat)
		{
			if(sr.contains("["))
			{
				protos.add(new Prototype(sr));
			}
		}
		return protos;
	}
	public void search(String file)
	{
		for(int i = 0; i < protos.size(); i++)
		{
			Prototype mini = protos.get(i);
			mini.count = 0;
			mini.weight = 0;
		}
		ArrayList<String> strings = inputRead(file);
		for(String whole: strings)
		{
			int chr = 0;
			int locus1 = 0; int locus2 = 0;
			double weight = 0;
			//System.out.println(whole);
			if(whole.contains("\t") && whole.contains(":")&&whole.contains("chr")&&whole.contains("-"))
			{
				weighting = true;
				chr  = Integer.parseInt(whole.substring(whole.indexOf("chr")+3,whole.indexOf(":")));
				locus1 = Integer.parseInt(whole.substring(whole.indexOf(":")+1,whole.indexOf("-")));
				locus2 = Integer.parseInt(whole.substring(whole.indexOf("-")+1,whole.indexOf("\t")));
				int locus = (locus1 +locus2)/2;
				weight = Double.parseDouble(whole.substring(whole.indexOf("\t")+1,whole.length()));
				for(int i = 0; i< protos.size(); i++)
				{
					for(int j = 0; j< protos.get(i).loci.length; j++)
					{
						if(protos.get(i).loci[j][0] == chr && protos.get(i).loci[j][1] <= locus && protos.get(i).loci[j][2]>=locus)
						{
							protos.get(i).weight += weight;
						}
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
				for(int i = 0; i< protos.size(); i++)
				{
					for(int j = 0; j< protos.get(i).loci.length; j++)
					{
						if(protos.get(i).loci[j][0] == chr && protos.get(i).loci[j][1] <= locus && protos.get(i).loci[j][2]>=locus)
						{
							protos.get(i).count ++;
						}
					}
				}
			}
			
		}
		heatMapping();
		repaint();
		///needs a heat mapping method + colorbar;
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
	public void findGraphSpace()
	{
		for(int i = 0; i<protos.size(); i++)
		{
			//System.out.println(protos.get(i).name.substring(protos.get(i).name.indexOf("r")+1, protos.get(i).name.indexOf(":")) + "("+protos.get(i).m[0]+ ","+protos.get(i).m[1]+ ")");
			//x
			if(protos.get(i).m[0]<minX)
			{
				minX = protos.get(i).m[0];
			}
			if(protos.get(i).m[0]>maxX)
			{
				maxX = protos.get(i).m[0];
			}
			//y
			if(protos.get(i).m[1]<minY)
			{
				minY = protos.get(i).m[1];
			}
			if(protos.get(i).m[1]>maxY)
			{
				maxY = protos.get(i).m[1];
			}	
		}
		if(-1*minX > scalarX)
			scalarX = -1*minX;
		if(-1*minY> scalarY)
			scalarY = -1*minY;
		if(maxY>scalarY)
			scalarY = maxY;
		if(maxX> scalarX)
			scalarX = maxX;
		//System.out.println(scalar);
	}
	public void colorCoder()
	{
		colors = new ArrayList<Color>();
		for(int j = 0; j<16; j++)
		{
			Color neww  = new Color(((int)(Math.random()*255)),((int)(Math.random()*255)),((int)(Math.random()*255)));
			colors.add(neww);
		}
	}
	//Need to add a relative scale to grid
	//add chrome selector is in SOM
	public void paintComponent(Graphics g)
	 {
	    super.paintComponent(g);
	    colorBar(g);
	    setBackground(Color.BLACK);
	    int relativeX = getWidth()/2;
	    int relativeY = getHeight()/2;
	    for(int i = 1; i<=16; i++)
	    {
	    	g.setColor(Color.BLACK);
	    	g.drawString("Chrome " + i, getWidth()- 80, i*30);
	    	g.setColor(colors.get(i-1));
	    	g.fillRect(getWidth() - 100, i*30, 20, 20);
	    }
	    for(int i = 0; i<protos.size(); i++)
	    {
	    	Prototype phunk = protos.get(i);
	    	int centerX = ((int)((phunk.m[0]/scalarX)*relativeX))+relativeX;
	    	int centerY = ((int)((phunk.m[1]/scalarY)*relativeY))+relativeY;
	    	//System.out.println(protos.get(i).name + "    Goo" + i);
	    	int ip = (int)Integer.parseInt(protos.get(i).name.substring(protos.get(i).name.indexOf("r")+1, protos.get(i).name.indexOf(":")))-1;
    		g.setColor(colors.get(ip));
    		g.fillOval(centerX-6, centerY-6, 12, 12);
	    }
	    if(coded)
	    {
	    	for(int i = 0; i<protos.size(); i++)
		    {
		    	Prototype phunk = protos.get(i);
		    	int centerX = ((int)((phunk.m[0]/scalarX)*relativeX))+relativeX;
		    	int centerY = ((int)((phunk.m[1]/scalarY)*relativeY))+relativeY;
		    	
		    	int ip = 0;
		    	for(int j = 0; j < protos.get(i).names.length; j++)
		    	{
		    		ip = (int)Integer.parseInt(protos.get(i).names[j].substring(protos.get(i).names[j].indexOf("r")+1, protos.get(i).names[j].indexOf(":")));
		    		if(ip == coder)
			    		break;
		    	}
		    	if(ip == coder)
		    	{
		    		g.setColor(look);
		    		g.fillOval(centerX-8, centerY-8, 16, 16);
		    	}
		    }
	    }
	}
	
	public void countingDPS(int chr) 
	{
		if(chr == 0)
		{
			colorCoder();
			coded = false;
		}
		else if(chr == -1)
		{
			for(int i = 0; i < colors.size(); i++)
			{
				colors.set(i, Color.BLUE);
			}
		}
		else
		{
			coder = chr;
			coded = true;
			for(int i = 0; i < colors.size(); i++)
			{
				if(i+1 == chr)
					colors.set(i, look);
				else
					colors.set(i, Color.LIGHT_GRAY);
			}
		}
		repaint();
	}
	public void heatMapping()
	{
		
	}
	private void colorBar(Graphics g) 
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
}
