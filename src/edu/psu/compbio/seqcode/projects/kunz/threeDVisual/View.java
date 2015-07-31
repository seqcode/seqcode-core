package edu.psu.compbio.seqcode.projects.kunz.threeDVisual;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Scanner;

import javax.swing.JPanel;
import javax.swing.Timer;

public class View extends JPanel implements ActionListener
{
	public double[][] coords;
	public Timer timer;
	public int dir;
	public double minX, maxX, minY, maxY, minZ, maxZ, scalarX, scalarY, scalarZ;
	public boolean rotating;
	public View(String s)
	{
		dir = 0;
		minX = 0; maxX =0; minY = 0; maxY = 0; minZ = 0; maxZ = 0; scalarX = 0; scalarY = 0; scalarZ = 0; 
		coords = reader(s);
		findGraphSpace();
		rotating = false;
		timer = new Timer(100,this);
		timer.start();
	}
	public void findGraphSpace()
	{
		minX = 0; maxX =0; minY = 0; maxY = 0; minZ = 0; maxZ = 0; scalarX = 0; scalarY = 0; scalarZ = 0; 
		for(int i = 0; i<coords.length; i++)
		{
			//System.out.println(protos.get(i).name.substring(protos.get(i).name.indexOf("r")+1, protos.get(i).name.indexOf(":")) + "("+protos.get(i).m[0]+ ","+protos.get(i).m[1]+ ")");
			//x
			if(coords[i][0]<minX)
			{
				minX = coords[i][0];
			}
			if(coords[i][0]>maxX)
			{
				maxX = coords[i][0];
			}
			//y
			if(coords[i][1]<minY)
			{
				minY = coords[i][1];
			}
			if(coords[i][1]>maxY)
			{
				maxY = coords[i][1];
			}	
			//z
			if(coords[i][2]<minZ)
			{
				minZ = coords[i][2];
			}
			if(coords[i][2]>maxZ)
			{
				maxZ = coords[i][2];
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
		if(maxZ>scalarZ)
			scalarZ = maxZ;
		if(-1*minZ> scalarZ)
			scalarZ = -1*minZ;
		//System.out.println(scalar);
	}
	public double[][] reader(String s)
	{
		ArrayList<String> StringMat = new ArrayList<String>();
		try 
		{
			Scanner in = new Scanner(new FileReader(s));
			//System.out.println(xo + "  x  "+ yo);
			while(in.hasNextLine())
			{
				StringMat.add(in.nextLine());
			}
		} catch (FileNotFoundException e) {e.printStackTrace();}
		return pointParse(StringMat);
	}
	public double[][] pointParse(ArrayList<String> s)
	{
		double[][] goof = new double[s.size()][3];
		for(int i = 0; i < s.size(); i++)
		{
			String r = s.get(i);
			goof[i][0] = Double.parseDouble(r.substring(0,r.indexOf(",")));
			r = r.substring(r.indexOf(",")+1);
			goof[i][1] = Double.parseDouble(r.substring(0,r.indexOf(",")));
			r = r.substring(r.indexOf(",")+1);
			goof[i][2] = Double.parseDouble(r.substring(0,r.length()));
		}
		
		return goof;
	}
	public void paintComponent(Graphics g)
	{
		super.paintComponent(g);
		rotate(5,dir);
		findGraphSpace();
	    setBackground(Color.GRAY);
	    g.setColor(Color.WHITE);
	    //set a constant scale to maintain "size"
    	
	    for(int i = 0; i < coords.length; i++)
	    {
	    	double shifter = .20 - (.15 * coords[i][2]/scalarZ);
	 	    double yScale = ((1.0-(2*shifter)) * getHeight())/2;
	     	double xScale = ((1.0-(2*shifter)) * getWidth())/2;
	     	double xShifter = shifter * getWidth();  
	     	double yShifter = shifter * getHeight();
	    	double c = (coords[i][2]/scalarZ);
	    	int col = (int) (127 + (c*127));
	    	
	    	g.setColor(new Color(0,col,0));
	    	int x = (int) (((coords[i][0]/scalarX) * xScale) + xScale+xShifter);
	    	int y = (int) (((coords[i][1]/scalarY) * yScale) + yScale+yShifter);
	    	g.fillOval(x-6, y-6, 12, 12);
	    	if(i<coords.length-1)
	    	{
	    		double shifter2 = .20 - (.15 * coords[i+1][2]/scalarZ);
	    		double yScale2 = ((1.0-(2*shifter2)) * getHeight())/2;
		     	double xScale2 = ((1.0-(2*shifter2)) * getWidth())/2;
		     	double xShifter2 = shifter2 * getWidth();  
		     	double yShifter2 = shifter2 * getHeight();
	    		int x1 = (int) (((coords[i][0]/scalarX) * xScale) + xScale+xShifter);
		    	int y1 = (int) (((coords[i][1]/scalarY) * yScale) + yScale+yShifter);
		    	int x2 = (int) (((coords[i+1][0]/scalarX) * xScale2) + xScale2+xShifter2);
		    	int y2 = (int) (((coords[i+1][1]/scalarY) * yScale2) + yScale2+yShifter2);
		    	g.drawLine(x1, y1, x2, y2);
	    	}
	    }
	}
	public void rotate(double degree, int way)
	{
		double rads = degree * .0174532925;
		/*	for rotating around x axis
		 *  |1     0           0| |x|   |        x        |   |x'|
    		|0   cos θ    -sin θ| |y| = |y cos θ - z sin θ| = |y'|
    		|0   sin θ     cos θ| |z|   |y sin θ + z cos θ|   |z'|
		 */
		if(way == 0)
		{
			for(int i = 0; i < coords.length; i++)
		    {
				coords[i][0] *= 1;
				coords[i][1] = (coords[i][1]* Math.cos(rads)) - (coords[i][2] * Math.sin(rads));
				coords[i][2] = (coords[i][1]* Math.sin(rads)) + (coords[i][2] * Math.cos(rads));
		    }
		}
		/* for rotating around Y
		    | cos θ    0   sin θ| |x|   | x cos θ + z sin θ|   |x'|
		    |   0      1       0| |y| = |         y        | = |y'|
		    |-sin θ    0   cos θ| |z|   |-x sin θ + z cos θ|   |z'|
		 */
		if(way == 1)
		{
			for(int i = 0; i < coords.length; i++)
		    {
				coords[i][0] = (coords[i][0]* Math.cos(rads)) + (coords[i][2] * Math.sin(rads));;
				coords[i][1] *= 1;
				coords[i][2] = (-1*coords[i][0]* Math.sin(rads)) + (coords[i][2] * Math.cos(rads));
		    }
		}
		/* for rotating around Z
		    |cos θ   -sin θ   0| |x|   |x cos θ - y sin θ|   |x'|
    		|sin θ    cos θ   0| |y| = |x sin θ + y cos θ| = |y'|
    		|  0       0      1| |z|   |        z        |   |z'|
		 */
		if(way == 2)
		{
			for(int i = 0; i < coords.length; i++)
		    {
				coords[i][0] = (coords[i][0]* Math.cos(rads)) - (coords[i][1] * Math.sin(rads));
				coords[i][1] = (coords[i][0]* Math.sin(rads)) + (coords[i][1] * Math.cos(rads));
				coords[i][2] *=1;
		    }
		}
	}
	public void actionPerformed(ActionEvent e) 
	{
		if(e.getSource()==timer)
		{
			if(rotating)
			{
				repaint();
			}
		}
	}		
	
}
