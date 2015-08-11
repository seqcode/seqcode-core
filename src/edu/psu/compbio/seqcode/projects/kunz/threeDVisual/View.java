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
	public double[][] coords, cool;
	public String file;
	public Timer timer;
	public int dir, count;
	public double minX, maxX, minY, maxY, minZ, maxZ, scalarX, scalarY, scalarZ, scalar;
	public boolean rotating;
	public View(String s)
	{
		file = s;
		dir = 0;
		minX = 0; maxX =0; minY = 0; maxY = 0; minZ = 0; maxZ = 0; scalarX = 0; scalarY = 0; scalarZ = 0; scalar = 0; 
		coords = reader(s);
		cool = coords;
		findGraphSpace();
		count = 0;
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
				minX = cool[i][0];
			}
			if(cool[i][0]>maxX)
			{
				maxX = cool[i][0];
			}
			//y
			if(cool[i][1]<minY)
			{
				minY = cool[i][1];
			}
			if(cool[i][1]>maxY)
			{
				maxY = cool[i][1];
			}	
			//z
			if(cool[i][2]<minZ)
			{
				minZ = cool[i][2];
			}
			if(cool[i][2]>maxZ)
			{
				maxZ = cool[i][2];
			}	
		}
		if(Math.abs(minX) > scalarX)
			scalarX = Math.abs(minX);
		if(maxX> scalarX)
			scalarX = maxX;
		if(Math.abs(minY)> scalarY)
			scalarY = Math.abs(minY);
		if(maxY>scalarY)
			scalarY = maxY;
		if(maxZ>scalarZ)
			scalarZ = maxZ;
		if(Math.abs(minZ)> scalarZ)
			scalarZ = Math.abs(minZ);
		
		scalar = scalarX;
		if(scalar < scalarZ)
			scalar = scalarZ;
		if(scalar < scalarY)
			scalar = scalarY;
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
		findGraphSpace();
	    setBackground(Color.GRAY);
	    g.setColor(Color.WHITE);
	    for(int i = 0; i < cool.length; i++)
	    {
	    	double shifter = .20 - (.15 * cool[i][2]/scalar);
	 	    double yScale = ((1.0-(2*shifter)) * getHeight())/2;
	     	double xScale = ((1.0-(2*shifter)) * getWidth())/2;
	     	double xShifter = shifter * getWidth();  
	     	double yShifter = shifter * getHeight();
	    	double c = (cool[i][2]/scalar);
	    	int col = (int) (127 + (c*127));
	    	g.setColor(new Color(0,col,0));
	    	int x = (int) (((cool[i][0]/scalar) * xScale) + xScale+xShifter);
	    	int y = (int) (((cool[i][1]/scalar) * yScale) + yScale+yShifter);
	    	int size = 6 + ((int) (2 * (cool[i][2]/scalar)));
	    	g.fillOval(x-size, y-size, size*2, size*2);
	    	if(i<cool.length-1)
	    	{
	    		double shifter2 = .20 - (.15 * cool[i+1][2]/scalar);
	    		double yScale2 = ((1.0-(2*shifter2)) * getHeight())/2;
		     	double xScale2 = ((1.0-(2*shifter2)) * getWidth())/2;
		     	double xShifter2 = shifter2 * getWidth();  
		     	double yShifter2 = shifter2 * getHeight();
	    		int x1 = (int) (((cool[i][0]/scalar) * xScale) + xScale+xShifter);
		    	int y1 = (int) (((cool[i][1]/scalar) * yScale) + yScale+yShifter);
		    	int x2 = (int) (((cool[i+1][0]/scalar) * xScale2) + xScale2+xShifter2);
		    	int y2 = (int) (((cool[i+1][1]/scalar) * yScale2) + yScale2+yShifter2);
		    	g.drawLine(x1, y1, x2, y2);
	    	}
	    }
	    cool = reader(file);
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
			for(int i = 0; i < cool.length; i++)
		    {
				//cool[i][0] *= 1;
				cool[i][1] = (cool[i][1]* Math.cos(rads)) - (cool[i][2] * Math.sin(rads));
				cool[i][2] = (cool[i][1]* Math.sin(rads)) + (cool[i][2] * Math.cos(rads));
		    }
		}
		/* for rotating around Y
		    | cos θ    0   sin θ| |x|   | x cos θ + z sin θ|   |x'|
		    |   0      1       0| |y| = |         y        | = |y'|
		    |-sin θ    0   cos θ| |z|   |-x sin θ + z cos θ|   |z'|
		 */
		if(way == 1)
		{
			for(int i = 0; i < cool.length; i++)
		    {
				cool[i][0] = (cool[i][0]* Math.cos(rads)) + (cool[i][2] * Math.sin(rads));
				//cool[i][1] *= 1;
				cool[i][2] = (-1*cool[i][0]* Math.sin(rads)) + (cool[i][2] * Math.cos(rads));
		    }
		}
		/* for rotating around Z
		    |cos θ   -sin θ   0| |x|   |x cos θ - y sin θ|   |x'|
    		|sin θ    cos θ   0| |y| = |x sin θ + y cos θ| = |y'|
    		|  0       0      1| |z|   |        z        |   |z'|
		 */
		if(way == 2)
		{
			for(int i = 0; i < cool.length; i++)
		    {
				cool[i][0] = (cool[i][0]* Math.cos(rads)) - (cool[i][1] * Math.sin(rads));
				cool[i][1] = (cool[i][0]* Math.sin(rads)) + (cool[i][1] * Math.cos(rads));
				//cool[i][2] *=1;
		    }
		}
	}
	public void actionPerformed(ActionEvent e) 
	{
		if(e.getSource()==timer)
		{
			if(rotating)
			{
				count++;
				double degree = 5;
				rotate(count * degree,dir);

				repaint();
				if(count*degree>=360)
					count = 0;
			}
		}
	}		
	
}
