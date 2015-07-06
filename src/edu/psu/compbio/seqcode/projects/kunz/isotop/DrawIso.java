package edu.psu.compbio.seqcode.projects.kunz.isotop;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.util.ArrayList;

import javax.swing.JPanel;


public class DrawIso extends JPanel
{
	public ArrayList<Prototype> protos; 
	public ArrayList<Color> colors;
	public double minX, maxX, minY, maxY, scalar;
	public DrawIso(ArrayList<Prototype> p)
	{
		protos = p;
		minX = 0; maxX =0; minY = 0; maxY = 0; 
	}
	public void findGraphSpace()
	{
		for(int i = 0; i<protos.size(); i++)
		{
			System.out.println(protos.get(i).name.substring(protos.get(i).name.indexOf("r")+1, protos.get(i).name.indexOf(":")) + "("+protos.get(i).m[0]+ ","+protos.get(i).m[1]+ ")");
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
		if(-1*minX > scalar)
			scalar = -1*minX;
		if(-1*minY> scalar)
			scalar = -1*minY;
		if(maxY>scalar)
			scalar = maxY;
		if(maxX> scalar)
			scalar = maxX;
		System.out.println(scalar);
	}
	public void colorCoder()
	{
		colors = new ArrayList<Color>();
		for(int j = 1; j<=16; j++)
		{
			Color neww  = new Color(((int)(Math.random()*255)),((int)Math.random()*255),((int)Math.random()*255));
			colors.add(neww);
		}
	}
	public void paintComponent(Graphics g)
	 {
	    super.paintComponent(g);
	    colorCoder();
	    setBackground(Color.WHITE);
	    int relativeX = getWidth()/2;
	    int relativeY = getHeight()/2;
	    for(int i = 1; i<16; i++)
	    {
	    	g.drawString("Chrome " + i + ":" , getWidth()- 50, i*10);
	    	g.drawRect(getWidth() - 70, i*10, 15, 15);
	    }
	    for(int i = 0; i<protos.size(); i++)
	    {
	    	Prototype phunk = protos.get(i);
	    	int centerX = ((int)((phunk.m[0]/scalar)*relativeX))+relativeX;
	    	int centerY = ((int)((phunk.m[1]/scalar)*relativeY))+relativeY;
	    	g.setColor(Color.BLACK);
	    	g.drawOval(centerX-2, centerY-2, 4, 4);
	    	for(Prototype pop: phunk.neighbors)
	    	{
	    		//g.drawLine(centerX, centerY, ((int)((pop.m[0]/scalar)*relativeX))+relativeX, ((int)((pop.m[1]/scalar)*relativeY))+relativeY);
	    	}

	    	g.setColor(colors.get((int)Integer.parseInt(protos.get(i).name.substring(protos.get(i).name.indexOf("r")+1, protos.get(i).name.indexOf(":")))-1));
	    	g.fillOval(centerX-3, centerY-3, 6, 6);
	    	
	    }
	}
}
