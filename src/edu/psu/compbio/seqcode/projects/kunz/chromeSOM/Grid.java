package edu.psu.compbio.seqcode.projects.kunz.chromeSOM;
import java.awt.BorderLayout;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;
import java.util.Scanner;

import javax.swing.*;


public class Grid extends JFrame
{
    static Grid window2;
    int winW, winH;
	public BatchMap m;
	public DrawHex d;
	private static final long serialVersionUID = 1L;
	public ArrayList<String> s;
	public Grid(int x, int y)
	  {
		 
		//title the window
	    super("Hexagonal Grid");
	    
	    BorderLayout gui = new BorderLayout();
	    setLayout(gui);
	    winW = x; winH = y;
	    
	    s = new ArrayList<String>();
	    }
	 public void drawGrid()
	 {
		d = new DrawHex();
		add(d, BorderLayout.CENTER);
		d.nodeBuild(700,700);
    	d.colors();
    	d.heatMapping();
	 }
	 public void batchMap(int x, int y)
	 {
		 System.out.println(" \'Veiw\' or \'Train\'");
		 Scanner sc = new Scanner(System.in);
		 String cmd = sc.next();
		 if(cmd.equalsIgnoreCase("Train"))
		 {
			 m = new BatchMap(x, y);
		 	 m.go();
		 }
		 window2.drawGrid();
	 }
	 public static void main(String[] args)
	  {
	    window2 = new Grid(700,700);

	    JMenuBar menubar = new JMenuBar();
	    window2.setJMenuBar(menubar);
	    
	    window2.batchMap(10,10);
	    window2.setBounds(0,0,700,700);
	    window2.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	    window2.setVisible(true);
	    window2.setResizable(true);	
	    
	    JMenu location = new JMenu("Location");
	    menubar.add(location);
	  }

}
