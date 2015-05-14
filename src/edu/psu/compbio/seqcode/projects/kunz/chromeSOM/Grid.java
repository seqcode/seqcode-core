package edu.psu.compbio.seqcode.projects.kunz.chromeSOM;
import java.awt.BorderLayout;
import java.awt.event.*;
import java.util.ArrayList;

import javax.swing.*;


public class Grid extends JFrame
{
    static Grid window2;
    int winW, winH;
	public BatchMap m;
	public DrawHex d;
	private static final long serialVersionUID = 1L;
	public ArrayList<String> s;
	public String find;
	public Grid(int x, int y)
	  {
		 
		//title the window
	    super("Super SOM");
	    
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
		 m = new BatchMap(x, y, 175);
		 m.go();
		 window2.drawGrid();
	 }
	 public void setUpMenu(JMenuBar menubar)
	 {
		 JButton button = new JButton("Seach");
	        menubar.add(button);
	        final JTextField texter = new JTextField("Chrome...", 20);
		    menubar.add(texter);
		    menubar.add(Box.createHorizontalGlue());
		    
		    button.addActionListener(new ActionListener() {
		    	 
	            public void actionPerformed(ActionEvent e)
	            {
	                //Execute when button is pressed
	                find = texter.getText();
	                int chr = Integer.parseInt(find);
	                System.out.print(chr);
	                d.countingDPS(chr);
	            }
	        });   
	 }
	 public static void main(String[] args)
	  {
		 //to view: args = "view"
		 //to train: args = "train" 'int x' 'int y' 'int sigma'
		if (args.length > 0)
		{
	    	if(args[0].equalsIgnoreCase(("train")))
	    	{
	    		int xArg = 0;
		    	int yArg = 0;
		    	int sigmaArg = 175;
	    	    
		    	try {
	    	        xArg = Integer.parseInt(args[1]);
	    	    } catch (NumberFormatException e) {
	    	        System.err.println("Argument" + args[1]);
	    	    }
	    	    try {
	    	        yArg = Integer.parseInt(args[2]);
	    	    } catch (NumberFormatException e) {
	    	        System.err.println("Argument" + args[2]);
	    	    }
	    	    if(args.length>3)
		    	{
	    	    	try {
		    	        sigmaArg = Integer.parseInt(args[3]);
		    	    } catch (NumberFormatException e) {
		    	        System.err.println("Argument" + args[3] + " is not an integer. Sigma set to 175 by default");
		    	    }
	    	    }
	    	    else
	    	    {
	    	    	System.err.println("No sigma argument entered. Sigma set to 175 by default");
	    	    }
		    	BatchMap o = new BatchMap(xArg,yArg,sigmaArg);
		    	
	   		 	o.go();
	    	}
	    	if(args[0].equalsIgnoreCase("view"))
	    	{
	    		window2 = new Grid(700,700);
	    	    JMenuBar menubar = new JMenuBar();
	    	    window2.setJMenuBar(menubar);
	    	    window2.drawGrid();
	    	         
	    	    window2.setUpMenu(menubar);
	    	    
	    	    window2.setBounds(0,0,700,700);
	    	    window2.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	    	    window2.setVisible(true);
	    	    window2.setResizable(true);	
	    	}
		}
	}
}
