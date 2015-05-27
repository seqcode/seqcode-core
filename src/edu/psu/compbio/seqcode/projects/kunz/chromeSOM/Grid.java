package edu.psu.compbio.seqcode.projects.kunz.chromeSOM;
import java.awt.BorderLayout;
import java.awt.event.*;
import java.io.File;
import java.util.ArrayList;

import javax.swing.*;
import javax.swing.filechooser.FileNameExtensionFilter;


public class Grid extends JFrame
{
    static Grid window2;
    int winW, winH;
	public BatchMap m;
	public DrawHex d;
	private static final long serialVersionUID = 1L;
	public ArrayList<String> s;
	public String find;
	public Grid(int x, int y, String p)
	  {
		 
		//title the window
	    super(p);
	    
	    BorderLayout gui = new BorderLayout();
	    setLayout(gui);
	    winW = x; winH = y;
	    
	    s = new ArrayList<String>();
	  }
	
	 //For viewing
	 public void drawGrid(String s)
	 {
		d = new DrawHex(s);
		add(d, BorderLayout.CENTER);
		d.nodeBuild(700,700);
    	d.colors();
    	d.heatMapping();
	 }
	 //Finds and shows chromes
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
	                int chr = 0;
	                try{chr = Integer.parseInt(find);
	                }catch (NumberFormatException u){chr = 0;}
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
		    	int sigmaArg = 3;
		    	int iter = 5000;
	    	    
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
	    	    try {
	    	        sigmaArg = Integer.parseInt(args[3]);
	    	    } catch (NumberFormatException e) {
	    	        System.err.println("Argument" + args[3] + " is not an integer. Sigma set to 3 by default");
	    	    }
	    	    try {
	    	        iter = Integer.parseInt(args[4]);
	    	    } catch (NumberFormatException e) {
	    	        System.err.println("Argument" + args[4] + " is not an integer. Iterations set to 5000 by default");
	    	    }
		    	BatchMap o = new BatchMap(xArg,yArg,sigmaArg,iter);
		    	
	   		 	o.go();
				}
	    	}
	    	if(args[0].equalsIgnoreCase("view"))
	    	{

	    	    String trained = "";
	    		JFileChooser chooser = new JFileChooser();
	    	    FileNameExtensionFilter filter = new FileNameExtensionFilter(
	    	        "TXT files", "txt");
	    	    chooser.setFileFilter(filter);
	    	    chooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
	    	    int returnVal = chooser.showOpenDialog(window2);
	    	    if(returnVal == JFileChooser.APPROVE_OPTION) {
	    	    	trained = chooser.getSelectedFile().getName();
	    	    	System.out.println("You chose to open this file: " +
	    	            chooser.getSelectedFile().getName());
	    	    }
	    		window2 = new Grid(700,700, trained);
	    	    JMenuBar menubar = new JMenuBar();
	    	    window2.setJMenuBar(menubar);
	    	    window2.drawGrid(trained);
	    	         
	    	    window2.setUpMenu(menubar);
	    	    
	    	    window2.setBounds(0,0,700,700);
	    	    window2.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	    	    window2.setVisible(true);
	    	    window2.setResizable(true);	
	    	}
		}
	}

