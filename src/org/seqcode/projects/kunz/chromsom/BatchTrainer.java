package org.seqcode.projects.kunz.chromsom;
import java.awt.BorderLayout;
import java.awt.event.*;
import java.io.File;
import java.util.ArrayList;

import javax.swing.*;
import javax.swing.filechooser.FileNameExtensionFilter;


public class BatchTrainer extends JFrame
{
    static BatchTrainer window2;
    int winW, winH;
	public BatchMap m;
	public DrawHex d;
	private static final long serialVersionUID = 1L;
	public ArrayList<String> s;
	public String find;
	public BatchTrainer(int x, int y, String p)
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
		 JButton search = new JButton("Search");
	     menubar.add(search);
	     JButton view = new JButton("View Swap");
	     menubar.add(view);
	     JButton button = new JButton("Chrome");
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
	    search.addActionListener(new ActionListener() {
	    	 
            public void actionPerformed(ActionEvent e)
            {
                //Execute when search is pressed
            	String trained = "";
 	    		JFileChooser chooser = new JFileChooser();
 	    	    FileNameExtensionFilter filter = new FileNameExtensionFilter(
 	    	        "TXT files", "txt");
 	    	    chooser.setFileFilter(filter);
 	    	    chooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
 	    	    int returnVal = chooser.showOpenDialog(window2);
 	    	    if(returnVal == JFileChooser.APPROVE_OPTION) {
 	    	    	trained = chooser.getSelectedFile().getName();
 	    	    	//System.out.println("You chose to open this file: " +chooser.getSelectedFile().getName());
 	    	    }
 	    	if(d != null)
 	    		d.search(trained);
            }
        }); 
	    view.addActionListener(new ActionListener() {
	    	 
            public void actionPerformed(ActionEvent e)
            {
                //Execute when button is pressed
                d.swap();
            }
        });  
	 }
	 public static void main(String[] args)
	  {
		 //to view: args = "view"
		 //to train: args = "train" 'int x' 'int y' 'int sigma'
		if (args.length > 0)
		{
			if(args[0].equalsIgnoreCase("tester"))
			{
				BatchMap o = new BatchMap(6,6,2,5000); o.go();
				o = new BatchMap(8,8,2,1000); o.go();
				o = new BatchMap(10,10,2,1000); o.go();
				o = new BatchMap(12,12,2,1000); o.go();
				
				o = new BatchMap(20,20,3,1000); o.go();
				o = new BatchMap(8,8,3,1000); o.go();
				o = new BatchMap(10,10,3,1000); o.go();
				o = new BatchMap(12,12,3,1000); o.go();
				o = new BatchMap(20,20,3,5000); o.go();
				
				o = new BatchMap(20,20,5,1000); o.go();
				o = new BatchMap(8,8,5,1000); o.go();
				o = new BatchMap(10,10,5,1000); o.go();
				o = new BatchMap(12,12,5,1000); o.go();
				
				o = new BatchMap(20,20,10,1000); o.go();
				o = new BatchMap(8,8,10,1000); o.go();
				o = new BatchMap(10,10,10,1000); o.go();
				o = new BatchMap(12,12,10,1000); o.go();
			}
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
	    	    	//System.out.println("You chose to open this file: " +chooser.getSelectedFile().getName());
	    	    }
	    		window2 = new BatchTrainer(700,700, trained);
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

