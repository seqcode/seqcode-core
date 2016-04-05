package edu.psu.compbio.seqcode.projects.kunz.isotop;
import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.File;
import java.util.ArrayList;
import java.util.Scanner;

import javax.swing.*;
import javax.swing.filechooser.FileNameExtensionFilter;

import edu.psu.compbio.seqcode.projects.kunz.chromeSOM.BatchMap;
import edu.psu.compbio.seqcode.projects.kunz.chromeSOM.BatchTrainer;


public class IsoGrid extends JFrame
{
    public int winW, winH;
    public Isochrome trainer;
    public static IsoGrid window2;
    public DrawIso d;
	private static final long serialVersionUID = 1L;
	public ArrayList<String> s;
	public ArrayList<Prototype> protos;
	public IsoGrid(int x, int y)
	  {
		 
		//title the window
	    super("IsoChrome");
	    
	    BorderLayout gui = new BorderLayout();
	    setLayout(gui);
	    winW = x; winH = y;
	    
	    s = new ArrayList<String>();
	    }
	 public void drawGrid(String file)
	 {
		d = new DrawIso(file);
		d.findGraphSpace();
		add(d, BorderLayout.CENTER);
		d.repaint();
		requestFocus();
		validate();
		repaint();
	 }
	 public void setUpMenu(JMenuBar menubar)
	 {
	        JButton search = new JButton("Search");
		     menubar.add(search);

			 JButton button = new JButton("Chrome");
		        menubar.add(button);
		     final JTextField texter = new JTextField("Chrome...", 20);
		    menubar.add(texter);
		    menubar.add(Box.createHorizontalGlue());
		    
		    
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
	            }
	        }); 
		    button.addActionListener(new ActionListener() {
		    	 
	            public void actionPerformed(ActionEvent e)
	            {
	                //Execute when button is pressed
	                String find = texter.getText();
	                int chr = 0;
	                try{chr = Integer.parseInt(find);
	                }catch (NumberFormatException u){chr = 0;}
	                d.countingDPS(chr);
	            }
	        });   
		   /* search.addActionListener(new ActionListener() {
		    	 
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
	        }); */
	 }
	 public void isochrome()
	 {
		 trainer = new Isochrome();
		 trainer.go();
		 System.exit(0);
	 }
	 public static void main(String[] args)
	  {
	    window2 = new IsoGrid(700,700);
	    
	    //window2.drawGrid();
	    if (args.length > 0)
		{
			
	    	if(args[0].equalsIgnoreCase(("train")))
	    	{
			    window2.isochrome();
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
	    	    if(returnVal == JFileChooser.APPROVE_OPTION) 
	    	    {
	    	    	trained = chooser.getSelectedFile().getName();
	    	    	System.out.println("You chose to open this file: " +
	    	            chooser.getSelectedFile().getName());
	    	    }
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
}
