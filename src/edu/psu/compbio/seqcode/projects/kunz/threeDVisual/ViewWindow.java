package edu.psu.compbio.seqcode.projects.kunz.threeDVisual;
import java.awt.BorderLayout;
import java.awt.event.*;
import java.io.File;
import java.util.ArrayList;

import javax.swing.*;
import javax.swing.filechooser.FileNameExtensionFilter;

//Rotate Z does not work at all
public class ViewWindow extends JFrame
{
    static ViewWindow window2;
    int winW, winH;
    public View v;
	private static final long serialVersionUID = 1L;
	public ArrayList<String> s;
	public String find;
	public ViewWindow(int x, int y)
	{
		 
		//title the window
	    super("Fun 3d");
	    
	    BorderLayout gui = new BorderLayout();
	    setLayout(gui);
	    winW = x; winH = y;
	 }
	public void drawGrid(String s)
	{
		v = new View(s);
		add(v, BorderLayout.CENTER);
		v.repaint();
		requestFocus();
		validate();
		repaint();
	}
	public void setUpMenu(JMenuBar menubar)
	 {
        JButton x = new JButton("Rotate X");
	     menubar.add(x);
	    x.addActionListener(new ActionListener() {	 
            public void actionPerformed(ActionEvent e)
            {
	 	    	if(v != null)
	 	    	{
	 	    		v.rotating = true;
	 	    		v.dir = 0;
	 	    		v.count = 0;
	            }
        }}); 
	    JButton y = new JButton("Rotate Y");
	     menubar.add(y);
	    y.addActionListener(new ActionListener() {	 
           public void actionPerformed(ActionEvent e)
           {
	 	    	if(v != null)
	 	    	{
	 	    		v.rotating = true;
	 	    		v.dir = 1;
	 	    		v.count = 0;
	            }
       }}); 
	    JButton z = new JButton("Rotate Z");
	     menubar.add(z);
	    z.addActionListener(new ActionListener() {	 
           public void actionPerformed(ActionEvent e)
           {
	 	    	if(v != null)
	 	    	{
	 	    		v.rotating = true;
	 	    		v.dir = 2;
	 	    		v.count = 0;
	            }
       }});
	    JButton stop = new JButton("Stop");
	     menubar.add(stop);
	    stop.addActionListener(new ActionListener() {	 
           public void actionPerformed(ActionEvent e)
           {
	 	    	if(v != null)
	 	    	{
	 	    		v.rotating = false;
	            }
       }}); 
		    
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
	public static void main(String[] args)
	{
		window2 = new ViewWindow(700,700);
		String trained = "";
		if (args.length > 0 )
		{
			trained = args[0];
		}
		else
	    {
    		JFileChooser chooser = new JFileChooser();
    	    FileNameExtensionFilter filter = new FileNameExtensionFilter(
    	        "CSV files", "csv");
    	    chooser.setFileFilter(filter);
    	    chooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
    	    int returnVal = chooser.showOpenDialog(window2);
    	    if(returnVal == JFileChooser.APPROVE_OPTION) {
    	    	trained = System.getProperty("user.dir")+ "/" + chooser.getSelectedFile().getName();
    	    	//System.out.println("You chose to open this file: " +chooser.getSelectedFile().getName());
    	    }
	    }
	    JMenuBar menubar = new JMenuBar();
	    window2.setJMenuBar(menubar);
	    window2.setUpMenu(menubar);
	    
	    window2.setBounds(0,0,700,700);
	    window2.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	    window2.setVisible(true);
	    window2.setResizable(true);	

	    window2.drawGrid(trained);
	}
}

