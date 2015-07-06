package edu.psu.compbio.seqcode.projects.kunz.isotop;
import java.awt.BorderLayout;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;
import java.util.Scanner;

import javax.swing.*;


public class IsoGrid extends JFrame
{
    public int winW, winH;
    public Isochrome trainer;
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
	 public void drawGrid()
	 {
		d = new DrawIso(protos);
		d.findGraphSpace();
		add(d, BorderLayout.CENTER);
	 }
	 public void isochrome()
	 {
		 trainer = new Isochrome();
		 protos = trainer.go();
	 }
	 public static void main(String[] args)
	  {
	    IsoGrid window2 = new IsoGrid(700,700);
	    window2.isochrome();
	    window2.setBounds(0,0,700,700);
	    window2.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	    window2.setVisible(true);
	    window2.setResizable(true);
	    window2.drawGrid();
	  }
}