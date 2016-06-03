package org.seqcode.projects.shaun.viz;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import javax.swing.JFrame;

import org.seqcode.gse.utils.ArgParser;
import org.seqcode.gse.viz.paintable.PaintableFrame;


public class MakeColorBarPlots{

	private ColorBarPlotPaintable painter;
	private PaintableFrame plotter;
	private int screenSizeX=1000, screenSizeY=900;
	private String inputFile;

	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("in")) { 
            System.err.println("Usage:\n " +
                               "MakeColorBarPlots " +
                               "--in <file name> ");
            return;
        }
        String dfile = ap.getKeyValue("in");
		MakeColorBarPlots maker = new MakeColorBarPlots(dfile);		
	}
	
	public MakeColorBarPlots(String df){
		inputFile = df;
		BarPlotData bpd = new BarPlotData();		
		bpd.readData(inputFile);
		int rows = bpd.getRows();
		
		painter = new ColorBarPlotPaintable(bpd);
		screenSizeY=rows+painter.getYborders();
		painter.setImageHeight(screenSizeY);
		painter.setImageWidth(screenSizeX);
		System.out.println(screenSizeX+"\t"+screenSizeY);
		
		plotter = new PaintableFrame("Color Bar Plots", painter);
		plotter.getContentPane().setBackground(Color.white);
		plotter.setBackground(Color.white);
		plotter.setSize(screenSizeX, screenSizeY);
		plotter.setVisible(true);
		plotter.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
		
	
}
