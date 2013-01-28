package edu.psu.compbio.seqcode.projects.shaun;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.List;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.MenuElement;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.LogAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.DefaultXYDataset;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

import edu.psu.compbio.seqcode.projects.shaun.viz.ScatterData;


/** 
 * iMNScatterPlotter: Legacy code that hardcoded the generation of some figures for the iMN manuscript
 * 
 * ScatterPlotter: Adapted from Multiple Dataset Demo 1 in org.jfree.chart.demo 
 * 
 * 
 * ===========================================================
 * JFreeChart : a free chart library for the Java(tm) platform
 * ===========================================================
 *
 * (C) Copyright 2000-2004, by Object Refinery Limited and Contributors.
 *
 * Project Info:  http://www.jfree.org/jfreechart/index.html
 *
 * This library is free software; you can redistribute it and/or modify it under the terms
 * of the GNU Lesser General Public License as published by the Free Software Foundation;
 * either version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with this
 * library; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
 * Boston, MA 02111-1307, USA.
 *
 * [Java is a trademark or registered trademark of Sun Microsystems, Inc. 
 * in the United States and other countries.]
 *
 * --------------------------
 * SecondaryDatasetDemo1.java
 * --------------------------
 * (C) Copyright 2004, by Object Refinery Limited.
 *
 * Original Author:  David Gilbert (for Object Refinery Limited).
 * 30-Jan-2004 : Version 1 (DG);
 *
 */

/**
 * A demo showing the addition and removal of multiple datasets / renderers.
 */
public class iMNScatterPlotter extends ApplicationFrame implements ActionListener {

   private static final long serialVersionUID = 1L;
	protected JFreeChart chart;
    protected ChartPanel chartPanel;
    protected XYPlot plot;
    protected int sImageWidth=500, sImageHeight=500;
    protected ScatterData data;
    protected boolean batch = true;
    protected boolean logScale = false; 
    
    /**
     * Constructor: initialize the plot & chart
     * @param title
     */
    public iMNScatterPlotter (String title) {
        super(title);
    
        data = new ScatterData();
        chart = ChartFactory.createScatterPlot(title,
                "X",
                "Y",
                new DefaultXYDataset(), PlotOrientation.VERTICAL,
                false,false,false);
        chart.setBackgroundPaint(Color.white);
        chart.setBorderVisible(false);

        plot = chart.getXYPlot();
        plot.setBackgroundPaint(Color.white);
        plot.setOutlineVisible(false);
        plot.setDomainGridlinesVisible(false);
        plot.setRangeGridlinesVisible(false);
        plot.setDatasetRenderingOrder(DatasetRenderingOrder.FORWARD);
        if(logScale){
        	plot.setDomainAxis(new LogAxis());
        	plot.setRangeAxis(new LogAxis());
        }
        ValueAxis daxis = this.plot.getDomainAxis();
        daxis.setAutoRange(true);
        daxis.setAxisLinePaint(Color.black);
        daxis.setLabelPaint(Color.black);
        daxis.setTickLabelPaint(Color.BLACK);
        daxis.setTickMarkPaint(Color.BLACK);
        ValueAxis raxis = this.plot.getRangeAxis();
        raxis.setAutoRange(true);
        raxis.setAxisLinePaint(Color.black);
        raxis.setLabelPaint(Color.black);
        raxis.setTickLabelPaint(Color.BLACK);
        raxis.setTickMarkPaint(Color.BLACK);
        
        //NumberAxis rangeAxis2 = new NumberAxis("Range Axis 2");
        //rangeAxis2.setAutoRangeIncludesZero(false);
        
        JPanel content = new JPanel(new BorderLayout());
        
        chartPanel = new ChartPanel(chart, true, true, true, true, true);
        MenuElement[] menuElement = chartPanel.getPopupMenu().getSubElements();
        int svgButtonPos = 2;
        for (int k = 0; k < menuElement.length; k++) {
        	JMenuItem menuItem = (JMenuItem)menuElement[k];
        	String label = menuItem.getText();
        	if(label.equals("Save as..."))
        		svgButtonPos=k+1;
        }
        chartPanel.getPopupMenu().insert(getSaveSVGAction(), svgButtonPos);
        content.add(chartPanel);
        
        final JButton button1 = new JButton("Random Dataset");
        button1.setActionCommand("RANDOM_DATASET");
        button1.addActionListener(this);
        
        final JButton button2 = new JButton("Add Dataset");
        button2.setActionCommand("LOAD_DATASET");
        button2.addActionListener(this);
        
        final JButton button3 = new JButton("Remove Dataset");
        button3.setActionCommand("REMOVE_DATASET");
        button3.addActionListener(this);

        final JPanel buttonPanel = new JPanel(new FlowLayout());
        buttonPanel.add(button1);
        buttonPanel.add(button2);
        buttonPanel.add(button3);
        
        content.add(buttonPanel, BorderLayout.SOUTH);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 800));
        setContentPane(content);
    }

    public boolean isBatch(){return batch;}
       
    /**
     * Handles various actions
     *
     * @param e  the action event.
     */
    public void actionPerformed(final ActionEvent e) {
    	if (e.getActionCommand().equals("LOAD_DATASET")) {
    		String pwdName = System.getProperty("user.dir");
            JFileChooser chooser;
            if(pwdName != null) 
                chooser = new JFileChooser(new File(pwdName));
            else
                chooser = new JFileChooser();
            
            int v = chooser.showOpenDialog(null);
            if(v == JFileChooser.APPROVE_OPTION) { 
                File f = chooser.getSelectedFile();
                data.loadFileDataset(f, "S"+(data.getDatasetIndex()+1), 1);
                this.plot.setDataset(
                		data.getDatasetIndex(), data.getDataset(data.getDatasetIndex())
                );
                Color c = new Color(75,75,75,75);
                data.editRenderer(data.getDatasetIndex(), 3, c);
                this.plot.setRenderer(data.getDatasetIndex(), data.getRenderer(data.getDatasetIndex()));
            }            
        }else if (e.getActionCommand().equals("RANDOM_DATASET")) {
            if (data.getDatasetIndex() < 20) {
                data.loadRandomDataset("Random"+(data.getDatasetIndex()+1));
                this.plot.setDataset(
                		data.getDatasetIndex(), data.getDataset(data.getDatasetIndex())
                );
                Color c = new Color(75,75,75,60);
                data.editRenderer(data.getDatasetIndex(), 8, c);
                this.plot.setRenderer(data.getDatasetIndex(), data.getRenderer(data.getDatasetIndex()));                
            }
        }else if (e.getActionCommand().equals("REMOVE_DATASET")) {
            if (data.getDatasetIndex() >= 0) {
                plot.setDataset(data.getDatasetIndex(), null);
                plot.setRenderer(data.getDatasetIndex(), null);
                data.removeLastDataset();
            }
        }
    }
    
    /** 
     * Action for saving an SVG
     * @return
     */
    public Action getSaveSVGAction() {
        return new AbstractAction("Save SVG...") { 
            /**
             * Comment for <code>serialVersionUID</code>
             */
            private static final long serialVersionUID = 1L;

            public void actionPerformed(ActionEvent e) { 
                String pwdName = System.getProperty("user.dir");
                JFileChooser chooser;
                if(pwdName != null) { 
                    chooser = new JFileChooser(new File(pwdName));
                } else {
                    chooser = new JFileChooser();
                }
                
                int v = 
                    chooser.showSaveDialog(null);
                if(v == JFileChooser.APPROVE_OPTION) { 
                    File f = chooser.getSelectedFile();
                    try {
                        saveSVG(f, sImageWidth, sImageHeight);
                        //System.out.println("Saved Image [" + sImageWidth + " by " + sImageHeight +  "]");
                    } catch(IOException ie) {
                        ie.printStackTrace(System.err);
                    }
                }
                
            }
        };
    }
    
    /**
     * Save SVG image
     */
    public void saveSVG(File f, int w, int h) 
    throws IOException { 
        DOMImplementation domImpl =
            GenericDOMImplementation.getDOMImplementation();
        // Create an instance of org.w3c.dom.Document
        Document document = domImpl.createDocument(null, "svg", null);
        // Create an instance of the SVG Generator
        SVGGraphics2D svgGenerator = new SVGGraphics2D(document);
        svgGenerator.setSVGCanvasSize(new Dimension(w,h));
        // Ask the test to render into the SVG Graphics2D implementation
        svgGenerator.setColor(Color.white);        
        svgGenerator.fillRect(0,0,w,h);
        chart.draw(svgGenerator, new Rectangle(0,0,w,h));

        // Finally, stream out SVG to the standard output using UTF-8
        // character to byte encoding
        boolean useCSS = true; // we want to use CSS style attribute
        FileOutputStream outStream = new FileOutputStream(f);
        Writer out = new OutputStreamWriter(outStream, "UTF-8");
        svgGenerator.stream(out, useCSS);
        outStream.flush();
        outStream.close();
    }


    /**
     * Hardcoded test method for generating iMN paper figs
     */
    public void iMNpaperFigs(){
    	//WD: C:\Shaun\ES2MN\iMN\docs\exprs_scatters
    	
    	ValueAxis daxis = this.plot.getDomainAxis();
    	ValueAxis raxis = this.plot.getRangeAxis();
    	
    	//Set 1: ctrl vs iMN with MN TFs
    	/*File f1 = new File("ctrl_vs_iMN.all.txt");
    	data.loadFileDataset(f1, "S1", 1);
        this.plot.setDataset(data.getDatasetIndex(), data.getDataset(data.getDatasetIndex()));
        Color c = new Color(75,75,75,60);
        data.editRenderer(data.getDatasetIndex(), 2, c);
        this.plot.setRenderer(data.getDatasetIndex(), data.getRenderer(data.getDatasetIndex()));
        File f2 = new File("ctrl_vs_iMN.goi1.txt");
        data.loadFileDataset(f2, "S2", 1);
        this.plot.setDataset(data.getDatasetIndex(), data.getDataset(data.getDatasetIndex()));
        Color dColor = Color.red;
        data.editRenderer(data.getDatasetIndex(), 5, dColor);
        this.plot.setRenderer(data.getDatasetIndex(), data.getRenderer(data.getDatasetIndex()));
        List<XYTextAnnotation> annots = data.getAnnotations(data.getDatasetIndex(), 0.0, -0.25, 18, Font.BOLD, dColor);
        for(XYTextAnnotation a : annots)
        	plot.addAnnotation(a);
        daxis.setLabel("Control expression");
        raxis.setLabel("NIL iMN expression");
    	*/        
        
        //Set 2: ctrl vs iMN with MN Term
    	/*File f1 = new File("ctrl_vs_iMN.all.txt");
    	data.loadFileDataset(f1, "S1", 1);
        this.plot.setDataset(data.getDatasetIndex(), data.getDataset(data.getDatasetIndex()));
        Color c = new Color(75,75,75,60);
        data.editRenderer(data.getDatasetIndex(), 2, c);
        this.plot.setRenderer(data.getDatasetIndex(), data.getRenderer(data.getDatasetIndex()));
        File f2 = new File("ctrl_vs_iMN.goi3.txt");
        data.loadFileDataset(f2, "S2", 1);
        this.plot.setDataset(data.getDatasetIndex(), data.getDataset(data.getDatasetIndex()));
        Color dColor = new Color(0,204,0);
        data.editRenderer(data.getDatasetIndex(), 5, dColor);
        this.plot.setRenderer(data.getDatasetIndex(), data.getRenderer(data.getDatasetIndex()));
        List<XYTextAnnotation> annots = data.getAnnotations(data.getDatasetIndex(), 0.0, -0.25, 18, Font.BOLD, dColor);
        for(XYTextAnnotation a : annots)
        	plot.addAnnotation(a);
        daxis.setLabel("Control expression");
        raxis.setLabel("NIL iMN expression");
        */
    	
    	//Set 3: MN vs iMN with MN TFs, term, labeled RC 
    	/*File f1 = new File("MN_vs_iMN.all.txt");
    	data.loadFileDataset(f1, "S1", 1);
        this.plot.setDataset(data.getDatasetIndex(), data.getDataset(data.getDatasetIndex()));
        Color c = new Color(75,75,75,60);
        data.editRenderer(data.getDatasetIndex(), 2, c);
        this.plot.setRenderer(data.getDatasetIndex(), data.getRenderer(data.getDatasetIndex()));
        File f2 = new File("MN_vs_iMN.goi3.txt");
        data.loadFileDataset(f2, "S2", 1);
        this.plot.setDataset(data.getDatasetIndex(), data.getDataset(data.getDatasetIndex()));
        Color dColor = new Color(0,204,0);
        data.editRenderer(data.getDatasetIndex(), 5, dColor);
        this.plot.setRenderer(data.getDatasetIndex(), data.getRenderer(data.getDatasetIndex()));
        File f3 = new File("MN_vs_iMN.goi1.txt");
        data.loadFileDataset(f3, "S3", 1);
        this.plot.setDataset(data.getDatasetIndex(), data.getDataset(data.getDatasetIndex()));
        dColor = Color.red;
        data.editRenderer(data.getDatasetIndex(), 5, dColor);
        this.plot.setRenderer(data.getDatasetIndex(), data.getRenderer(data.getDatasetIndex()));
        File f4 = new File("MN_vs_iMN.goi2.txt");
        data.loadFileDataset(f4, "S4", 1);
        this.plot.setDataset(data.getDatasetIndex(), data.getDataset(data.getDatasetIndex()));
        dColor = Color.blue;
        data.editRenderer(data.getDatasetIndex(), 5, dColor);
        this.plot.setRenderer(data.getDatasetIndex(), data.getRenderer(data.getDatasetIndex()));
        List<XYTextAnnotation> annots = data.getAnnotations(data.getDatasetIndex(), 0.0, -0.25, 18, Font.BOLD, dColor);
        for(XYTextAnnotation a : annots)
        	plot.addAnnotation(a);
        daxis.setLabel("RA/Hh MN expression");
        raxis.setLabel("NIL iMN expression");
        */
    	
    	//Set 3.1: MN vs iMN with MN TFs, term, no RC genes
    	File f1 = new File("MN_vs_iMN.all.txt");
    	data.loadFileDataset(f1, "S1", 1);
        this.plot.setDataset(data.getDatasetIndex(), data.getDataset(data.getDatasetIndex()));
        Color c = new Color(75,75,75,60);
        data.editRenderer(data.getDatasetIndex(), 2, c);
        this.plot.setRenderer(data.getDatasetIndex(), data.getRenderer(data.getDatasetIndex()));
        File f2 = new File("MN_vs_iMN.goi3.txt");
        data.loadFileDataset(f2, "S2", 1);
        this.plot.setDataset(data.getDatasetIndex(), data.getDataset(data.getDatasetIndex()));
        Color dColor = new Color(0,204,0);
        data.editRenderer(data.getDatasetIndex(), 9, dColor);
        this.plot.setRenderer(data.getDatasetIndex(), data.getRenderer(data.getDatasetIndex()));
        File f3 = new File("MN_vs_iMN.goi1.txt");
        data.loadFileDataset(f3, "S3", 1);
        this.plot.setDataset(data.getDatasetIndex(), data.getDataset(data.getDatasetIndex()));
        dColor = Color.red;
        data.editRenderer(data.getDatasetIndex(), 9, dColor);
        this.plot.setRenderer(data.getDatasetIndex(), data.getRenderer(data.getDatasetIndex()));
        //File f4 = new File("MN_vs_iMN.goi2.txt");
        //data.loadFileDataset(f4, "S4", 1);
        //this.plot.setDataset(data.getDatasetIndex(), data.getDataset(data.getDatasetIndex()));
        //dColor = Color.blue;
        //data.editRenderer(data.getDatasetIndex(), 5, dColor);
        //this.plot.setRenderer(data.getDatasetIndex(), data.getRenderer(data.getDatasetIndex()));
        //List<XYTextAnnotation> annots = data.getAnnotations(data.getDatasetIndex(), 0.0, -0.25, 18, Font.BOLD, dColor);
        //for(XYTextAnnotation a : annots)
        //	plot.addAnnotation(a);
        daxis.setLabel("RA/Hh MN expression");
        raxis.setLabel("NIL iMN expression");
        
        
    	//Set 4: MN vs iMN+RA with MN TFs, term, labeled RC 
    	/*File f1 = new File("MN_vs_iMN+RA.all.txt");
    	data.loadFileDataset(f1, "S1", 1);
        this.plot.setDataset(data.getDatasetIndex(), data.getDataset(data.getDatasetIndex()));
        Color c = new Color(75,75,75,60);
        data.editRenderer(data.getDatasetIndex(), 2, c);
        this.plot.setRenderer(data.getDatasetIndex(), data.getRenderer(data.getDatasetIndex()));
        File f2 = new File("MN_vs_iMN+RA.goi3.txt");
        data.loadFileDataset(f2, "S2", 1);
        this.plot.setDataset(data.getDatasetIndex(), data.getDataset(data.getDatasetIndex()));
        Color dColor = new Color(0,204,0);
        data.editRenderer(data.getDatasetIndex(), 5, dColor);
        this.plot.setRenderer(data.getDatasetIndex(), data.getRenderer(data.getDatasetIndex()));
        File f3 = new File("MN_vs_iMN+RA.goi1.txt");
        data.loadFileDataset(f3, "S3", 1);
        this.plot.setDataset(data.getDatasetIndex(), data.getDataset(data.getDatasetIndex()));
        dColor = Color.red;
        data.editRenderer(data.getDatasetIndex(), 5, dColor);
        this.plot.setRenderer(data.getDatasetIndex(), data.getRenderer(data.getDatasetIndex()));
        File f4 = new File("MN_vs_iMN+RA.goi2.txt");
        data.loadFileDataset(f4, "S4", 1);
        this.plot.setDataset(data.getDatasetIndex(), data.getDataset(data.getDatasetIndex()));
        dColor = Color.blue;
        data.editRenderer(data.getDatasetIndex(), 5, dColor);
        this.plot.setRenderer(data.getDatasetIndex(), data.getRenderer(data.getDatasetIndex()));
        daxis.setLabel("RA/Hh MN expression");
        raxis.setLabel("NIL+RA iMN expression");
        */
    	
    	/*
    	//Isl1/2 vs iLhx3 ChIP-seq
    	File f1 = new File("../../chip_compare/scatters/isl_vs_lhx.unscaled");
    	data.loadFileDataset(f1, "S1", 1);
        this.plot.setDataset(data.getDatasetIndex(), data.getDataset(data.getDatasetIndex()));
        Color c = new Color(75,75,75,60);
        data.editRenderer(data.getDatasetIndex(), 2, c);
        this.plot.setRenderer(data.getDatasetIndex(), data.getRenderer(data.getDatasetIndex()));
        File f2 = new File("../../chip_compare/scatters/isl_vs_lhx.outliers.unscaled");
        data.loadFileDataset(f2, "S2", 1);
        this.plot.setDataset(data.getDatasetIndex(), data.getDataset(data.getDatasetIndex()));
        Color dColor = new Color(0,0,255,100);
        data.editRenderer(data.getDatasetIndex(), 3, dColor);
        this.plot.setRenderer(data.getDatasetIndex(), data.getRenderer(data.getDatasetIndex()));
        daxis.setLabel("Isl1/2 ChIP-seq read count");
        raxis.setLabel("iLhx3 ChIP-seq read count");
        */
        
        
        daxis.setLowerBound(1.0);
        daxis.setUpperBound(1000.0);
        raxis.setLowerBound(1.0);
        raxis.setUpperBound(1000.0);
        if(daxis instanceof org.jfree.chart.axis.NumberAxis)
        	((NumberAxis)daxis).setTickUnit(new NumberTickUnit(2.0));
        if(daxis instanceof org.jfree.chart.axis.LogAxis)
        	((LogAxis)daxis).setTickUnit(new NumberTickUnit(1.0));
        if(raxis instanceof org.jfree.chart.axis.NumberAxis)
        	((NumberAxis)raxis).setTickUnit(new NumberTickUnit(2.0));
        if(raxis instanceof org.jfree.chart.axis.LogAxis)
        	((LogAxis)raxis).setTickUnit(new NumberTickUnit(1.0));
         
        daxis.setLabelFont(new Font("Tahoma", Font.BOLD, 28));
        raxis.setLabelFont(new Font("Tahoma", Font.BOLD, 28));
        raxis.setTickLabelFont(new Font("Tahoma", Font.PLAIN, 18));
        daxis.setTickLabelFont(new Font("Tahoma", Font.PLAIN, 18));
    }
    
    /**
     * Starting point for the demonstration application.
     *
     * @param args 
     */
    public static void main(final String[] args) {

        iMNScatterPlotter scatter = new iMNScatterPlotter("");

        if(scatter.isBatch()){
        	scatter.iMNpaperFigs();
        	scatter.pack();
        	RefineryUtilities.centerFrameOnScreen(scatter);
        	scatter.setVisible(true);
        	try {
				scatter.saveSVG(new File("out.svg"), 500, 500);
			} catch (IOException e) {
				e.printStackTrace();
			}
        }else{
        	scatter.pack();
        	RefineryUtilities.centerFrameOnScreen(scatter);
        	scatter.setVisible(true);
        }

    }

}

