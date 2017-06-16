package org.seqcode.viz.scatter;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.List;

import javax.imageio.ImageIO;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.LogAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.DatasetRenderingOrder;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.DefaultXYDataset;
import org.jfree.ui.TextAnchor;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

import Jama.Matrix;

/**
 * ScatterPlot: Adapted from Multiple Dataset Demo 1 in org.jfree.chart.demo
 * This class contains the chart, plot, data, and the rest of the main drawing
 * items for a scatter plot.
 * 
 * Extended by Shaun Mahony
 * 
 * =========================================================== JFreeChart : a
 * free chart library for the Java(tm) platform
 * ===========================================================
 *
 * (C) Copyright 2000-2004, by Object Refinery Limited and Contributors.
 *
 * Project Info: http://www.jfree.org/jfreechart/index.html
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 2.1 of the License, or (at your option)
 * any later version.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
 *
 * [Java is a trademark or registered trademark of Sun Microsystems, Inc. in the
 * United States and other countries.]
 *
 * -------------------------- SecondaryDatasetDemo1.java
 * -------------------------- (C) Copyright 2004, by Object Refinery Limited.
 *
 * Original Author: David Gilbert (for Object Refinery Limited). 30-Jan-2004 :
 * Version 1 (DG);
 *
 */

public class ScatterPlot {
	protected JFreeChart chart;
	protected XYPlot plot;
	protected ValueAxis daxis = null;
	protected ValueAxis raxis = null;
	protected int width = 800, height = 800;
	protected ScatterData data;
	protected boolean batch = false;
	protected boolean xLogScale = false;
	protected boolean yLogScale = false;
	protected String xAxisLabel = "";
	protected String yAxisLabel = "";
	protected int defaultDotSize = 3;
	protected Color defaultDotColor = new Color(75, 75, 75, 50);
	protected int dataLabelFontSize = 24;
	protected int xAxisLabelFontSize = 28;
	protected int yAxisLabelFontSize = 28;
	protected int xAxisTickLabelFontSize = 14;
	protected int yAxisTickLabelFontSize = 14;

	/**
	 * Constructor: initialize the plot & chart
	 * 
	 * @param title
	 */
	public ScatterPlot(String title) {
		data = new ScatterData();
		chart = ChartFactory.createScatterPlot(title, "X", "Y", new DefaultXYDataset(), PlotOrientation.VERTICAL, false,
				false, false);
		chart.setBackgroundPaint(Color.white);
		chart.setBorderVisible(false);

		plot = chart.getXYPlot();
		plot.setBackgroundPaint(Color.white);
		plot.setOutlineVisible(false);
		plot.setDomainGridlinesVisible(false);
		plot.setRangeGridlinesVisible(false);
		plot.setDomainZeroBaselineVisible(true);
		plot.setRangeZeroBaselineVisible(true);
		plot.setDatasetRenderingOrder(DatasetRenderingOrder.FORWARD);

		this.configAxes();
	}

	/**
	 * Configure axes to log or linear scales
	 */
	protected void configAxes() {
		if (xLogScale)
			plot.setDomainAxis(new LogAxis());
		else
			plot.setDomainAxis(new NumberAxis());
		if (yLogScale)
			plot.setRangeAxis(new LogAxis());
		else
			plot.setRangeAxis(new NumberAxis());
		daxis = this.plot.getDomainAxis();
		daxis.setAxisLinePaint(Color.black);
		daxis.setLabelPaint(Color.black);
		daxis.setTickLabelPaint(Color.BLACK);
		daxis.setTickMarkPaint(Color.BLACK);
		raxis = this.plot.getRangeAxis();
		raxis.setAxisLinePaint(Color.black);
		raxis.setLabelPaint(Color.black);
		raxis.setTickLabelPaint(Color.BLACK);
		raxis.setTickMarkPaint(Color.BLACK);

		if (xLogScale) {
			daxis.setLowerBound(1.0);
			daxis.setUpperBound(1000000.0);
			daxis.setAutoRange(false);
			if (daxis instanceof org.jfree.chart.axis.LogAxis)
				((LogAxis) daxis).setTickUnit(new NumberTickUnit(1.0));
		} else {
			daxis.setAutoRange(true);
		}
		if (yLogScale) {
			raxis.setLowerBound(1.0);
			raxis.setUpperBound(1000000.0);
			raxis.setAutoRange(false);
			if (daxis instanceof org.jfree.chart.axis.LogAxis)
				((LogAxis) raxis).setTickUnit(new NumberTickUnit(1.0));
		} else {
			raxis.setAutoRange(true);
		}
		daxis.setLabel(xAxisLabel);
		raxis.setLabel(yAxisLabel);

		daxis.setLabelFont(new Font("Tahoma", Font.BOLD, xAxisLabelFontSize));
		raxis.setLabelFont(new Font("Tahoma", Font.BOLD, yAxisLabelFontSize));
		daxis.setTickLabelFont(new Font("Tahoma", Font.PLAIN, xAxisTickLabelFontSize));
		raxis.setTickLabelFont(new Font("Tahoma", Font.PLAIN, yAxisTickLabelFontSize));
	}

	/**
	 * Save image If raster is false, save a SVG, otherwise save a PNG
	 * 
	 */
	public void saveImage(File f, int w, int h, boolean raster) throws IOException {
		if (raster) {
			BufferedImage im = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
			Graphics g = im.getGraphics();
			Graphics2D g2 = (Graphics2D) g;
			g2.setRenderingHints(
					new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
			chart.draw(g2, new Rectangle(0, 0, w, h));
			ImageIO.write(im, "png", f);
		} else {
			DOMImplementation domImpl = GenericDOMImplementation.getDOMImplementation();
			// Create an instance of org.w3c.dom.Document
			Document document = domImpl.createDocument(null, "svg", null);
			// Create an instance of the SVG Generator
			SVGGraphics2D svgGenerator = new SVGGraphics2D(document);
			svgGenerator.setSVGCanvasSize(new Dimension(w, h));
			// Ask the test to render into the SVG Graphics2D implementation
			svgGenerator.setColor(Color.white);
			svgGenerator.fillRect(0, 0, w, h);
			chart.draw(svgGenerator, new Rectangle(0, 0, w, h));

			// Finally, stream out SVG to the standard output using UTF-8
			// character to byte encoding
			boolean useCSS = true; // we want to use CSS style attribute
			FileOutputStream outStream = new FileOutputStream(f);
			Writer out = new OutputStreamWriter(outStream, "UTF-8");
			svgGenerator.stream(out, useCSS);
			outStream.flush();
			outStream.close();
		}
	}

	public JFreeChart getChart() {
		return chart;
	}

	public int getWidth() {
		return width;
	}

	public int getHeight() {
		return height;
	}

	public boolean isBatch() {
		return batch;
	}

	public boolean getXLogScale() {
		return xLogScale;
	}

	public boolean getYLogScale() {
		return yLogScale;
	}

	// Modifiers
	public void setWidth(int w) {
		width = w;
	}

	public void setHeight(int h) {
		height = h;
	}

	public void setBatch(boolean b) {
		batch = b;
	}

	public void setGridlinesVisible(boolean v) {
		plot.setDomainGridlinesVisible(false);
		plot.setRangeGridlinesVisible(false);
	}

	public void setXLogScale(boolean l) {
		this.xLogScale = l;
		this.configAxes();
	}

	public void setYLogScale(boolean l) {
		this.yLogScale = l;
		this.configAxes();
	}

	public void setXAxisLabel(String s) {
		this.xAxisLabel = s;
		daxis.setLabel(xAxisLabel);
	}

	public void setYAxisLabel(String s) {
		this.yAxisLabel = s;
		raxis.setLabel(yAxisLabel);
	}

	public void setDataLabelFontSize(int d) {
		dataLabelFontSize = d;
		chart.fireChartChanged();
	}

	public void setXAxisLabelFontSize(int d) {
		xAxisLabelFontSize = d;
		daxis.setLabelFont(new Font("Tahoma", Font.BOLD, xAxisLabelFontSize));
	}

	public void setYAxisLabelFontSize(int d) {
		yAxisLabelFontSize = d;
		raxis.setLabelFont(new Font("Tahoma", Font.BOLD, yAxisLabelFontSize));
	}

	public void setXAxisTickLabelFontSize(int d) {
		xAxisTickLabelFontSize = d;
		daxis.setLabelFont(new Font("Tahoma", Font.BOLD, xAxisTickLabelFontSize));
	}

	public void setYAxisTickLabelFontSize(int d) {
		yAxisTickLabelFontSize = d;
		raxis.setLabelFont(new Font("Tahoma", Font.BOLD, yAxisTickLabelFontSize));
	}

	public void setXAutoRange(boolean a) {
		daxis.setAutoRange(a);
	}

	public void setYAutoRange(boolean a) {
		raxis.setAutoRange(a);
	}

	public void setXRange(double min, double max) {
		daxis.setLowerBound(min);
		daxis.setUpperBound(max);
	}

	public void setYRange(double min, double max) {
		raxis.setLowerBound(min);
		raxis.setUpperBound(max);
	}

	public void setXRangeFromData() {
		daxis.setLowerBound(data.getMin(0));
		daxis.setUpperBound(data.getMax(0) + 1);
	}

	public void setYRangeFromData() {
		raxis.setLowerBound(data.getMin(1));
		raxis.setUpperBound(data.getMax(1) + 1);
	}

	/**
	 * Add a dataset to the plot
	 * 
	 * @param name
	 *            String
	 * @param datasetMatrix
	 *            2D matrix
	 * @param c
	 *            Color
	 * @param dotSize
	 *            int
	 * @param drawAnnotations
	 *            boolean
	 */
	public void addDataset(String name, Matrix datasetMatrix, Color c, int dotSize, boolean drawAnnotations,
			boolean drawConnectingLines, boolean drawDots) {
		data.loadDataset(name, datasetMatrix);
		this.plot.setDataset(data.getDatasetIndex(), data.getDataset(data.getDatasetIndex()));
		data.editRenderer(data.getDatasetIndex(), dotSize, c, drawConnectingLines, drawDots);
		this.plot.setRenderer(data.getDatasetIndex(), data.getRenderer(data.getDatasetIndex()));
		if (drawAnnotations) {
			List<XYTextAnnotation> annots = data.getAnnotations(data.getDatasetIndex(), 0.0, 0.0,
					TextAnchor.BOTTOM_RIGHT, dataLabelFontSize, Font.PLAIN, c);
			for (XYTextAnnotation a : annots)
				plot.addAnnotation(a);
		}
	}

	public void addDataset(String name, Matrix datasetMatrix, Color c, int dotSize) {
		addDataset(name, datasetMatrix, c, dotSize, false, false, true);
	}

	public void addDataset(String name, Matrix datasetMatrix, Color c) {
		addDataset(name, datasetMatrix, c, defaultDotSize, false, false, true);
	}

	public void addDataset(String name, Matrix datasetMatrix) {
		addDataset(name, datasetMatrix, defaultDotColor, defaultDotSize, false, false, true);
	}

	/**
	 * Add a dataset from a file
	 * 
	 * @param f
	 *            File
	 * @param c
	 *            Color
	 * @param dotSize
	 *            integer
	 */
	public void addFileDataset(File f, Color c, int dotSize, boolean drawAnnotations) {
		data.loadFileDataset(f, "S" + (data.getDatasetIndex() + 1), 1);
		this.plot.setDataset(data.getDatasetIndex(), data.getDataset(data.getDatasetIndex()));
		data.editRenderer(data.getDatasetIndex(), dotSize, c);
		this.plot.setRenderer(data.getDatasetIndex(), data.getRenderer(data.getDatasetIndex()));
		if (drawAnnotations) {
			List<XYTextAnnotation> annots = data.getAnnotations(data.getDatasetIndex(), 0.0, 0.0,
					TextAnchor.BOTTOM_RIGHT, dataLabelFontSize, Font.PLAIN, c);
			for (XYTextAnnotation a : annots)
				plot.addAnnotation(a);
		}
		configAxes();
	}

	/**
	 * Add a random dataset (for testing)
	 */
	public void addRandomDataset(Color c, int dotSize) {
		data.loadRandomDataset("Random" + (data.getDatasetIndex() + 1));
		this.plot.setDataset(data.getDatasetIndex(), data.getDataset(data.getDatasetIndex()));
		data.editRenderer(data.getDatasetIndex(), dotSize, c);
		this.plot.setRenderer(data.getDatasetIndex(), data.getRenderer(data.getDatasetIndex()));
		configAxes();
	}

	/**
	 * Add a domain marker to the plot
	 */
	public void addDomainMarker(double position) {
		ValueMarker marker = new ValueMarker(position); // position is the value
														// on the axis
		marker.setPaint(Color.red);
		plot.addDomainMarker(marker);
	}

	/**
	 * Add a range marker to the plot
	 */
	public void addRangeMarker(double position) {
		ValueMarker marker = new ValueMarker(position); // position is the value
														// on the axis
		marker.setPaint(Color.red);
		plot.addRangeMarker(marker);
	}

	/**
	 * Remove the last dataset on the pile
	 */
	public void removeLastDataset() {
		if (data.getDatasetIndex() >= 0) {
			plot.setDataset(data.getDatasetIndex(), null);
			plot.setRenderer(data.getDatasetIndex(), null);
			data.removeLastDataset();
		}
	}
}
