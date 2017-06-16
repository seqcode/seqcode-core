package org.seqcode.viz.scatter;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.MenuElement;

import org.jfree.chart.ChartPanel;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;

/**
 * ScatterPlotter: Adapted from Multiple Dataset Demo 1 in org.jfree.chart.demo
 * This class codes the application frame that drives a GUI version of the
 * scatter plot generator.
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
public class ScatterPlotter extends ApplicationFrame implements ActionListener {

	private static final long serialVersionUID = 1L;
	protected ChartPanel chartPanel;
	protected ScatterPlot scatter;
	protected int addedDatasets = 0;

	/**
	 * Constructor: initialize the plot & chart
	 * 
	 * @param title
	 */
	public ScatterPlotter(String title) {
		super(title);
		scatter = new ScatterPlot(title);

		JPanel content = new JPanel(new BorderLayout());

		chartPanel = new ChartPanel(scatter.getChart(), true, true, true, true, true);
		MenuElement[] menuElement = chartPanel.getPopupMenu().getSubElements();
		int svgButtonPos = 2;
		for (int k = 0; k < menuElement.length; k++) {
			JMenuItem menuItem = (JMenuItem) menuElement[k];
			String label = menuItem.getText();
			if (label.equals("Save as..."))
				svgButtonPos = k + 1;
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

		final JButton button4 = new JButton("Toggle Log Scale");
		button4.setActionCommand("LOG_TOGGLE");
		button4.addActionListener(this);

		final JPanel buttonPanel = new JPanel(new FlowLayout());
		buttonPanel.add(button1);
		buttonPanel.add(button2);
		buttonPanel.add(button3);
		buttonPanel.add(button4);

		content.add(buttonPanel, BorderLayout.SOUTH);
		chartPanel.setPreferredSize(new java.awt.Dimension(scatter.getWidth(), scatter.getHeight()));
		setContentPane(content);

		pack();
	}

	/**
	 * Handles various actions
	 *
	 * @param e
	 *            the action event.
	 */
	public void actionPerformed(final ActionEvent e) {
		if (e.getActionCommand().equals("LOAD_DATASET")) {
			String pwdName = System.getProperty("user.dir");
			JFileChooser chooser;
			if (pwdName != null)
				chooser = new JFileChooser(new File(pwdName));
			else
				chooser = new JFileChooser();

			int v = chooser.showOpenDialog(null);
			if (v == JFileChooser.APPROVE_OPTION) {
				File f = chooser.getSelectedFile();
				// Color and dotSize should be options
				Color c = new Color(75, 75, 75, 50);
				// hack
				if (addedDatasets > 0)
					c = new Color(0, 0, 255, 50);
				int dotSize = 3;
				boolean drawAnnotations = false;
				scatter.addFileDataset(f, c, dotSize, drawAnnotations);
				addedDatasets++;
			}
		} else if (e.getActionCommand().equals("RANDOM_DATASET")) {
			if (addedDatasets < 20) {
				// Color and dotSize should be options
				Color c = new Color(75, 75, 75, 50);
				int dotSize = 8;
				scatter.addRandomDataset(c, dotSize);
				addedDatasets++;
			}
		} else if (e.getActionCommand().equals("REMOVE_DATASET")) {
			if (addedDatasets >= 0) {
				scatter.removeLastDataset();
				addedDatasets--;
			}
		} else if (e.getActionCommand().equals("LOG_TOGGLE")) {
			scatter.setXLogScale(!scatter.getXLogScale());
			scatter.setYLogScale(!scatter.getYLogScale());
		}
	}

	/**
	 * Action for saving an SVG
	 * 
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
				if (pwdName != null) {
					chooser = new JFileChooser(new File(pwdName));
				} else {
					chooser = new JFileChooser();
				}

				int v = chooser.showSaveDialog(null);
				if (v == JFileChooser.APPROVE_OPTION) {
					File f = chooser.getSelectedFile();
					try {
						scatter.saveImage(f, scatter.getWidth(), scatter.getHeight(), false);
					} catch (IOException ie) {
						ie.printStackTrace(System.err);
					}
				}

			}
		};
	}

	// Modifiers
	public void setWidth(int w) {
		scatter.setWidth(w);
		chartPanel.setPreferredSize(new java.awt.Dimension(scatter.getWidth(), scatter.getHeight()));
	}

	public void setHeight(int h) {
		scatter.setHeight(h);
		chartPanel.setPreferredSize(new java.awt.Dimension(scatter.getWidth(), scatter.getHeight()));
	}

	public boolean isBatch() {
		return scatter.isBatch();
	}

	public ScatterPlot getScatterPlot() {
		return scatter;
	}

	/**
	 * Starting point for the demonstration application.
	 *
	 * @param args
	 */
	public static void main(final String[] args) {

		ScatterPlotter sp = new ScatterPlotter("");

		if (sp.getScatterPlot().isBatch()) {
			sp.pack();
			RefineryUtilities.centerFrameOnScreen(sp);
			sp.setVisible(false);
			// try {
			// scatter.saveSVG(new File("out.svg"), 500, 500);
			// } catch (IOException e) {
			// e.printStackTrace();
			// }
		} else {
			sp.pack();
			RefineryUtilities.centerFrameOnScreen(sp);
			sp.setVisible(true);
		}

	}

}
