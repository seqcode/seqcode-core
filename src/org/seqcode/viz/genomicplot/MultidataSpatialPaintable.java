package org.seqcode.viz.genomicplot;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.seqcode.deepseq.experiments.Sample;
import org.seqcode.genome.location.ExonicGene;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.gseutils.Pair;
import org.seqcode.viz.paintable.AbstractPaintable;

public class MultidataSpatialPaintable extends AbstractPaintable {

	private boolean reverseIt = true;
	private int screenSizeX, screenSizeY;
	private final int geneHeight = 15;
	private final int boxHeight = 30, boxWidth = 18, siteWidth = 10, litWidth = 4, motifWidth = 5, motifHeight = 10;
	private final int exptTrackHeight = 120;
	private final int colorbarHeight = 14, colorbarWidth = 200;
	private final int rTrack = 75;
	private final double maxExp = 6.0, minExp = -3.0;
	private final double maxSiteZ = 200, minSiteZ = 60;
	private Color geneColor = new Color(232, 232, 232);
	private Color motifColor = new Color(51, 102, 0);
	private Color siteMaxColor = new Color(92, 0, 184);
	private Color siteMinColor = new Color(184, 0, 184);
	private Color rareColDR5 = new Color(0, 153, 0);
	private Color rareColDR2 = new Color(153, 153, 255);
	private final int topBorder = 50, bottomBorder = 50;
	private final int leftBorder = 25, rightBorder = 25;
	private int topBound, bottomBound, leftBound, rightBound, baseLine;
	private List<String> timeLabels;
	private List<Sample> times;
	private List<ExonicGene> genes;
	private Map<String, List<Double>> expression;
	private Map<String, List<Region>> sites = null;
	private List<Pair<String, StrandedRegion>> motifs = new ArrayList<Pair<String, StrandedRegion>>();
	private List<Region> literature = new ArrayList<Region>();
	private Region gRegion;
	private int rstart, rstop, rwidth;
	private String chr;
	private Map<Sample, ThinOverlapPaintable> exptPainters = new HashMap<Sample, ThinOverlapPaintable>();
	private List<ThinOverlapPaintable> thinpaints = new ArrayList<ThinOverlapPaintable>();

	public MultidataSpatialPaintable(List<String> timepointLabels, List<Sample> timepoints, Region genomeRegion,
			List<ExonicGene> genes, Map<String, List<Double>> expression, Map<String, List<Region>> sites,
			List<Pair<String, StrandedRegion>> mHits, List<Region> lits) {
		timeLabels = timepointLabels;
		times = timepoints;
		rstart = genomeRegion.getStart();
		rstop = genomeRegion.getEnd();
		chr = genomeRegion.getChrom();
		gRegion = genomeRegion;
		rwidth = rstop - rstart;
		this.genes = genes;
		this.expression = expression;
		this.sites = sites;
		motifs = mHits;
		literature = lits;

		for (Sample t : times) {
			ThinOverlapPaintable tp = new ThinOverlapPaintable(gRegion, sites.get(t), t, 200, false);
			tp.setReverse(reverseIt);
			tp.setMaxOverlap(220);
			exptPainters.put(t, tp);
			thinpaints.add(tp);
		}
	}

	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		Graphics2D g2d = (Graphics2D) g;
		FontMetrics metrics = g2d.getFontMetrics();
		screenSizeX = x2 - x1;
		screenSizeY = y2 - y1;
		g2d.setColor(Color.white);
		g2d.fillRect(0, 0, screenSizeX, screenSizeY);
		topBound = topBorder;
		bottomBound = screenSizeY - bottomBorder;
		leftBound = leftBorder + 30;
		rightBound = screenSizeX - rightBorder;
		baseLine = (topBound + bottomBound) / 2;

		// AlignedMotifPaintable timMotifPaint = new
		// AlignedMotifPaintable(gRegion, motifs, thinpaints.get(1),
		// thinpaints.get(0), reverseIt, motifWidth);

		// Draw the genes
		boolean drawGeneNum = true;
		int offset = 0, time = 0;
		for (Sample t : times) {
			ThinOverlapPaintable p = exptPainters.get(t);
			// p.paintItem(g, leftBound,
			// baseLine-offset-motifHeight-exptTrackHeight, rightBound,
			// baseLine-offset-motifHeight-1);
			p.paintItem(g, leftBound, baseLine - offset - exptTrackHeight, rightBound, baseLine - offset);
			// Draw the RAREs
			// drawRAREs(g2d, baseLine-offset-motifHeight);
			// Draw the axis
			g2d.setColor(Color.lightGray);
			g2d.setStroke(new BasicStroke(1.0f));
			// g2d.drawLine(leftBound, baseLine-offset-motifHeight, rightBound,
			// baseLine-offset-motifHeight);

			// Draw the timepoint name
			g2d.setColor(Color.darkGray);
			g2d.setFont(new Font("Ariel", Font.BOLD, 14));
			metrics = g2d.getFontMetrics();
			String currT = timeLabels.get(time);
			int ypos = baseLine - offset - (metrics.getHeight() / 2);

			offset += (exptTrackHeight + motifHeight + 10);
			time++;
		}
		offset = 0;
		time = 0;
		for (Sample t : times) {
			drawGenes(g2d, x1, baseLine - offset, x2, time, drawGeneNum);
			offset += (exptTrackHeight + motifHeight + (geneHeight + 14));
			time++;
			drawGeneNum = false;
		}

		// Draw some coordinates
		g2d.setColor(Color.black);
		g2d.setFont(new Font("Ariel", Font.PLAIN, 14));
		metrics = g2d.getFontMetrics();
		AffineTransform oldtrans = g2d.getTransform();
		AffineTransform newtrans = new AffineTransform();
		String text = reverseIt ? new String("chr" + chr + ":" + rstop) : new String("chr" + chr + ":" + rstart);
		newtrans.translate(leftBound, baseLine + 30 + metrics.stringWidth(text));
		newtrans.rotate(Math.toRadians(-90));
		g2d.setTransform(newtrans);
		g2d.drawString(text, 0, 0);
		g2d.setTransform(oldtrans);
		newtrans = new AffineTransform();
		text = reverseIt ? new String("chr" + chr + ":" + rstart) : new String("chr" + chr + ":" + rstop);
		newtrans.translate(rightBound + (metrics.getHeight() / 2), baseLine + 30 + metrics.stringWidth(text));
		newtrans.rotate(Math.toRadians(-90));
		g2d.setTransform(newtrans);
		g2d.drawString(text, 0, 0);
		g2d.setTransform(oldtrans);

		// draw the colorbars
		// drawExpColorBar(g2d, (leftBound+rightBound)/2-(colorbarWidth/2),
		// topBound);

		/*
		 * //Special g2d.setFont(new Font("Ariel", Font.BOLD, 18));
		 * g2d.setColor(rareColDR2); g2d.drawString("DR2", leftBound+10,
		 * currLine+(boxHeight/2)+(metrics.getHeight()/2));
		 * g2d.setColor(Color.black); g2d.drawString(" & ",
		 * leftBound+10+metrics.stringWidth("DR2"),
		 * currLine+(boxHeight/2)+(metrics.getHeight()/2));
		 * g2d.setColor(rareColDR5); g2d.drawString("DR5",
		 * leftBound+10+metrics.stringWidth("DR2 & "),
		 * currLine+(boxHeight/2)+(metrics.getHeight()/2));
		 * 
		 * //Draw the literature sites currLine = baseLine+rTrack+boxHeight;
		 * for(Region r : literature){ int gx1 =
		 * xcoord(r.getStart())-(litWidth/2); g2d.setColor(motifColor);
		 * g2d.fillRect(gx1, currLine, litWidth, boxHeight); }//Special
		 * g2d.setColor(motifColor); g2d.drawString("Literature", leftBound+10,
		 * currLine+(boxHeight/2)+(metrics.getHeight()/2));
		 */
		// drawLegend(g2d, (leftBound+rightBound)/2+colorbarWidth, topBound);
	}

	private int xcoord(int coord) {
		double frac = ((double) (coord - rstart) / (double) (rstop - rstart));
		int pos;
		if (reverseIt)
			pos = leftBound + (int) ((double) (rightBound - leftBound) * (1 - frac));
		else
			pos = leftBound + (int) ((double) (rightBound - leftBound) * frac);
		return (pos);
	}

	private Color expColor(double v) {
		Color c;
		if (v > 0) {
			Color maxColor = Color.yellow;
			Color minColor = Color.black;

			double sVal = v > maxExp ? 1 : (v) / (maxExp);
			int red = (int) (maxColor.getRed() * sVal + minColor.getRed() * (1 - sVal));
			int green = (int) (maxColor.getGreen() * sVal + minColor.getGreen() * (1 - sVal));
			int blue = (int) (maxColor.getBlue() * sVal + minColor.getBlue() * (1 - sVal));
			c = new Color(red, green, blue);
		} else {
			Color maxColor = Color.black;
			Color minColor = Color.blue;
			double sVal = v < minExp ? 1 : (minExp - v) / (minExp);
			int red = (int) (maxColor.getRed() * sVal + minColor.getRed() * (1 - sVal));
			int green = (int) (maxColor.getGreen() * sVal + minColor.getGreen() * (1 - sVal));
			int blue = (int) (maxColor.getBlue() * sVal + minColor.getBlue() * (1 - sVal));
			c = new Color(red, green, blue);
		}
		return (c);
	}

	private Color getSiteColor(double v) {
		Color c;
		Color maxColor = siteMaxColor;
		Color minColor = siteMinColor;
		Color zeroColor = Color.white;
		if (v > maxSiteZ)
			v = maxSiteZ;
		if (v < minSiteZ)
			v = minSiteZ;

		double sVal = (v - minSiteZ) / (maxSiteZ - minSiteZ);
		if (v == 0) {
			return (zeroColor);
		}

		int red = (int) (maxColor.getRed() * sVal + minColor.getRed() * (1 - sVal));
		int green = (int) (maxColor.getGreen() * sVal + minColor.getGreen() * (1 - sVal));
		int blue = (int) (maxColor.getBlue() * sVal + minColor.getBlue() * (1 - sVal));
		c = new Color(red, green, blue);
		return (c);
	}

	private void drawGenes(Graphics2D g2d, int x1, int y, int x2, int exprTime, boolean drawNames) {
		int[] a = new int[7];
		int[] b = new int[7];
		FontMetrics metrics = g2d.getFontMetrics();
		for (ExonicGene gene : genes) {
			int trackHeight = 200;
			int gy1 = y;
			int halfGeneHeight = geneHeight / 2;
			int gmy = y + halfGeneHeight;

			int geneStart = gene.getStart(), geneEnd = gene.getEnd();
			boolean strand = reverseIt ? (gene.getStrand() == '-') : (gene.getStrand() == '+');
			int gx1 = reverseIt ? xcoord(geneEnd) : xcoord(geneStart);
			int gx2 = reverseIt ? xcoord(geneStart) : xcoord(geneEnd);

			g2d.setColor(Color.black);
			g2d.drawLine(gx1, gmy, gx2, gmy);
			arrangeArrow(a, b, strand, trackHeight, gx1, gx2, gmy);
			g2d.drawPolyline(a, b, 7);

			Iterator<Region> exons = gene.getExons();
			while (exons.hasNext()) {
				Region exon = exons.next();
				int ex1 = reverseIt ? xcoord(exon.getEnd()) : xcoord(exon.getStart());
				int ex2 = reverseIt ? xcoord(exon.getStart()) : xcoord(exon.getEnd());
				int eleft = Math.max(x1, ex1);
				int eright = Math.min(x2, ex2);

				int rectwidth = eright - eleft + 1;

				// Expression filled rectangles (if they exist)
				g2d.setColor(geneColor);
				if (expression.containsKey(gene.getName())) {
					Double val = expression.get(gene.getName()).get(exprTime);
					// Convert color
					Color eCol = expColor(val.doubleValue());
					g2d.setColor(eCol);
				}
				g2d.fillRect(eleft, gy1, rectwidth, geneHeight);
				g2d.setColor(Color.black);
				g2d.drawRect(eleft, gy1, rectwidth, geneHeight);
			}
			// Gene name
			if (drawNames) {
				g2d.setFont(new Font("Ariel", Font.BOLD, 24));
				metrics = g2d.getFontMetrics();
				int nx = strand ? gx1 : gx2;
				int ny = gmy + (3 * geneHeight);
				g2d.setColor(Color.black);
				Rectangle2D textrect = metrics.getStringBounds(gene.getName(), g2d);
				int diff = metrics.getHeight() / 2;
				String cleanName = gene.getName();
				cleanName = cleanName.replaceAll("Hoxa", "");
				cleanName = cleanName.replaceAll("Hoxb", "");
				cleanName = cleanName.replaceAll("Hoxc", "");
				cleanName = cleanName.replaceAll("Hoxd", "");
				g2d.drawString(cleanName, nx - (metrics.stringWidth(cleanName)) / 2, ny);
			}
			// Draw the axis
			// g2d.setColor(Color.lightGray);
			// g2d.setStroke(new BasicStroke(1.0f));
			// g2d.drawLine(leftBound, y, rightBound, y);
		}
	}

	private void drawRAREs(Graphics2D g2d, int currLine) {
		HashMap<String, Color> motifsProcessed = new HashMap<String, Color>();
		for (Pair<String, StrandedRegion> pair : motifs) {
			if (!motifsProcessed.containsKey(pair.car())) {
				if (pair.car().equals("DR2"))
					motifsProcessed.put(pair.car(), rareColDR2);
				else if (pair.car().equals("DR5"))
					motifsProcessed.put(pair.car(), rareColDR5);
				else
					motifsProcessed.put(pair.car(), new Color((int) (Math.random() * 256), (int) (Math.random() * 256),
							(int) (Math.random() * 256)));
			}
			Color currCol = motifsProcessed.get(pair.car());
			StrandedRegion sr = pair.cdr();
			int start = sr.getStart();
			if (sr.getStrand() == '-') {
				start = sr.getEnd();
			}

			int gx1 = xcoord(start) - (motifWidth / 2);
			g2d.setColor(currCol);
			if ((sr.getStrand() == '+' && !reverseIt) || (sr.getStrand() == '-' && reverseIt)) {
				int[] xPoints = { gx1, gx1, gx1 + motifWidth };
				int[] yPoints = { currLine, currLine + motifHeight, currLine + (motifHeight / 2) };
				g2d.fillPolygon(xPoints, yPoints, 3);
			} else {
				gx1 = xcoord(sr.getEnd()) + (motifWidth / 2);
				int[] xPoints = { gx1, gx1, gx1 - motifWidth };
				int[] yPoints = { currLine, currLine + motifHeight, currLine + (motifHeight / 2) };
				g2d.fillPolygon(xPoints, yPoints, 3);
			}
		}
	}

	private void arrangeArrow(int[] a, int[] b, boolean strand, int geneHeight, int gx1, int gx2, int my) {
		double arrowHt = 0.005 * geneHeight;
		double arrowWd = 1;
		int a1, a2, a3;

		// forward arrow
		if (strand) {
			int startX = gx1;
			a1 = startX;
			a2 = (int) Math.round(startX + (arrowWd * 6));
			a3 = (int) Math.round(startX + (arrowWd * 10));

		} else {
			// backward arrow
			int startX = gx2 + 1;
			a1 = startX;
			a2 = (int) Math.round(startX - (arrowWd * 6));
			a3 = (int) Math.round(startX - (arrowWd * 10));
		}

		a[0] = a1;
		a[1] = a1;
		a[2] = a2;
		a[3] = a2;
		a[4] = a3;
		a[5] = a2;
		a[6] = a2;

		int b1 = (int) Math.round(my);
		int b2 = (int) Math.round(my + (arrowHt * 13));
		int b3 = (int) Math.round(my + (arrowHt * 10));
		int b4 = (int) Math.round(my + (arrowHt * 16));

		b[0] = b1;
		b[1] = b2;
		b[2] = b2;
		b[3] = b3;
		b[4] = b2;
		b[5] = b4;
		b[6] = b2;
	}

	private void drawLegend(Graphics2D g2d, int x, int y) {
		g2d.setColor(rareColDR5);
		int[] x1Points = { x, x, x + 16 };
		int[] y1Points = { y, y + 16, y + 8 };
		g2d.fillPolygon(x1Points, y1Points, 3);
		g2d.setFont(new Font("Ariel", Font.BOLD, 16));
		FontMetrics metrics = g2d.getFontMetrics();
		int textY = y + 4 + (metrics.getHeight() / 2);
		g2d.drawString("DR5 RARE Motif", x + 15, textY);
		g2d.setColor(rareColDR2);
		int[] x2Points = { x, x, x + 16 };
		int[] y2Points = { y + 16, y + 32, y + 24 };
		g2d.fillPolygon(x2Points, y2Points, 3);
		textY = y + 20 + (metrics.getHeight() / 2);
		g2d.drawString("DR2 RARE Motif", x + 15, textY);
	}

	private void drawExpColorBar(Graphics2D g2d, int x, int y) {
		// Draw colors
		GradientPaint colorbar = new GradientPaint(x, y, Color.blue, x + colorbarWidth / 2, y, Color.black, false);
		g2d.setPaint(colorbar);
		g2d.fillRect(x, y, colorbarWidth / 2, colorbarHeight);
		colorbar = new GradientPaint(x + colorbarWidth / 2, y, Color.black, x + colorbarWidth, y, Color.yellow, false);
		g2d.setPaint(colorbar);
		g2d.fillRect(x + (colorbarWidth / 2), y, colorbarWidth / 2, colorbarHeight);

		// Draw border
		g2d.setPaint(Color.black);
		g2d.setColor(Color.black);
		g2d.setStroke(new BasicStroke(1.0f));
		g2d.drawRect(x, y, colorbarWidth, colorbarHeight);

		// Legend
		g2d.setFont(new Font("Ariel", Font.PLAIN, 14));
		FontMetrics metrics = g2d.getFontMetrics();
		int textY = y + colorbarHeight + (metrics.getHeight());
		g2d.drawString("0", x + (colorbarWidth / 2) - (metrics.stringWidth("0") / 2), textY);
		g2d.drawString(String.format("%.1f", minExp), x - (metrics.stringWidth(String.format(".1f", minExp)) / 2),
				textY);
		g2d.drawString(String.format("%.1f", maxExp),
				x + colorbarWidth - (metrics.stringWidth(String.format(".1f", maxExp)) / 2), textY);

		// Title
		g2d.setFont(new Font("Ariel", Font.ITALIC, 14));
		metrics = g2d.getFontMetrics();
		g2d.drawString("log2-foldchange", x + (colorbarWidth / 2) - (metrics.stringWidth("log2-foldchange") / 2),
				y - (metrics.getHeight()) / 2);
	}

	private void drawSiteColorBar(Graphics2D g2d, int x, int y) {
		// Draw colors
		GradientPaint colorbar = new GradientPaint(x, y, siteMinColor, x + colorbarWidth, y, siteMaxColor, false);
		g2d.setPaint(colorbar);
		g2d.fillRect(x, y, colorbarWidth, colorbarHeight);

		// Draw border
		g2d.setPaint(Color.black);
		g2d.setColor(Color.black);
		g2d.setStroke(new BasicStroke(3.0f));
		g2d.drawRect(x, y, colorbarWidth, colorbarHeight);

		// Legend
		g2d.setFont(new Font("Ariel", Font.PLAIN, 16));
		FontMetrics metrics = g2d.getFontMetrics();
		int textY = y + colorbarHeight + (metrics.getHeight());
		g2d.drawString(String.format("%.1f", minSiteZ), x - (metrics.stringWidth(String.format("%.1f", minSiteZ))),
				textY);
		g2d.drawString(String.format("%.1f", maxSiteZ),
				x + colorbarWidth - (metrics.stringWidth(String.format(".1f", maxSiteZ)) / 2), textY);

		// Title
		g2d.setFont(new Font("Ariel", Font.ITALIC, 16));
		metrics = g2d.getFontMetrics();
		g2d.drawString("RAR Peak (-10 log(p-value))",
				x + (colorbarWidth / 2) - (metrics.stringWidth("RAR Peak (-10 log(p-value))") / 2),
				y - (metrics.getHeight()) / 2);
	}
}
