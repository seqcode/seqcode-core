package org.seqcode.viz.genomicplot;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.seqcode.data.motifdb.WeightMatrix;
import org.seqcode.deepseq.experiments.Sample;
import org.seqcode.genome.location.Gene;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.genome.location.StrandedRegion;
import org.seqcode.gsebricks.verbs.location.RefGeneGenerator;
import org.seqcode.gsebricks.verbs.motifs.WeightMatrixScoreProfile;
import org.seqcode.gsebricks.verbs.motifs.WeightMatrixScorer;
import org.seqcode.gseutils.Pair;
import org.seqcode.viz.paintable.AbstractPaintable;

public class ChipSeqFigurePaintable extends FigurePaintable {

	private int topBorder = 50, bottomBorder = 50;
	private int leftBorder = 25, rightBorder = 25;
	private List<String> exptNames;
	private Map<String, AbstractPaintable> exptPainters = new HashMap<String, AbstractPaintable>();
	private RefGeneGenerator<Region> geneGen;
	protected List<List<Gene>> geneSets = new ArrayList<List<Gene>>();
	protected List<String> geneSetNames = new ArrayList<String>();

	public ChipSeqFigurePaintable(FigureOptions opts) {
		options = opts;
		reverseIt = options.reverseOrder;
		chr = options.gRegion.getChrom();
		rstart = options.gRegion.getStart();
		rstop = options.gRegion.getEnd();
		exptNames = options.exptNames;

		// Initialize genes
		if (options.useDBGenes) {
			geneSetNames.add("Reference");
			ArrayList<Gene> geneSet = new ArrayList<Gene>();
			geneGen = new RefGeneGenerator<Region>(options.genome, "refGene");
			geneGen.retrieveExons(true);
			geneGen.setWantAlias(true);
			Iterator<Gene> gi = geneGen.execute(options.gRegion);
			while (gi.hasNext()) {
				geneSet.add(gi.next());
			}
			geneSets.add(geneSet);
		}
		if (options.transcriptGTF != null) {// Load GTF
			geneSetNames.add(options.transcriptGTF.getName());
			ArrayList<Gene> geneSet = new ArrayList<Gene>();
			geneSet.addAll(loadGenes(options.transcriptGTF, options.gRegion));
			geneSets.add(geneSet);
		}
		// layout.setRegions(genes);

		// Initialize thin painters
		for (String t : exptNames) {
			List<Pair<Point, Point>> inters = options.experiments.get(t).inters;
			Sample sample = options.experiments.get(t).exptSample;
			ArrayList<Region> sites = options.experiments.get(t).peaks;
			if (options.experiments.get(t).preFormattedDataFile != null) {
				// Initialize a DiffSeqExptPaintable, even for non-diff
				// experiments
				DiffSeqExptPaintable dp = new DiffSeqExptPaintable(options.gconfig, options.gRegion,
						options.experiments.get(t).preFormattedDataFile, options.experiments.get(t).diffWinWidth,
						options.experiments.get(t).diffWinStep, options.experiments.get(t).scaling);
				dp.setMaxLogFold(options.experiments.get(t).yMax);
				if (options.experiments.get(t).isDiff)
					dp.setMinLogFold(options.experiments.get(t).yMax);
				else
					dp.setMinLogFold(0);
				dp.setReverse(options.reverseOrder);
				dp.setNegColor(options.diffExptNegColor);
				dp.setPosColor(options.diffExptPosColor);
				exptPainters.put(t, dp);
			} else if (sample != null) {
				if (options.experiments.get(t).isDiff && options.experiments.get(t).baseExpt != null) {
					Sample otherSample = options.experiments.get(t).baseExptSample;
					// Initialize a DiffSeqExptPaintable
					DiffSeqExptPaintable dp = new DiffSeqExptPaintable(options.gconfig, options.gRegion, sample,
							otherSample, options.experiments.get(t).diffWinWidth,
							options.experiments.get(t).diffWinStep, options.experiments.get(t).scaling);
					dp.setMaxLogFold(options.experiments.get(t).yMax);
					dp.setMinLogFold(options.experiments.get(t).yMax);
					dp.setReverse(options.reverseOrder);
					dp.setNegColor(options.diffExptNegColor);
					dp.setPosColor(options.diffExptPosColor);
					exptPainters.put(t, dp);
				} else {
					// Initialize the ThinOverlapPaintable
					ThinOverlapPaintable tp = new ThinOverlapPaintable(options.gRegion, sites, sample, options.readExt,
							options.experiments.get(t).pairedReads, options.experiments.get(t).inters);
					tp.setReverse(options.reverseOrder);
					tp.setMaxOverlap(options.experiments.get(t).yMax);
					tp.setLoopColor(options.loopColor);
					tp.setInterColor(options.interColor);
					tp.setBgColor(options.experiments.get(t).exptBgColor);
					tp.setBgThick(options.exptBgThick);
					tp.setHighlightColor(options.experiments.get(t).exptPeakColor);
					tp.setHighlightThick(options.exptPeakThick);
					tp.setFilledColumns(options.experiments.get(t).filledColumns);
					exptPainters.put(t, tp);
				}
			}
		}

	}

	public void paintItem(Graphics g, int x1, int y1, int x2, int y2) {
		Graphics2D g2d = (Graphics2D) g;
		FontMetrics metrics = g2d.getFontMetrics();
		int screenSizeX = x2 - x1;
		int screenSizeY = y2 - y1;
		g2d.setColor(Color.white);
		g2d.fillRect(0, 0, screenSizeX, screenSizeY);
		topBound = topBorder;
		bottomBound = screenSizeY - bottomBorder;
		leftBound = leftBorder + 30;
		rightBound = screenSizeX - rightBorder;
		baseLine = bottomBound - 200;

		int offset = 0;
		for (String t : exptNames) {
			AbstractPaintable p = exptPainters.get(t);
			p.paintItem(g, leftBound, baseLine - offset - options.experiments.get(t).exptTrackHeight, rightBound,
					baseLine - offset);

			// Draw motifs
			if (options.experiments.get(t).motif != null) {
				drawMotifs(g2d, baseLine - offset, options.experiments.get(t).motif,
						options.experiments.get(t).motifThres);
			}

			// Draw the axis
			g2d.setColor(Color.lightGray);
			g2d.setStroke(new BasicStroke(1.0f));
			// g2d.drawLine(leftBound, baseLine-offset-motifHeight, rightBound,
			// baseLine-offset-motifHeight);

			// Experiment label
			if (options.drawExptLabels) {
				g2d.setColor(Color.gray);
				g2d.setFont(new Font("Ariel", Font.BOLD, options.labelFontSize));
				metrics = g2d.getFontMetrics();
				AffineTransform oldtrans = g2d.getTransform();
				AffineTransform newtrans = new AffineTransform();
				newtrans.translate(leftBound - (metrics.getHeight()),
						baseLine - offset + metrics.getHeight() - (options.experiments.get(t).exptTrackHeight / 2));
				newtrans.rotate(Math.toRadians(-90));
				g2d.setTransform(newtrans);
				g2d.drawString(t, (-1 * metrics.stringWidth(t)) / 2, 0);
				g2d.setTransform(oldtrans);
			}
			offset += (options.experiments.get(t).exptTrackHeight + options.motifHeight + 10);
		}
		// Gene tracks
		offset = options.motifHeight;
		for (int i = 0; i < geneSetNames.size(); i++) {
			offset += drawGenes(g2d, x1, baseLine + offset + options.geneHeight, x2, geneSets.get(i),
					options.drawGeneLabels, geneSetNames.get(i));
		}

		// Draw some coordinates
		g2d.setColor(Color.black);
		g2d.setFont(new Font("Ariel", Font.PLAIN, options.fontSize));
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
	}

	private void drawMotifs(Graphics2D g2d, int currLine, WeightMatrix wm, double thres) {
		ArrayList<StrandedRegion> motifLocs = new ArrayList<StrandedRegion>();
		WeightMatrixScorer scorer = new WeightMatrixScorer(wm);
		WeightMatrixScoreProfile profiler = scorer.execute(options.gRegion);
		for (int z = 0; z < options.gRegion.getWidth(); z++) {
			double currScore = profiler.getMaxScore(z);
			if (currScore >= thres) {
				motifLocs.add(new StrandedRegion(options.genome, chr, rstart + z, rstart + z + wm.length(),
						profiler.getMaxStrand(z)));
			}
		}

		for (StrandedRegion sr : motifLocs) {
			Color currCol = options.motifColor;
			int start = sr.getStart();
			if (sr.getStrand() == '-') {
				start = sr.getEnd();
			}

			int gx1 = xcoord(start) - (options.motifWidth / 2);
			g2d.setColor(currCol);
			if ((sr.getStrand() == '+' && !reverseIt) || (sr.getStrand() == '-' && reverseIt)) {
				int[] xPoints = { gx1, gx1, gx1 + options.motifWidth };
				int[] yPoints = { currLine, currLine + options.motifHeight, currLine + (options.motifHeight / 2) };
				g2d.fillPolygon(xPoints, yPoints, 3);
			} else {
				gx1 = xcoord(sr.getEnd()) + (options.motifWidth / 2);
				int[] xPoints = { gx1, gx1, gx1 - options.motifWidth };
				int[] yPoints = { currLine, currLine + options.motifHeight, currLine + (options.motifHeight / 2) };
				g2d.fillPolygon(xPoints, yPoints, 3);
			}
		}
	}

}
