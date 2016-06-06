package org.seqcode.motifs;

import java.util.*;

import javax.imageio.ImageIO;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.*;

import org.seqcode.data.motifdb.*;
import org.seqcode.genome.Genome;
import org.seqcode.utils.*;


/** Generates a .png file with a set of motif logos in it.
 * Needs --species and --out (output file name) on the command line.
 * Reads motifs on STDIN.  Fields are tab or semicolon separated and 
 * are name, version, and an optional additional label
 *
 */

public class DrawMotifs {

    public static void main(String args[]) throws Exception {
        MarkovBackgroundModel bgModel = null;
        String bgmodelname = Args.parseString(args,"bgmodel","whole genome zero order");
        BackgroundModelMetadata md = BackgroundModelLoader.getBackgroundModel(bgmodelname,
                                                                              1,
                                                                              "MARKOV",
                                                                              Args.parseGenome(args).cdr().getDBID());
        if (md != null) {
            bgModel = BackgroundModelLoader.getMarkovModel(md);
        } else {
            System.err.println("Couldn't get metadata for " + bgmodelname);
        }
        String outfilename = Args.parseString(args,"out","motifs.png");
        Genome genome = Args.parseGenome(args).cdr();
        int species = genome.getSpeciesDBID();
        BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
        String line = null;
        WeightMatrixLoader loader = new WeightMatrixLoader();
        List<WeightMatrix> matrices = new ArrayList<WeightMatrix>();
        while ((line = reader.readLine()) != null) {
            String pieces[] = line.split("[\\t;]");
            for (WeightMatrix m : loader.query(pieces[0], pieces[1], null)) {
                if (pieces.length >= 3) {
                    m.version = m.version + " (" + pieces[2] + ")";
                }
                matrices.add(m);
            }
        }
        if (bgModel == null) {
            for (WeightMatrix m : matrices) {
                m.toLogOdds();
            }            
        } else {
            for (WeightMatrix m : matrices) {
                m.toLogOdds(bgModel);
            }            
        }
        ClusterMotifs.drawCluster(matrices,
                                  outfilename,
                                  800,
                                  200,
                                  2);



    }
    
	/**
	 * Print a motif logo
	 * @param wm
	 * @param f
	 * @param pixheight
	 */
	public static void printMotifLogo(WeightMatrix wm, File f, int pixheight, String label, boolean drawAxes){
		int pixwidth = (pixheight-WeightMatrixPainter.Y_MARGIN*3-WeightMatrixPainter.YLABEL_SIZE) * wm.length() /2 +WeightMatrixPainter.X_MARGIN*2;
		System.setProperty("java.awt.headless", "true");
		BufferedImage im = new BufferedImage(pixwidth, pixheight,BufferedImage.TYPE_INT_ARGB);
        Graphics g = im.getGraphics();
        Graphics2D g2 = (Graphics2D)g;
        g2.setRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
        WeightMatrixPainter wmp = new WeightMatrixPainter();
        g2.setColor(Color.WHITE);
        g2.fillRect(0,0,pixwidth, pixheight);
        wmp.paint(wm,g2,0,0,pixwidth,pixheight, label, drawAxes);
        try {
            ImageIO.write(im,"png",f);
        }  catch (IOException ex) {
            ex.printStackTrace();
        }
	}
	public static void printMotifLogo(WeightMatrix wm, File f, int pixheight, String label){printMotifLogo(wm, f, pixheight, label, false);}

}