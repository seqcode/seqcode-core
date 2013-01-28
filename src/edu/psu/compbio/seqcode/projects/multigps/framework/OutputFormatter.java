package edu.psu.compbio.seqcode.projects.multigps.framework;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import javax.imageio.ImageIO;

import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.multigps.experiments.ControlledExperiment;

public class OutputFormatter {

	protected Config config;
	
	public OutputFormatter(Config c){ 
		config = c;	
	}
	
    /**
     * Plot a sets of read distributions, indexed in the Map by replicates 
     */
    public void plotAllReadDistributions(Map<ControlledExperiment, List<BindingModel>> models){
		Color[] colors = {Color.black, Color.red, Color.blue, Color.green, Color.cyan, Color.orange};
		for(ControlledExperiment rep : models.keySet()){
			String replicateName = rep.getCondName()+"-"+rep.getRepName();
			String filename = config.getOutputImagesDir()+File.separator+config.getOutBase()+"_"+replicateName + "_Read_Distributions.png";
			File f = new File(filename);
			int w = 1000;
			int h = 600;
			int margin= 50;
		    BufferedImage im = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
		    Graphics g = im.getGraphics();
		    Graphics2D g2 = (Graphics2D)g;
		    g2.setRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
		    g2.setColor(Color.white);
		    g2.fillRect(0, 0, w, h);	
		    g2.setColor(Color.gray);
		    g2.drawLine(20, h-margin, w-20, h-margin);		// x-axis
		    g2.drawLine(w/2, margin, w/2, h-margin);	// y-axis    
		    g.setFont(new Font("Arial",Font.PLAIN,16));
		    for (int p=-2;p<=2;p++){
		    	g2.drawLine(w/2+p*200, h-margin-10, w/2+p*200, h-margin);	// tick  
		    	g2.drawString(""+p*200, w/2+p*200-5, h-margin+22);			// tick label
		    }
		    
		    double maxProb = 0;
		    List<BindingModel> currModels = models.get(rep);
		    for (BindingModel bm:currModels){
		    	int summit = bm.getSummit();
		    	maxProb = Math.max(maxProb, bm.probability(summit));
		    }
		    
		    for (int i=0;i<currModels.size();i++){
		    	BindingModel m = currModels.get(i);
		    	List<Pair<Integer, Double>> points = m.getEmpiricalDistribution();
			    g2.setColor(colors[i % colors.length]);
			    g2.setStroke(new BasicStroke(4));
			    for (int p=0;p<points.size()-1;p++){
			    	int x1=points.get(p).car()+w/2;
			    	int y1=(int) (h-points.get(p).cdr()/maxProb*(h-margin*2)*0.8)-margin;
			    	int x2=points.get(p+1).car()+w/2;
			    	int y2=(int) (h-points.get(p+1).cdr()/maxProb*(h-margin*2)*0.8)-margin;
			    	g2.drawLine(x1, y1, x2, y2);	    
			    }
			    g.setFont(new Font("Arial",Font.PLAIN,20));
			    g2.drawString(String.format("%d", i), w-300, i*25+margin+25);
		    }
	
		    try{
		    	ImageIO.write(im, "png", f);
		    }
		    catch(IOException e){
		    	System.err.println("Error in printing file "+filename);
		    }
		}
	}

}
