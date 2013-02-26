package edu.psu.compbio.seqcode.gse.seqview.paintable;

import java.awt.*;
import java.util.*;

import edu.psu.compbio.seqcode.gse.datasets.*;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.ChipChipData;
import edu.psu.compbio.seqcode.gse.datasets.chipchip.GenericExperiment;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.ewok.nouns.*;
import edu.psu.compbio.seqcode.gse.seqview.model.ChipChipDataModel;
import edu.psu.compbio.seqcode.gse.seqview.model.ChipChipScaleModel;
import edu.psu.compbio.seqcode.gse.seqview.model.Model;
import edu.psu.compbio.seqcode.gse.utils.*;

public class ChipChipIntensitiesPainter extends ChipChipAdvPainter {
    private ChipChipDataModel model;
    public ChipChipIntensitiesPainter (ChipChipData data, ChipChipDataModel model) {
        super(data,model);
        this.model = model;
    }
    
    public void paintItem(Graphics2D g, 
                          int x1, int y1, 
                          int x2, int y2) {
        //        super.paintItem(g,x1,y1,x2,y2);        
        if (model.getGenericExperiment() instanceof ChipChipData) {
            double maxval = 65535;
            Region region = model.getRegion();
            int rs = region.getStart(), re = region.getEnd();
            ChipChipData data = (ChipChipData)model.getGenericExperiment();
            int circleradius = 1;
            int circlediameter = 1;
            for (int i = 0; i < data.getCount(); i++) {
                int px = getXPos(data.getPos(i),
                                 rs,re,x1,x2);
                for (int j = 0; j < data.getReplicates(i); j++) {
                    int ippy = getYPos(data.getIP(i,j),
                                       0,maxval,
                                       y1,y2,getProperties().IntensitiesOnLogScale);
                    int wcepy = getYPos(data.getWCE(i,j),
                                        0,maxval,
                                        y1,y2,getProperties().IntensitiesOnLogScale);
                    g.setColor(Color.GREEN);
                    g.fillOval(px-circleradius,ippy-circleradius,circlediameter,circlediameter);
                    g.setColor(Color.RED);
                    g.fillOval(px-circleradius,wcepy-circleradius,circlediameter,circlediameter);
                }
            }
        }
        if (getProperties().DrawTrackLabel) {
            g.setColor(Color.BLACK);
            g.drawString(getLabel(),x1,y1 + g.getFont().getSize() * 2);    	
        }
    }
}

