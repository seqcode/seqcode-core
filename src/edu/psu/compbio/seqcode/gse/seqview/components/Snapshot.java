package edu.psu.compbio.seqcode.gse.seqview.components;

import java.sql.SQLException;
import java.io.File;
import java.io.IOException;

import edu.psu.compbio.seqcode.gse.seqview.SeqViewOptions;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

public class Snapshot {

    
    public static void main(String args[]) {
        String picturename = null;
        int w = 1600, h = 1200;
        boolean exit = true;
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--picture")) {
                picturename = args[++i];
            }
            if (args[i].equals("--width")) {
                w = Integer.parseInt(args[++i]);
            }
            if (args[i].equals("--height")) {
                h = Integer.parseInt(args[++i]);
            } 
            if (args[i].equals("--noexit")) {
                exit = false;
            }

        }
        if (picturename == null) {return;}
        try {
            SeqViewOptions opts = SeqViewOptions.parseCL(args);
            SeqViewStatus s = new SeqViewStatus();
            RegionPanel panel = new RegionPanel(opts, s, new File(System.getProperty("user.home")));
            File file = new File(picturename);
            panel.computeLayout(0,0,w,h);
            while (!panel.allCanPaint()) {
                Thread.yield();
            }
            panel.saveImage(file,w,h,true);
            panel.close();
		} catch (NotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
            e.printStackTrace();
        } catch (SQLException e) {
            e.printStackTrace();
        }
        if (exit) {
            System.exit(0);
        }
    }


}