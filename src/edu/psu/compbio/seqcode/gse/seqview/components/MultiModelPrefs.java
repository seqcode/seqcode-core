package edu.psu.compbio.seqcode.gse.seqview.components;

/** totally ripped off from edu.psu.compbio.seqcode.gse.viz.eye.ModelPrefs but I want
 * to include multiple panels with different types without making Tim scream
 * when he sees what I did to his code...
 *
 * Also, this version doesn't do listeners since that's more complicated.
 */

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.*;

import edu.psu.compbio.seqcode.gse.utils.models.*;
import edu.psu.compbio.seqcode.gse.viz.eye.*;

public class MultiModelPrefs extends JFrame {
	
	private Collection<? extends Model> models;
	private ArrayList<PrefsPanel<Model>> panels;
	private JButton ok, cancel;
    private RegionPanel regionpanel;
    private boolean batch=false; //allows prefs to be broadcast to all appropriate tracks in region panel
	
	public MultiModelPrefs(Collection<? extends Model> models, RegionPanel rp, boolean batchUpdate) { 
		super(batchUpdate ? "Batch Update Preferences" : "Update Track Preferences");
        this.models = models;
		regionpanel = rp;
        batch = batchUpdate;
        Container c = (Container)getContentPane();
		c.setLayout(new BorderLayout());
        panels = new ArrayList<PrefsPanel<Model>>();
		
        JPanel mainpanel = new JPanel();
        GridBagLayout gridbag = new GridBagLayout();
        mainpanel.setLayout(gridbag);
        GridBagConstraints constraints = new GridBagConstraints();        
        constraints.weightx = 1.0;
        constraints.fill = GridBagConstraints.BOTH;
        constraints.gridwidth = GridBagConstraints.REMAINDER;
        int height = 0;
        for (Model m : models) {
            if (m.getFields().isEmpty()) {
                continue;
            }
            PrefsPanel panel = new PrefsPanel<Model>(m);
            JLabel label = new JLabel(m.getClass().toString().replaceAll("^.*\\.",""));            
            height += (int)panel.getPreferredSize().getHeight() + 50;            
            gridbag.setConstraints(label,constraints);
            mainpanel.add(label);
            gridbag.setConstraints(panel,constraints);
            mainpanel.add(panel);
            panels.add(panel);
        }
        if (!panels.isEmpty()) {
            JScrollPane pane = new JScrollPane(mainpanel);
            pane.setPreferredSize(new Dimension(700,height));
            c.add(pane, BorderLayout.CENTER);
            JPanel buttons = new JPanel();
            buttons.setLayout(new FlowLayout());
            buttons.add(ok = new JButton(createOkAction()));
            buttons.add(cancel = new JButton(createCancelAction()));
            gridbag.setConstraints(buttons,constraints);
            c.add(buttons, BorderLayout.SOUTH);
            
            setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        } else {
            this.dispose();
        }
	}
	
	public Action createOkAction() { 
		return new AbstractAction("OK") { 
			public void actionPerformed(ActionEvent e) { 
				ok();
			}
		};
	}
	
	public Action createCancelAction() { 
		return new AbstractAction("Cancel") { 
			public void actionPerformed(ActionEvent e) { 
				cancel();
			}
		};
	}
	
	public void ok() { 
        for (PrefsPanel p : panels) {
            p.saveToModel();            
        }
        dispose();
        if(batch){
        	regionpanel.batchUpdateModels(models);
        }else{
        	wakeWaiters();
        }
        regionpanel.forceModelUpdate();
        regionpanel.repaint();
	}
	
	private synchronized void wakeWaiters() { 
		notifyAll();
	}
	
	public void cancel() { 
		dispose();
		wakeWaiters();
	}
	
	public void display() { 
		SwingUtilities.invokeLater(new Runnable() { 
			public void run() { 
				setLocation(100, 100);
				setVisible(true);
				pack();
			}
		});
	}   
}


