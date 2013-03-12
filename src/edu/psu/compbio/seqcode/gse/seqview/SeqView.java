package edu.psu.compbio.seqcode.gse.seqview;


import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;

import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JTextField;

import edu.psu.compbio.seqcode.gse.datasets.binding.BindingEvent;
import edu.psu.compbio.seqcode.gse.datasets.general.Region;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.seqview.components.ArrayDesignSelectFrame;
import edu.psu.compbio.seqcode.gse.seqview.components.BindingEventAnnotationPanel;
import edu.psu.compbio.seqcode.gse.seqview.components.BindingScanSelectFrame;
import edu.psu.compbio.seqcode.gse.seqview.components.RegionAnnotationPanel;
import edu.psu.compbio.seqcode.gse.seqview.components.RegionListPanel;
import edu.psu.compbio.seqcode.gse.seqview.components.RegionPanel;
import edu.psu.compbio.seqcode.gse.seqview.components.SaveRegionsAsFasta;
import edu.psu.compbio.seqcode.gse.seqview.components.SeqViewOptionsFrame;
import edu.psu.compbio.seqcode.gse.seqview.components.SeqViewToolsMenu;
import edu.psu.compbio.seqcode.gse.seqview.components.SpeciesAlignFrame;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;
import edu.psu.compbio.seqcode.gse.viz.DynamicAttribute;

/**
 * SeqView is a driver class that derives a little from WarpDrive's RegionFrame. 
 * 
 * It aims to be the top-level class, running the others, whereas WarpOptionsFrame
 * was the top-level class in WarpDrive.
 * This class assumes you will view one genome version at a time (RegionFrame did not).
 * 
 * @author mahony
 *
 */
public class SeqView extends JFrame {
	protected SeqViewOptions options;
	protected SeqViewOptionsFrame optionsFrame;
	protected Genome selectedGenome; 
	protected RegionPanel panel=null;
	protected boolean imageraster;
	protected int imageheight = 1200, imagewidth = 1600;

    public SeqView(SeqViewOptions opt) {
    	options = opt;
    	selectedGenome = options;
        setTitle(options.species + " " + options.genome);
        
        //panel = new RegionPanel(opts);
        //getContentPane().add(new ImageCachingPanel(panel));
        //getContentPane().add(panel);
        setJMenuBar(createDefaultJMenuBar(options));
        setSize(600,400);
        setLocation(50,50);
        setVisible(true);
        imageraster = true;

        this.addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent arg0) {
                System.out.println("WindowClosing: " + arg0.toString());
                if(panel==null)
                	panel.handleWindowClosing();
                System.exit(0);
            }
        });
    }
    
    public RegionPanel getRegionPanel() { return panel; }
    
    private JMenuBar createDefaultJMenuBar(SeqViewOptions opts) { 
        JMenuBar jmb = new JMenuBar();
        JMenu filemenu, imagemenu, navigationmenu, displaymenu, toolsmenu; 
        JMenuItem item;
        final SeqView thisframe = this;
        final RegionPanel thispanel = panel;
        jmb.add((filemenu = new JMenu("File")));
        filemenu.add((item = new JMenuItem("Configure Tracks")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {                    
                    SeqViewOptions current = panel.getCurrentOptions();
                    SeqViewOptionsFrame frame;
                    try {
                        frame = new SeqViewOptionsFrame(panel.getCurrentOptions());
                        frame.addPainterContainer(panel);
                    } catch (NotFoundException e1) {
                        e1.printStackTrace();
                    }
                }
            });
        filemenu.add((item = new JMenuItem("Add a track from a file")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {                    
                    thispanel.addTrackFromFile();
                }
            });
        filemenu.add((item = new JMenuItem("Save Region as FASTA")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    new SaveRegionsAsFasta(thispanel.getRegion());
                }
            });
        filemenu.addSeparator();
        filemenu.add((item = new JMenuItem("Close")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    thisframe.dispose();
                }
            });
        filemenu.add((item = new JMenuItem("Exit")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) { 
                    thispanel.close();
                    System.exit(0);
                }
            });        

        jmb.add((imagemenu = new JMenu("Image")));
        imagemenu.add((item = new JMenuItem("Save Image")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    JFileChooser chooser;
                    chooser = new JFileChooser(new File(System.getProperty("user.dir")));
                    
                    int v = chooser.showSaveDialog(null);
                    if(v == JFileChooser.APPROVE_OPTION) { 
                        File f = chooser.getSelectedFile();
                        boolean raster = imageraster;
                        if (f.getAbsolutePath().matches(".svg")) {
                            raster = false;
                        }
                        try {
                            panel.saveImage(f,imagewidth,imageheight,raster);
                        } catch (IOException ex) {
                            ex.printStackTrace();
                        }
                    }
                } 
            });
        imagemenu.add((item = new JMenuItem("Image Settings")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    new ImageConfigurationFrame(thisframe);
                }
            });
        jmb.add((navigationmenu = new JMenu("Navigation")));
        final JCheckBoxMenuItem linkeditem;
        navigationmenu.add((linkeditem = new JCheckBoxMenuItem("Link to Alignments")));
        linkeditem.setSelected(false);
        linkeditem.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    if (linkeditem.isSelected()) {
                        SpeciesAlignFrame.addRegionPanel(panel);
                        panel.setRegion(panel.getRegion());
                    } else {
                        SpeciesAlignFrame.removeRegionPanel(panel);
                    }
                }
            });


        navigationmenu.add((item = new JMenuItem("Binding Scan Annotation")));
        item.addActionListener(new ActionListener()  {
                public void actionPerformed(ActionEvent e) {
                    BindingEventAnnotationPanel beap = new BindingEventAnnotationPanel(thispanel, new ArrayList<BindingEvent>());
                    JFrame f = new BindingEventAnnotationPanel.Frame(beap);
                    f.pack();
                    new BindingScanSelectFrame(thispanel.getGenome(),beap);
                }
            });

        navigationmenu.add((item = new JMenuItem("Binding Event List")));
        item.addActionListener(new ActionListener()  {
                public void actionPerformed(ActionEvent e) {
                    RegionListPanel rlp = new RegionListPanel(thispanel,null);
                    RegionListPanel.makeFrame(rlp);
                    new BindingScanSelectFrame(thispanel.getGenome(),rlp);
                }
            });               

        navigationmenu.add((item = new JMenuItem("Array Tiled Regions")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    RegionListPanel rlp = new RegionListPanel(thispanel,null);
                    RegionListPanel.makeFrame(rlp);
                    new ArrayDesignSelectFrame(thispanel.getGenome(),rlp);
                }
            });

        navigationmenu.add((item = new JMenuItem("Array Tiled Regions and Genes")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    RegionAnnotationPanel beap = new RegionAnnotationPanel(thispanel,null);
                    JFrame f = new RegionAnnotationPanel.Frame(beap);
                    f.pack();
                    new ArrayDesignSelectFrame(thispanel.getGenome(),beap);
                }
            });

        navigationmenu.add((item = new JMenuItem("Open Region List")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    JFileChooser chooser;
                    chooser = new JFileChooser(new File(System.getProperty("user.dir")));
                    
                    int v = chooser.showOpenDialog(null);
                    if(v == JFileChooser.APPROVE_OPTION) { 
                        File f = chooser.getSelectedFile();
                        java.util.List<Region> regions = RegionPanel.readRegionsFromFile(panel.getGenome(),f.getAbsolutePath());
                        RegionListPanel p = new RegionListPanel(panel,
                                                                regions);
                        RegionListPanel.makeFrame(p);
                    }
                    
                }
            });
        navigationmenu.add((item = new JMenuItem("New Region List")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    RegionListPanel p = new RegionListPanel(panel,
                                                            new ArrayList<Region>());
                    RegionListPanel.makeFrame(p);
                }
            });
        ButtonGroup group = new ButtonGroup();
        jmb.add(displaymenu = new JMenu("Display"));
        displaymenu.add((item = new JRadioButtonMenuItem("Screen")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    DynamicAttribute.getGlobalAttributes().setType(DynamicAttribute.SCREEN);
                    thisframe.repaint();
                }
            });
        item.setSelected(true);
        group.add(item);
        displaymenu.add((item = new JRadioButtonMenuItem("Display")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    DynamicAttribute.getGlobalAttributes().setType(DynamicAttribute.DISPLAY);
                    thisframe.repaint();
                }
            });
        group.add(item);
        displaymenu.add((item = new JRadioButtonMenuItem("Print")));
        item.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    DynamicAttribute.getGlobalAttributes().setType(DynamicAttribute.PRINT);
                    thisframe.repaint();
                }
            });
        group.add(item);
        jmb.add(new SeqViewToolsMenu(panel));

        return jmb;
    }

    class ImageConfigurationFrame extends JFrame implements ActionListener {
        private SeqView parent;
        private JCheckBox rasterbox;
        private JTextField widthfield, heightfield;
        private JButton okbutton, cancelbutton;

        public ImageConfigurationFrame(SeqView p) {
            parent = p;
            JLabel boxlabel = new JLabel("Configure Save-as-image");
            rasterbox = new JCheckBox("Raster Image?",parent.imageraster);
            JLabel widthlabel = new JLabel("Width");
            JLabel heightlabel = new JLabel("Height");
            widthfield = new JTextField(Integer.toString(parent.imagewidth));
            heightfield = new JTextField(Integer.toString(parent.imageheight));
            okbutton = new JButton("OK");
            cancelbutton = new JButton("Cancel");
            okbutton.addActionListener(this);
            cancelbutton.addActionListener(this);
            
            JPanel toppanel = new JPanel();
            toppanel.setLayout(new BorderLayout());
            
            JPanel buttonpanel = new JPanel();
            buttonpanel.add(okbutton);
            buttonpanel.add(cancelbutton);
            toppanel.add(buttonpanel,BorderLayout.SOUTH);

            JPanel infopanel = new JPanel();
            infopanel.setLayout(new BorderLayout());
            infopanel.add(rasterbox,BorderLayout.NORTH);
            JPanel textpanel = new JPanel();
            textpanel.setLayout(new GridLayout(2,2));
            textpanel.add(widthlabel);
            textpanel.add(widthfield);
            textpanel.add(heightlabel);
            textpanel.add(heightfield);
            infopanel.add(textpanel,BorderLayout.CENTER);
            
            toppanel.add(infopanel,BorderLayout.CENTER);

            getContentPane().add(toppanel);
            setMinimumSize(new Dimension(150,150));
            setSize(getPreferredSize());
            pack();
            setVisible(true);
        }

        public void actionPerformed (ActionEvent e) {
            if (e.getSource() == okbutton) {
                parent.imageraster = rasterbox.isSelected();
                try {
                    parent.imageheight = Integer.parseInt(heightfield.getText());
                } catch (NumberFormatException ex) {
                }
                try {
                    parent.imagewidth = Integer.parseInt(widthfield.getText());
                } catch (NumberFormatException ex) {
                }


                this.dispose();
            } else if (e.getSource() == cancelbutton) {
                this.dispose();
            }
        }     
   
    }

	
	public static void main(String args[]) throws NotFoundException, SQLException, IOException {
        SeqViewOptions opts = SeqViewOptions.parseCL(args);
        SeqView frame = new SeqView(opts);
    }
}
