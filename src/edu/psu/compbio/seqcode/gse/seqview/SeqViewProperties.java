package edu.psu.compbio.seqcode.gse.seqview;

import java.util.*;
import java.io.*;
import java.net.URL;
import javax.swing.*;

import edu.psu.compbio.seqcode.gse.seqview.components.MultiModelPrefs;
import edu.psu.compbio.seqcode.gse.seqview.components.RegionPanel;
import edu.psu.compbio.seqcode.gse.utils.json.*;
import edu.psu.compbio.seqcode.gse.utils.models.*;
import edu.psu.compbio.seqcode.gse.viz.components.RegexFileFilter;

/** Base class for SeqView (paintable and model) properties.
 * Uses Tim's JSON model stuff to handle serializing and deserializing.
 */
public abstract class SeqViewProperties extends Model {

    // controls whether some exception stack traces are printed
    private boolean debugging = true;
    /* true when a window is currently open to configure this SeqViewProperties.
       Used to prevent multiple windows from working on the same properties at once.
    */
    private Boolean configuring;

    public SeqViewProperties() {
        super();
        configuring = Boolean.FALSE;
    }

    static class DoneConfiguring implements Runnable {
        private Collection<SeqViewProperties> p;
        private MultiModelPrefs mmp;
        public DoneConfiguring(Collection<SeqViewProperties> p, MultiModelPrefs mmp ) {this.p = p;this.mmp = mmp;} 
        public void run() {
            synchronized (mmp) {
                boolean done = false;
                while (!done) {
                    try {
                        mmp.wait(); 
                        for (SeqViewProperties prop : p) {
                            prop.configuring = false;   
                        }
                        done = true;
                    } catch (InterruptedException e) {
                        // ignore it and go back to sleeping
                    }
                }
            }
        }            
    }
    public static void configure(Collection<? extends SeqViewProperties> props, RegionPanel regionpanel, boolean batchUpdate) {
        Collection<SeqViewProperties> touse = new ArrayList<SeqViewProperties>();
        for (SeqViewProperties p : props) {
            synchronized(p) {
                if (p.configuring) {
                    continue;
                }
                p.configuring = true;            
            }
            touse.add(p);
        }
        
        MultiModelPrefs frame = new MultiModelPrefs(touse,regionpanel, batchUpdate);
        (new Thread(new DoneConfiguring(touse,frame))).start();
        frame.setSize(frame.getPreferredSize());
        frame.setVisible(true);
        frame.pack();
        frame.setLocationRelativeTo(null);
    }
    /** 
     * Since both PaintableProperties and ModelProperties inherit from here, we don't want to
     * confuse the two types of properties.  So each subclass needs to define the suffix that its
     * files use.
     */
    public abstract String fileSuffix();
    public String defaultName() {
        return getClass().toString().replaceAll("^.*\\.","") + ".defaults." + fileSuffix();
    }
    public File currentFile() {
        return new File(System.getProperty("user.dir") + System.getProperty("file.separator") + defaultName());
    }
    public File defaultFile() {
        return new File(System.getProperty("user.home") + System.getProperty("file.separator") + defaultName());
    }
    /** Displays a JFileChooser to select the name of the file to which you want to save the properties
     */
    public void saveToFile() {
        JFileChooser chooser = new JFileChooser();
        RegexFileFilter filter = new RegexFileFilter(".*\\." + fileSuffix() + "$", "SeqView Prefs",true);
        chooser.setFileFilter(filter);
        chooser.setSelectedFile(defaultFile());
        int returnVal = chooser.showSaveDialog(null);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            File f = chooser.getSelectedFile();
            saveToFile(f);
        }
    }
    /** Saves these properties to the specified file
     */
    public void saveToFile(File f) {
        JSONObject json = asJSON();
        try {
            PrintWriter pw = new PrintWriter(f);
            pw.println(json.toString());
            pw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    /** Displays a JFileChooser to select the file from which you want to load the properties
     */
    public void loadFromFile() {
        JFileChooser chooser = new JFileChooser();
        RegexFileFilter filter = new RegexFileFilter(".*\\." + fileSuffix() + "$", "SeqView Prefs",true);
        chooser.setFileFilter(filter);
        chooser.setSelectedFile(defaultFile());
        int returnVal = chooser.showOpenDialog(null);
        if (returnVal == JFileChooser.APPROVE_OPTION) {
            File f = chooser.getSelectedFile();
            try {
                loadFromFile(f);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }
    /** Load these properties from the specified file
     */
    public void loadFromFile(File f) throws FileNotFoundException {
        loadFromStream(new FileReader(f));
    }
    public void loadFromStream(InputStreamReader reader) {
        try {
            BufferedReader breader = new BufferedReader(reader);
            JSONObject json = new JSONObject(breader.readLine());
            breader.close();
            reader.close();
            setFromJSON(json);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    /** Load default values for properties.  This looks for a file of the form
     * FOO.defaults.svpp (seq view paintable properties) where FOO is the runtime class name in
     * - current directory
     * - home directory
     * - the whole classpath
     */
    public void loadDefaults() {
        File name = currentFile();
        try {
            if (name.exists()) {
                loadFromFile(name);
            } else {
                name = defaultFile();
                if (name.exists()) {
                    loadFromFile(name);
                } else {
                    ClassLoader cl = ClassLoader.getSystemClassLoader();
                    URL url = cl.getResource(defaultName());
                    if (url != null) {
                        loadFromStream(new InputStreamReader(url.openStream()));
                    }
                }
            }
        } catch (Exception e) {
            if (debugging) {
                e.printStackTrace();
            }
        }
    }

}