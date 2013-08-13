package edu.psu.compbio.seqcode.gse.seqview;

import java.util.*;
import java.util.prefs.*;
import java.io.*;
import java.sql.SQLException;

import edu.psu.compbio.seqcode.gse.datasets.chipchip.ExptNameVersion;
import edu.psu.compbio.seqcode.gse.datasets.expression.Experiment;
import edu.psu.compbio.seqcode.gse.datasets.motifs.*;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqAlignment;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqDataLoader;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqExpt;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocatorMatchedExpt;
import edu.psu.compbio.seqcode.gse.datasets.species.Genome;
import edu.psu.compbio.seqcode.gse.datasets.species.Organism;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

public class SeqViewOptions {

	/**
	 * Constants for accessing settings
	 */
	private static final String WINDOW_WIDTH = "WINDOW_WIDTH";
	private static final String WINDOW_HEIGHT = "WINDOW_HEIGHT";
	private static final String WINDOW_IS_CENTERED = "WINDOW_IS_CENTERED";
	private static final String WINDOW_TOP_LEFT_X = "WINDOW_X";
	private static final String WINDOW_TOP_LEFT_Y = "WINDOW_Y";
	
	/**
	 * Constants for default values and limits on settings
	 */
	public static final int DEFAULT_WINDOW_WIDTH = 900;
	public static final int DEFAULT_WINDOW_HEIGHT = 650;
	public static final int MIN_WINDOW_WIDTH = 400;
	public static final int MAX_WINDOW_WIDTH = 2400;
	public static final int MIN_WINDOW_HEIGHT = 400;
	public static final int MAX_WINDOW_HEIGHT = 1600;
	public static final boolean DEFAULT_WINDOW_IS_CENTERED = true;
	public static final int DEFAULT_TOP_LEFT_X = 50;
	public static final int DEFAULT_TOP_LEFT_Y = 50;
	public static final int MIN_TOP_LEFT_X = 0;
	public static final int MAX_TOP_LEFT_X = 1000;
	public static final int MIN_TOP_LEFT_Y = 0;
	public static final int MAX_TOP_LEFT_Y = 700;
	public static final int STATUS_BAR_HEIGHT=16;
	
	// General display settings
	private int preferredWindowWidth;
	private int preferredWindowHeight;
	private boolean isWindowCentered;
	private int preferredWindowTopLeftX;
	private int preferredWindowTopLeftY;
	
	// General connection info
    public Genome genome;

    // where to start the display.
    // Either use (chrom,start,stop), gene, position (which will be parsed
    // into either chrom/start/stop or gene), or a regionListFile
    public String chrom, gene, position, regionListFile;
    public int start, stop;

    // tracks to paint and their options
    public boolean hash=true, relative=false, seqletters=true, polya=false, gccontent=false, pyrpurcontent=false, cpg=false, regexmatcher=false;
    public ArrayList<String> genes, ncrnas, otherannots;
    public ArrayList<ExptNameVersion> agilentdata;
    public ArrayList<WeightMatrix> motifs;
    public ArrayList<Experiment> exprExperiments;
    public ArrayList<SeqLocatorMatchedExpt> seqExpts;
    // filename to label mappings.  These are loaded from a file and the data held statically
    public HashMap<String,String> regionTracks, regexes;
    
    // options for saving image to a file
    public boolean saveimage;
    public String filename;

    // startup-only options
    public boolean seqHistogramPainter = true;

    /* These constants correspond to the different input arrays.  They are used
       in WarpPaintable to store the type of input that created the paintable 
       (this must be set by the creator, currently done in RegionPanel).  This information
       is useful for removing a painter because it lets you figure out which field of
       the SeqViewOptions class should be modified */
    public static final int BINDINGSCAN = 1,
        GENES = 2,
        NCRNAS = 3,
        OTHERANNOTS = 4,
        AGILENTDATA = 5,
        MOTIFS = 11,
        PEAKS = 12,
        SEQLETTERS = 13,
        GCCONTENT = 14,
        EXPRESSION = 15,
        CPG = 16,
        MOTIFSCANS = 17,
        REGEXMATCHER = 18,
        REGIONTRACKS = 19,
        PYRPURCONTENT = 20,
        CHIPSEQANALYSES = 21,
    	SEQDATA = 22;
    
    
    public SeqViewOptions(Genome g){
    	init(g);
    }
    public SeqViewOptions(String gname) {
        try {
        	genome = Organism.findGenome(gname);
        } catch (NotFoundException e) {
            e.printStackTrace();
            throw new IllegalArgumentException(gname);
        }
        init(genome);
    }

    public SeqViewOptions() {
    	String genomeStr="mm9";
        try {        
            ResourceBundle res = ResourceBundle.getBundle("defaultgenome");
            genomeStr = res.getString("genome"); 
        } catch (MissingResourceException e) {
            // who cares, we're just getting defaults
        } catch (Exception e) {
            // ditto
        }
        
        try {
        	genome = Organism.findGenome(genomeStr);
        } catch (NotFoundException e) {
            e.printStackTrace();
            throw new IllegalArgumentException(genomeStr);
        }
        init(genome);
    }

    public void init(Genome g){
    	this.loadOptions();
    	genome = g;
        List<String> chroms = genome.getChromList();
        java.util.Collections.sort(chroms);
        chrom = chroms.get(0);
        start = 10000;
        stop = 20000;
        gene = null;
        hash = true;
        seqHistogramPainter = true;
        genes = new ArrayList<String>();
        ncrnas = new ArrayList<String>();
        otherannots = new ArrayList<String>();
        agilentdata= new ArrayList<ExptNameVersion>();
        seqExpts = new ArrayList<SeqLocatorMatchedExpt>();
        motifs = new ArrayList<WeightMatrix>();
        exprExperiments = new ArrayList<Experiment>();
        regionTracks = new HashMap<String,String>();
        regexes = new HashMap<String,String>();
    }
    /* adds options from this into union.  For lists, it generates the 
       union.  For other settings, this takes priority */
    public void mergeInto(SeqViewOptions union) throws IllegalArgumentException {
        if (union == null) {throw new NullPointerException("Must supply options to mergeInto");}
        if (genome == null) {throw new NullPointerException("Tried to call mergeInto when genome is null");}
        if (chrom != null) {
            union.chrom = chrom;
            union.start = start;
            union.stop = stop;
        }
        if (gene != null) {
            union.gene = gene;
        }
        union.hash = hash;
        union.relative = relative;
        union.gccontent = gccontent;
        union.pyrpurcontent = pyrpurcontent;
        union.cpg = cpg;
        union.gene = gene;
        union.regexmatcher = regexmatcher;
        union.seqletters = seqletters;
        union.seqHistogramPainter = (seqHistogramPainter && union.seqHistogramPainter);
        mergeInto(genes,union.genes);
        mergeInto(ncrnas,union.ncrnas);
        mergeInto(otherannots,union.otherannots);
        mergeInto(agilentdata,union.agilentdata);
        mergeInto(seqExpts,union.seqExpts);
        mergeInto(motifs,union.motifs);
        mergeInto(exprExperiments,union.exprExperiments);
        mergeInto(regionTracks,union.regionTracks);
        mergeInto(regexes, union.regexes);
    }

    public void mergeInto(ArrayList source, ArrayList target) {
        for (int i = 0; i < source.size(); i++) {
            if (!target.contains(source.get(i))) {
                target.add(source.get(i));
            }
        }        
    }
    public void mergeInto(Map source, Map target) {
        for (Object k : source.keySet()) {
            if (!target.containsKey(k)) {
                target.put(k,source.get(k));
            }
        }
    }

    /* deletes options from this that are also present in other.
       Doesn't mess with chrom, start, stop, or gene 
     */
    public void differenceOf(SeqViewOptions other) {
        if (genome != null &&
            other.genome != null &&
            !other.genome.equals(genome)) {
            throw new IllegalArgumentException("Genome must match in SeqViewOptions.mergeInto");
        }
        if (other.position != null && !other.gene.equals("")) {
            position = other.position;
        }
        if (other.regionListFile != null) {
            regionListFile = other.regionListFile;
        }
        if (other.gene != null && !other.gene.equals("")) {
            gene = other.gene;
        }
        hash = hash && (!other.hash);        
        relative = relative || other.relative;
        gccontent = gccontent && (!other.gccontent);
        pyrpurcontent = pyrpurcontent && (!other.pyrpurcontent);
        cpg = cpg && (!other.cpg);
        regexmatcher = regexmatcher && (!other.regexmatcher);        
        seqletters = seqletters || other.seqletters;
        seqHistogramPainter = (seqHistogramPainter && other.seqHistogramPainter);
        differenceOf(genes,other.genes);
        differenceOf(ncrnas,other.ncrnas);
        differenceOf(otherannots,other.otherannots);
        differenceOf(agilentdata,other.agilentdata);
        differenceOf(seqExpts,other.seqExpts);
        differenceOf(motifs,other.motifs);
        differenceOf(exprExperiments, other.exprExperiments);
        differenceOf(regionTracks,other.regionTracks);
        differenceOf(regexes,other.regexes);
    }

    public void differenceOf(ArrayList removeFrom, ArrayList other) {
        removeFrom.removeAll(other);
    }
    public void differenceOf(Map removeFrom, Map other) {
        for (Object k : other.keySet()) {
            removeFrom.remove(k);
        }
    }

    public SeqViewOptions clone() {
        SeqViewOptions o = new SeqViewOptions();
        o.genome = genome;
        o.chrom = chrom;
        o.gene = gene;
        o.position = position;
        o.regionListFile = regionListFile;
        o.start = start;
        o.stop = stop;
        o.hash = hash;
        o.gccontent = gccontent;
        o.pyrpurcontent = pyrpurcontent;
        o.cpg = cpg;
        o.relative = relative;
        o.seqletters = seqletters;
        o.regexmatcher = regexmatcher;
        o.seqHistogramPainter = seqHistogramPainter;
        o.genes = (ArrayList<String>) genes.clone();
        o.ncrnas = (ArrayList<String>) ncrnas.clone();
        o.otherannots = (ArrayList<String>) otherannots.clone();
        o.agilentdata = (ArrayList<ExptNameVersion>)agilentdata.clone();
        o.seqExpts = (ArrayList<SeqLocatorMatchedExpt>)seqExpts.clone();
        o.motifs = (ArrayList<WeightMatrix>)motifs.clone();
        o.exprExperiments = (ArrayList<Experiment>)exprExperiments.clone();
        o.regionTracks = (HashMap<String,String>)regionTracks.clone();
        o.regexes = (HashMap<String,String>)regexes.clone();
        return o;
    }

    /* Fills in a SeqViewOptions from command line arguments */
    public static SeqViewOptions parseCL(String[] args) throws NotFoundException, SQLException, IOException {
        SeqViewOptions opts = new SeqViewOptions();
        String genomeStr=null, speciesStr=null;
        
        //TODO: restore weight matrices
        //WeightMatrixLoader wmloader = new WeightMatrixLoader();
        
        SeqDataLoader seqloader = new SeqDataLoader();

        try {        
            ResourceBundle res = ResourceBundle.getBundle("defaultgenome");
            speciesStr = res.getString("species"); 
            genomeStr = res.getString("genome"); 
        } catch (MissingResourceException e) {
            // who cares, we're just getting defaults
        } catch (Exception e) {
            // ditto
        }
        
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--species")) {
                speciesStr = args[++i];
                if (speciesStr.indexOf(';') != -1) {
                    String[] pieces = speciesStr.split(";");
                    speciesStr = pieces[0];
                    genomeStr = pieces[1];
                }
            }        
            if (args[i].equals("--genome") || args[i].equals("--genomeversion")) {
                genomeStr = args[++i];
            }
            if (args[i].equals("--oldchipseq")) {
                opts.seqHistogramPainter = false;
                System.err.println("Will use old ChipSeq painters");
            }

        }
        try {
        	Organism organism = null;
            if (speciesStr != null && genomeStr != null) {
                organism = new Organism(speciesStr);
                opts.genome = organism.getGenome(genomeStr);
            }
            for (int i = 0; i < args.length; i++) {
                if (args[i].equals("--chrom") || args[i].equals("--region")) {
                    if (args[i+1].matches(".+:.+\\-.+")) {
                        int colon = args[i+1].indexOf(':');
                        int dash = args[i+1].indexOf('-');
                        opts.chrom = args[i+1].substring(0, colon);
                        String startstring = args[i+1].substring(colon+1,dash);
                        String stopstring = args[i+1].substring(dash+1);
                        opts.start = stringToNum(startstring);
                        opts.stop = stringToNum(stopstring);
                        i++;
                    } else {
                        opts.chrom = args[i+1];
                    }
                }
                if (args[i].equals("--saveimage")) {
                    opts.saveimage = true;
                    opts.filename = args[++i];
                }
                if (args[i].equals("--gene")) {
                    opts.gene = args[++i];
                }
                if (args[i].equals("--hash")) { 
                    opts.hash = true;
                }
                if (args[i].equals("--gccontent")) {
                    opts.gccontent = true;
                }
                if (args[i].equals("--pyrpurcontent")) {
                    opts.pyrpurcontent = true;
                }
                if (args[i].equals("--cpg")) {
                    opts.cpg = true;
                }
                if (args[i].equals("--repeatmasked")) {
                    opts.otherannots.add("RepeatMasker");
                }
                if (args[i].equals("--repeatmasker")) {
                    opts.otherannots.add("RepeatMasker");
                }
                if (args[i].equals("--cpgislands")) {
                    opts.otherannots.add("CpGIslands");
                }
                if (args[i].equals("--relative")) {
                    opts.relative = true;
                }
                if (args[i].equals("--genes")) {
                    opts.genes.add(args[++i]);;
                }            
                
                if (args[i].equals("--seqexpt")) {
                	//This is a really awkward (and possibly bug-filled) way to load experiments, 
                	//but hopefully loading from command-line is rare anyway
                    String pieces[] = args[++i].split(";");
                    SeqLocator loc=null;
                    if (pieces.length == 2) {
                        loc = new SeqLocator(pieces[0], pieces[1]);
                    } else if (pieces.length >= 3) {
                        Set<String> repnames = new HashSet<String>();
                        for (int j = 1; j < pieces.length - 1; j++) {
                            repnames.add(pieces[j]);
                        }
                        loc = new SeqLocator(pieces[0], repnames, pieces[pieces.length-1]);
                    } else {
                        System.err.println("Couldn't parse --seqexpt " + args[i]);
                    }
                    if(loc!=null){
                    	SeqDataLoader loader = new SeqDataLoader();
                    	Collection<SeqAlignment> aligns = loader.loadAlignments(loc, opts.genome);
                    	SeqExpt expt = null;
                    	for(SeqAlignment a : aligns){
                    		if(expt==null)
                    			expt = a.getExpt();
                    		else if(!expt.equals(a.getExpt()))
                    			System.err.println("SeqLocator "+loc.toString()+" returned multiple SeqExpts.\tUsing: "+expt.getDBID()+" but something needs to be fixed!");
                    	}
                    	List<SeqExpt> expts = new ArrayList<SeqExpt>();
                    	expts.add(expt);
                    	opts.seqExpts.add(new SeqLocatorMatchedExpt(expts, loc));
                    }
                }
                                
                if (args[i].equals("--sgdOther")) {
                    opts.otherannots.add("sgdOther");
                }
                if (args[i].equals("--regionList")) {
                    opts.regionListFile = args[++i];
                }
                if (args[i].equals("--otherannot")) {
                    opts.otherannots.add(args[++i]);
                }
                /*
                if (args[i].equals("--wm")) {
                    String[] pieces = args[++i].split(";");
                    Organism thisorg = organism;
                    if (pieces.length == 3) {
                        thisorg = new Organism(pieces[2]);
                    }                
                    WeightMatrix matrix = wmloader.query(thisorg.getDBID(),pieces[0],pieces[1]);
                    opts.motifs.add(matrix);
                }*/
                if (args[i].equals("--regex")) {
                    opts.regexmatcher = true;
                    String[] pieces = args[++i].split(";");
                    if (pieces.length == 1) {
                        opts.regexes.put("regex" + i, pieces[0]);
                    } else if (pieces.length >= 2) {
                        opts.regexes.put(pieces[1],pieces[0]);
                    }
                }
                if (args[i].equals("--fileTrack")) {
                    String pieces[] = args[++i].split(";");
                    if (pieces.length == 1) {
                        opts.regionTracks.put(pieces[0],pieces[0]);
                    } else {
                        opts.regionTracks.put(pieces[0],pieces[1]);
                    }
                }


            }
        } finally {
            //wmloader.close();
            seqloader.close();
        }

        return opts;
    }


    /**
     * Import the options from an input stream
     * @param is the input stream from which to read the options
     * @throws IOException if reading from the specified output stream results in an IOException.
     * @throws InvalidPreferencesFormatException Data on input stream does not constitute a valid XML document with the mandated document type.
     */
    public void importOptions(InputStream is) throws IOException, SeqViewException {
    	try {
			Preferences.importPreferences(is);
		} 
    	catch (InvalidPreferencesFormatException ipfex) {
    		throw new SeqViewException(ipfex);
		}
    	this.loadOptions();
    	
    	//TODO the WarpOptionsPane needs to be updated after these options are
    	//imported
    }
    
    
    /**
     * Export the options to an output stream
     * @param os the output stream to which to write the options
     * @throws IOException if writing to the specified output stream results in an IOException.
     * @throws BackingStoreException if preference data cannot be read from backing store.
     */
    public void exportOptions(OutputStream os) throws IOException, SeqViewException {
    	this.saveOptions();
    	try {
    		Preferences prefs = Preferences.userNodeForPackage(this.getClass());
    		prefs.exportSubtree(os);
    	}
    	catch (BackingStoreException bsex) {
    		throw new SeqViewException(bsex);
    	}
    }
    
    
    /**
     * Save these options to the system
     *
     */
    public void saveOptions() {
    	Preferences prefs = Preferences.userNodeForPackage(this.getClass());
    	prefs.putInt(WINDOW_WIDTH, preferredWindowWidth);
    	prefs.putInt(WINDOW_HEIGHT, preferredWindowHeight);
    	prefs.putBoolean(WINDOW_IS_CENTERED, isWindowCentered);
    	prefs.putInt(WINDOW_TOP_LEFT_X, preferredWindowTopLeftX);
    	prefs.putInt(WINDOW_TOP_LEFT_Y, preferredWindowTopLeftY);
    }
    
    
    /**
     * Load the options from the system
     *
     */
    public void loadOptions() {
    	Preferences prefs = Preferences.userNodeForPackage(this.getClass());
    	preferredWindowWidth = prefs.getInt(WINDOW_WIDTH, DEFAULT_WINDOW_WIDTH);
    	preferredWindowHeight = prefs.getInt(WINDOW_HEIGHT, DEFAULT_WINDOW_HEIGHT);
    	isWindowCentered = prefs.getBoolean(WINDOW_IS_CENTERED, DEFAULT_WINDOW_IS_CENTERED);
    	preferredWindowTopLeftX = prefs.getInt(WINDOW_TOP_LEFT_X, DEFAULT_TOP_LEFT_X);
    	preferredWindowTopLeftY = prefs.getInt(WINDOW_TOP_LEFT_Y, DEFAULT_TOP_LEFT_Y);
    }
    
    
    /**
     * Checks whether the specified width is within the allowable range for
     * the SeqView Main Frame Width
     * @param testWidth the width to test
     */
    public boolean checkPreferredWindowWidth(int testWidth) {
    	return ((testWidth <= MAX_WINDOW_WIDTH) && (testWidth >= MIN_WINDOW_WIDTH));
    }
    
    
    /**
     * Checks whether the specified height is within the allowable range for
     * the SeqView Main Frame Height
     * @param testHeight the height to test
     */
    public boolean checkPreferredWindowHeight(int testHeight) {
    	return ((testHeight <= MAX_WINDOW_HEIGHT) && (testHeight >= MIN_WINDOW_HEIGHT));	
    }
    
    
    /**
     * Checks whether the specified X coordinate is within the allowable range
     * for the location of the top left corner of the SeqView Main Frame
     * @param testX the x-coord to test
     */
    public boolean checkPreferredTopLeftX(int testX) {
    	return ((testX <= MAX_TOP_LEFT_X) && (testX >= MIN_TOP_LEFT_X));
    }
    
    
    /**
     * Checks whether the specified Y coordinate is within the allowable range
     * for the location of the top left corner of the SeqView Main Frame
     * @param the y-coord to test
     */
    public boolean checkPreferredTopLeftY(int testY) {
    	return ((testY <= MAX_TOP_LEFT_Y) && (testY >= MIN_TOP_LEFT_Y));
    }
 
    
    /**
     * Returns the preferred width for the SeqView Main Frame
     * @return the preferred width for the SeqView Main Frame
     */
    public int getPreferredWindowWidth() {
    	return preferredWindowWidth;
    }
    
    
    /**
     * Sets the preferred width for the SeqView Main Frame
     * @param newWidth the preferred width for the SeqView Main Frame
     */
    public boolean setPreferredWindowWidth(int newWidth) {
    	if (this.checkPreferredWindowWidth(newWidth)) {
    		this.preferredWindowWidth = newWidth;
    		return true;
    	}
    	else {
    		return false;
    	}
    }
    
    
    /**
     * Returns the preferred height for the SeqView Main Frame
     * @return the preferred height for the SeqView Main Frame
     */
    public int getPreferredWindowHeight() {
    	return preferredWindowHeight;
    }
    
    
    /**
     * Sets the preferred height for the SeqView Main Frame
     * @param newHeight the preferred height for the SeqView Main Frame
     */
    public boolean setPreferredWindowHeight(int newHeight) {
    	if (this.checkPreferredWindowHeight(newHeight)) {
    		this.preferredWindowHeight = newHeight;
    		return true;
    	}
    	else {
    		return false;
    	}
    	
    }

    
    /**
     * Returns whether the SeqView Main Frame should be centered on the screen
     * @return true if the SeqView Main Frame should be centered on the screen
     */
    public boolean isWindowCentered() {
    	return isWindowCentered;
    }
    
    
    /**
     * Sets whether the SeqView Main Frame should be centered on the screen
     * @param isCentered true causes the SeqView Main Frame to be centered on the screen
     */
    public void setWindowCentered(boolean isCentered) {
    	this.isWindowCentered = isCentered;
    }
    
    
    /**
     * Returns the x-coord of the top left corner of the preferred location for 
     * the SeqView Main Frame
     * @return the x-coord of the top left corner of the preferred location 
     */
    public int getPreferredTopLeftX() {
    	return preferredWindowTopLeftX;
    }
    
    
    /**
     * Sets the x-coord of the top left corner of the preferred location for 
     * the SeqView Main Frame
     * @param newX the x-coord of the top left corner of the preferred location
     */
    public boolean setPreferredTopLeftX(int newX) {
    	if (this.checkPreferredTopLeftX(newX)) {
    		this.preferredWindowTopLeftX = newX;
    		return true;
    	}
    	else {
    		return false;
    	}
    }
    
   
    /**
     * Returns the y-coord of the top left corner of the preferred location for 
     * the SeqView Main Frame
     * @return the y-coord of the top left corner of the preferred location 
     */
    public int getPreferredTopLeftY() {
    	return preferredWindowTopLeftY;
    }
    
    
    /**
     * Sets the y-coord of the top left corner of the preferred location for 
     * the SeqView Main Frame
     * @param newY the y-coord of the top left corner of the preferred location 
     */
    public boolean setPreferredTopLeftY(int newY) {
    	if (this.checkPreferredTopLeftY(newY)) {
    		this.preferredWindowTopLeftY = newY;
    		return true;
    	}
    	else {
    		return false;
    	}
    }
 
    
    
    
    public static int stringToNum(String s) {
        if (s.matches(".*[kK]$")) {
            return 1000 * Integer.parseInt(s.substring(0,s.length()-1));
        }
        if (s.matches(".*[mM]$")) {
            return 1000000 * Integer.parseInt(s.substring(0,s.length()-1));
        } 
        return Integer.parseInt(s);
    }
}
