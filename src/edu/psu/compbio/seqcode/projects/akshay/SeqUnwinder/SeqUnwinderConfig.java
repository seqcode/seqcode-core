package edu.psu.compbio.seqcode.projects.akshay.SeqUnwinder;

import java.io.File;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;

/**
 * @author akshaykakumanu
 * @twitter ikaka89
 * @email auk262@psu.edu
 */
public class SeqUnwinderConfig {
	
	// General options
	protected GenomeConfig gcon;
	protected SequenceGenerator<Region> seqgen = null;
	
	// Feature options (for making the arff file and beyond)
	protected int minK=4;
	protected int maxK=5;
	/** Number K-mers in the model */
	protected int numK;
	protected int win=150;
	protected String outArffFileName="out.arff";
	protected String designFileName="SeqUnwinder.design";
	
	// Classifier options
	protected double m_Ridge=10;
	/** Augmented Lagrangian parameter rho */
	protected double m_ADMM_pho=1.7;
	protected double m_ADMM_numThreads=5;
	protected int m_SeqUnwinder_MaxIts=20;
	protected int m_ADMM_MaxIts=500;
	protected int sm_NumLayers=2;
	/** Relaxation parameter (to help faster convergence) */
	protected static final double ADMM_ALPHA=1.9;
	/** Absolute feasibility tolerance for the primal and dual feasibility conditions */
	protected final double ADMM_ABSTOL = 1E-2; 
	/** Relative  feasibility tolerance for the primal and dual feasibility conditions */
	protected final double ADMM_RELTOL = 1E-2;
	/** Tolerence for internal Nodes convergence */
	protected final double NODES_tol = 1E-2;
	/** The maximum allowed value for pho */
	protected double ADMM_pho_max = 1000;
	protected double ADMM_pho_fold = 1.0;
	
	// Discrim options
	/** Minimum length to consider for motif finding */
	protected int minM=6;
	/** Maximum length for motif finding */
	protected int maxM=10;
	/** The base name of the output directory */
	protected String outbase;
	/** The output directory file */
	protected File outdir;
	/** The minimum value of model scan score to consider form motif finding */
	protected double thresold_hills = 0.1;
	/** The number of hills to consider for clustering and motif finding for a given K-mer model */
	public final int numHills = 1500;
	/** No of iterations of clustering */
	protected int its_CLUS=10;
	/** Number of cluster; ideally each cluster corresponds to a motif */
	protected int numClus_CLUS=3;
	// Meme parameters
	protected String MEMEpath;
	protected String MEMEargs = " -dna -mod zoops -revcomp -nostatus ";
	protected int MEMEminw = 6;
	protected int MEMEmaxw = 11;
	protected int MEMEnmotifs = 3;
	protected int MEMEwin = 16;
	
	
	

}
