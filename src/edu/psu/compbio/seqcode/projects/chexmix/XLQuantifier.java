package edu.psu.compbio.seqcode.projects.chexmix;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentCondition;
import edu.psu.compbio.seqcode.deepseq.experiments.ExperimentManager;
import edu.psu.compbio.seqcode.deepseq.experiments.ExptConfig;
import edu.psu.compbio.seqcode.genome.GenomeConfig;
import edu.psu.compbio.seqcode.genome.location.Point;
import edu.psu.compbio.seqcode.genome.location.Region;
import edu.psu.compbio.seqcode.genome.location.StrandedPoint;
import edu.psu.compbio.seqcode.genome.location.StrandedRegion;
import edu.psu.compbio.seqcode.gse.datasets.motifs.WeightMatrix;
import edu.psu.compbio.seqcode.gse.gsebricks.verbs.sequence.SequenceGenerator;
import edu.psu.compbio.seqcode.gse.utils.sequence.SequenceUtils;

public class XLQuantifier {

	protected ExperimentManager manager;
	protected GenomeConfig gconfig;
	protected ExptConfig econfig;
	protected ChExMixConfig cconfig;
	protected ProteinDNAInteractionModel model;
	protected CompositeModelMixture mixtureModel;
	protected CompositeTagDistribution signalComposite;
	protected CompositeTagDistribution controlComposite;
	protected List<StrandedPoint> compositePoints;
	
	
	public XLQuantifier(GenomeConfig gcon, ExptConfig econ, ChExMixConfig ccon){
		gconfig = gcon;
		econfig = econ;
		cconfig = ccon;
		cconfig.makeChExMixOutputDirs(true);
		manager = new ExperimentManager(econfig);
	}
	
	public void execute(){
		//Load appropriate data
		model = ProteinDNAInteractionModel.loadFromFile(cconfig, new File(cconfig.getModelFilename()));
		compositePoints = cconfig.getCompositePoints();
		int winSize = cconfig.getCompositeWinSize();
		
		//Build the composite distribution(s)
		signalComposite = new CompositeTagDistribution(compositePoints, manager, winSize, true);
		//controlComposite = new CompositeTagDistribution(points, manager, winSize, false);
		controlComposite =null;
		
		//Initialize the mixture model 
		mixtureModel = new CompositeModelMixture(signalComposite, controlComposite, gconfig, econfig, cconfig, manager);
		
		//Set the loaded model
		mixtureModel.setModel(model);
		
		//ML assignment
		System.err.println("ML assignment");
		mixtureModel.assignML(true);
		//Get the per-site assignments
		List<CompositeModelSiteAssignment> assignments = mixtureModel.getSiteAssignments();
		Map<StrandedPoint, String> pointSeqs = retrieveSequences(assignments, signalComposite.getWinSize());
		 
		//Build per-XL-point weighted PWMs
		Map<CompositeModelComponent, WeightMatrix> componentMotifs = getPerComponentWeightMatrices(assignments, pointSeqs, 20);
		
		
		//Report
		for(ExperimentCondition cond : manager.getConditions()){
			String compositeFileName = cconfig.getOutputParentDir()+File.separator+cconfig.getOutBase()
					+"_composite."+cond.getName()+".txt";
			signalComposite.printProbsToFile(cond, compositeFileName);
			
			String perSiteRespFileName = cconfig.getOutputParentDir()+File.separator+cconfig.getOutBase()
					+"_site-component-ML."+cond.getName()+".txt";
			mixtureModel.printPerSiteComponentResponsibilitiesToFile(cond, perSiteRespFileName);
		}
		
	}
	
	/**
	 * Get sequences corresponding to the stranded points in the site assignments
	 * @param siteAssignments
	 * @param seqWin
	 * @return
	 */
	protected Map<StrandedPoint, String> retrieveSequences(List<CompositeModelSiteAssignment> siteAssignments, int seqWin){
		Map<StrandedPoint, String> pointSequences = new HashMap<StrandedPoint, String>();
		SequenceGenerator<Region> seqgen = gconfig.getSequenceGenerator();
		
		for(CompositeModelSiteAssignment sa : siteAssignments){
			Region reg = sa.getPoint().expand(seqWin/2);
			String seq = seqgen.execute(reg);
			if(sa.getPoint().getStrand()=='-')
				seq = SequenceUtils.reverseComplement(seq);
			pointSequences.put(sa.getPoint(), seq);
		}
		return pointSequences;
	}
	
	/**
	 * Calculate weighted motifs for each non-zero XL component
	 * @return
	 */
	protected Map<CompositeModelComponent, WeightMatrix> getPerComponentWeightMatrices(List<CompositeModelSiteAssignment> siteAssignments, Map<StrandedPoint, String> pointSeqs, int motifWidth){
		Map<CompositeModelComponent, WeightMatrix> xlComponentMotifs = new HashMap<CompositeModelComponent, WeightMatrix>();
		
		for(CompositeModelComponent xlComp : model.getXLComponents()){
			if(xlComp.isNonZero()){
				int compOffset = xlComp.getPosition();
				int motifStartPos = compOffset - (motifWidth/2);
				double[][] freq =  new double[motifWidth][WeightMatrix.MAXLETTERVAL];
				double[] sums =  new double[motifWidth];
				for(int x=0; x<motifWidth; x++){ sums[x]=0; for(int y=0; y<WeightMatrix.MAXLETTERVAL; y++){freq[x][y]=0;}}
				
				//Weighted freq matrix
				for(CompositeModelSiteAssignment sa : siteAssignments){
					StrandedPoint pt = sa.getPoint();
					String currSeq = pointSeqs.get(pt);
					
					for(ExperimentCondition cond : manager.getConditions()){
						double currWeight = sa.getCompResponsibility(cond, xlComp.getIndex());
						for(int x=0; x<motifWidth; x++){ 
							char letter = currSeq.charAt(motifStartPos+x);
							freq[x][letter]+=currWeight;
							sums[x]+=currWeight;
						}
					}
				}
				//Normalize freq matrix
				for(int x=0; x<motifWidth; x++){ for(int y=0; y<WeightMatrix.MAXLETTERVAL; y++){
					if(sums[x]>0)
						freq[x][y]/=sums[x];
				}}
				
				//Temp print
				System.out.println("XLComponent motif: "+xlComp.getPosition());
				for(int x=0; x<motifWidth; x++){
					System.out.println(x+"\t"+freq[x]['A']+"\t"+freq[x]['C']+"\t"+freq[x]['G']+"\t"+freq[x]['T']);
				}System.out.println("");
				
				//Log-odds over standard motif
				
				//Save
			}
		}		
		
		return xlComponentMotifs;
	}
	
	public static void main(String[] args){
		System.setProperty("java.awt.headless", "true");
		System.err.println("ChExMix version: "+ChExMixConfig.version);
		GenomeConfig gcon = new GenomeConfig(args);
		ExptConfig econ = new ExptConfig(gcon.getGenome(), args);
		//econ.setLoadRead2(false);//Enforce for chip-exo
		ChExMixConfig ccon = new ChExMixConfig(gcon, args);
		if(ccon.helpWanted()){
			System.err.println(ccon.getArgsList());
		}else if(ccon.getModelFilename()==null){
			System.err.println("Error: no ChExMix model provided. Use --model");
		}else{
			XLQuantifier xlQuant = new XLQuantifier(gcon, econ, ccon);
			xlQuant.execute();
		}
		
	}
}
