package edu.psu.compbio.seqcode.gse.viz.metaprofile;

import java.awt.Color;
import java.util.List;

import edu.psu.compbio.seqcode.gse.tools.utils.Args;

public class MetaConfig {

	public boolean printHelp = false;
	public boolean batchRun= false;
	public boolean cluster = false;
	public int winLen = 10000;
	public int bins = 100;
	public int readExt = 0;
	public double lineMin =0;
	public double lineMax = 100;
	public int lineThick = 1;
	public double pbMax = 100;
	public char strand = '.';
	public boolean drawColorBar = true;
	public boolean saveSVG = false;
	public boolean transparent = false;
	public boolean printMatrix = false;
	public String profilerType = "fiveprime";	
	public String profileStyle = "Line";	
	public List<String> peakFiles = null;
	public String outName = "meta";
	public Color color = Color.blue;
	public char baseLimit='.'; //Only draw meta-plots for tags with this character at baseLimitRelPosition relative to 5' end (. = draw all tags)
	public int baseLimitRelPosition=0;
	public int fivePrimeShift=0;
	
	public MetaConfig(String [] args){
		if(args.length < 2){
			printHelp = true;
		}else{
			winLen = Args.parseInteger(args,"win", 10000);
			bins = Args.parseInteger(args,"bins", 100);
			readExt = Args.parseInteger(args,"readext", 0);
			fivePrimeShift = Args.parseInteger(args,"readshift", 0);
			lineMin = Args.parseDouble(args,"linemin", 0);
			lineMax = Args.parseDouble(args,"linemax", 100);
			lineThick = Args.parseInteger(args,"linethick", 1);
			pbMax = Args.parseDouble(args,"pbmax", 100);
			strand = Args.parseString(args, "strand", ".").charAt(0);
			drawColorBar = !Args.parseFlags(args).contains("nocolorbar");
			saveSVG = Args.parseFlags(args).contains("svg");
			transparent = Args.parseFlags(args).contains("transparent");
			printMatrix = Args.parseFlags(args).contains("printMatrix");
			profilerType = Args.parseString(args, "profiler", "simplechipseq");	
			profileStyle = Args.parseString(args, "style", "Line");	
			peakFiles = (List<String>)Args.parseStrings(args, "peaks");
			outName = Args.parseString(args, "out", "meta");
			batchRun = Args.parseFlags(args).contains("batch");
			cluster = Args.parseFlags(args).contains("cluster");
			color = Color.blue;
			String newCol = Args.parseString(args, "color", "blue");
			if(newCol.equals("red"))
				color=Color.red;
			if(newCol.equals("green"))
				color=Color.green;
			for(int s=0; s<args.length; s++){
				if(args[s].equals("--color4")){
					Integer R = new Integer(args[s+1]);
					Integer G = new Integer(args[s+2]);
					Integer B = new Integer(args[s+3]);
					Integer A = new Integer(args[s+4]);
					color= new Color(R,G,B,A);
				}
			}
			baseLimit = Args.parseString(args, "baselimit", ".").charAt(0);
			baseLimitRelPosition = Args.parseInteger(args, "baselimitposition", 0);
		}
	}
	
	public boolean helpWanted(){return printHelp;}
	
	/**
	 * Returns a string describing the arguments handled by this config parser. 
	 * @return String
	 */
	public String getArgsList(){
		return(new String(""+
				"MetaConfig:\n"+
				"\t--profiler <simplechipseq/fiveprime/nucleosome> \n" +
				"\t--win <profile width> --bins <num bins> \n" +
				"\t--linemin <min>  --linemax <max> \n" +
				"\t--readext <extension>\n" +
				"\t--readshift <5' shift>\n" +
				"\t--peaks <peaks file name>\n" +
				"\t--out <output root name> \n" +
				"\t--color <red/green/blue> or --color4 <R G B A>\n" +
				"\t--strand <+-.>\n" +
				"\t--baselimit <./A/C/G/T: only draw tags with this base at below position>\n" +
				"\t--baselimitposition <only draw tags with above base at this position>\n" +
				"\t--printMatrix [flag to print the matrix of tags] \n"+
				"\t--cluster [flag to cluster in batch mode] \n" +
				"\t--batch [a flag to run without displaying the window]\n" +
				"\t--nocolorbar [flag to turn off colorbar in batch mode]\n" +
				"\t--transparent [flag for transparent background]\n" +
				"\t--style <Line/Histo>\n" +
				"\t--svg [flag to save an SVG image]\n"));
	}
}
