package edu.psu.compbio.seqcode.projects.akshay.Utils;

import edu.psu.compbio.seqcode.gse.tools.utils.Args;
import edu.psu.compbio.seqcode.gse.utils.ArgParser;

public class MemeER {
	
	protected String MEMEpath;
	protected String MEMEargs;
	protected Float pseudo = (float)0.001;
	public MemeER(String path, String args) {
		this.MEMEpath = path;
		this.MEMEargs = args;
		
	}
	public static void main(String[] args){
		ArgParser ap = new ArgParser(args);
		String memeargs = Args.parseString(args, "memeargs", " -dna -mod zoops -revcomp -nostatus ");
		int MEMEminw = Args.parseInteger(args, "mememinw", 6);
		//MEME maxw
		int MEMEmaxw = Args.parseInteger(args, "mememaxw", 18);
		//MEME nmotifs option
		int MEMEnmotifs = Args.parseInteger(args,"memenmotifs", 3);
		memeargs = memeargs + " -nmotifs "+MEMEnmotifs + " -minw "+MEMEminw+" -maxw "+MEMEmaxw;
		MemeER meme = new MemeER(Args.parseString(args, "memepath", ""), memeargs);
		
		
		
	}
	
}
