package org.seqcode.projects.shaun;

/**
 * ConsensusSequence: simple class to represent consensus DNA sequences
 * @author mahony
 *
 */
public class ConsensusSequence {
	private static int MAXLETTERVAL = Math.max(Math.max(Math.max('A','C'),Math.max('T','G')),
            Math.max(Math.max('a','c'),Math.max('t','g')))
            + 1;
	private static char[] allLetters = {'A','a','C','c','T','t','G','g'};
	private static char[] letters = {'A','C','T','G'};
	
	private boolean[][] consensusMatrix; //Represents whether the consensus allows each character at each position
	private String consensus;
	private int length;
	
	public ConsensusSequence(String cons){this(cons.toCharArray());}
	public ConsensusSequence(char[] cons){
		consensus = new String(cons);
		consensusMatrix = new boolean[cons.length][MAXLETTERVAL];
		length = cons.length;
		for(int i=0; i<cons.length; i++)
			for(int j=0; j<MAXLETTERVAL; j++)
				consensusMatrix[i][j]=false;
		
		for(int i=0; i<cons.length; i++){
			if(Character.toUpperCase(cons[i])=='A' || Character.toUpperCase(cons[i])=='C' || Character.toUpperCase(cons[i])=='G' || Character.toUpperCase(cons[i])=='T')
				consensusMatrix[i][Character.toUpperCase(cons[i])]=true;
			else if(Character.toUpperCase(cons[i])=='N'){
				consensusMatrix[i]['A']=true; consensusMatrix[i]['C']=true;
				consensusMatrix[i]['C']=true; consensusMatrix[i]['T']=true;
			}else if(Character.toUpperCase(cons[i])=='R'){
				consensusMatrix[i]['A']=true; consensusMatrix[i]['G']=true;
			}else if(Character.toUpperCase(cons[i])=='Y'){
				consensusMatrix[i]['C']=true; consensusMatrix[i]['T']=true;
			}else if(Character.toUpperCase(cons[i])=='M'){
				consensusMatrix[i]['A']=true; consensusMatrix[i]['C']=true;
			}else if(Character.toUpperCase(cons[i])=='K'){
				consensusMatrix[i]['G']=true; consensusMatrix[i]['T']=true;
			}else if(Character.toUpperCase(cons[i])=='S'){
				consensusMatrix[i]['C']=true; consensusMatrix[i]['G']=true;
			}else if(Character.toUpperCase(cons[i])=='W'){
				consensusMatrix[i]['A']=true; consensusMatrix[i]['T']=true;
			}else if(Character.toUpperCase(cons[i])=='H'){
				consensusMatrix[i]['A']=true; consensusMatrix[i]['C']=true; consensusMatrix[i]['T']=true;
			}else if(Character.toUpperCase(cons[i])=='B'){
				consensusMatrix[i]['C']=true; consensusMatrix[i]['G']=true; consensusMatrix[i]['T']=true;
			}else if(Character.toUpperCase(cons[i])=='V'){
				consensusMatrix[i]['A']=true; consensusMatrix[i]['C']=true; consensusMatrix[i]['G']=true;
			}else if(Character.toUpperCase(cons[i])=='D'){
				consensusMatrix[i]['A']=true; consensusMatrix[i]['G']=true; consensusMatrix[i]['T']=true;
			}
			consensusMatrix[i]['a'] = consensusMatrix[i]['A'];
			consensusMatrix[i]['c'] = consensusMatrix[i]['C'];
			consensusMatrix[i]['g'] = consensusMatrix[i]['G'];
			consensusMatrix[i]['t'] = consensusMatrix[i]['T'];
		}
	}
	
	//Accessors
	public String getSequence(){return consensus;}
	public int getMaxMismatch(){return length;}
	public int getLength(){return length;}
	public boolean[][] getMatrix(){ return consensusMatrix;}
	
}
