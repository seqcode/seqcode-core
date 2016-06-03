package org.seqcode.projects.shaun;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

import org.seqcode.genome.Genome;
import org.seqcode.genome.Species;
import org.seqcode.genome.location.Point;
import org.seqcode.genome.location.Region;
import org.seqcode.gse.datasets.seqdata.SeqLocator;
import org.seqcode.gse.gsebricks.verbs.location.PointParser;
import org.seqcode.gse.gsebricks.verbs.location.RegionParser;
import org.seqcode.gse.projects.gps.DeepSeqExpt;
import org.seqcode.gse.projects.gps.ReadHit;
import org.seqcode.gse.tools.utils.Args;
import org.seqcode.gse.utils.ArgParser;
import org.seqcode.gse.utils.NotFoundException;
import org.seqcode.gse.utils.Pair;

import cern.jet.random.Poisson;
import cern.jet.random.engine.DRand;

public class JointBindingEventSimulator {
	private Species org=null;
	private Genome gen =null;
	private ArrayList<Point> peaks=null;
	private boolean haveControl=false;
	private HashMap<Point,List<ReadHit>> peakIPReads = new HashMap<Point,List<ReadHit>>();
	private HashMap<Point,List<ReadHit>> peakCtrlReads = new HashMap<Point,List<ReadHit>>();
	private int window=1000;
	private int event_spacing=100000;
	private int readLen = 36;
	private int numEvents = 500;
	private int peakpeakDistance = 50;
	private int singleEventFactor = 20;
	protected DeepSeqExpt IP=null;
	protected DeepSeqExpt CTRL=null;
	protected boolean dbconnected=false;
	protected String outBase = "out";
	protected boolean fillControlWithRandom = true;
    private int maxChromSize = 200000000;
		
	public static void main(String[] args) {
		ArgParser ap = new ArgParser(args);
        if(!ap.hasKey("species") ||(!ap.hasKey("peaks"))) { 
            System.err.println("Usage:\n " +
                               "JointBindingEventSimulator \n" +
                               "Input a set of peaks from an experiment, outputs reads corresponding to simulated combinations of these peaks. \n" +
                               " Required: \n" +
                               "  --species <species;version> " +
                               "  --peaks <file containing coordinates of peaks> \n" +
                               "  --(rdb)expt <IP expt names> \n" +
                               "  --(rdb)ctrl <ctrl expt names> \n" +
                               " More Info: \n"+
                               "  --win <window of sequence around peaks> \n"+
                               " Options: \n" +
                               "  --out output filename\n" +
                               "  --readlen <rlen>\n" +
                               "  --ppd <peak-peak dist>\n" +
                               "  --numevents <number of joint events>\n" +
                               "  --singlefactor <number of single events as a factor of numevents>\n" +
                               "  --format <ELAND/NOVO/BOWTIE>\n" +                               
                               "");
            return;
        }
        try {
        	Pair<Species, Genome> pair = Args.parseGenome(args);
        	Species currorg = pair.car();
        	Genome currgen = pair.cdr();
			
	        String peaksFile = ap.hasKey("peaks") ? ap.getKeyValue("peaks") : null;
	    	int win = ap.hasKey("win") ? new Integer(ap.getKeyValue("win")).intValue():1000;
	    	int ppd = ap.hasKey("ppd") ? new Integer(ap.getKeyValue("ppd")).intValue():50;
	    	int numEvents = ap.hasKey("numevents") ? new Integer(ap.getKeyValue("numevents")).intValue():50;
	    	int singleFac = ap.hasKey("singlefactor") ? new Integer(ap.getKeyValue("singlefactor")).intValue():20;
	        int rL = ap.hasKey("readlen") ? new Integer(ap.getKeyValue("readlen")).intValue():32;
	        String outFile = ap.hasKey("out") ? ap.getKeyValue("out") : "out";
	        
        	/////////////////////////////////////////////////////////////////////////
	        ///////////// START 
        
			//initialize
	        JointBindingEventSimulator simulator = new JointBindingEventSimulator(currorg, currgen);
	        simulator.setWin(win);
	        simulator.setPPD(ppd);
	        simulator.setNumEvents(numEvents);
	        simulator.setSingleFac(singleFac);
	        simulator.setReadLen(rL);
	        simulator.setOut(outFile);
			
	        //Load data
			simulator.loadPeaks(peaksFile);
			simulator.loadExperiments(args);
			simulator.loadPeakReads();
			
			//Simulate reads
			//All simulated joint peaks will be placed in chromosome "A". 
			//All simulated single peaks (controls) will be placed in chromosome "B" onwards
			simulator.simulate();
			
			simulator.close();
        } catch (NotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public JointBindingEventSimulator(Species o, Genome g){
		org = o;
		gen=g;
	}
	///////////////////////////////////////////////////////////////////////
	//Options first
	///////////////////////////////////////////////////////////////////////

	
	
	///////////////////////////////////////////////////////////////////////
	
	public void setWin(int w){window=w;}
	public void setPPD(int p){peakpeakDistance = p;}
	public void setNumEvents(int n){numEvents=n;}
	public void setSingleFac(int sef){singleEventFactor=sef;}
	public void setReadLen(int rl){readLen=rl;}
	public void setOut(String o){outBase = o;}
	
	//load peaks
	public ArrayList<Point> loadPeaks(String fname){
		peaks = new ArrayList<Point>();
		try{
			File pFile = new File(fname);
			if(!pFile.isFile()){System.err.println("Invalid positive file name");System.exit(1);}
	        BufferedReader reader = new BufferedReader(new FileReader(pFile));
	        String line; //= reader.readLine();//Ignore first line
	        while ((line = reader.readLine()) != null) {
	            line = line.trim();
	            String[] words = line.split("\\s+");
	            
	            if(words.length>=3  && words[2].indexOf(":")>0){
	                PointParser pparser = new PointParser(gen);
	            	Point p = pparser.execute(words[2]);
	            	peaks.add(p);
		    }else if(words.length>=1 && words[0].indexOf(":")>0 && words[0].indexOf("-")>0){
	            	RegionParser rparser = new RegionParser(gen);
	            	Region r = rparser.execute(words[0]);
	            	if(r!=null){peaks.add(r.getMidpoint());}
		    }else if(words.length>=1 && words[0].indexOf(":")>0){
			    PointParser pparser = new PointParser(gen);
			    Point p = pparser.execute(words[0]);
			    if(p!=null){peaks.add(p);}
		    }
	        }reader.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return(peaks);
	}
	
	
	//Load peak reads
	public void loadPeakReads(){
		
		for(Point p : peaks){
			int rstart = p.getLocation()-(window/2)<1 ? 1:p.getLocation()-(window/2);
			int rend = p.getLocation()+(window/2)>gen.getChromLength(p.getChrom()) ? gen.getChromLength(p.getChrom()):p.getLocation()+(window/2)-1;
			Region r = new Region(p.getGenome(), p.getChrom(), rstart, rend);
			
			List<ReadHit> ihits = IP.loadHits(r);
			peakIPReads.put(p, ihits);
			if(haveControl){
				List<ReadHit> chits = CTRL.loadHits(r);
				peakCtrlReads.put(p, chits);
			}
		}
	}
	
	protected void loadExperiments(String [] args){
		List<SeqLocator> dbexpts = Args.parseSeqExpt(args,"dbexpt");
        List<SeqLocator> dbctrls = Args.parseSeqExpt(args,"dbctrl");
        List<SeqLocator> rdbexpts = Args.parseSeqExpt(args,"rdbexpt");
        List<SeqLocator> rdbctrls = Args.parseSeqExpt(args,"rdbctrl");
        List<File> expts = Args.parseFileHandles(args, "expt");
        List<File> ctrls = Args.parseFileHandles(args, "ctrl");
        String fileFormat = Args.parseString(args, "format", "ELAND");
        if(expts.size()>0 && dbexpts.size() == 0 && rdbexpts.size()==0){
        	IP = new DeepSeqExpt(gen, expts, false, fileFormat, readLen);
        }else if(dbexpts.size()>0 && expts.size() == 0){
        	IP = new DeepSeqExpt(gen, dbexpts, "db", readLen);
        	dbconnected=true;
        }else if(rdbexpts.size()>0 && expts.size() == 0){
        	IP = new DeepSeqExpt(gen, rdbexpts, "readdb", readLen);
        	dbconnected=true;
        }else{}
        if(ctrls.size()>0 && dbctrls.size() == 0){
        	CTRL = new DeepSeqExpt(gen, ctrls, false, fileFormat, readLen);
        }else if(dbctrls.size()>0 && ctrls.size() == 0){
        	CTRL = new DeepSeqExpt(gen, dbctrls, "db", readLen); 
        	dbconnected=true;
        }else if(rdbctrls.size()>0 && ctrls.size() == 0){
        	CTRL = new DeepSeqExpt(gen, rdbctrls, "readdb", readLen); 
        	dbconnected=true;
        }else{
        	if(dbctrls.size()>0 && ctrls.size()>0){
        		System.err.println("Cannot mix files and db loading yet...");System.exit(1);
        	}else{
        		CTRL=null;
        	}
        }
        if(CTRL==null){
        	haveControl=false;
        	fillControlWithRandom=true;
        }else{
        	haveControl=true;
        }
	}
	
	//Simulate:
	//Generate numEvents joint binding events that are peakpeakDistance apart and place them every 2*window bp in fictonal chromosome A
	//Generate 10 x numEvents single binding events and place them every 2*window bp in fictional chromosome B
	public void simulate(){
		int numPeaks = peaks.size();
		Random rng = new Random();
		int currCoord = window+1;
		int totalIPReads=0, totalCtrlReads=0;
		//Set up output files: IP reads, ctrl reads, peak positions
		try{
			FileWriter ipFile = new FileWriter(outBase+"_ip.bed");
			FileWriter ctrlFile = new FileWriter(outBase+"_ctrl.bed");
			FileWriter peakFile = new FileWriter(outBase+".peaks");
			
			char chromosome = 'A';
			//Joint events first
			for(int i=1; i<=numEvents; i++){
				//Pick two random peaks
				Point p1 = peaks.get(rng.nextInt(numPeaks));
				Point p2 = peaks.get(rng.nextInt(numPeaks));
				
				//Place peak1's reads around currCoord
				//IP
				for(ReadHit r : peakIPReads.get(p1)){
					int rstart = currCoord+(r.getStart()-p1.getLocation());
					int rend = currCoord+(r.getEnd()-p1.getLocation());
					totalIPReads++;
					ipFile.write("chr"+chromosome+"\t"+rstart+"\t"+rend+"\tU\t0\t"+r.getStrand()+"\n");				
				}
				if(haveControl){
					//CTRL
					for(ReadHit r : peakCtrlReads.get(p1)){
						int rstart = currCoord+(r.getStart()-p1.getLocation());
						int rend = currCoord+(r.getEnd()-p1.getLocation());
						totalCtrlReads++;
						ctrlFile.write("chr"+chromosome+"\t"+rstart+"\t"+rend+"\tU\t0\t"+r.getStrand()+"\n");				
					}
				}
				//Print peak position to peak file
				peakFile.write("chr"+chromosome+":"+currCoord+"\n");
				
				//Place peak2's reads around currCoord+peakpeakDistance
				//IP
				for(ReadHit r : peakIPReads.get(p2)){
					int rstart = currCoord+peakpeakDistance+(r.getStart()-p2.getLocation());
					int rend = currCoord+peakpeakDistance+(r.getEnd()-p2.getLocation());
					totalIPReads++;
					ipFile.write("chr"+chromosome+"\t"+rstart+"\t"+rend+"\tU\t0\t"+r.getStrand()+"\n");				
				}
				if(haveControl){
					//CTRL
					for(ReadHit r : peakCtrlReads.get(p2)){
						int rstart = currCoord+peakpeakDistance+(r.getStart()-p2.getLocation());
						int rend = currCoord+peakpeakDistance+(r.getEnd()-p2.getLocation());
						totalCtrlReads++;
						ctrlFile.write("chr"+chromosome+"\t"+rstart+"\t"+rend+"\tU\t0\t"+r.getStrand()+"\n");				
					}
				}
				//Print peak position to peak file
				peakFile.write("chr"+chromosome+":"+(currCoord+peakpeakDistance)+"\n");
				
				//Jump to next coordinate
				currCoord += event_spacing;
			}
			//More control reads...
			if(fillControlWithRandom){
				while(totalCtrlReads<totalIPReads){
					int rstart = rng.nextInt(currCoord-readLen-1);
					char strand = rng.nextInt(100)%2 ==0 ? '+' : '-';
					int rend = rstart+readLen;
					ctrlFile.write("chr"+chromosome+"\t"+rstart+"\t"+rend+"\tU\t0\t"+strand+"\n");
					totalCtrlReads++;
				}
			}
			
			System.out.println("chr"+chromosome+"\t"+currCoord);
			currCoord=window+1;
			
			
			//Now single events
			chromosome++;
			for(int i=1; i<=singleEventFactor*numEvents; i++){
				//Pick two random peaks
				Point p = peaks.get(rng.nextInt(numPeaks));
				
				//Place peak1's reads around currCoord
				//IP
				for(ReadHit r : peakIPReads.get(p)){
					int rstart = currCoord+(r.getStart()-p.getLocation());
					int rend = currCoord+(r.getEnd()-p.getLocation());
					totalIPReads++;
					//This should go to IP file
					ipFile.write("chr"+chromosome+"\t"+rstart+"\t"+rend+"\tU\t0\t"+r.getStrand()+"\n");				
				}
				if(haveControl){
					//CTRL
					for(ReadHit r : peakCtrlReads.get(p)){
						int rstart = currCoord+(r.getStart()-p.getLocation());
						int rend = currCoord+(r.getEnd()-p.getLocation());
						totalCtrlReads++;
						//This should go to IP file
						ctrlFile.write("chr"+chromosome+"\t"+rstart+"\t"+rend+"\tU\t0\t"+r.getStrand()+"\n");				
					}
				}
				//Print peak position to peak file
				peakFile.write("chr"+chromosome+":"+currCoord+"\n");
				
				//Jump to next coordinate
				currCoord += event_spacing;
				
				if(currCoord>maxChromSize){
					//Even more control reads...
					if(fillControlWithRandom){
						while(totalCtrlReads<totalIPReads){
							int rstart = rng.nextInt(currCoord-readLen-1);
							char strand = rng.nextInt(100)%2 ==0 ? '+' : '-';
							int rend = rstart+readLen;
							ctrlFile.write("chr"+chromosome+"\t"+rstart+"\t"+rend+"\tU\t0\t"+strand+"\n");
							totalCtrlReads++;
						}
					}			
					System.out.println("chr"+chromosome+"\t"+currCoord);
					
					chromosome++;
					currCoord=window+1;
				}
			}
			//Fill the last chromosome
			if(fillControlWithRandom){
				while(totalCtrlReads<totalIPReads){
					int rstart = rng.nextInt(currCoord-readLen-1);
					char strand = rng.nextInt(100)%2 ==0 ? '+' : '-';
					int rend = rstart+readLen;
					ctrlFile.write("chr"+chromosome+"\t"+rstart+"\t"+rend+"\tU\t0\t"+strand+"\n");
					totalCtrlReads++;
				}
			}			
			System.out.println("chr"+chromosome+"\t"+currCoord);
			
			
			ipFile.close();
			ctrlFile.close();
			peakFile.close();
			System.out.println("Files written with basename: "+outBase);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	//Poisson threshold for needle filter
	public int getPoissonThreshold(double threshold, double count, double hitLength, double seqLen, double mappable, double binWidth, double binOffset){
		int countThres=0;
		DRand re = new DRand();
		Poisson P = new Poisson(0, re);
		double pMean = (count*(hitLength/binOffset + binWidth/binOffset))/(seqLen*mappable/binOffset); 
		P.setMean(pMean);
		double l=1;
		for(int b=1; l>threshold; b++){
			l=1-P.cdf(b);
			countThres=b;
		}
		return(Math.max(1,countThres));
	}
	
	//Clean up the loaders
	public void close(){
		if(IP!=null)
			IP.closeLoaders();
		if(CTRL!=null)
			CTRL.closeLoaders();
	}
}

