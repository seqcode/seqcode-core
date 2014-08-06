package edu.psu.compbio.seqcode.projects.shaun.metaexperiment;

import java.io.File;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import edu.psu.compbio.seqcode.genome.Genome;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqExptHandler;
import edu.psu.compbio.seqcode.gse.datasets.seqdata.SeqLocator;
import edu.psu.compbio.seqcode.gse.projects.gps.DeepSeqExpt;
import edu.psu.compbio.seqcode.gse.utils.NotFoundException;

public class SeqExptInfo {
	private Genome g;
	private String name;
	private String type; //chip/ctrl
	private String ctrlName; //Control experiment for this experiment
	private double readLen=26;
	private double readShift=0;
	private double zThres=2.33;
	private boolean usingFile=false;
	private String fileFormat="ELAND";
	private boolean useNonUnique=false;
	private DeepSeqExpt dse;
	private double totalReads=0;
	private List<File> exptFiles =null;
	private List<SeqLocator> exptLocs =null;
	
	//Constructor 
	public SeqExptInfo(Genome gen, String eName, String eType, String eCtrl, double eReadLen, double eReadShift, List<SeqLocator> expts){
		usingFile=false;
		useNonUnique=false;
		exptLocs = expts;
		g=gen;
		name=eName;
		ctrlName = eCtrl;
		readLen = eReadLen;
		readShift = eReadShift;
		if(eType.equals("chip") || eType.equals("ctrl")){
			type = eType;
		}else{System.err.println("Invalid experiment type; must be chip or ctrl");}
		
		//dse = new DeepSeqExpt(gen, expts);
		//initialize(gen, eName, eType, eCtrl, eReadLen, eReadShift);
	}
	public SeqExptInfo(Genome gen, String eName, String eType, String eCtrl, double eReadLen, double eReadShift, List<File> files, boolean useNonUnique, String format){
		usingFile=true;
		this.useNonUnique = useNonUnique;
		fileFormat = format;
		exptFiles = files;
		g=gen;
		name=eName;
		ctrlName = eCtrl;
		readLen = eReadLen;
		readShift = eReadShift;
		if(eType.equals("chip") || eType.equals("ctrl")){
			type = eType;
		}else{System.err.println("Invalid experiment type; must be chip or ctrl");}
		
		//dse = new DeepSeqExpt(gen, files, useNonUnique, format);
		//initialize(gen, eName, eType, eCtrl, eReadLen, eReadShift);
	}
	public void initialize(){
		if(usingFile)
			dse = new DeepSeqExpt(g, exptFiles, useNonUnique, fileFormat, (int)readLen);
		else
			dse = new DeepSeqExpt(g, exptLocs, "db", (int)readLen);
		dse.setShift((int)readShift);
		totalReads=dse.getHitCount();
		System.err.print(String.format(name+" %.0f reads loaded\n", totalReads));
	}
	
	//Accessors
	public String getName(){return name;}
	public String getType(){return type;}
	public String getCtrl(){return ctrlName;}
	public double getReadLen(){return readLen;}
	public double getReadShift(){return readShift;}
	public double getZThres(){return zThres;}
	public DeepSeqExpt getDeepSeqHandle(){return dse;}
	public double getTotalReads(){return totalReads;}
	public void setName(String n){name = n;}
	public void setType(String t){type = t;}
	public void setCtrl(String c){ctrlName = c;}
	public void setReadLen(int l){readLen = l;}
	public void setReadShift(int s){readShift = s;}
	public void setZThres(double z){zThres = z;}
	
}
