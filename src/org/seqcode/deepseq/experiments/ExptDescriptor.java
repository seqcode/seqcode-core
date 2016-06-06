package org.seqcode.deepseq.experiments;

import java.util.ArrayList;
import java.util.List;

import org.seqcode.utils.Pair;


/**
 * ExptDescriptor: simple class for describing an experiment to be loaded
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class ExptDescriptor {
	
	public String target="";
	public String expttype = "";
	public String condition = "";
	public String replicate = "";
	public boolean signal; //true=signal, false=control
	public List<Pair<String,String>> sources = new ArrayList<Pair<String,String>>(); //name & type/format
	public float perBaseMaxReads = -1;
	
	/**
	 * ExptDescriptor constructor
	 * @param type : (optional) Experiment type (CHIPSEQ/CHIPEXO/INPUT/DNASESEQ/RNASEQ/etc) [String]
	 * @param targ : (optional) Target of the assay; e.g. protein name for ChIP experiments [String]
	 * @param cond : Experiment condition name [String]
	 * @param rep : Experiment replicate name [String]
	 * @param sig : Sample contains signal as opposed to control [boolean]
	 * @param src : List of data source and type/format pairs [List<Pair<String, String>>]
	 * @param perBPMax : Per base hit count maximum (-1 = no limit; >0 = fixed limit; P = poisson determined limit)
	 */
	public ExptDescriptor(String type, String targ, String cond, String rep, boolean sig, List<Pair<String, String>> src, float perBPMax) {
		expttype = type.equals("") ? "SEQEXPT" : type;
		target = targ.equals("") ? "NA" : targ;
		condition = cond;
		replicate = rep;
		signal = sig;
		sources.addAll(src);
		perBaseMaxReads = perBPMax;
	}
	
	/**
	 * ExptDescriptor constructor
	 * @param type : (optional) Experiment type (CHIPSEQ/CHIPEXO/INPUT/DNASESEQ/RNASEQ/etc) [String]
	 * @param targ : (optional) Target of the assay; e.g. protein name for ChIP experiments [String]
	 * @param cond : Experiment condition name [String]
	 * @param rep : Experiment replicate name [String]
	 * @param sig : Sample contains signal as opposed to control [boolean]
	 * @param src : Data source and type/format pair [Pair<String, String>]
	 * @param perBPMax : Per base hit count maximum (-1 = no limit; >0 = fixed limit; P = poisson determined limit)
	 */
	public ExptDescriptor(String type, String targ, String cond, String rep, boolean sig, Pair<String, String> src, float perBPMax) {
		expttype = type.equals("") ? "SEQEXPT" : type;
		target = targ.equals("") ? "NA" : targ;
		condition = cond;
		replicate = rep;
		signal = sig;
		sources.add(src);
		perBaseMaxReads = perBPMax;
	}
	
}
