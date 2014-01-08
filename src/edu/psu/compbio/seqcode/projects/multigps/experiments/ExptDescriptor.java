package edu.psu.compbio.seqcode.projects.multigps.experiments;

import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.multigps.framework.BindingModel;

/**
 * ExptDescriptor: simple class for describing an experiment to be loaded
 * @author Shaun Mahony
 * @version	%I%, %G%
 */
public class ExptDescriptor {
	
	public String condition = "";
	public String replicate = "";
	public boolean signal; //true=signal, false=control
	public List<Pair<String,String>> sources = new ArrayList<Pair<String,String>>(); //name & type/format
	public BindingModel bindingModel=null;
	public float perBaseMaxReads = -1;
	
	public ExptDescriptor(String cond, String rep, boolean sig, List<Pair<String, String>> src) {this(cond, rep, sig, src, null, -1);}
	public ExptDescriptor(String cond, String rep, boolean sig, List<Pair<String, String>> src, BindingModel mod, float perBPMax) {
		condition = cond;
		replicate = rep;
		signal = sig;
		sources.addAll(src);
		bindingModel=mod;
		perBaseMaxReads = perBPMax;
	}
	
	public ExptDescriptor(String cond, String rep, boolean sig, Pair<String, String> src) {this(cond, rep, sig, src, null, -1);}
	public ExptDescriptor(String cond, String rep, boolean sig, Pair<String, String> src, BindingModel mod, float perBPMax) {
		condition = cond;
		replicate = rep;
		signal = sig;
		sources.add(src);
		bindingModel=mod;
		perBaseMaxReads = perBPMax;
	}
	
	// added by Akshay, to define ExptDescriptor class in the absence of binding model and still having a non default perBPMax
	public ExptDescriptor(String cond, String rep, boolean sig, Pair<String, String> src, float perBPMax){this(cond, rep, sig, src, null, perBPMax);} 
	
}
