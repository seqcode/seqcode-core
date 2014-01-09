package edu.psu.compbio.seqcode.projects.akshay.bayesments.experiments;

import java.util.ArrayList;
import java.util.List;

import edu.psu.compbio.seqcode.gse.utils.Pair;
import edu.psu.compbio.seqcode.projects.multigps.framework.BindingModel;

public class ExptDescriptor {
	
	public String condition;
	public String replicate;
	public String feature;
	public boolean signal;
	public List<Pair<String,String>> sources; // src and type pairs (eg ("/akshay/atf3.bed", "BED"))
	public float perBaseReadLimit;
	public BindingModel bindingmodel = null; // always null: had to include for compatability with multigps code
	
	
	
	public ExptDescriptor(String condition, String replicate, String feature, boolean signal, Pair<String, String> src) {
		this(condition, replicate, feature, signal, src, -1);
	}
	public ExptDescriptor(String condition, String replicate, String feature, boolean signal, Pair<String, String> src, float perBaseReadLimit ) {
		this.condition = condition;
		this.replicate = replicate;
		this.feature = feature;
		this.signal = signal;
		this.sources = new ArrayList<Pair<String, String>>();
		this.sources.add(src);
	}
	
	
	

}
