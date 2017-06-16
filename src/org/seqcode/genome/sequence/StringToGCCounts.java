/*
 * Created on Apr 18, 2007
 */
package org.seqcode.genome.sequence;

import java.util.Map;

import org.seqcode.gsebricks.types.*;
import org.seqcode.gsebricks.verbs.Mapper;
import org.seqcode.gseutils.Pair;

public class StringToGCCounts implements Mapper<String, Pair<Integer, Integer>>, SelfDescribingVerb {

	public Pair<Integer, Integer> execute(String a) {
		int length = a.length();
		int count = 0;
		for (int i = 0; i < a.length(); i++) {
			char c = Character.toUpperCase(a.charAt(i));
			if (c == 'G' || c == 'C') {
				count += 1;
			}
		}
		return new Pair<Integer, Integer>(count, length);
	}

	private static final String[] inputNames = { "Strings" };
	private static final EchoType[] inputTypes = { new ClassType(String.class) };
	private static final EchoType outputType = new PairType(new ClassType(Integer.class), new ClassType(Integer.class));

	public EchoType[] getInputClasses() {
		return inputTypes;
	}

	public String[] getInputNames() {
		return inputNames;
	}

	public EchoType getOutputClass() {
		return outputType;
	}

	public EchoType[] getParameterClasses() {
		return null;
	}

	public String[] getParameterNames() {
		return null;
	}

	public void init(Map<String, Object> params) {
	}

}
