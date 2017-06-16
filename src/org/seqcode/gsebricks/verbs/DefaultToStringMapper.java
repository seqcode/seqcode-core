package org.seqcode.gsebricks.verbs;

public class DefaultToStringMapper implements Mapper<Object, String> {

	public String execute(Object o) {
		return o.toString();
	}

}