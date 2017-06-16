package org.seqcode.data.connections;

public class UnknownRoleException extends RuntimeException {
	public UnknownRoleException(String s) {
		super(s);
	}

	public UnknownRoleException(String s, Exception e) {
		super(s, e);
	}

}
