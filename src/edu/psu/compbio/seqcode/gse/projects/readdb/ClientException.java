package edu.psu.compbio.seqcode.gse.projects.readdb;

public class ClientException extends Exception {
	private static final long serialVersionUID = 1L;

	public ClientException (String reason) {
        super(reason);
    }
}