package org.seqcode.viz.paintable;

public class PaintableChangedEvent {

	private Paintable fSource;

	public PaintableChangedEvent(Paintable p) {
		fSource = p;
	}

	public Paintable getSource() {
		return fSource;
	}
}
