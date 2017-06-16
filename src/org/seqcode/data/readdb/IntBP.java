package org.seqcode.data.readdb;

import java.nio.*;

/**
 * contains an IntBuffer ib that is derived from the ByteBuffer bb
 */
public class IntBP extends ByteBP {
	private IntBuffer ib;

	public IntBP(ByteBuffer b) {
		super(b);
		setib(b.asIntBuffer());
	}

	public IntBP(int size) {
		super(ByteBuffer.allocate(size * 4));
		setib(bb.asIntBuffer());
	}

	public IntBP slice(int start, int length) {
		ByteBuffer b;
		synchronized (bb) {
			bb.position(start * 4);
			b = bb.slice();
		}
		b.limit(length * 4);
		return new IntBP(b);
	}

	public int get(int i) {
		return getib().get(i);
	}

	public void put(int index, int val) {
		getib().put(index, val);
	}

	public int limit() {
		return getib().limit();
	}

	public int size() {
		return getib().limit();
	}

	public IntBuffer getib() {
		return ib;
	}

	public void setib(IntBuffer ib) {
		this.ib = ib;
	}
}
