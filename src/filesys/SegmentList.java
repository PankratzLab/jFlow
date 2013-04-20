package filesys;

import java.io.Serializable;

import common.Files;

public class SegmentList implements Serializable {
	public static final long serialVersionUID = 1L;

	private Segment[] list;

	public SegmentList(Segment[] list) {
		this.list = list;
	}

	public Segment[] getList() {
		return list;
	}

	public void serialize(String filename) {
		Files.writeSerial(this, filename);
	}

	public static SegmentList load(String filename, boolean jar) {
		return (SegmentList)Files.readSerial(filename, jar, true);
	}
}
