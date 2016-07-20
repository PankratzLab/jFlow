package org.genvisis.filesys;

import java.io.PrintWriter;
import java.io.Serializable;

import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.filesys.PlainTextExport;

public class SegmentList implements Serializable, PlainTextExport {
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
	
	@Override
	public void exportToText(String outputFile, Logger log) {
        PrintWriter writer;
        
        writer = Files.getAppropriateWriter(outputFile);
        writer.println("Chr\tStart\tStop");
        for (Segment seg : list) {
            writer.println(seg.toAnalysisString());
        }
        writer.flush();
        writer.close();
	}
	
}