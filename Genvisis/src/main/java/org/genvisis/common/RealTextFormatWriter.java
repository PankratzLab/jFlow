package org.genvisis.common;

import java.io.PrintWriter;

public class RealTextFormatWriter {
	private PrintWriter writer;
	private boolean rtfOutput; 

	static final String NEW_LINE = "{\\line }";
	static final String NEW_PARAGRAPH = "{\\par }";

	public RealTextFormatWriter(String filename, boolean useRtfOutput) {
		this.writer = Files.getWriter(filename);
		this.rtfOutput = useRtfOutput;
		
		if (rtfOutput) {
			writer.println("{\\rtf1 ");
		}
	}
	
	public void print(String str) {
		str = str.replaceAll("<super>", rtfOutput?"}{\\\\super ":"");
		str = str.replaceAll("<i>", rtfOutput?"}{\\\\i ":"");
		str = str.replaceAll("<b>", rtfOutput?"}{\\\\b ":"");
		str = str.replaceAll("<u>", rtfOutput?"}{\\\\ul ":"");
		str = str.replaceAll("</>", rtfOutput?"}{":"");
		
		writer.print(str);
	}

	public void println(String str) {
		if (rtfOutput) {
			print("{");
		}
		print(str);
		if (rtfOutput) {
			print("}");
		}
		print("\n");
		if (rtfOutput) {
			writer.println(NEW_LINE);
		}
	}
	
	public void close() {
		if (rtfOutput) {
			writer.println("}");
		}
		writer.close();
	}

	public void newParagraph() {
		if (rtfOutput) {
			writer.println("\\par");
		}
	}

}
