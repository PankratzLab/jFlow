package org.genvisis.widgets;

import java.io.IOException;

import org.genvisis.common.ext;

public class QuoteMe {
	public static boolean NEWLINE = false;

	public static void main(String[] args) throws IOException {
		try {
			String str = ext.getClipboard();
			String out = "";
			String[] line, subline;
			boolean matrix;

			str = ext.replaceAllWith(str, "\r\n", "\n");

			while (str.endsWith("\n")) {
				str = str.substring(0, str.length() - 1);
			}

			line = str.split("\n", -1);
			matrix = line.length > 1;

			for (String element : line) {
				subline = element.split("\t", -1);
				out += matrix ? "{" : "";
				for (int j = 0; j < subline.length; j++) {
					out += (j == 0 ? "" : ", " + (NEWLINE ? "\n" : ""))	+ "\""
									+ ext.replaceQuotesWithSlashQuotes(subline[j]) + "\"";
				}
				out += matrix ? "},\n" : "";
			}

			ext.setClipboard(out);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
