// Each choice should be available through Cuts/
package org.genvisis.widgets;

import java.io.IOException;

import org.genvisis.common.ext;

public class OutMe {
	public static final int SYSTEM_OUT = 0;
	public static final int WRITER = 1;
	public static final int QUOTES = 2;

	public static void main(String[] args) throws IOException {
		try {
			String str = ext.getClipboard();
			String out = "";
			String[] line;
			int choice = WRITER;

			if (args.length > 0) {
				choice = Integer.parseInt(args[0]);
			}

			str = ext.replaceAllWith(str, "\r\n", "\n");
			str = ext.replaceAllWith(str, "\"", "~!@#$%^");
			str = ext.replaceAllWith(str, "~!@#$%^", "\\\"");

			while (str.endsWith("\n")) {
				str = str.substring(0, str.length() - 1);
			}

			line = str.split("\n", -1);

			for (String element : line) {
				switch (choice) {
					case SYSTEM_OUT:
						out += "System.out.println(\"" + element + "\");\n";
						break;
					case WRITER:
						out += "writer.println(\"" + element + "\");\n";
						break;
					case QUOTES:
						out += "\"" + element + "\"+\"\\n\"+\n";
						break;
					default:
						break;
				}
			}

			ext.setClipboard(out);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
