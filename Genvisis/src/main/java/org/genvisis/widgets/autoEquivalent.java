// want to autogenerate code for getVar and setVar?
package org.genvisis.widgets;

import org.genvisis.common.PSF;
import org.genvisis.common.ext;

public class autoEquivalent {
	public static void equivocate() {
		String[] line, clip;
		String trav, trans;

		clip = ext.getClipboard().split("\\n");
		trans = "";
		for (String element : clip) {
			line = element.split(PSF.Regex.GREEDY_WHITESPACE);
			trav = line[line.length - 1];
			trav = trav.substring(0, trav.length() - 1);
			trans += "this." + trav + " = " + trav + ";\n";
		}
		ext.setClipboard(trans);
	}

	public static void main(String[] argv) {
		equivocate();
	}
}
