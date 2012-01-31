// want to autogenerate code for getVar and setVar?
package widgets;

import common.ext;

public class autoEquivalent {
	public static void equivocate() {
		String[] line, clip;
		String trav, trans;

		clip = ext.getClipboard().split("\\n");
		trans = "";
		for (int i = 0; i<clip.length; i++) {
			line = clip[i].split("[\\s]+");
			trav = line[line.length-1];
			trav = trav.substring(0, trav.length()-1);
			trans += "this."+trav+" = "+trav+";\n";
		}
		ext.setClipboard(trans);
	}

	public static void main(String[] argv) {
		equivocate();
	}
}
