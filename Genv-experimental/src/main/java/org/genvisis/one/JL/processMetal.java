package org.genvisis.one.JL;

import org.genvisis.common.Files;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

/**
 * @author Kitty One off for processing .metal files to .beta files
 */
public class processMetal {

	public static void main(String[] args) {
		String[] conv = new String[] { "rsID", "ref", "alt", "beta", "p" };
		String[] required = new String[] { "MarkerName", "Allele1", "Allele2", "Effect", "P.value" };
		String dir = "/Volumes/Beta/data/wbcGwasMetal/oldWbcMeta/";
		String[] metals = Files.listFullPaths(dir, ".metal", false);
		Logger log = new Logger();
		for (String metal : metals) {
			String[] header = Files.getHeaderOfFile(metal, log);
			int[] indices = ext.indexFactors(required, header, true, true);
			String[][] file = HashVec.loadFileToStringMatrix(metal, false, indices, false);
			file[0] = conv;
			Files.writeMatrix(file, ext.rootOf(metal, false) + ".beta", "\t");
		}

	}

}
