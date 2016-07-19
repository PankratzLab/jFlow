package filesys;

import common.*;

public class VariantList {
	public static final String[][] HEADER = {{"chr"}, {"position"}, {"ref"}, {"alt"}};
			
	public static String[][] parseList(String filename, Logger log) {
		String[] header;
		int[] indices;
		
		header = Files.getHeaderOfFile(filename, log);
		indices = ext.indexFactors(HEADER, header, false, true, true, log, true);
		
		return HashVec.loadFileToStringMatrix(filename, true, new int[] {indices[0], indices[1], indices[2], indices[3]}, false);
	}
}
