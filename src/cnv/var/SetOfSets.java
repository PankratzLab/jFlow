package cnv.var;

import java.util.*;

public class SetOfSets {
	private ArrayList<int[]> sets;
	private ArrayList<String> filenames;

	public SetOfSets() {
		sets = new ArrayList<int[]>();
		filenames = new ArrayList<String>();
	}

	public void add(int[] set, String filename) {
		sets.add(set);
		filenames.add(filename);
	}

	public int[][] getSets() {
		return sets.toArray(new int[sets.size()][]);
	}

	public String[] getFilenames() {
		return filenames.toArray(new String[filenames.size()]);
	}
}
