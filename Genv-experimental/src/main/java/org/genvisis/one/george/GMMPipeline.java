package org.genvisis.one.george;

import java.io.File;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class GMMPipeline {
	String fcs_dirpath;
	
	public GMMPipeline(String fcs_dirpath) {
		this.fcs_dirpath = fcs_dirpath;
	}
	// Fits GMM on fcs datasets converted to csvs
	public void run() {
		StringBuilder out = new StringBuilder();
		
		// get list of csv data sets
		Tuple<String[],String[]> tup = FileManipulator.filesWithExt(".csv", fcs_dirpath);
		String[] fileListExt = tup.x; // just filenames
		String[] pathListExt = tup.y; // file paths
		
		// run gmm model
		
		// plot
	}

	public static void main(String[] args) {
		GMMPipeline pipeline = new GMMPipeline("./fcs");
		pipeline.run();
	}
}
