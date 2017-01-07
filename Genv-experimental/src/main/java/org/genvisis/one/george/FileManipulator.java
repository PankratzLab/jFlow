package org.genvisis.one.george;

import java.io.File;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.genvisis.one.ben.fcs.FCSDataLoader;

// TEST THIS CODE
public class FileManipulator {
	/** Returns list of file names ending with extension @param ext (e.g. ".csv") 
	 * @param dir_path - path to fcs directory, e.g. "./fcs" */
	public static Tuple<String[], String[]> filesWithExt(String ext, String dir_path) {
		Pattern p = Pattern.compile(".*\\" + ext);
		File dir = new File(dir_path);
		File[] fileList = dir.listFiles();
		String[] fileListExt = new String[fileList.length];
		String[] pathListExt = new String[fileList.length];
		
		int j = 0;
		for (int i = 0; i < fileList.length; i++) {
			String filename = fileList[i].getName();
			Matcher m = p.matcher(filename);
			if (m.matches()) {
				fileListExt[j] = filename;
				pathListExt[j] = dir_path + "/" + filename;
				j++;
			}
		}
		// j is the length of filled array
		String[] fileListExtFilled = new String[j];
		String[] pathListExtFilled = new String[j];
		
		for (int i = 0; i < j; i++) {
			fileListExtFilled[i] = fileListExt[i];
			pathListExtFilled[i] = pathListExt[i];
		}
		return new Tuple<>(fileListExtFilled, pathListExtFilled);	
	}
	/** @param fcs_dirpath - e.g. "./fcs" (no trailing forward slash) */
	public static void exportDir(String fcs_dirpath) {
		Tuple<String[], String[]> tup = filesWithExt(".fcs", fcs_dirpath);
		String[] fileListExt = tup.x;
		String[] pathListExt = tup.y;	
		
		for (int i = 0; i < fileListExt.length; i++) {
			if (fileListExt[i] != null) {
				FCSDataLoader loader = new FCSDataLoader();
				
				try {loader.loadData(pathListExt[i]);}
				catch (IOException e) {};
				loader.waitForData();
				
				String filename = fileListExt[i];
				int extIndex = filename.indexOf(".fcs");
				String nameNoExt = filename.substring(0, extIndex);
				loader.exportData(fcs_dirpath + "/" + nameNoExt + ".csv", FCSDataLoader.DATA_SET.ALL);
			}
		}
	}
	public static void exportHumanLabels(String xln_dirpath) {
		
		Tuple<String[], String[]> tup = filesWithExt(".xln", xln_dirpath);
		String[] fileListExt = tup.x;
		String[] pathListExt = tup.y;	
		
		for (int i = 0; i < fileListExt.length; i++) {
			if (fileListExt[i] != null) {
				HumanLabels hl = new HumanLabels(pathListExt[i]);
				String out_filepath = pathListExt[i].substring(0, pathListExt[i].indexOf(".xln")) + "_mtrx.xln"; 
				hl.exportData(out_filepath);
			}
		}
	}
}
