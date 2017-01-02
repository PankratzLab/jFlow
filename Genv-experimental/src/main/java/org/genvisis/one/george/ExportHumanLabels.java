package org.genvisis.one.george;

import java.io.File;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ExportHumanLabels {
	public static void main(String[] args) {
		if (args.length < 1) {
			System.out.println("err: requires directory name "
					+ "(containing .xln files) as first command line arg");
			return;
		}
		else {
			Pattern p = Pattern.compile(".*xln");

			String dirname = args[0];
			File dir = new File(dirname);
			File[] fileList = dir.listFiles();
			String[] fileListExt = new String[fileList.length];
			String[] pathListExt = new String[fileList.length];
		
			int j = 0;
			for (int i = 0; i < fileList.length; i++) {
				String filename = fileList[i].getName();
				Matcher m = p.matcher(filename);
				if (m.matches()) {
					fileListExt[j] = filename;
					pathListExt[j] = dirname + "/" + filename;
					j++;
				}
			}
			
			for (int i = 0; i < fileListExt.length; i++) {
				if (fileListExt[i] != null) {
					HumanLabels hl = new HumanLabels(pathListExt[i]);
					String out_filepath = pathListExt[i].substring(0, pathListExt[i].indexOf(".xln")) + "_mtrx.xln"; 
					hl.exportData(out_filepath);
				}
			}
		}
	}
}
