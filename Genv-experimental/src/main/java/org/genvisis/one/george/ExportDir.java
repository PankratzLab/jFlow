package org.genvisis.one.george;
import java.io.File;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.io.IOException;

import org.genvisis.one.ben.fcs.FCSDataLoader;
import org.genvisis.one.ben.fcs.FCSDataLoader.LOAD_STATE;

// deprecated, see FileManipulator
public class ExportDir {
	// exports all fcs in directory to csv files
	public static void main(String[] args) {
		if (args.length < 1) {
			System.out.println("err: requires directory name "
					+ "(containing .fcs files) as first command line arg");
			return;
		}
		else {
			Pattern p = Pattern.compile(".*fcs");
			
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
					FCSDataLoader loader = new FCSDataLoader();
					
					try {loader.loadData(pathListExt[i]);}
					catch (IOException e) {};
					loader.waitForData();
					
					String filename = fileListExt[i];
					int extIndex = filename.indexOf(".");
					String nameNoExt = filename.substring(0, extIndex);
					loader.exportData(dirname + "/" + nameNoExt + ".csv", FCSDataLoader.DATA_SET.ALL);
				}
			}
		}
	}
}
