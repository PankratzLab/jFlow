package one.JL;

import java.io.IOException;
import java.util.ArrayList;

import one.ScriptExecutor;
import common.Array;
import common.CmdLine;
import common.Files;
import common.HashVec;
import common.Logger;
import common.ext;

/**
 * One timer but could be more, for getting specific file types from amazon buckets, I think there is an actual s3cmd to do this but oh well
 *
 */
public class S3cmdDL {

	private static void getEm(String dir, String bucket, String exten, int numThreads) {
		String sh = dir + "ls.sh";
		String txt = dir + "ls.txt";
		String get = dir + "get.sh";

		Logger log = new Logger();
		ArrayList<String> cmd1 = new ArrayList<String>();
		cmd1.add("s3cmd");
		cmd1.add("ls");
		cmd1.add(bucket);
		cmd1.add(">");
		cmd1.add(txt);
		Files.write(Array.toStr(Array.toStringArray(cmd1), " "), sh);
		Files.chmod(sh);
		CmdLine.runCommandWithFileChecks(new String[] { sh }, "", null, null, true, true, false, log);

		String[] files = HashVec.loadFileToStringArray(txt, false, null, false);
		ArrayList<String> filesGet = new ArrayList<String>();
		for (int i = 0; i < files.length; i++) {
			out: if (files[i].endsWith(exten)) {
				String actual = files[i].substring(files[i].indexOf(" s3:"));
				String dl = dir + ext.removeDirectoryInfo(actual);
				// String dlVerify = dl+".verify";
				if (!Files.exists(dl)) {
					// ||!Files.exists(dlVerify)
					int count = 0;
					for (int j = 0; j < actual.length(); j++) {
						if ((actual.charAt(j) + "").equals("/")) {
							count++;
						}
						if (count > 3) {
							break out;
						}
					}
					System.out.println(actual);
					filesGet.add("s3cmd get " + actual + "\n");
					// "; touch " + dlVerify +
				}
			}
		}

		Files.write(Array.toStr(Array.toStringArray(filesGet), ""), get);
		Files.chmod(get);
		log.reportTimeInfo("Beginning retrieve : " + filesGet.size() + " files of type " + exten);
		ScriptExecutor executor = new ScriptExecutor(numThreads);
		try {
			executor.run(get, null);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		// CmdLine.runCommandWithFileChecks(new String[] { get }, "", null, null, true, true, false, log);

	}

	public static void main(String[] args) {
		// TODO, all this if needed.
		String dir = "/scratch.global/lane0212/Project_Spector_Project_014/";
		String bucket = "s3://Project_Spector_Project_014/";
		String it = "gvcf";
		String it2 = "gvcf.idx";
		int numThreads = 24;
		getEm(dir, bucket, it2, numThreads);
		getEm(dir, bucket, it, numThreads);

	}

}
