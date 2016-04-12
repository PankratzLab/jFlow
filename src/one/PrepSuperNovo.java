package one;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Hashtable;

import common.Files;
import common.Logger;
import common.ext;

/**
 * Currently specific to spector data
 *
 */
public class PrepSuperNovo {
	private static String[] TRIO_ENDINGS = new String[] { "C", "D", "M" };

	public static void prepDir(String dir, String extension, String outputDir, Logger log) {
		String[] bams = Files.list(dir, extension, false);
		Hashtable<String, Hashtable<String, String>> trios = new Hashtable<String, Hashtable<String, String>>();
		for (int i = 0; i < bams.length; i++) {
			String id = bams[i].split("_")[0];
			String trioId = id.substring(0, id.length() - 1);
			char trioMember = id.charAt(id.length() - 1);

			int index = ext.indexOfStr(trioMember + "", TRIO_ENDINGS);
			if (index < 0) {
				System.out.println("ERROR not trio\t" + trioMember + "\t" + id);
				return;
			} else if (!trios.containsKey(trioId)) {
				trios.put(trioId, new Hashtable<String, String>());
			}
			trios.get(trioId).put(trioMember + "", ext.removeDirectoryInfo(bams[i]));

		}
		String trioListFile = outputDir + "trios.trio";
		new File(outputDir).mkdirs();
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(trioListFile));

			for (String trio : trios.keySet()) {
				String trioentry = trio;
				boolean write = true;
				for (int i = 0; i < TRIO_ENDINGS.length; i++) {
					if (trios.get(trio).containsKey(TRIO_ENDINGS[i])) {
						trioentry += "\t" + trios.get(trio).get(TRIO_ENDINGS[i]);
					} else {
						write = false;
					}
				}
				if(write){
				writer.println(trioentry);
				}else{
					log.reportTimeError("Skipping un filled trio for " +trioentry);
				}
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + trioListFile);
			log.reportException(e);
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String dir = "/lustre/lanej0/Project_Spector_Project_014/140516_SN261_0548_AC4GVNACXX/sam/";
		String extension = ".sorted.dedup.realigned.recal.bam";
		String outputDir = "/lustre/lanej0/Project_Spector_Project_014/140516_SN261_0548_AC4GVNACXX/superNovo/";

		prepDir(dir, extension, outputDir, new Logger());
	}
}
