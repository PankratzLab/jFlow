package org.genvisis.one.JL;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;

import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.VCFOps.VcfPopulation;

public class GenVCFTrio {

	public static void genList(String list, Logger log) {
		Hashtable<String, HashSet<String>> trios = new Hashtable<String, HashSet<String>>();

		try {
			BufferedReader reader = Files.getAppropriateReader(list);
			while (reader.ready()) {
				String[] line = reader.readLine().trim().split("[\\s]+");
				String fam = line[0].substring(0, line[0].length() - 1);
				if (!trios.containsKey(fam)) {
					trios.put(fam, new HashSet<String>());
				}
				trios.get(fam).add(line[0]);
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + list + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + list + "\"");
			return;
		}
		String ouput = ext.addToRoot(list, ".vcfPop");
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(ouput));
			writer.println(ArrayUtils.toStr(VcfPopulation.HEADER));
			Set<String> fams = trios.keySet();
			for (String fam : fams) {
				HashSet<String> tmp = trios.get(fam);
				for (String ind : tmp) {
					if (ind.endsWith("C")) {
						writer.println(ind + "\t" + fam + "\tOFFSPRING");
					} else {
						writer.println(ind + "\t" + fam + "\tCONTROL");
					}
				}
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + ouput);
			log.reportException(e);
		}
	}

	public static void main(String[] args) {
		String osList = "D:/data/logan/OSv2_seq/RegNovo/OsSamps.txt";
		genList(osList, new Logger());
	}

}
