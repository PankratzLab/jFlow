package org.genvisis.one.JL;

import java.io.FileWriter;
import java.io.PrintWriter;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.common.ArrayUtils;
import org.genvisis.filesys.GeneTrack;

public class ProblemGenes {

	public static void dumpProblems(Project proj, String output, String[] startWithPatters) {
		GeneTrack geneTrack = GeneTrack.load(proj.getGeneTrackFilename(false), false);
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(output));
			for (int i = 0; i < geneTrack.getGenes().length; i++) {
				for (int j = 0; j < geneTrack.getGenes()[i].length; j++) {
					String geneName = geneTrack.getGenes()[i][j].getGeneName();
					for (String startWithPatter : startWithPatters) {
						if (geneName.startsWith(startWithPatter)) {
							writer.println(geneTrack.getGenes()[i][j].getUCSClocation() + "\t" + geneName);
						}
					}
				}
			}
			writer.close();
		} catch (Exception e) {
			proj.getLog().reportError("Error writing to " + output);
			proj.getLog().reportException(e);
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		String[] startWithPatters = new String[] {"HLA", "ZNF", "MUC", "OR"};
		String output = "badGenes.txt";
		String usage = "\n" + "one.JL.ProblemGenes requires 0-1 arguments\n";
		usage += "   (1) project filename (i.e. proj=" + filename + " (default))\n" + "";
		usage += "   (2) gene patters to remove, comma delimited (i.e. remove="
							+ ArrayUtils.toStr(startWithPatters, ",") + " (default))\n" + "";
		usage += "   (3) gene patters to remove, comma delimited (i.e. out="	+ output + " (default))\n"
							+ "";

		for (String arg : args) {
			if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (arg.startsWith("proj=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else if (arg.startsWith("remove=")) {
				filename = arg.split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + arg);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			Project proj = new Project(filename, false);
			dumpProblems(proj, proj.PROJECT_DIRECTORY.getValue() + output, startWithPatters);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
