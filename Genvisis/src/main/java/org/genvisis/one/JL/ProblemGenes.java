package one.JL;

import java.io.FileWriter;
import java.io.PrintWriter;

import common.Array;

import cnv.filesys.Project;
import filesys.GeneTrack;

public class ProblemGenes {

	public static void dumpProblems(Project proj, String output, String[] startWithPatters) {
		GeneTrack geneTrack = GeneTrack.load(proj.getGeneTrackFilename(false), false);
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(output));
			for (int i = 0; i < geneTrack.getGenes().length; i++) {
				for (int j = 0; j < geneTrack.getGenes()[i].length; j++) {
					String geneName = geneTrack.getGenes()[i][j].getGeneName();
					for (int k = 0; k < startWithPatters.length; k++) {
						if (geneName.startsWith(startWithPatters[k])) {
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
		String[] startWithPatters = new String[] { "HLA", "ZNF", "MUC","OR" };
		String output = "badGenes.txt";
		String usage = "\n" + "one.JL.ProblemGenes requires 0-1 arguments\n";
		usage += "   (1) project filename (i.e. proj=" + filename + " (default))\n" + "";
		usage += "   (2) gene patters to remove, comma delimited (i.e. remove=" + Array.toStr(startWithPatters, ",") + " (default))\n" + "";
		usage += "   (3) gene patters to remove, comma delimited (i.e. out=" + output + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("remove=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
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
