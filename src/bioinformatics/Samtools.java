package bioinformatics;

import java.io.*;
import java.util.*;
import java.util.Map.Entry;

import one.SuperNovo;
import common.*;

public class Samtools {

	public static void extractRegions(String outputFilename, String fullPathToSuperNovoPhase2Output, String fullPathToTrioNameList, int windowInBp, Logger log) {
		Hashtable<String, Vector<String>> hash = new Hashtable<String, Vector<String>>();
		String[][] bamFilenamesByTrios;

		hash = loadFromFile(fullPathToSuperNovoPhase2Output, log);
		bamFilenamesByTrios = SuperNovo.loadNamesFromList(fullPathToTrioNameList);
		writerToFile(outputFilename, null, hash, bamFilenamesByTrios, windowInBp, log);
	}
	
	public static Hashtable<String, Vector<String>> loadFromFile(String filename, Logger log) {
		BufferedReader reader;
		String[] line;
		Hashtable<String, Vector<String>> hash;
		Vector<String> v;

		hash = new Hashtable<String, Vector<String>>();
		try {
			reader = new BufferedReader(new FileReader(filename));
			while (reader.ready()) {
				line = reader.readLine().split("\t");
				if(hash.containsKey(line[2])) {
					v = hash.get(line[2]);
				} else {
					v = new Vector<String>();
				}
				v.add(line[4]);
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return hash;
	}

	public static void writerToFile(String filename, String xmlDir, Hashtable<String, Vector<String>> hash, String[][] bamFilenamesByTrios, int windowInBp, Logger log) {
		PrintWriter writer;
		Vector<String> v;
		int index;
		Set<String> keySet;
		String chr;
		int pos;

		keySet = hash.keySet();
		try {
			for (String key : keySet) {
				chr = key.split(":")[0];
				pos = Integer.parseInt(key.split(":")[1]);
				v = hash.get(key);
				for (int i = 0; i < v.size(); i++) {
					index = getIndex(v.elementAt(i), bamFilenamesByTrios, log);

					writer = new PrintWriter(new FileWriter(xmlDir + ".bamScript"));
					writer.println("samtools view " + bamFilenamesByTrios[index][1] + "-b " + chr + ":" + (pos - windowInBp) + "-" + (pos + windowInBp) + " > " + bamFilenamesByTrios[index][0] + "_chr" + chr + "_" + pos+"_C");
					writer.close();

					writer = new PrintWriter(new FileWriter(filename + ".bamScript"));
					writer.println("samtools view " + bamFilenamesByTrios[index][2] + "-b " + chr + ":" + (pos - windowInBp) + "-" + (pos + windowInBp) + " > " + bamFilenamesByTrios[index][0] + "_chr" + chr + "_" + pos+"_D");
					writer.close();

					writer = new PrintWriter(new FileWriter(filename + ".bamScript"));
					writer.println("samtools view " + bamFilenamesByTrios[index][3] + "-b " + chr + ":" + (pos - windowInBp) + "-" + (pos + windowInBp) + " > " + bamFilenamesByTrios[index][0] + "_chr" + chr + "_" + pos+"_M");
					writer.close();

					writer = new PrintWriter(new FileWriter(filename + ".IGV_Script"));
					writer.println(getIgvScript(xmlDir, bamFilenamesByTrios[index]));
					writer.close();
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static int getIndex(String trioId, String[][] bamFilenamesByTrios, Logger log) {
		int index;

		index = -1;
		for (int i = 0; i < bamFilenamesByTrios.length; i++) {
			if (bamFilenamesByTrios[i][0].equals(trioId)) {
				index = i;
				break;
			}
		}

		return index;
	}
	
	public static String getIgvScript(String dir, String[] bamFilenamesByTrios) {
		return "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>"
				+ "\n<Session genome=\"hg19\" hasGeneTrack=\"true\" hasSequenceTrack=\"true\" locus=\"chr10:104574306-104574443\" version=\"8\">"
					+ "\n<Resources>"
						+ "\n<Resource path=\"" + dir + bamFilenamesByTrios[1] + "\"/>"
						+ "\n<Resource path=\"" + dir + bamFilenamesByTrios[2] + "\"/>"
						+ "\n<Resource path=\"" + dir + bamFilenamesByTrios[3] + "\"/>"
					+ "\n</Resources>"
				+ "</Session>";
	}

	/**
	 * 
	 * @param filename
	 * @param trioIdChrPos		trioId + chr + pos
	 * @param bamFilenamesByTrios
	 * @param windowInBp
	 * @param log
	 */
	public static void writerToFile(String filename, Vector<String> trioIdChrPos, String[][] bamFilenamesByTrios, int windowInBp, Logger log) {
		PrintWriter writer;
		Vector<String> v;
		int index;
		String[] line;
		String chr;
		int pos;

		try {
			writer = new PrintWriter(new FileWriter(filename));
			for (String key : trioIdChrPos) {
				line = key.split("\t");
				chr = line[1];
				pos = Integer.parseInt(line[2]);
				index = getIndex(line[0], bamFilenamesByTrios, log);

				writer.print("samtools view " + bamFilenamesByTrios[index][1] + " -b chr" + chr + ":" + Math.max(0, pos - windowInBp) + "-" + (pos + windowInBp) + " > " + bamFilenamesByTrios[index][0] + "_chr" + chr + "_" + pos+"_C.bam\n"
							+ "samtools index " + bamFilenamesByTrios[index][0] + "_chr" + chr + "_" + pos+"_C.bam\n"
							+ "samtools view " + bamFilenamesByTrios[index][2] + " -b chr" + chr + ":" + Math.max(0, pos - windowInBp) + "-" + (pos + windowInBp) + " > " + bamFilenamesByTrios[index][0] + "_chr" + chr + "_" + pos+"_D.bam\n"
							+ "samtools index " + bamFilenamesByTrios[index][0] + "_chr" + chr + "_" + pos+"_D.bam\n"
							+ "samtools view " + bamFilenamesByTrios[index][3] + " -b chr" + chr + ":" + Math.max(0, pos - windowInBp) + "-" + (pos + windowInBp) + " > " + bamFilenamesByTrios[index][0] + "_chr" + chr + "_" + pos+"_M.bam\n"
							+ "samtools index " + bamFilenamesByTrios[index][0] + "_chr" + chr + "_" + pos+"_M.bam\n"
							);
			}
			
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	public static HashSet<String> listFiles(String miniBamDir, Logger log) {
		String[] filenames;
		boolean[] isVisited;
		String[] line1;
		String[] line2;
		HashSet<String> result;
		int isAllThreeFilesExist;

		result = new HashSet<String>();
		filenames = Files.list(miniBamDir, ".bam", false);
		isVisited = new boolean[filenames.length];
		for (int i = 0; i < filenames.length; i++) {
			if (! isVisited[i]) {
				isAllThreeFilesExist = 1;
				line1 = filenames[i].split("_");
				for (int j = 0; j < filenames.length; j++) {
					line2 = filenames[j].split("_");
					if (! isVisited[j] && j != i && line2[0].equals(line1[0]) && line2[1].equals(line1[1]) && line2[2].equals(line1[2])) {
						isVisited[j] = true;
						isAllThreeFilesExist ++;
					}
				}
				if (isAllThreeFilesExist == 3) {
					result.add(line1[0] + "_" + line1[1] + "_" + line1[2]);
					if (! new File(miniBamDir + "/xmls/" + line1[1] + "_" + line1[2] + ".xml").exists()) {
						Files.write(getIgvScript(miniBamDir, new String[0]), miniBamDir + "/xmls/" + line1[1] + "_" + line1[2] + ".xml");
					}
				}
			}
		}

		return result;
	}


	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "Samtools.dat";
		String fullPathToSuperNovoPhase2Output = "N:/statgen/OS_Logan/SuperNovo/rawOutput/SuperNovo_summary.txt";
		String fullPathToTrioNameList = "Samtools.dat";
		int windowInBp = 5000;
		String logfile = null;
		Logger log;

		String usage = "\n" + "bioinformatics.Samtools requires 0-1 arguments\n" + "   (1) filename (i.e. file=" + filename + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				filename = ext.parseStringArg(args[i], null);
				numArgs--;
			} else if (args[i].startsWith("window=")) {
				windowInBp = ext.parseIntArg(args[i]);
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = ext.parseStringArg(args[i], null);
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
			log = new Logger(logfile);
			extractRegions(filename, fullPathToSuperNovoPhase2Output, fullPathToTrioNameList, windowInBp, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
