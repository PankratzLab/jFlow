package dead;

import java.io.*;

import common.*;

public class ParseChisquare {
	public static final String DEFAULT_DIR = "C:\\Documents and Settings\\npankrat\\My Documents\\PD-singleton\\";

	public static final String DEFAULT_FILE = "09_Susceptibility\\03_chisquares\\plink.model";

	public static final String[] HEADER = {"CHR", "SNP", "TEST", "AFF", "UNAFF", "CHISQ", "DF", "P"};

	public static final String[] MODELS = {"GENO", "TREND", "ALLELIC", "DOM", "REC"};

	public static final String GENES = "genes.xls";

	// public static final String MAP = "pd.recode.map";
	public static final String MAP = "pd.recode.map-composite.out";

	public static final String FREQS = "pd.frq";

	public static final int N_TOP_HITS = 50;

	public static final boolean IGNORE_SEX = true;

	public static void parseChisquare(String dir, String filename) {
		BufferedReader reader = null;
		PrintWriter writer = null;
		String[] line;
		String trav;
		double best;
		int bestIndex;

		try {
			reader = new BufferedReader(new FileReader(dir+filename));
			writer = new PrintWriter(new FileWriter(dir+filename+".out"));
			ext.checkHeader(reader.readLine().trim().split("[\\s]+"), HEADER, true);
			writer.println("SNP\t"+Array.toStr(MODELS)+"\tBest Model\tp-value");
			while (reader.ready()) {
				trav = "";
				best = 1;
				bestIndex = -1;
				for (int i = 0; i<MODELS.length; i++) {
					line = reader.readLine().trim().split("[\\s]+");
					if (i==0) {
						writer.print(line[1]);
						trav = line[1];
					} else if (!line[1].equals(trav)) {
						System.err.println("Error mixup at the intersection of "+trav+" and "+line[1]);
					}
					if (!line[2].equals(MODELS[i])) {
						System.err.println("Models out of whack at "+line[1]);
					}
					writer.print("\t"+line[7]);
					if (!line[7].equals("NA")&&Double.parseDouble(line[7])<best) {
						best = Double.parseDouble(line[7]);
						bestIndex = i;
					}
				}
				writer.println(bestIndex==-1?"\t.\t.":"\t"+MODELS[bestIndex]+"\t"+best);
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+dir+filename+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+dir+filename+"\"");
			System.exit(2);
		}

	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String dir = DEFAULT_DIR;
		String filename = DEFAULT_FILE;

		String usage = "\n"+"park.parseChisquare requires 0-1 arguments\n"+"   (1) filename (i.e. "+filename+" (default))\n"+"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else {
				filename = args[i];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			parseChisquare(dir, filename);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
