package link;

import java.io.*;
import common.*;

public class Recode implements Runnable {
	private String dir;

	private int chr;

	public Recode(String dir, int chr) {
		this.dir = dir;
		this.chr = chr;
	}

	public void run() {
		BufferedReader reader;
		PrintWriter writer;
		String[] line, markerNames;
		String chrome = ext.chrome(chr);
		CountVector[] cvs;
		int numCovariates = 0;
		String[][] orderedAlleles;
		LinkageMap map;

		map = new LinkageMap(dir, chr);
		markerNames = map.getMarkerNames();
		cvs = CountVector.initArray(markerNames.length);

		try {
			reader = new BufferedReader(new FileReader(dir+"chrom"+chrome+".pre"));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (line[0].equals("")) {
					continue;
				}
				if (line.length%2==1) {
					numCovariates = 1;
					System.out.println("Assuming there is a covaraite in chrom"+chrome+".pre");
				}
				for (int i = 0; i<markerNames.length; i++) {
					cvs[i].add(line[6+numCovariates+i*2+0]);
					cvs[i].add(line[6+numCovariates+i*2+1]);
				}
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+dir+"chrom"+chrome+".pre"+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading file \""+dir+"chrom"+chrome+".pre"+"\"");
			System.exit(2);
		}

		orderedAlleles = map.recode(cvs);
		map.updateFile();

		try {
			reader = new BufferedReader(new FileReader(dir+"chrom"+chrome+".pre"));
			writer = new PrintWriter(new FileWriter(dir+"re_chrom"+chrome+".pre"));
			while (reader.ready()) {
				line = reader.readLine().trim().split("[\\s]+");
				if (line[0].equals("")) {
					continue;
				}
				writer.print(line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\t"+line[4]+"\t"+line[5]);
				for (int i = 0; i<numCovariates; i++) {
					writer.print("\t"+line[6+i]);
				}
				for (int i = 0; i<markerNames.length; i++) {
					writer.print("\t"+(1+ext.indexOfStr(line[6+numCovariates+i*2+0], orderedAlleles[i]))+"\t"+(1+ext.indexOfStr(line[6+numCovariates+i*2+1], orderedAlleles[i])));
				}
				writer.println();
			}
			reader.close();
			writer.close();
		} catch (FileNotFoundException fnfe) {
			System.err.println("Error: file \""+dir+"chrom"+chrome+".pre"+"\" not found in current directory");
			System.exit(1);
		} catch (IOException ioe) {
			System.err.println("Error reading \""+dir+"chrom"+chrome+".pre"+"\" or writing \""+"re_chrom"+chrome+".pre"+"\"");
			System.exit(2);
		}
	}

	public static void recode(String dir, int start, int stop, int numThreads) {
		Thread threads[] = new Thread[numThreads];
		boolean done = false;
		
		for (int chr = start; chr<=stop || !done;) {
			done = true;
			for (int i = 0; i<threads.length; i++) {
				if ((threads[i]==null || !threads[i].isAlive()) && chr<=stop) {
					threads[i] = new Thread(new Recode(dir, chr));
					threads[i].start();
					chr++;
					done = false;
				} else if (threads[i].isAlive()) {
					done = false;
				}

			}

			try {
				Thread.sleep(500L);
			} catch (InterruptedException ex) {}
		}

	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String dir = "";
		// String dir = "C:\\Documents and Settings\\npankrat\\My
		// Documents\\hearing\\";
		int start = -1;
		int stop = -1;
		int threads = 1;

		String usage = "\n"+"link.recode requires 0-3 arguments\n"+"   recode all chromosomes: recode\n"+"   recode in specific directory: recode dir=check/\n"+"   recode chromosome 6: recode 6\n"+"   recode chromosomes 1 through 5: recode 1 5\n"+"   recode all using 4 threads: recode threads=4\n";

		for (int i = 0; i<args.length; i++) {
			try {
				if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
					System.err.println(usage);
					System.exit(1);
				} else if (args[i].startsWith("dir=")) {
					dir = args[i].split("=")[1];
					numArgs--;
				} else if (args[i].startsWith("threads=")) {
					threads = Integer.parseInt(args[i].split("=")[1]);
					numArgs--;
				} else if (args[i].substring(args[i].lastIndexOf("\\")+1).startsWith("chrom")&&args[i].endsWith(".pre")) {
					start = Integer.parseInt(args[i].substring(args[i].lastIndexOf("\\")+1).substring(5, 7));
					numArgs--;
				} else if (start==-1) {
					start = Integer.parseInt(args[i]);
					numArgs--;
				} else {
					stop = Integer.parseInt(args[i]);
					numArgs--;
				}
			} catch (NumberFormatException nfe) {
				System.err.println("Error - '"+args[i]+"' is an invalid argument");
			}
		}

		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}

		if (start==-1) {
			start = 1;
			stop = 23;
		} else if (stop==-1) {
			stop = start;
		}

		try {
			recode(dir, start, stop, threads);
		} catch (RuntimeException re) {
			System.err.println("Error: "+re.getMessage());
			if (re.getMessage().length()<10)
				re.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
