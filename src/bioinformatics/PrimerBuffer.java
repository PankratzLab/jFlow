package bioinformatics;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;

import seq.manage.ReferenceGenome;
import common.Array;
import common.Files;
import common.Logger;
import common.Positions;
import common.ext;
import filesys.Segment;

public class PrimerBuffer {
	private static final String[] HEADER = new String[] { "CHR", "START", "STOP", "TARET_SEQUENCE" };
	private static final String[] HEADER_OUT_ADD = new String[] { "BUFFER_LOC", "BUFFER_SEQUENCE_TOTAL_LENGTH", "SEQUENCE" };

	private static void extractVariants(String queryFile, String referenceGenomeFast, int bpBuffer, Logger log) {
		ReferenceGenome referenceGenome = new ReferenceGenome(referenceGenomeFast, log);
		String output = ext.addToRoot(queryFile, ".query");
		ArrayList<ReferenceAlleleQuery> rAlleleQueries = new ArrayList<ReferenceAlleleQuery>();
		try {
			BufferedReader reader = Files.getAppropriateReader(queryFile);
			int[] header = ext.indexFactors(reader.readLine().trim().split("\t"), HEADER, true, true);
			while (reader.ready()) {
				String[] line = reader.readLine().trim().split("\t");
				byte chr = Positions.chromosomeNumber(line[header[0]]);
				int start = Integer.parseInt(line[header[1]]);
				int stop = Integer.parseInt(line[header[2]]);
				Segment seg = new Segment(chr, start, stop);
				String[] q = new String[seg.getSize()];
				Arrays.fill(q, ReferenceAlleleQuery.BLANK);
				if (line.length > 3) {
					String tmp = line[header[3]];
					q = new String[tmp.length()];
					for (int i = 0; i < q.length; i++) {
						q[i] = tmp.charAt(i) + "";
					}
				}
				rAlleleQueries.add(new ReferenceAlleleQuery(seg, q));
			}
			reader.close();
		} catch (FileNotFoundException fnfe) {
			log.reportError("Error: file \"" + queryFile + "\" not found in current directory");
			return;
		} catch (IOException ioe) {
			log.reportError("Error reading file \"" + queryFile + "\"");
			return;
		}
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(output));
			writer.println(Array.toStr(HEADER) + "\t" + Array.toStr(HEADER_OUT_ADD));
			for (int i = 0; i < rAlleleQueries.size(); i++) {
				rAlleleQueries.get(i).populateQuery(referenceGenome, bpBuffer);
				writer.println(Array.toStr(rAlleleQueries.get(i).getResult()));
			}
			writer.close();
		} catch (Exception e) {
			log.reportError("Error writing to " + output);
			log.reportException(e);
		}
	}

	private static class ReferenceAlleleQuery {
		private static final String BLANK = "N";
		private Segment seg;
		private String[] sequence;
		private String[] surroundingSeguence;
		private Segment buffered;
		private int targetStart;
		private int targetStop;

		private ReferenceAlleleQuery(Segment seg, String[] sequence) {
			super();
			this.seg = seg;
			this.sequence = sequence;
			if (seg.getSize() != sequence.length && Array.countIf(sequence, BLANK) != sequence.length) {
				throw new IllegalArgumentException("Segment size must equal sequence size");
			}
		}

		public String[] getResult() {
			//	private static final String[] HEADER = new String[] { "CHR", "START", "STOP", "TARET_SEQUENCE" };
			// private static final String[] HEADER_OUT_ADD = new String[] { "BUFFER_LOC", "BUFFER_SEQUENCE_TOTAL_LENGTH", "SEQUENCE" };

			ArrayList<String> result = new ArrayList<String>();
			result.add(seg.getChr() + "");
			result.add(seg.getStart() + "");
			result.add(seg.getStop() + "");
			result.add(Array.toStr(sequence, ""));
			result.add(buffered.getUCSClocation());
			result.add(buffered.getSize() + "");
			result.add(Array.toStr(surroundingSeguence, ""));

			return Array.toStringArray(result);
			
		}

		public String[] getSurroundingSeguence() {
			return surroundingSeguence;
		}

		private void populateQuery(ReferenceGenome referenceGenome, int bpBuffer) {
			this.buffered = seg.getBufferedSegment(bpBuffer);
			this.surroundingSeguence = referenceGenome.getSequenceFor(buffered);
			this.targetStart = seg.getStart() - buffered.getStart();
			this.targetStop = targetStart+seg.getSize()-1;
			String[] targets = Array.subArray(surroundingSeguence, targetStart, targetStop + 1);
			if (targets.length != sequence.length) {
				System.out.println(Array.toStr(targets));
				System.out.println(Array.toStr(sequence));

				throw new IllegalStateException("Internal error, extracted mismatching bases");
			} else {
				for (int i = 0; i < targets.length; i++) {
					if (!sequence[i].equals(BLANK) && !targets[i].equals(sequence[i])) {
						throw new IllegalStateException("Mismatched target sequence extraction, verify correct reference");

					}
				}
			}
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String referenceGenome = "ref.fasta";
		int bpBuffer = 250;
		String queryFile = null;
		String usage = "\n" + "bioinformatics.PrimerBuffer requires 0-1 arguments\n";
		usage += "   (1) filename to query (i.e. file=" + queryFile + " (default))\n" + "";
		usage += "   (2) reference genome to query (i.e. ref=" + referenceGenome + " (default))\n" + "";
		usage += "   (3) bp buffer on each side (i.e. buf=" + bpBuffer + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("file=")) {
				queryFile = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("ref=")) {
				referenceGenome = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("buf=")) {
				bpBuffer = ext.parseIntArg(args[i]);
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
			Logger log = new Logger(ext.rootOf(queryFile) + ".log");

			extractVariants(queryFile, referenceGenome, bpBuffer, log);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
