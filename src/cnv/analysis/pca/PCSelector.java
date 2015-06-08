package cnv.analysis.pca;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;

import common.Array;
import common.Sort;
import common.ext;
import link.sepByFam;
import stats.StatsCrossTabs.STAT_TYPE;
import stats.StatsCrossTabs.StatsCrossTabRank;
import stats.StatsCrossTabs.VALUE_TYPE;
import cnv.filesys.Project;
import cnv.qc.LrrSd;
import cnv.qc.SampleQC;

/**
 * @author lane0212 <br>
 *         Try to find most significant PCs associated with our quality metrics from {@link LrrSd}
 */
public class PCSelector implements Iterator<StatsCrossTabRank> {
	private Project proj;
	private PrincipalComponentsResiduals pResiduals;
	private SampleQC sampleQC;
	private STAT_TYPE sType;

	private boolean valid;
	private int index;

	public PCSelector(Project proj, STAT_TYPE sType) {
		super();
		this.proj = proj;
		this.sType = sType;
		load();
	}

	public PrincipalComponentsResiduals getpResiduals() {
		return pResiduals;
	}

	private void load() {
		this.pResiduals = proj.loadPcResids();
		valid = pResiduals != null;
		if (valid) {
			this.sampleQC = SampleQC.loadSampleQC(proj);
			valid = sampleQC != null;
			this.index = 0;
		}
	}

	public boolean isValid() {
		return valid;
	}

	@Override
	public boolean hasNext() {
		return index < LrrSd.NUMERIC_COLUMNS.length;
	}

	@Override
	public StatsCrossTabRank next() {
		String currentQC = LrrSd.NUMERIC_COLUMNS[index];
		proj.getLog().reportTimeInfo("Analyzing QC metric " + currentQC);
		StatsCrossTabRank sRank = pResiduals.getStatRankFor(sampleQC.getDataFor(LrrSd.NUMERIC_COLUMNS[index]), null, null, currentQC, sType, VALUE_TYPE.STAT, proj.getLog());
		index++;
		return sRank;
	}

	@Override
	public void remove() {
	}

	public static int[] select(Project proj, double absStatMin) {
		PCSelector selector = new PCSelector(proj, STAT_TYPE.PEARSON_CORREL);
		ArrayList<Integer> sigPCs = new ArrayList<Integer>();

		if (selector.isValid()) {
			String output = ext.addToRoot(proj.INTENSITY_PC_FILENAME.getValue(), ".significantPCs");
			try {
				HashSet<Integer> has = new HashSet<Integer>();
				PrintWriter writer = new PrintWriter(new FileWriter(output));
				writer.println("TYPE\tQC_METRIC\t" + Array.toStr(selector.getpResiduals().getPcTitles()));
				while (selector.hasNext()) {
					StatsCrossTabRank sRank = selector.next();
					writer.println("SIG\t" + sRank.getRankedTo() + "\t" + Array.toStr(sRank.getSigs()));
					writer.println("STAT\t" + sRank.getRankedTo() + "\t" + Array.toStr(sRank.getStats()));
					for (int i = 0; i < sRank.getStats().length; i++) {
						if (has.add(i + 1) && Math.abs(sRank.getStats()[i]) > absStatMin) {
							sigPCs.add(i + 1);
						}
					}
				}
				writer.close();
			} catch (Exception e) {
				proj.getLog().reportError("Error writing to " + output);
				proj.getLog().reportException(e);
			}
		} else {
			proj.getLog().reportTimeError("Could not select QC associated PCs...");
		}
		return Array.toIntArray(sigPCs);
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = null;
		double absStatMin = 0.10;

		String usage = "\n" + "cnv.analysis.pca.PCSelector requires 0-1 arguments\n";
		usage += "   (1) project filename (i.e. proj=" + filename + " ( no default))\n" + "";
		usage += "   (2) the minimum (absolute value of the test statistic) across all qc metrics (i.e. statMin=" + absStatMin + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("statMin=")) {
				absStatMin = ext.parseDoubleArg(args[i]);
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		Project proj = new Project(filename, false);
		select(proj, absStatMin);
		try {
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
