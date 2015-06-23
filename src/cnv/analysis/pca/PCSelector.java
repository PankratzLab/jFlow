package cnv.analysis.pca;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;

import common.Array;
import common.ext;
import stats.SimpleM;
import stats.Rscript.PLOT_DEVICE;
import stats.Rscript.RScatter;
import stats.Rscript.RScatters;
import stats.Rscript.SCATTER_TYPE;
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

	public SampleQC getSampleQC() {
		return sampleQC;
	}

	private void load() {
		this.pResiduals = proj.loadPcResids();
		pResiduals.fillInMissing();
		valid = pResiduals != null;
		if (valid) {
			this.sampleQC = SampleQC.loadSampleQC(proj);
			valid = sampleQC != null;
			this.index = 0;
		}
	}

	public int determineEffectiveNumberOfTests() {
		SimpleM.Builder builder = new SimpleM.Builder();
		builder.verbose(false);
		SimpleM simpleMQcVar = builder.build(sampleQC.getQcMatrix(), proj.getLog());
		int effQc = simpleMQcVar.determineM();
		SimpleM simpleMPCVar = builder.build(pResiduals.getPcBasis(), proj.getLog());
		int effPC = simpleMPCVar.determineM();
		int numTests = effPC * effQc;
		proj.getLog().reportTimeInfo("Effective num QC metrics tested = " + effQc);
		proj.getLog().reportTimeInfo("Effective num Pcs tested = " + effPC);
		proj.getLog().reportTimeInfo("Effective number of tests = " + effQc + " X " + effPC + " = " + numTests);
		return numTests;
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

	public enum SELECTION_TYPE {
		/**
		 * filtered by the {@link STAT_TYPE } actual stat
		 */
		STAT, /**
		 * filtered by the {@link STAT_TYPE } p-value
		 */
		P_VAL,
		/**
		 * pval cutoff after effective M correction, see {@link SimpleM}
		 */
		EFFECTIVE_M_CORRECTED;
	}

	public static SelectionResult select(Project proj, double filterValue, STAT_TYPE sType, SELECTION_TYPE selType) {
		PCSelector selector = new PCSelector(proj, sType);
		ArrayList<Integer> sigPCs = new ArrayList<Integer>();
		SelectionResult rankResult = null;
		StatsCrossTabRank[] ranks = new StatsCrossTabRank[LrrSd.NUMERIC_COLUMNS.length];

		if (selector.isValid()) {
			proj.getLog().reportTimeInfo("Stat type: " + sType);
			proj.getLog().reportTimeInfo("Selection type: " + selType);

			int index = 0;
			while (selector.hasNext()) {
				ranks[index] = selector.next();
				index++;
			}
			String title = "QC_Association: ";

			switch (selType) {
			case EFFECTIVE_M_CORRECTED:
				int numTests = selector.determineEffectiveNumberOfTests();
				proj.getLog().reportTimeInfo("Controling type I error at " + filterValue + "; " + filterValue + "/" + numTests + " = " + ((double) filterValue / numTests));
			//	String originalP = filterValue + "";
				filterValue = (double) filterValue / numTests;
				title += "p < " + ext.prettyP(filterValue) + " (" + numTests + " tests)";
				proj.getLog().reportTimeInfo("Filter value : " + filterValue);
				break;
			case P_VAL:
				proj.getLog().reportTimeInfo("Filter value : " + filterValue);
				title += "p < " + ext.prettyP(filterValue);
				break;
			case STAT:
				proj.getLog().reportTimeInfo("Filter value : " + filterValue);
				title += "abs(r) > " + filterValue;
				break;
			default:
				proj.getLog().reportTimeError("Invalid selection type " + selType);
				break;

			}
			Hashtable<String, Integer> has = new Hashtable<String, Integer>();

			for (int i = 0; i < ranks.length; i++) {
				StatsCrossTabRank sRank = ranks[i];
				for (int j = 0; j < sRank.getStats().length; j++) {
					boolean add = false;
					if (!has.containsKey((j + 1) + "")) {
						switch (selType) {
						case EFFECTIVE_M_CORRECTED:
							if (sRank.getSigs()[j] < filterValue) {
								add = true;
							}
							break;
						case P_VAL:
							if (sRank.getSigs()[j] < filterValue) {
								add = true;
							}
							break;
						case STAT:
							if (Math.abs(sRank.getStats()[j]) > filterValue) {
								add = true;
							}
							break;
						default:
							proj.getLog().reportTimeError("Invalid selection type " + selType);
							break;
						}
					}
					if (add) {
						sigPCs.add(j + 1);
						has.put((j + 1) + "", (j + 1));
					}
				}
			}
			
			int[] finalSelection = Array.toIntArray(sigPCs);
			proj.getLog().reportTimeInfo("Found " + sigPCs.size() + " pcs passing threshold of " + filterValue);
			
			String dir = proj.PROJECT_DIRECTORY.getValue()+"PC_QC/";
			new File(dir).mkdirs();
			String base = ext.removeDirectoryInfo(proj.INTENSITY_PC_FILENAME.getValue());
			
			String outputAll = dir + ext.addToRoot(base, ".QC_assoc.all");
			String outputSelect = dir + ext.addToRoot(base, ".QC_assoc.selected");
			String outputBoth = dir + ext.addToRoot(base, ".QC_assoc");

			String[] titles = summarize(ranks, selector.getpResiduals().getPcTitles(), null, outputAll);
			summarize(ranks, selector.getpResiduals().getPcTitles(), finalSelection, outputSelect);

			RScatter rScatterAll = plot(proj, sType, " n=" + selector.getpResiduals().getPcTitles().length + " total PCs", selector.getpResiduals().getPcTitles().length + 1, outputAll, titles);
			RScatter rScatterSelect = plot(proj, sType, title + "; n=" + sigPCs.size() + " QC PCs", selector.getpResiduals().getPcTitles().length + 1, outputSelect, titles);

			RScatters rScatters = new RScatters(new RScatter[] { rScatterAll, rScatterSelect }, outputBoth + ".rscript", outputBoth + ".pdf", null, PLOT_DEVICE.PDF, proj.getLog());
			rScatters.execute();
			rankResult = new SelectionResult(Array.toIntArray(sigPCs), rScatterAll, rScatterSelect);

		}
		return rankResult;
	}

	private static RScatter plot(Project proj, STAT_TYPE sType, String title,int xmax, String output, String[] titles) {
		RScatter rScatter = new RScatter(output, output + ".rscript", ext.rootOf(output), output + ".pdf", "PC", titles, SCATTER_TYPE.POINT, proj.getLog());
		rScatter.setyLabel(sType.toString());
		rScatter.setOverWriteExisting(true);
		rScatter.setxLabel("PC");
		rScatter.setTitle(title);
		rScatter.setxRange(new double[] { 0, xmax  });
		return rScatter;
	}

	private static String[] summarize(StatsCrossTabRank[] ranks, String[] titles, int[] sigPCs, String output) {
		ArrayList<String> titleSummary = new ArrayList<String>();

		try {
			PrintWriter writer = new PrintWriter(new FileWriter(output));
			writer.print("PCTitle\tPC");
			for (int i = 0; i < ranks.length; i++) {
				titleSummary.add(ranks[i].getRankedTo());
				writer.print("\t" + ranks[i].getRankedTo() + "\t" + ranks[i].getRankedTo() + "_p");
			}
			writer.println();
			for (int i = 0; i < titles.length; i++) {
				if (sigPCs == null || ext.indexOfInt((i + 1), sigPCs) >= 0) {
					writer.print(titles[i] + "\t" + (i + 1));
					for (int j = 0; j < ranks.length; j++) {
						writer.print("\t" + ranks[j].getStats()[i] + "\t" + ranks[j].getSigs()[i]);
					}
					writer.println();
				}
			}
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return titleSummary.toArray(new String[titleSummary.size()]);
	}



	public static class SelectionResult {
		private int[] order;
		private RScatter rScatterAll, rScatterSelect;

		public SelectionResult(int[] order, RScatter rScatterAll, RScatter rScatterSelect) {
			super();
			this.order = order;
			this.rScatterAll = rScatterAll;
			this.rScatterSelect = rScatterSelect;
		}

		public int[] getOrder() {
			return order;
		}

		public RScatter getrScatter() {
			return rScatterAll;
		}

		public RScatter getrScatterSelect() {
			return rScatterSelect;
		}

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
		//select(proj, absStatMin, STAT_TYPE.SPEARMAN_CORREL, SELECTION_TYPE.STAT);
		try {
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}

//public static SelectionResult select(Project proj, double absStatMin, STAT_TYPE sType, SELECTION_TYPE selType) {
//
//	PCSelector selector = new PCSelector(proj, sType);
//	ArrayList<Integer> sigPCs = new ArrayList<Integer>();
//	SelectionResult rankResult = null;
//	if (selector.isValid()) {
//		double filterValue = Double.NaN;
//
//		String output = ext.addToRoot(proj.INTENSITY_PC_FILENAME.getValue(), ".significantPCs");
//		try {
//			ArrayList<StatsCrossTabRank> ranks = new ArrayList<StatsCrossTabRank>();
//			Hashtable<String, Integer> has = new Hashtable<String, Integer>();
//			PrintWriter writer = new PrintWriter(new FileWriter(output));
//
//			writer.println("TYPE\tQC_METRIC\t" + Array.toStr(selector.getpResiduals().getPcTitles()));
//			while (selector.hasNext()) {
//				StatsCrossTabRank sRank = selector.next();
//				ranks.add(sRank);
//				writer.println("SIG\t" + sRank.getRankedTo() + "\t" + Array.toStr(sRank.getSigs()));
//				writer.println("STAT\t" + sRank.getRankedTo() + "\t" + Array.toStr(sRank.getStats()));
//				double bonf = (double) absStatMin / LrrSd.NUMERIC_COLUMNS.length * selector.getpResiduals().getPcTitles().length;
//
//				for (int j = 0; j < sRank.getStats().length; j++) {
//					boolean add = false;
//					if (!has.containsKey((j + 1) + "")) {
//						switch (selType) {
//						case EFFECTIVE_M_CORRECTED:
//
//							break;
//						case P_VAL:
//							break;
//						case STAT:
//							if (Math.abs(sRank.getStats()[j]) > absStatMin) {
//								add = true;
//
//							}
//							break;
//						default:
//							proj.getLog().reportTimeError("Invalid selection type " + selType);
//							break;
//						}
//					}
//					if (add) {
//						sigPCs.add(j + 1);
//						has.put((j + 1) + "", (j + 1));
//					}
//				}
//			}
//			writer.close();
//			proj.getLog().reportTimeInfo("Found " + sigPCs.size() + " pcs passing threshold of " + absStatMin);
//			String[] minMax = new String[] { "Min_" + sType, "Max_" + sType };
//			String outputT = ext.addToRoot(output, ".transposedStat");
//			writer = new PrintWriter(new FileWriter(outputT));
//			writer.print("PCTitle\tPC\t" + Array.toStr(minMax));
//			for (int i = 0; i < ranks.size(); i++) {
//				writer.print("\t" + ranks.get(i).getRankedTo());
//			}
//			writer.println();
//			String[] titles = selector.getpResiduals().getPcTitles();
//			for (int i = 0; i < titles.length; i++) {
//				writer.print(titles[i] + "\t" + (i + 1) + "\t" + (-1 * absStatMin) + "\t" + absStatMin);
//				for (int j = 0; j < ranks.size(); j++) {
//					writer.print("\t" + ranks.get(j).getStats()[i]);
//				}
//				writer.println();
//			}
//
//			writer.close();
//			String title = "QC_Association: n=" + sigPCs.size() + " PCs at abs(r) > " + absStatMin;
//			RScatter rScatter = new RScatter(outputT, outputT + ".rscript", ext.rootOf(outputT), outputT + ".pdf", "PC", Array.concatAll(minMax, LrrSd.NUMERIC_COLUMNS), SCATTER_TYPE.POINT, proj.getLog());
//			rScatter.setyLabel(sType.toString());
//			rScatter.setOverWriteExisting(true);
//			rScatter.setxLabel("PC");
//			rScatter.setTitle(title);
//			rScatter.execute();
//			//rankResult = new SelectionResult(Array.toIntArray(sigPCs), rScatter);
//		} catch (Exception e) {
//			proj.getLog().reportError("Error writing to " + output);
//			proj.getLog().reportException(e);
//		}
//
//	} else {
//		proj.getLog().reportTimeError("Could not select QC associated PCs...");
//	}
//	return rankResult;
//
//}
