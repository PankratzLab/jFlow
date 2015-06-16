package cnv.analysis.pca;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;

import common.Array;
import common.ext;
import stats.SimpleM;
import stats.Rscript.RScatter;
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

	public static SelectionResult v(Project proj, double filterValue, STAT_TYPE sType, SELECTION_TYPE selType) {
		PCSelector selector = new PCSelector(proj, sType);
		ArrayList<Integer> sigPCs = new ArrayList<Integer>();
		SelectionResult rankResult = null;
		ArrayList<StatsCrossTabRank> ranks = new ArrayList<StatsCrossTabRank>();

		if (selector.isValid()) {
			proj.getLog().reportTimeInfo("Stat type: " + sType);
			proj.getLog().reportTimeInfo("Selection type: " + selType);

			while (selector.hasNext()) {
				ranks.add(selector.next());
			}
			switch (selType) {
			case EFFECTIVE_M_CORRECTED:
				int numTests = selector.determineEffectiveNumberOfTests();
				proj.getLog().reportTimeInfo("Controling type I error at " + filterValue + "; " + filterValue + "/" + numTests + " = " + ((double) filterValue / numTests));
				filterValue = (double) filterValue / numTests;
				proj.getLog().reportTimeInfo("Filter value : " + filterValue);
				break;
			case P_VAL:
				proj.getLog().reportTimeInfo("Filter value : " + filterValue);
				break;
			case STAT:
				proj.getLog().reportTimeInfo("Filter value : " + filterValue);
				break;
			default:
				proj.getLog().reportTimeError("Invalid selection type " + selType);
				break;

			}
		}

		return rankResult;

	}

	public static SelectionResult select(Project proj, double absStatMin, STAT_TYPE sType, SELECTION_TYPE selType) {

		PCSelector selector = new PCSelector(proj, sType);
		ArrayList<Integer> sigPCs = new ArrayList<Integer>();
		SelectionResult rankResult = null;
		if (selector.isValid()) {
			double filterValue = Double.NaN;

			String output = ext.addToRoot(proj.INTENSITY_PC_FILENAME.getValue(), ".significantPCs");
			try {
				ArrayList<StatsCrossTabRank> ranks = new ArrayList<StatsCrossTabRank>();
				Hashtable<String, Integer> has = new Hashtable<String, Integer>();
				PrintWriter writer = new PrintWriter(new FileWriter(output));

				writer.println("TYPE\tQC_METRIC\t" + Array.toStr(selector.getpResiduals().getPcTitles()));
				while (selector.hasNext()) {
					StatsCrossTabRank sRank = selector.next();
					ranks.add(sRank);
					writer.println("SIG\t" + sRank.getRankedTo() + "\t" + Array.toStr(sRank.getSigs()));
					writer.println("STAT\t" + sRank.getRankedTo() + "\t" + Array.toStr(sRank.getStats()));
					double bonf = (double) absStatMin / LrrSd.NUMERIC_COLUMNS.length * selector.getpResiduals().getPcTitles().length;

					for (int i = 0; i < sRank.getStats().length; i++) {
						boolean add = false;
						if (!has.containsKey((i + 1) + "")) {
							switch (selType) {
							case EFFECTIVE_M_CORRECTED:

								break;
							case P_VAL:
								break;
							case STAT:
								if (Math.abs(sRank.getStats()[i]) > absStatMin) {
									add = true;

								}
								break;
							default:
								proj.getLog().reportTimeError("Invalid selection type " + selType);
								break;
							}
						}
						if (add) {
							sigPCs.add(i + 1);
							has.put((i + 1) + "", (i + 1));
						}
					}
				}
				writer.close();
				proj.getLog().reportTimeInfo("Found " + sigPCs.size() + " pcs passing threshold of " + absStatMin);
				String[] minMax = new String[] { "Min_" + sType, "Max_" + sType };
				String outputT = ext.addToRoot(output, ".transposedStat");
				writer = new PrintWriter(new FileWriter(outputT));
				writer.print("PCTitle\tPC\t" + Array.toStr(minMax));
				for (int i = 0; i < ranks.size(); i++) {
					writer.print("\t" + ranks.get(i).getRankedTo());
				}
				writer.println();
				String[] titles = selector.getpResiduals().getPcTitles();
				for (int i = 0; i < titles.length; i++) {
					writer.print(titles[i] + "\t" + (i + 1) + "\t" + (-1 * absStatMin) + "\t" + absStatMin);
					for (int j = 0; j < ranks.size(); j++) {
						writer.print("\t" + ranks.get(j).getStats()[i]);
					}
					writer.println();
				}

				writer.close();
				String title = "QC_Association: n=" + sigPCs.size() + " PCs at abs(r) > " + absStatMin;
				RScatter rScatter = new RScatter(outputT, outputT + ".rscript", ext.rootOf(outputT), outputT + ".pdf", "PC", Array.concatAll(minMax, LrrSd.NUMERIC_COLUMNS), SCATTER_TYPE.POINT, proj.getLog());
				rScatter.setyLabel(sType.toString());
				rScatter.setOverWriteExisting(true);
				rScatter.setxLabel("PC");
				rScatter.setTitle(title);
				rScatter.execute();
				rankResult = new SelectionResult(Array.toIntArray(sigPCs), rScatter);
			} catch (Exception e) {
				proj.getLog().reportError("Error writing to " + output);
				proj.getLog().reportException(e);
			}

		} else {
			proj.getLog().reportTimeError("Could not select QC associated PCs...");
		}
		return rankResult;

	}

	public static class SelectionResult {
		private int[] order;
		private RScatter rScatter;

		public SelectionResult(int[] order, RScatter rScatter) {
			super();
			this.order = order;
			this.rScatter = rScatter;
		}

		public int[] getOrder() {
			return order;
		}

		public RScatter getrScatter() {
			return rScatter;
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
		select(proj, absStatMin, STAT_TYPE.SPEARMAN_CORREL, SELECTION_TYPE.STAT);
		try {
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
