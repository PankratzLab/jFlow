package org.genvisis.cnv.plots;

import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.JFrame;
import javax.swing.WindowConstants;

import org.genvisis.CLI.Arg;
import org.genvisis.cnv.plots.data.AbstractPipe.DropIfMatchAnyPipe;
import org.genvisis.cnv.plots.data.DataFile;
import org.genvisis.cnv.plots.data.DataListener;
import org.genvisis.cnv.plots.data.DataPipe;
import org.genvisis.cnv.plots.data.FilterPipe;
import org.genvisis.cnv.plots.data.TransformPipe;
import org.genvisis.common.Aliases;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.stats.Maths.OPERATOR;

public class ManhattanPlot extends JFrame {

	private static final int SNP_LINKER = 0;
	private static final int CHR_LINKER = 1;
	private static final int POS_LINKER = 2;
	private static final int PVAL_LINKER = 3;

	public static final String[][] LINKERS = {Aliases.MARKER_NAMES, Aliases.CHRS, Aliases.POSITIONS,
																						Aliases.PVALUES};
	private static final boolean[] REQ_LINKS = {false, true, true, true};

	static final int LIN_CHR_BUFFER = 10000000;
	static final int LIN_LOC_START = Integer.MIN_VALUE + LIN_CHR_BUFFER;

	/**
	 * TODO possible virtual columns: <br />
	 * - FID + IID <br />
	 * - CHR:POS <br />
	 * - linearized CHR:POS <br />
	 * - transformed values of any sort<br />
	 */

	ManhattanPanel manPan;
	Logger log;
	DataFile data;
	HashMap<String, DataPipe> dataPipes;

	public static void main(String[] args) {
		ManhattanPlot mp = new ManhattanPlot();
		mp.loadFile("D:/data/ny_choanal/omni2.5v1.2/plink/plink.assoc");
	}

	public ManhattanPlot() {
		super("Genvisis - ManhattanPlot");
		setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);

		manPan = new ManhattanPanel(this);
		dataPipes = new HashMap<>();

		addTestDataPipes();

		getContentPane().add(manPan);
		setBounds(100, 100, 800, 600);

		setVisible(true);
		log = new Logger();
	}

	private void addTestDataPipes() {
		DataPipe pipe = new DataPipe();
		pipe.addPipe(new TransformPipe() {
			@Override
			public String transformValue(String value) {
				Double v = Double.parseDouble(value);
				return "" + (-Math.log10(v));
			}
		});
		dataPipes.put("P", pipe);
	}

	public ArrayList<ManhattanDataPoint> createData(String[] markerNames, int[] chrs, int[] pos,
																									double[] pvals) {
		return new ArrayList<>(); // TODO
	}

	public void setData(ArrayList<ManhattanDataPoint> dataPoints) {
		this.cachedData = dataPoints;
		// TODO new object "manualData", use cachedData for transforms / filters
	}

	public void loadFile(String filename) {
		if (this.data != null && this.data.getFilename().equals(filename)) {
			return;
		}
		if (!Files.exists(filename)) {
			log.reportError("File not found: " + filename);
			return;
		}

		DataFile dataFile = new DataFile(filename, log, LINKERS, REQ_LINKS);
		if (!dataFile.hasRequiredData()) {
			// TODO error, cancel
		}
		this.data = dataFile;

		// add DataPipes to remove missing data
		HashMap<String, DataPipe> selTrans = new HashMap<>();
		DataPipe nonMissPipe = new DataPipe();
		nonMissPipe.addPipe(new DropIfMatchAnyPipe(ext.MISSING_VALUES, true));
		nonMissPipe.addPipe(new FilterPipe<Integer>(OPERATOR.GREATER_THAN, 0, Integer.class));
		selTrans.put(dataFile.getLinkedColumnName(CHR_LINKER), nonMissPipe);
		selTrans.put(dataFile.getLinkedColumnName(POS_LINKER), nonMissPipe);

		// add DataPipe to remove missing data and remove pVals < 0.05
		DataPipe pValPipe = new DataPipe();
		pValPipe.addPipe(new DropIfMatchAnyPipe(ext.MISSING_VALUES, true));
		pValPipe.addPipe(new FilterPipe<Double>(OPERATOR.LESS_THAN_OR_EQUAL, 0.05, Double.class));
		selTrans.put(dataFile.getLinkedColumnName(PVAL_LINKER), pValPipe);

		// No data pipe, but indicate that we'd like to load SNP data also
		// but check to ensure we have SNP data first
		if (dataFile.hasLinkedColumn(SNP_LINKER)) {
			selTrans.put(dataFile.getLinkedColumnName(SNP_LINKER), null);
		}

		// TODO allow column selection? Load other info?

		dataFile.setSortColumns(new String[] {dataFile.getLinkedColumnName(CHR_LINKER),
																					dataFile.getLinkedColumnName(POS_LINKER)},
														new Arg[] {Arg.NUMBER, Arg.NUMBER});


		DataListener pinger = new DataListener() {
			@Override
			public void ping(DataFile dataFile) {
				ManhattanPlot.this.manPan.paintAgain();
			}
		};
		dataFile.pingWhenReady(pinger);
		dataFile.loadData(selTrans, false);
		this.manPan.paintAgain();
	}

	public ArrayList<ManhattanDataPoint> getData() {
		ArrayList<ManhattanDataPoint> pts = getCachedData();
		if (pts == null) {
			if (data == null) {
				return null;
			}
			pts = new ArrayList<>();
			if (!data.isLoaded()) {
				return pts;
			}
			String chrCol = data.getLinkedColumnName(CHR_LINKER);
			String posCol = data.getLinkedColumnName(POS_LINKER);
			String pValCol = data.getLinkedColumnName(PVAL_LINKER);
			int chrInd = data.getIndexOfData(chrCol);
			int posInd = data.getIndexOfData(posCol);
			int pValInd = data.getIndexOfData(pValCol);
			int mkrInd = data.getIndexOfData(data.getLinkedColumnName(SNP_LINKER));
			// TODO add other data / info / pipes

			ArrayList<String[]> dataRows = data.getAllData();
			for (int i = 0, cnt = dataRows.size(); i < cnt; i++) {
				String[] row = dataRows.get(i);
				String chr = row[chrInd];
				String pos = row[posInd];
				String pVal = row[pValInd];
				if (dataPipes.containsKey(chrCol)) {
					chr = dataPipes.get(chrCol).pipe(chr);
					if (chr == null)
						continue;
				}
				if (dataPipes.containsKey(posCol)) {
					pos = dataPipes.get(posCol).pipe(pos);
					if (pos == null)
						continue;
				}
				if (dataPipes.containsKey(pValCol)) {
					pVal = dataPipes.get(pValCol).pipe(pVal);
					if (pVal == null)
						continue;
				}
				// TODO other columns and pipes
				int chrI = Integer.parseInt(chr);
				int posI = Integer.parseInt(pos);
				int linLoc = LIN_LOC_START;
				if (pts.size() > 0) {
					ManhattanDataPoint prev = pts.get(pts.size() - 1);
					if (chrI > prev.chr) {
						linLoc = prev.linearLoc + LIN_CHR_BUFFER;
					} else {
						linLoc = prev.linearLoc + (Math.min(posI - prev.pos, LIN_CHR_BUFFER));
					}
				}
				double val = Double.parseDouble(pVal);
				pts.add(new ManhattanDataPoint(mkrInd == -1 ? chr + ":" + pos : row[mkrInd], chrI, posI,
																			 linLoc, val));
				// TODO add other data here
			}
			cachedData = pts;
		}
		return pts;
	}

	ArrayList<ManhattanDataPoint> cachedData;

	private ArrayList<ManhattanDataPoint> getCachedData() {
		return cachedData; // TODO caching based on x/y transforms
	}

	class ManhattanDataPoint {
		public ManhattanDataPoint(String mkr, int chrI, int posI, int linLoc, double val) {
			this.mkr = mkr;
			this.chr = chrI;
			this.pos = posI;
			this.linearLoc = linLoc;
			this.transformedPVal = val;
			this.otherData = new HashMap<>();
		}

		String mkr;
		int chr, pos;
		int linearLoc;
		double transformedPVal;
		HashMap<String, String> otherData;
	}

}
