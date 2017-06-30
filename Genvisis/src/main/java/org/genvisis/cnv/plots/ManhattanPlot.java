package org.genvisis.cnv.plots;

import java.awt.event.ActionEvent;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.AbstractAction;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.WindowConstants;

import org.genvisis.CLI.Arg;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.plots.data.AbstractPipe;
import org.genvisis.cnv.plots.data.AbstractPipe.DropIfMatchAnyPipe;
import org.genvisis.cnv.plots.data.DataFile;
import org.genvisis.cnv.plots.data.DataListener;
import org.genvisis.cnv.plots.data.DataPipe;
import org.genvisis.cnv.plots.data.FilterPipe;
import org.genvisis.common.Aliases;
import org.genvisis.common.ArrayUtils;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.stats.Maths.OPERATOR;

public class ManhattanPlot extends JFrame {

	static final int SNP_LINKER = 0;
	static final int CHR_LINKER = 1;
	static final int POS_LINKER = 2;
	static final int PVAL_LINKER = 3;

	static final String[][] LINKERS = {Aliases.MARKER_NAMES, Aliases.CHRS, Aliases.POSITIONS,
																		 Aliases.PVALUES, Aliases.EFFECTS, Aliases.STD_ERRS};
	static final boolean[] REQ_LINKS = {false, true, true, true, false, false};

	static final int LIN_CHR_BUFFER = 10000000;
	static final int LIN_LOC_START = Integer.MIN_VALUE + LIN_CHR_BUFFER;

	ManhattanPanel manPan;
	Project proj;
	Logger log;
	DataFile data;
	HashMap<String, DataPipe> dataPipes;

	public ManhattanPlot(Project proj) {
		super("Genvisis - ManhattanPlot");
		setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);

		createMenubar();

		this.proj = proj;
		manPan = new ManhattanPanel(this);
		dataPipes = new HashMap<>();

		getContentPane().add(manPan);
		setBounds(100, 100, 800, 600);

		log = proj == null ? new Logger() : proj.getLog();
	}

	private void createMenubar() {
		JMenuBar menuBar = new JMenuBar();

		JMenu fileM = new JMenu("File");
		fileM.setMnemonic('F');
		menuBar.add(fileM);

		JMenuItem openFileItem = new JMenuItem();
		openFileItem.setAction(new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent arg0) {
				loadFile();
			}
		});
		openFileItem.setText("Open File");
		openFileItem.setMnemonic('O');
		fileM.add(openFileItem);


		this.setJMenuBar(menuBar);
	}

	public Project getProject() {
		return proj;
	}

	/**
	 * Automatically loads a file based on defaults. Only loads CHR/POS/PVAL/SNP columns and filters
	 * PVAL for <= 0.05.
	 * 
	 * @param filename File to load
	 */
	public void loadFileAuto(String filename) {
		if ("".equals(filename)) {
			return;
		}
		if (!Files.exists(filename)) {
			log.reportTime("ERROR - file not found: " + filename);
			return;
		}

		DataFile dataFile = new DataFile(filename, log, LINKERS, REQ_LINKS);
		if (!dataFile.hasRequiredData()) {
			boolean missChr = dataFile.getLinkedColumnName(CHR_LINKER) == null;
			boolean missPos = dataFile.getLinkedColumnName(POS_LINKER) == null;
			boolean missPVal = dataFile.getLinkedColumnName(PVAL_LINKER) == null;
			String miss = "";
			if (missChr) {
				miss += "Chr";
			}
			if (missPos) {
				if (miss.length() > 0) {
					miss += ", ";
				}
				miss += "Pos";
			}
			if (missPVal) {
				if (miss.length() > 0) {
					miss += ", ";
				}
				miss += "P-val";
			}
			log.reportTime("ERROR - missing " + miss + " data!");
			return;
		}

		this.data = dataFile;
		String[] cols = new String[0];
		if (dataFile.getLinkIndices()[SNP_LINKER] >= 0) {
			cols = new String[] {dataFile.getLinkedColumnName(SNP_LINKER)};
		}
		boolean[] chrs = ArrayUtils.booleanArray(27, true);
		chrs[0] = false;
		double filt = 0.05;

		readyData(chrs, filt, cols);
	}

	/**
	 * 
	 * @param markerNames
	 * @param chrs
	 * @param pos
	 * @param pvals Raw P-values, will be transformed with a -Log10 function.
	 * @return List of data to be used with setData()
	 */
	public ArrayList<ManhattanDataPoint> createData(String[] markerNames, int[] chrs, int[] pos,
																									double[] pvals) {
		int len1 = markerNames.length;
		int len2 = chrs.length;
		int len3 = pos.length;
		int len4 = pvals.length;
		if (len1 == len2 && len2 == len3 && len3 == len4) {
			ArrayList<ManhattanDataPoint> data = new ArrayList<>();
			for (int i = 0; i < len1; i++) {
				int linLoc = LIN_LOC_START;
				if (data.size() > 0) {
					int chrI = chrs[i];
					int posI = pos[i];
					ManhattanDataPoint prev = data.get(data.size() - 1);
					if (chrI > prev.chr) {
						linLoc = prev.linearLoc + LIN_CHR_BUFFER;
					} else {
						linLoc = prev.linearLoc + (Math.min(posI - prev.pos, LIN_CHR_BUFFER));
					}
				}
				double pVal = -Math.log10(pvals[i]);
				data.add(new ManhattanDataPoint(markerNames[i], chrs[i], pos[i], linLoc, pVal, pvals[i]));
			}
			return data;
		} else {
			throw new IllegalArgumentException("All variable arrays must be the same length!");
		}
	}

	public void setData(ArrayList<ManhattanDataPoint> dataPoints) {
		this.cachedData = dataPoints;
		// TODO new object "manualData", use cachedData for transforms / filters
	}


	public void loadFile() {
		ManhattanLoadGUI mlg = new ManhattanLoadGUI();
		mlg.setVisible(true);

		if (mlg.getCloseCode() != JFileChooser.APPROVE_OPTION) {
			return;
		}

		String[] cols = mlg.getSelectedColumns();
		double filt = mlg.getPValueThreshold();
		boolean[] chrs = mlg.getSelectedChrs();
		this.data = mlg.getDataFile();

		readyData(chrs, filt, cols);
	}

	private void readyData(boolean[] chrsToLoad, double pvalThresh, String[] extraColsToLoad) {
		HashMap<String, DataPipe> selTrans = new HashMap<>();

		DataPipe chrPipe = new DataPipe();
		chrPipe.addPipe(new DropIfMatchAnyPipe(ext.MISSING_VALUES, true));
		int numChrs = ArrayUtils.booleanArraySum(chrsToLoad);
		if (numChrs != chrsToLoad.length) {
			if (numChrs == 1) {
				int chr = ArrayUtils.booleanArrayToIndices(chrsToLoad)[0];
				chrPipe.addPipe(new FilterPipe<Integer>(OPERATOR.EQUAL, chr, Integer.class));
			} else {
				String[] dropIfMatch = ArrayUtils.toStringArray(ArrayUtils.booleanArrayToIndices(ArrayUtils.booleanNegative(chrsToLoad)));
				chrPipe.addPipe(new DropIfMatchAnyPipe(dropIfMatch, true));
			}
		}
		selTrans.put(data.getLinkedColumnName(CHR_LINKER), chrPipe);

		DataPipe posPipe = new DataPipe();
		posPipe.addPipe(new DropIfMatchAnyPipe(ext.MISSING_VALUES, true));
		selTrans.put(data.getLinkedColumnName(POS_LINKER), posPipe);

		DataPipe pValPipe = new DataPipe();
		pValPipe.addPipe(new DropIfMatchAnyPipe(ext.MISSING_VALUES, true));
		if (!Double.isNaN(pvalThresh)) {
			pValPipe.addPipe(new FilterPipe<Double>(OPERATOR.LESS_THAN_OR_EQUAL, pvalThresh, Double.class));
		}
		selTrans.put(data.getLinkedColumnName(PVAL_LINKER), pValPipe);

		for (String s : extraColsToLoad) {
			if (!selTrans.containsKey(s)) {
				selTrans.put(s, null);
			}
		}

		data.setSortColumns(new String[] {data.getLinkedColumnName(CHR_LINKER),
																			data.getLinkedColumnName(POS_LINKER)},
												new Arg[] {Arg.NUMBER, Arg.NUMBER});

		DataListener pinger = new DataListener() {
			@Override
			public void ping(DataFile dataFile) {
				dataPipes.clear();
				DataPipe pipe = new DataPipe();
				pipe.addPipe(new AbstractPipe() {
					@Override
					public String pipeValue(String value) {
						Double v = Double.parseDouble(value);
						return "" + (-Math.log10(v));
					}
				});
				dataPipes.put(data.getLinkedColumnName(PVAL_LINKER), pipe);
				ManhattanPlot.this.manPan.paintAgain();
			}
		};
		data.pingWhenReady(pinger);
		data.loadData(selTrans, false);
		this.manPan.paintAgain();
	}

	public boolean isDataLoaded() {
		return (data != null && data.isLoaded()) /* || manualData != null */;
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

			String[] loadedRest = data.getLoadedColumns();
			loadedRest = ArrayUtils.removeFromArray(loadedRest,
																							new String[] {
																														data.getLinkedColumnName(SNP_LINKER),
																														chrCol,
																														posCol,
																														pValCol
																							});
			HashMap<String, Integer> otherColInds = new HashMap<>();
			for (String s : loadedRest) {
				otherColInds.put(s, data.getIndexOfData(s));
			}

			ArrayList<String[]> dataRows = data.getAllData();
			for (int i = 0, cnt = dataRows.size(); i < cnt; i++) {
				String[] row = dataRows.get(i);
				String chr = row[chrInd];
				String pos = row[posInd];
				String pVal = row[pValInd];
				String transPVal = pVal;
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
					transPVal = dataPipes.get(pValCol).pipe(pVal);
					if (transPVal == null)
						continue;
				}
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
				double val = Double.parseDouble(transPVal);
				double val2 = Double.parseDouble(pVal);
				ManhattanDataPoint mdp = new ManhattanDataPoint(mkrInd == -1 ? chr + ":" + pos
																																		: row[mkrInd], chrI, posI,
																												linLoc, val, val2);
				for (String s : loadedRest) {
					String v = row[otherColInds.get(s)];
					mdp.otherData.put(s, dataPipes.containsKey(s) ? dataPipes.get(s).pipe(v) : v);
				}
				pts.add(mdp);
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
		public ManhattanDataPoint(String mkr, int chrI, int posI, int linLoc, double val, double valOrig) {
			this.mkr = mkr;
			this.chr = chrI;
			this.pos = posI;
			this.linearLoc = linLoc;
			this.transformedPVal = val;
			this.originalPVal = valOrig;
			this.otherData = new HashMap<>();
		}

		String mkr;
		int chr, pos;
		int linearLoc;
		double transformedPVal;
		double originalPVal;
		HashMap<String, String> otherData;
	}

	public static void main(String[] args) {
		Project proj = null;
		// proj = new Project("D:/projects/ny_choanal_v12.properties", false);
		ManhattanPlot mp = new ManhattanPlot(proj);
		mp.setVisible(true);
		// mp.loadFile();// "D:/data/ny_choanal/omni2.5v1.2/plink/plink.assoc");
	}

}
