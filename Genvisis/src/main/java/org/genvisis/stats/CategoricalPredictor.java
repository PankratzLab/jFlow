package org.genvisis.stats;

import java.util.ArrayList;

import org.genvisis.common.Array;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

/**
 * @author lane0212 Pilfering from http://www.psychstat.missouristate.edu/multibook/mlt08m.html Class to do the processing of a categorical variable for use in regression analysis
 */
public class CategoricalPredictor {
	private String[] catVariables;
	private String[] masks;
	private Logger log;

	/**
	 * @param catVariables
	 *            we use a string for generalization
	 * @param log
	 */
	public CategoricalPredictor(String[] catVariables, String[] masks, Logger log) {
		super();
		this.catVariables = catVariables;
		this.masks = masks;
		this.log = log;
	}

	public Logger getLog() {
		return log;
	}

	public DummyCoding createDummyCoding(boolean debug) {
		String[] uniq = Array.unique(catVariables);
		ArrayList<String> finalUniq = new ArrayList<String>();
		for (int i = 0; i < uniq.length; i++) {
			if (ext.indexOfStr(uniq[i], masks) < 0) {
				finalUniq.add(uniq[i]);
			} else {
				log.reportTimeWarning("Category " + uniq[i] + " is being masked");
			}
		}
		uniq = Array.toStringArray(finalUniq);
		log.reportTimeInfo(uniq.length + " categories");
		double[][] dummyData = new double[uniq.length - 1][catVariables.length];
		boolean[] dummyBoolean = Array.booleanArray(catVariables.length, true);
		for (int i = 0; i < catVariables.length; i++) {
			if (ext.indexOfStr(catVariables[i], masks) < 0) {
				int index = ext.indexOfStr(catVariables[i], uniq);
				if (index < uniq.length - 1) {
					dummyData[index][i] = 1;
				}
			} else {
				dummyBoolean[i] = false;
			}

		}
		if (debug) {
			for (int i = 0; i < dummyBoolean.length; i++) {
				String debugS = catVariables[i] + "\tUSE = " + dummyBoolean[i];
				for (int j2 = 0; j2 < dummyData.length; j2++) {
					debugS += "\t" + dummyData[j2][i];
				}

				log.reportTimeInfo(debugS);
				try {
					Thread.sleep(10);
				} catch (InterruptedException ie) {
				}
			}
		}
		return new DummyCoding(dummyData, dummyBoolean, Array.subArray(uniq, 0, uniq.length - 1));
	}

	public static class DummyCoding {
		private double[][] dummyData;// data recoded
		private boolean[] dummyBoolean;// data represented
		private String[] titles;

		public DummyCoding(double[][] dummyData, boolean[] dummyBoolean, String[] titles) {
			super();
			this.dummyData = dummyData;
			this.dummyBoolean = dummyBoolean;
			this.titles = titles;
		}

		public double[][] getDummyData() {
			return dummyData;
		}

		public boolean[] getDummyBoolean() {
			return dummyBoolean;
		}

		public String[] getTitles() {
			return titles;
		}

	}
}
