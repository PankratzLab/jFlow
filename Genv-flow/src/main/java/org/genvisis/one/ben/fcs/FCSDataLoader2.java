package org.genvisis.one.ben.fcs;

import java.io.File;
import java.io.IOException;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.List;

import org.genvisis.common.ArrayUtils;
import org.genvisis.jfcs.FCSKeywords;
import org.genvisis.jfcs.FCSReader;
import org.genvisis.one.ben.fcs.AbstractPanel2.AXIS_SCALE;
import org.genvisis.one.ben.fcs.AbstractPanel2.AxisTransform;

public class FCSDataLoader2 {

	private static final String COMPENSATED_PREPEND = "Comp-";
	private static final int COMP_LEN = COMPENSATED_PREPEND.length();
	static final String GATING_KEY = "GATING";

	public static enum LOAD_STATE {
		LOADED, LOADING, PARTIALLY_LOADED, LOADING_REMAINDER, UNLOADED;
	}

	private volatile LOAD_STATE state = LOAD_STATE.UNLOADED;

	ArrayList<String> paramNamesInOrder;
	ArrayList<String> paramShortNamesInOrder;
	ArrayList<String> paramLongNamesInOrder;
	private HashMap<String, AXIS_SCALE> paramScales;
	public HashMap<String, AxisTransform> paramTransforms;
	private final ArrayList<Integer> ranges;
	LinkedHashSet<String> compensatedNames;
	HashMap<String, Integer> compensatedIndices;
	int eventCount = -1;
	int loadedCount = 0;
	String[] presetGating;
	private long startLoadTime;
	private String loadedFile = null;
	int paramsCount = -1;
	Date lastModified;
	Date runDate;
	Date beginTime;
	Date endTime;
	double[] paramScaling;
	private FCSReader reader;

	public enum DATA_SET {
		ALL, COMPENSATED, UNCOMPENSATED;
	}

	public FCSDataLoader2() {
		paramNamesInOrder = new ArrayList<String>();
		paramShortNamesInOrder = new ArrayList<String>();
		paramLongNamesInOrder = new ArrayList<String>();
		paramScales = new HashMap<>();
		paramTransforms = new HashMap<String, AxisTransform>();
		ranges = new ArrayList<Integer>();
		compensatedNames = new LinkedHashSet<String>();
		compensatedIndices = new HashMap<String, Integer>();
	}

	public void emptyAndReset() {
		setState(LOAD_STATE.UNLOADED);
		paramNamesInOrder = new ArrayList<String>();
		paramShortNamesInOrder = new ArrayList<String>();
		paramLongNamesInOrder = new ArrayList<String>();
		paramScales = new HashMap<>();
		paramTransforms = new HashMap<String, AxisTransform>();
		compensatedNames = new LinkedHashSet<String>();
		compensatedIndices = new HashMap<String, Integer>();
		eventCount = -1;
		loadedCount = 0;
		presetGating = null;
		loadedFile = null;
		paramsCount = -1;
		paramScaling = null;
		runDate = null;
		beginTime = null;
		endTime = null;
		reader.dispose();
		reader = null;
		System.gc();
	}

	public String getLoadedFile() {
		return loadedFile;
	}

	private synchronized void setState(LOAD_STATE state) {
		this.state = state;
	}

	public synchronized LOAD_STATE getLoadState() {
		return state;
	}

	public int[] getLoadedStatus() {
		if (getLoadState() == LOAD_STATE.LOADED) {
			return new int[] {1, 1}; // complete
		} else if (eventCount >= 0) {
			return new int[] {loadedCount, eventCount}; // in progress
		} else {
			return new int[] {-1, -1}; // indeterminate
		}
	}

	public int getCount() {
		LOAD_STATE currState = getLoadState();
		if (currState == LOAD_STATE.LOADED) {
			return eventCount;
		} else if (currState == LOAD_STATE.UNLOADED || currState == LOAD_STATE.LOADING) {
			return 0;
		} else {
			return loadedCount;
		}
	}

	public Date getLastModified() {
		return lastModified;
	}

	public Date getRunDate() {
		return runDate;
	}

	public Date getBeginTime() {
		return beginTime;
	}

	public Date getEndTime() {
		return endTime;
	}

	private void loadTimeData(FCSKeywords keys) {
		String dKey = "$DATE";
		String bKey = "$BTIM";
		String eKey = "$ETIM";

		String dVal = null, eVal = null, bVal = null;
		try {
			dVal = keys.getKeyword(dKey);
		} catch (Exception e) {
		}
		try {
			bVal = keys.getKeyword(bKey);
		} catch (Exception e) {
		}
		try {
			eVal = keys.getKeyword(eKey);
		} catch (Exception e) {
		}

		SimpleDateFormat sdfDate = new SimpleDateFormat("dd-MMM-yyyy");
		SimpleDateFormat sdfTime = new SimpleDateFormat("HH:mm:ss");
		if (dVal != null) {
			try {
				runDate = sdfDate.parse(dVal);
			} catch (ParseException e1) {
				System.out.println("Unable to parse date: " + dVal);
			}
		}
		if (bVal != null) {
			try {
				beginTime = sdfTime.parse(bVal);
			} catch (ParseException e1) {
				System.out.println("Unable to parse time: " + bVal);
			}
		}
		if (eVal != null) {
			try {
				endTime = sdfTime.parse(eVal);
			} catch (ParseException e1) {
				System.out.println("Unable to parse time: " + eVal);
			}
		}
	}

	public void loadData(String fcsFilename) throws IOException {
		// synchronized(this) {
		if (getLoadState() != LOAD_STATE.UNLOADED) {
			return;
		}
		setState(LOAD_STATE.LOADING);
		// }
		loadedFile = fcsFilename;

		startLoadTime = System.nanoTime();
		FCSReader reader = FCSReader.open(fcsFilename);

		FCSKeywords keys = reader.getKeywords();
		eventCount = keys.getEventCount();
		paramsCount = keys.getParameterCount();

		loadTimeData(keys);

		lastModified = null;// keys.getLastModified();
		if (lastModified == null) {
			lastModified = new Date(new File(fcsFilename).lastModified());
		}

		String gating = null;
		try {
			gating = keys.getKeyword(FCSDataLoader2.GATING_KEY);
		} catch (Exception e) {
			// System.err.println("Info - no precreated gating info available.");
		}
		if (gating != null) {
			String[] sp = gating.split(",");
			presetGating = sp;
		} else {
			presetGating = null;
		}

		HashMap<String, String> names = new HashMap<String, String>();
		for (int i = 0; i < paramsCount; i++) {
			String name = null;
			String shortName = null;
			try {
				name = keys.getParameterLongName(i);
			} catch (Exception e) {
			}
			try {
				shortName = keys.getParameterShortName(i);
			} catch (Exception e) {
			}
			// if (name == null) {
			// // TODO Error, no name set (included in spillover matrix?) -- will this scenario ever
			// happen?
			// name = "P" + i;
			// }
			String actName = shortName;
			if (name != null) {
				actName = shortName + " (" + name + ")";
			}
			paramNamesInOrder.add(actName);
			paramShortNamesInOrder.add(shortName);
			if (name != null) {
				paramLongNamesInOrder.add(name);
			}
			names.put(shortName, actName);

			String axisKeywork = "P" + (i + 1) + "DISPLAY";
			AXIS_SCALE scale = AXIS_SCALE.LIN;
			try {
				scale = AXIS_SCALE.valueOf(keys.getKeyword(axisKeywork));
			} catch (Exception e) {
				// System.err.println("Warning - no axis scale set for parameter " +
				// paramNamesInOrder.get(i)
				// + "; assuming a linear scale.");
			} ;
			if (scale == AXIS_SCALE.LOG) {
				scale = AXIS_SCALE.BIEX;
			}
			paramScales.put(actName, scale);
			paramTransforms.put(actName, getDefaultTransform(scale));

			String rangeKeyword = "$P" + (i + 1) + "R";
			int rng = -1;
			try {
				rng = keys.getKeywordInt(rangeKeyword);
			} catch (Exception e) {
				System.err.println("Warning - no parameter range value for parameter "
													 + paramNamesInOrder.get(i) + "; assuming standard of 262144");
				rng = 262144;
			}
			ranges.add(rng);
		}

		String[] arr = reader.getSpillover().getParameterNames();
		for (int i = 0, count = arr.length; i < count; i++) {
			compensatedNames.add(names.get(arr[i]));
			compensatedIndices.put(names.get(arr[i]), i);
			paramScales.put(COMPENSATED_PREPEND + names.get(arr[i]), getScaleForParam(names.get(arr[i])));
			paramTransforms.put(COMPENSATED_PREPEND + names.get(arr[i]),
													getDefaultTransform(getScaleForParam(names.get(arr[i]))));
		}

		this.reader = reader;
		setState(LOAD_STATE.LOADED);
	}

	private AxisTransform getDefaultTransform(AXIS_SCALE scale) {
		switch (scale) {
			case LOG:
				System.err.println("Error - NO DEFINED AXIS TRANSFORM FOR LOG AXES!");
			case BIEX:
				return AxisTransform.createBiexTransform();
			case LIN:
			default:
				return AxisTransform.createLinearTransform(0, 262144, 1);
		}
	}

	public AxisTransform getParamTransform(String param) {
		return paramTransforms.get(getInternalParamName(param));
	}

	public void setTransformMap(HashMap<String, AxisTransform> map) {
		HashMap<String, AxisTransform> newMap = new HashMap<>();
		for (String key : map.keySet()) {
			String nm = getInternalParamName(key);
			newMap.put(nm, map.get(key));
		}
		for (String key : this.paramTransforms.keySet()) {
			if (newMap.containsKey(key)) {
				continue;
			}
			newMap.put(key, this.paramTransforms.get(key));
		}
		this.paramTransforms = newMap;
	}

	public double[] getData(String colName, boolean waitIfNecessary) {
		boolean compensated = false;
		String columnName = (compensated = colName.startsWith(COMPENSATED_PREPEND))
																																								? COMPENSATED_PREPEND
																																									+ getInternalParamName(colName.substring(COMP_LEN))
																																								: getInternalParamName(colName);
		LOAD_STATE currState = getLoadState();
		double[] data;
		if (currState == LOAD_STATE.LOADED) {
			data = reader.getParamAsDoubles(columnName, compensated);
		} else {
			if (currState != LOAD_STATE.UNLOADED && waitIfNecessary) {
				while ((currState = getLoadState()) != LOAD_STATE.LOADED) {
					Thread.yield();
				}
				data = reader.getParamAsDoubles(columnName, compensated);
			} else {
				int len = eventCount == -1 ? 0 : eventCount;
				return ArrayUtils.doubleArray(len, Double.NaN);
			}
		}
		return data;
	}

	public String getInternalParamName(String name) {
		String nm = name;
		boolean prepend = nm.startsWith(COMPENSATED_PREPEND);
		if (prepend) {
			nm = nm.substring(COMP_LEN);
		}
		int ind = paramNamesInOrder.indexOf(nm);
		if (ind != -1) {
			return prepend ? COMPENSATED_PREPEND + nm : nm;
		} else if ((ind = paramShortNamesInOrder.indexOf(nm)) != -1) {
			return prepend ? COMPENSATED_PREPEND + paramNamesInOrder.get(ind)
										 : paramNamesInOrder.get(ind);
		} else if ((ind = paramLongNamesInOrder.indexOf(nm)) != -1) {
			return prepend ? COMPENSATED_PREPEND + paramNamesInOrder.get(ind)
										 : paramNamesInOrder.get(ind);
		}
		return prepend ? COMPENSATED_PREPEND + nm : nm;
	}

	public String getPresetGateAssignment(int line) {
		if (presetGating == null || presetGating.length <= line || line < 0) {
			return "";
		}
		return presetGating[line];
	}

	public double[] getDataLine(List<String> params, int ind) {
		double[] line = new double[params.size()];
		double[] data = reader.getEventAsDoubles(ind, false);
		double[] compData = reader.getEventAsDoubles(ind, true);

		LOAD_STATE currState = getLoadState();
		if (currState != LOAD_STATE.UNLOADED) {
			while ((currState = getLoadState()) != LOAD_STATE.LOADED) {
				Thread.yield();
			}
		}
		for (int i = 0; i < params.size(); i++) {
			String columnName = getInternalParamName(params.get(i));
			if (columnName.startsWith(COMPENSATED_PREPEND)) {
				line[i] = compData[compensatedIndices.get(columnName.substring(COMP_LEN))];
			} else {
				line[i] = data[paramNamesInOrder.indexOf(columnName)];
			}
		}
		return line;
	}

	public boolean containsParam(String paramName) {
		String nm = paramName;
		if (paramName.startsWith(COMPENSATED_PREPEND)) {
			nm = paramName.substring(COMP_LEN);
		}
		return (compensatedNames.contains(nm) || paramNamesInOrder.contains(nm))
					 || paramShortNamesInOrder.contains(nm) || paramLongNamesInOrder.contains(nm);
	}

	public ArrayList<String> getAllDisplayableNames(DATA_SET set) {
		ArrayList<String> names = new ArrayList<String>();
		switch (set) {
			case ALL:
				names = new ArrayList<String>(paramNamesInOrder);
				for (int i = 0; i < names.size(); i++) {
					if (compensatedNames.contains(names.get(i))) {
						names.add(COMPENSATED_PREPEND + names.get(i));
					}
				}
				break;
			case COMPENSATED:
				names = new ArrayList<String>(paramNamesInOrder);
				for (int i = 0; i < names.size(); i++) {
					if (compensatedNames.contains(names.get(i))) {
						names.set(i, COMPENSATED_PREPEND + names.get(i));
					}
				}
				break;
			case UNCOMPENSATED:
				names = new ArrayList<String>(paramNamesInOrder);
				break;
			default:
				break;
		}
		return names;
	}

	public void exportData(String outputFilename, DATA_SET set) {}

	public AXIS_SCALE getScaleForParam(String string) {
		AXIS_SCALE scale = paramScales.get(getInternalParamName(string));
		return scale == null ? AXIS_SCALE.LIN : scale;
	}

	public void setScaleForParam(String dataName, AXIS_SCALE scale) {
		paramScales.put(getInternalParamName(dataName), scale);
		paramTransforms.put(getInternalParamName(dataName), getDefaultTransform(scale));
	}

}