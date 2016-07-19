package org.genvisis.db;

import org.genvisis.common.Array;
import org.genvisis.common.Elision;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class Filter {
	public static final String[] OPERATORS = {">=", ">", "<=", "<", "!=", "="};
	
	private String filter;
	private String label;
	private String[] variableNames;
	private String[] operators;
	private String[] thresholds;
	private String[][][] recodes;
	private double[] dThresholds;
	
	private int[] variableIndices;

	public Filter(String filter, String label) throws Elision {
		String[] elements, recodeElements;
		int operator, recodeIndex;
		
		this.filter = filter;
		this.label = label;
		
		// currently only uses AND; user might also want to use OR
		elements = filter.split("&");
		
		variableNames = new String[elements.length];
		operators = new String[elements.length];
		thresholds = new String[elements.length];
		dThresholds = new double[elements.length];
		recodes = new String[elements.length][][];
		
		try {
			for (int i = 0; i < elements.length; i++) {
				if (elements[i].contains(";")) {
					recodeElements = elements[i].substring(elements[i].indexOf(";")+1).split(";");
					elements[i] = elements[i].substring(0, elements[i].indexOf(";"));
					recodes[i] = new String[recodeElements.length][2];
					for (int j = 0; j < recodeElements.length; j++) {
						recodeIndex = recodeElements[j].indexOf("->");
						if (recodeIndex == -1) {
							throw new Elision("Recode needs to be split up by a \"->\"  as in NA->0 (currently listed as "+recodeElements[j]+")");
						}
						recodes[i][j][0] = recodeElements[j].substring(0, recodeIndex);
						recodes[i][j][1] = recodeElements[j].substring(recodeIndex+2);
					}
				} else {
					recodes[i] = null;
				}

				operator = -1;
				for (int j = 0; operator == -1 && j < OPERATORS.length; j++) {
					if (elements[i].indexOf(OPERATORS[j]) >= 0) {
						operator = j;
					}
				}
				if (operator == -1) {
					throw new Elision("Element '"+elements[i]+"' does not contain a valid operator ( "+ext.listWithCommas(OPERATORS, true)+" )");
				}
				
				variableNames[i] = elements[i].substring(0, elements[i].indexOf(OPERATORS[operator]));
				operators[i] = OPERATORS[operator];
				thresholds[i] = elements[i].substring(elements[i].indexOf(OPERATORS[operator])+OPERATORS[operator].length());
				try {
					dThresholds[i] = Double.parseDouble(thresholds[i]);
				} catch (NumberFormatException nfe) {
					if (ext.indexOfStr(operators[i], OPERATORS) < 4) {
							throw new Elision("Threshold "+thresholds[i]+" is not a valid number and cannot be used with "+operators[i]+" as listed in "+filter);
					}
					dThresholds[i] = Double.MAX_VALUE;
				}
			}
		
		} catch (Exception e) {
			System.err.println("Error - problem with: "+filter);
			e.printStackTrace();
		}
	}
	
	public String getLabel(String[] line) {
		String str;

		str = "";
		for (int i=0; i<variableNames.length; i++) {
			str += (i==0?"":", ")+(line[variableIndices[i]].length()>5?line[variableIndices[i]].substring(0,5):line[variableIndices[i]]);
		}

		if (line != null) {
			return label+" ("+str+")";
		} else {
			return label;
		}
	}

	public String getFailedValues(String[] line, Logger log) {
		String str;

		str = "";
		for (int i=0; i<variableNames.length; i++) {
			str += (i==0?"":" ")+(line[variableIndices[i]].length()>5?line[variableIndices[i]].substring(0,5):line[variableIndices[i]]);
		}
		
		return str;
	}

	public void determineIndices(String[] header) throws Elision {
		try {
			variableIndices = ext.indexFactors(variableNames, header, false, false);
			if (Array.min(variableIndices) == -1) {
				throw new Elision("Filter variable name was not found in the file to be filtered");
			}
		} catch (Exception e) {
			System.err.println("Error - problem with: "+filter);
			e.printStackTrace();
		}
	}
	
	public boolean meetsCriteria(String[] line, Logger log) {
		boolean meetsAllCriteria;
		double value;
		
		meetsAllCriteria = true;
		
		for (int i=0; i<variableNames.length; i++) {
			if (dThresholds[i] == Double.MAX_VALUE) {
				if (operators[i].equals("!=")) {
					if(line[variableIndices[i]].equals(thresholds[i])) {
						meetsAllCriteria = false;
					}
				} else if (operators[i].equals("=")) {
					if(!line[variableIndices[i]].equals(thresholds[i])) {
						meetsAllCriteria = false;
					}
				} else {
					log.reportError("This error should not be able to be reached if caught in the constructor");
				}
			} else {
				value = Double.MAX_VALUE;
				for (int j = 0; recodes[i] != null && j < recodes[i].length; j++) {
					if (line[variableIndices[i]].equals(recodes[i][j][0])) {
						value = Double.parseDouble(recodes[i][j][1]);
					}
				}
				if (ext.isMissingValue(line[variableIndices[i]]) && value == Double.MAX_VALUE) {
					log.reportError("Error - a missing value ("+line[variableIndices[i]]+") in column "+(variableIndices[i]+1)+" was not properly recoded; could cause a filter to malfunction");
					return false;
				}
				
				if (value == Double.MAX_VALUE) {
					try {
						value = Double.parseDouble(line[variableIndices[i]]);
					} catch (NumberFormatException nfe) {
						log.reportError("Error - could not filter on value '"+line[variableIndices[i]]+"'");
						return false;
					}
				}

				if (Double.isNaN(value)) {
				} else if (operators[i].equals(">=")) {
					if(value < dThresholds[i]) {
						meetsAllCriteria = false;
					}
				} else if (operators[i].equals(">")) {
					if(value <= dThresholds[i]) {
						meetsAllCriteria = false;
					}
				} else if (operators[i].equals("<=")) {
					if(value > dThresholds[i]) {
						meetsAllCriteria = false;
					}
				} else if (operators[i].equals("<")) {
					if(value >= dThresholds[i]) {
						meetsAllCriteria = false;
					}
				} else if (operators[i].equals("!=")) {
					if(value == dThresholds[i]) {
						meetsAllCriteria = false;
					}
				} else if (operators[i].equals("=")) {
					if(value != dThresholds[i]) {
						meetsAllCriteria = false;
					}
				} else {
					log.reportError("This error should not be able to be reached if caught in the constructor");
				}
			}
		}
		
		return meetsAllCriteria;
	}
}
