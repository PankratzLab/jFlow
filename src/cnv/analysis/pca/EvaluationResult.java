package cnv.analysis.pca;

import java.io.Serializable;
import java.util.ArrayList;

import cnv.analysis.pca.CorrectionIterator.ITERATION_TYPE;
import cnv.analysis.pca.CorrectionIterator.MODEL_BUILDER_TYPE;
import cnv.analysis.pca.CorrectionIterator.ORDER_TYPE;
import common.Files;
import common.Logger;
import stats.ICC;

class EvaluationResult implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private String title;
	private ORDER_TYPE orType;
	private ITERATION_TYPE itType;
	private MODEL_BUILDER_TYPE bType;
	private double[] estimateData;
	private double rSquared;
	private ArrayList<ICC> iccs;
	private ArrayList<String> iccTitles;
	private ArrayList<double[]> pearsonCorrels;
	private ArrayList<double[]> spearmanCorrel;
	private ArrayList<String> correlTitles;

	public EvaluationResult(String title, double[] estimateData, double rSquared) {
		super();
		this.title = title;
		this.rSquared = rSquared;
		this.estimateData = estimateData;
		this.iccs = new ArrayList<ICC>();
		this.iccTitles = new ArrayList<String>();
		this.pearsonCorrels = new ArrayList<double[]>();
		this.spearmanCorrel = new ArrayList<double[]>();
		this.correlTitles = new ArrayList<String>();
	}

	public void setOrType(ORDER_TYPE orType) {
		this.orType = orType;
	}

	public void setbType(MODEL_BUILDER_TYPE bType) {
		this.bType = bType;
	}

	public void setItType(ITERATION_TYPE itType) {
		this.itType = itType;
	}

	public String[] getHeader() {
		ArrayList<String> tmp = new ArrayList<String>();
		tmp.add("Evaluated");
		tmp.add("IterationType");
		tmp.add("OrderType");
		tmp.add("Model Building type");

		tmp.add("Rsquare_correction");
		for (int i = 0; i < iccTitles.size(); i++) {
			tmp.add("ICC_" + iccTitles.get(i));
		}

		for (int i = 0; i < correlTitles.size(); i++) {
			tmp.add("PEARSON_CORREL_" + correlTitles.get(i));
			tmp.add("PEARSON_P_" + correlTitles.get(i));
			tmp.add("SPEARMAN_CORREL_" + correlTitles.get(i));
			tmp.add("SPEARMAN_P_" + correlTitles.get(i));
		}
		return tmp.toArray(new String[tmp.size()]);
	}

	public double[] getEstimateData() {
		return estimateData;
	}

	public void shrink() {
		for (int i = 0; i < iccs.size(); i++) {
			iccs.get(i).shrink();
		}
	}

	public String[] getData() {
		ArrayList<String> tmp = new ArrayList<String>();
		tmp.add(title);
		tmp.add(itType == null ? "NA" : itType.toString());
		tmp.add(orType == null ? "NA" : orType.toString());
		tmp.add(bType == null ? "NA" : bType.toString());

		tmp.add(rSquared + "");
		for (int i = 0; i < iccs.size(); i++) {
			tmp.add(iccs.get(i).getICC() + "");
		}
		for (int i = 0; i < pearsonCorrels.size(); i++) {
			for (int j = 0; j < pearsonCorrels.get(i).length; j++) {
				tmp.add(pearsonCorrels.get(i)[j] + "");
			}
			for (int j = 0; j < spearmanCorrel.get(i).length; j++) {
				tmp.add(spearmanCorrel.get(i)[j] + "");
			}
		}
		return tmp.toArray(new String[tmp.size()]);
	}

	public String getTitle() {
		return title;
	}

	public ArrayList<ICC> getIccs() {
		return iccs;
	}

	public ArrayList<String> getIccTitles() {
		return iccTitles;
	}

	public ArrayList<double[]> getPearsonCorrels() {
		return pearsonCorrels;
	}

	public ArrayList<double[]> getSpearmanCorrel() {
		return spearmanCorrel;
	}

	public ArrayList<String> getCorrelTitles() {
		return correlTitles;
	}

	public static void serialize(EvaluationResult[] results, String fileName) {
		Files.writeSerial(results, fileName, true);
	}

	public static EvaluationResult[] readSerial(String fileName, Logger log) {
		return (EvaluationResult[]) Files.readSerial(fileName, false, log, false, true);
	}

	// @Override
	// public String[] getIndexKeys() {
	// return new String[] { title };
	// }

}