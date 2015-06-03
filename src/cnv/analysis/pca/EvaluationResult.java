package cnv.analysis.pca;

import java.util.ArrayList;

import stats.ICC;

class EvaluationResult {
	private String title;
	private ArrayList<ICC> iccs;
	private ArrayList<String> iccTitles;
	private ArrayList<double[]> pearsonCorrels;
	private ArrayList<double[]> spearmanCorrel;
	private ArrayList<String> correlTitles;

	public EvaluationResult(String title) {
		super();
		this.title = title;
		this.iccs = new ArrayList<ICC>();
		this.iccTitles = new ArrayList<String>();
		this.pearsonCorrels = new ArrayList<double[]>();
		this.spearmanCorrel = new ArrayList<double[]>();
		this.correlTitles = new ArrayList<String>();
	}

	public String[] getHeader() {
		ArrayList<String> tmp = new ArrayList<String>();
		tmp.add("Evaluated");

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

	public String[] getData() {
		ArrayList<String> tmp = new ArrayList<String>();
		tmp.add(title);

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

}