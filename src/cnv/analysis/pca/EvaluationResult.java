package cnv.analysis.pca;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Hashtable;

import link.Heritability;
import cnv.analysis.pca.CorrectionIterator.ITERATION_TYPE;
import cnv.analysis.pca.CorrectionIterator.MODEL_BUILDER_TYPE;
import cnv.analysis.pca.CorrectionIterator.ORDER_TYPE;
import cnv.filesys.Project;
import common.Array;
import common.Files;
import common.HashVec;
import common.Logger;
import common.ext;
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

	public static EvalHeritabilityResult prepareHeritability(Project proj, String ped, boolean[] samplesToEvaluate, String serFile) {
		Logger log = proj.getLog();
		EvaluationResult[] evaluationResults = readSerial(serFile, log);
		log.reportTimeInfo("Loaded " + evaluationResults.length + " evaluation results");
		String db = ext.rootOf(serFile, false) + ".heritability.dat";
		String crf = ext.rootOf(serFile, false) + ".heritability.crf";
//		System.out.println(db);
//		System.out.println(crf);
//		System.out.println(ext.rootOf(crf, false)+"_summary.xln");
//		System.exit(1);
		generateHeritabilityDb(proj, evaluationResults, samplesToEvaluate, db, ped, crf, log);
		Heritability.fromParameters(crf, log);
		EvalHeritabilityResult evalHeritabilityResult = new EvalHeritabilityResult(ped, db, crf);
		return evalHeritabilityResult;
	}

	public static class EvalHeritabilityResult {
		private String ped;
		private String db;
		private String crf;

		public EvalHeritabilityResult(String ped, String db, String crf) {
			super();
			this.ped = ped;
			this.db = db;
			this.crf = crf;
		}

		public String getPed() {
			return ped;
		}

		public String getDb() {
			return db;
		}

		public String getCrf() {
			return crf;
		}
	}

	private static void generateHeritabilityDb(Project proj, EvaluationResult[] results, boolean[] samplesToEvaluate, String output, String ped, String crf, Logger log) {
		log.reportTimeWarning("Assuming stored estimate results are in project order to create heritability db " + output);
		log.reportTimeWarning("Assuming ped file has DNA listed in the last column of  " + output);
		if (samplesToEvaluate != null) {
			log.reportTimeInfo("Using " + Array.booleanArraySum(samplesToEvaluate) + " samples that are not excluded");
		}
		try {
			Hashtable<String, String> pedHash = HashVec.loadFileToHashString(ped, 6, new int[] { 0, 1, 2, 3, 4, 5 }, "\t", false);
			String[] samples = proj.getSamples();
			String[] titles = new String[results.length];
			PrintWriter writer = new PrintWriter(new FileWriter(output));
			writer.print("IID\tFID");
			for (int i = 0; i < results.length; i++) {
				writer.print("\t" + i);
				titles[i] = i + "";
			}
			writer.println();

			ArrayList<String> sampsNotSeen = new ArrayList<String>();
			for (int i = 0; i < samples.length; i++) {
				if (samplesToEvaluate == null || samplesToEvaluate[i]) {
					if (pedHash.containsKey(samples[i])) {
						String[] fidIid = Array.subArray(pedHash.get(samples[i]).split("\t"), 0, 2);
						writer.print(fidIid[1] + "\t" + fidIid[0]);
						for (int j = 0; j < results.length; j++) {
							writer.print("\t" + results[j].getEstimateData()[i]);
						}
						writer.println();

					} else {
						sampsNotSeen.add(samples[i]);
					}
				}
			}

			writer.close();
			if (sampsNotSeen.size() > 0) {
				String missing = ext.addToRoot(ped, ".missing");
				log.reportTimeWarning(sampsNotSeen.size() + " samples were not found in the ped file , writing to " + missing);
				Files.writeList(sampsNotSeen.toArray(new String[sampsNotSeen.size()]), missing);
			}
			Heritability.developCrf(ped, output, crf, ext.rootOf(output), titles, log);
		} catch (Exception e) {
			log.reportError("Error writing to " + output);
			log.reportException(e);
		}
	}

	// @Override
	// public String[] getIndexKeys() {
	// return new String[] { title };
	// }

}
