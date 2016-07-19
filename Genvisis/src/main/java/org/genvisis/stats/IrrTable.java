package org.genvisis.stats;

import java.util.HashSet;
import java.util.Hashtable;

import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.qc.SampleQC;
import org.genvisis.common.Array;
import org.genvisis.common.Logger;

/**
 * Forms the basis of several Inter-rater reliability stats (Percent Agreement,Cohen�s Kappa,Fleiss' kappa, etc) <br>
 * Currently supports: <br>
 * Cohen's kappa (http://en.wikipedia.org/wiki/Cohen%27s_kappa) <br>
 * Percent Agreement<br>
 */
public class IrrTable {
	private static final String JUDGE = "judge";
	private static final String SUBJECT = "subject";

	private int[][] ratings;// judge,ratings of judge
	private Hashtable<Integer, Integer>[] judgeTotals;// stores the number of a particular rating given by a judge
	private HashSet<Integer> uniqRatings;
	private Hashtable<String, String> rated;// stores subjects judged already
	private int totalRated;
	private int[] judgedAgreementsBySample;// sample,agreement,1=yes,2=no
	private int[] judgedAgreementsByCategory;// ratingType,agreement, agreement total

	private Hashtable<Integer, Integer> agreementIndex;
	private double[] byChance;
	private boolean verbose;
	private Logger log;

	public IrrTable(int numJudges, int numSubjects, boolean verbose, Logger log) {
		this.log = log;
		this.verbose = verbose;
		this.ratings = new int[numJudges][numSubjects];
		this.rated = new Hashtable<String, String>();
		this.judgeTotals = initJudgeTotals(numJudges);
		this.uniqRatings = new HashSet<Integer>();
		this.totalRated = 0;
		this.judgedAgreementsBySample = null;
		this.judgedAgreementsByCategory = null;
		this.byChance = null;

	}

	public int[] getUniqRatings() {
		int[] uniqs = new int[uniqRatings.size()];
		int index = 0;
		for (int cat : uniqRatings) {
			uniqs[index] = cat;
			index++;
		}
		return uniqs;
	}

	/**
	 * @param ratingCategory
	 *            the actual category, not its index
	 * @return
	 */
	public double getPercentAgreementFor(int ratingCategory) {
		if (!agreementIndex.containsKey(ratingCategory)) {
			log.reportTimeError("Could not find rating category " + ratingCategory + " ");
			return Double.NaN;
		}
		double percentAgreement = 0;
		int totalRatings = 0;
		for (int i = 0; i < ratings.length; i++) {
			if (judgeTotals[i].containsKey(ratingCategory)) {
				totalRatings += judgeTotals[i].get(ratingCategory);
			}
		}
		double agreements = 0;
		agreements = judgedAgreementsByCategory[agreementIndex.get(ratingCategory)];
		if (totalRatings > 0) {
			percentAgreement = (double) (agreements * 2) / (double) totalRatings;
		}
		return percentAgreement;
	}

	public double getCohensKappa() {
		if (ratings.length > 2) {
			log.reportTimeWarning("Generally cohen's kappa is only used on two raters, currently computing with " + ratings.length + " raters ");
		}
		double chance = Array.sum(byChance);
		// System.out.println("CHANCE"+chance);

		double numerator = (double) Array.sum(judgedAgreementsBySample) / ratings[0].length;

		if (verbose) {
			log.reportTimeInfo("Pr(a) -agreement by judges: " + numerator);
			log.reportTimeInfo("Pr(e) -agreement by chance: " + chance);
		}

		// System.out.println("CHANCE"+chance+" TOP "+numerator);
		numerator -= chance;
		double denominator = (double) 1 - chance;
		return numerator / denominator;
	}

	/**
	 * Determines the agreement among the judges, tracks by sample and by category
	 */
	public boolean parseAgreement() {
		this.agreementIndex = initIndices(uniqRatings);
		log.reportTimeInfo("Finding agreement among " + ratings.length + " judges over " + uniqRatings.size() + " categories ");
		boolean parsed = true;
		judgedAgreementsBySample = new int[ratings[0].length];
		judgedAgreementsByCategory = new int[uniqRatings.size()];
		int shouldHave = ratings.length * ratings[0].length;

		if (totalRated != shouldHave) {
			log.reportTimeError("Found " + totalRated + " total ratings and was expecting " + shouldHave + " total ratings");
			parsed = false;
		} else {
			for (int i = 0; i < ratings[0].length; i++) {// for each sample
				boolean agree = true;
				int current = ratings[0][i];
				for (int j = 0; j < ratings.length; j++) {// check if all judges agree
					if (ratings[j][i] != current) {
						agree = false;
					}
				}
				if (agree) {
					if (!agreementIndex.containsKey(current)) {
						log.reportTimeError("Could not find rating " + current + " in data set");
					}
					judgedAgreementsByCategory[agreementIndex.get(current)]++;
					judgedAgreementsBySample[i] = 1;
				} else {
					judgedAgreementsBySample[i] = 0;
				}
			}
		}
		log.reportTimeInfo(ratings.length + " judges agreed on " + Array.sum(judgedAgreementsBySample) + " of " + judgedAgreementsBySample.length + " subjects");
		if (parsed) {
			parsed = parseByChance();
		}
		return parsed;
	}

	/**
	 * @param judge
	 * @param judgedRatings
	 *            must be the same length as initilized values for the number of samples
	 * @return if the value was added
	 */
	public boolean addRatings(int judge, int[] judgedRatings) {
		boolean added = true;
		if (judgedRatings.length != ratings[judge].length) {
			log.reportTimeError("Mismatched array size for judge " + judge + ", trying to add " + judgedRatings.length + " judgments and should have " + ratings[judge].length);
			added = false;
		} else {
			for (int i = 0; i < judgedRatings.length; i++) {
				added = addRating(judge, judgedRatings[i], i);
			}
		}
		return added;
	}

	/**
	 * Developes the probability model used in the kappa statistic
	 */
	private boolean parseByChance() {
		boolean parsed = true;
		int[] uniqs = getUniqRatings();
		int numJudges = ratings.length;
		int numSubjects = ratings[0].length;
		byChance = new double[uniqs.length];
		for (int i = 0; i < uniqs.length; i++) {
			double chance = 0;
			if (judgeTotals[0].containsKey(uniqs[i])) {
				chance = (double) judgeTotals[0].get(uniqs[i]) / numSubjects;
			} else if (verbose) {
				log.reportTimeWarning("judge 0 did not make any judgments for category " + uniqs[i]);
			}
			for (int j = 1; j < numJudges; j++) {
				if (judgeTotals[j].containsKey(uniqs[i])) {
					chance *= (double) judgeTotals[j].get(uniqs[i]) / numSubjects;
				} else {
					if (verbose) {
						log.reportTimeWarning("judge " + j + " did not make any judgments for category " + uniqs[i]);
					}
					chance = 0;
				}
			}
			byChance[i] = chance;
		}

		return parsed;
	}

	/**
	 * A given judge is allowed a single rating per subject
	 */
	private boolean addRating(int judge, int judgedRating, int subject) {
		boolean added = true;
		if (rated.containsKey(getJudgeSubjectKey(judge, subject))) {
			log.reportTimeError(JUDGE + " " + judge + " has already rated " + SUBJECT + " " + subject + ", skipping");
			added = false;
		} else {
			rated.put(getJudgeSubjectKey(judge, subject), getJudgeSubjectKey(judge, subject));
			ratings[judge][subject] = judgedRating;
			totalRated++;
			uniqRatings.add(judgedRating);
			if (!judgeTotals[judge].containsKey(judgedRating)) {
				judgeTotals[judge].put(judgedRating, 1);
			} else {
				int curVal = judgeTotals[judge].get(judgedRating);
				judgeTotals[judge].put(judgedRating, curVal + 1);
			}
		}
		return added;
	}

	private static String getJudgeSubjectKey(int judge, int subject) {
		return JUDGE + "_" + judge + ":" + SUBJECT + "_" + subject;
	}

	@SuppressWarnings("unchecked")
	private static Hashtable<Integer, Integer>[] initJudgeTotals(int numJudges) {
		Hashtable<Integer, Integer>[] judgeTotals = new Hashtable[numJudges];
		for (int i = 0; i < judgeTotals.length; i++) {
			judgeTotals[i] = new Hashtable<Integer, Integer>();
		}
		return judgeTotals;
	}

	private static Hashtable<Integer, Integer> initIndices(HashSet<Integer> uniq) {
		Hashtable<Integer, Integer> indices = new Hashtable<Integer, Integer>();
		int index = 0;
		for (int cat : uniq) {
			indices.put(cat, index);
			index++;
		}
		return indices;

	}

	public static void test() {
		Project proj = new Project(null, false);
		SampleQC sampleQC = SampleQC.loadSampleQC(proj);
		Quantiles[] quantiles = Quantiles.qetQuantilesFor(100, sampleQC.getQcMatrix(), sampleQC.getQctitles(), proj.getLog());
		IrrTable rIrrTable = new IrrTable(2, proj.getSamples().length, true, proj.getLog());
		rIrrTable.addRatings(0, quantiles[1].getQuantileMembershipAsRoundedInt());
		rIrrTable.addRatings(1, quantiles[1].getQuantileMembershipAsRoundedInt());
		rIrrTable.parseAgreement();
		System.out.println(rIrrTable.getPercentAgreementFor(10));
		System.out.println(rIrrTable.getCohensKappa());

	}

	public static void main(String[] args) {
		test();
	}

}
