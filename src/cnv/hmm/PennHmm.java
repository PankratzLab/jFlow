package cnv.hmm;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;

import cnv.analysis.PennCNV;
import cnv.filesys.MarkerSet;
import cnv.filesys.Project;
import cnv.filesys.Project.ARRAY;
import cnv.filesys.Sample;
import cnv.qc.GcAdjustor;
import cnv.qc.GcAdjustor.GcModel;
import cnv.var.CNVariant;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.OpdfGaussian;
import common.Array;
import common.Files;
import common.Logger;

/**
 * @author lane0212 <br>
 *         Mimics the hmm functionality used in in PennCNV (http://penncnv.openbioinformatics.org/)<br>
 *         Really, we mimic the kext C package <br>
 *         See below for an example .hmm file
 */
// M=6
// N=6
// A:
// 0.936719716 0.006332139 0.048770575 0.000000001 0.008177573 0.000000001
// 0.000801036 0.949230924 0.048770575 0.000000001 0.001168245 0.000029225
// 0.000004595 0.000047431 0.999912387 0.000000001 0.000034971 0.000000621
// 0.000049998 0.000049998 0.000049998 0.999750015 0.000049998 0.000049998
// 0.000916738 0.001359036 0.048770575 0.000000001 0.948953653 0.000000002
// 0.000000001 0.000000001 0.027257213 0.000000001 0.000000004 0.972742785
// B:
// 0.950000 0.000001 0.050000 0.000001 0.000001 0.000001
// 0.000001 0.950000 0.050000 0.000001 0.000001 0.000001
// 0.000001 0.000001 0.999995 0.000001 0.000001 0.000001
// 0.000001 0.000001 0.050000 0.950000 0.000001 0.000001
// 0.000001 0.000001 0.050000 0.000001 0.950000 0.000001
// 0.000001 0.000001 0.050000 0.000001 0.000001 0.950000
// pi:
// 0.000001 0.000500 0.999000 0.000001 0.000500 0.000001
// B1_mean:
// -3.527211 -0.664184 0.000000 100.000000 0.395621 0.678345
// B1_sd:
// 1.329152 0.284338 0.159645 0.211396 0.209089 0.191579
// B1_uf:
// 0.010000
// B2_mean:
// 0.000000 0.250000 0.333333 0.500000 0.500000
// B2_sd:
// 0.016372 0.042099 0.045126 0.034982 0.304243
// B2_uf:
// 0.010000
// B3_mean:
// -2.051407 -0.572210 0.000000 0.000000 0.361669 0.626711
// B3_sd:
// 2.132843 0.382025 0.184001 0.200297 0.253551 0.353183
// B3_uf:
// 0.010000

public class PennHmm {
	private static final double NOT_ZERO_PI = 0.000000001; // 1e-9
	private static final double STATE_CHANGE = 100000.0;
	private static final double VITHUGE = 100000000000.0;
	private int M;
	private int N;// number of states
	private double[] pi;
	private double[][] a;
	private BStatus B1;// for LRR measure from SNP markers
	private BStatus B2;// for BAF measure from SNP markers
	private BStatus B3;// for LRR measure from CN only markers

	private Logger log;

	public PennHmm(int m, int n, double[] pi, double[][] a, BStatus b1, BStatus b2, BStatus b3, Logger log) {
		super();
		this.M = m;
		this.N = n;
		this.pi = pi;
		this.a = a;
		this.B1 = b1;
		this.B2 = b2;
		this.B3 = b3;
		this.log = log;
	}

	public PennHmm(PennHmm pennHmm) {
		this.M = pennHmm.M;
		this.N = pennHmm.N;
		this.pi = pennHmm.pi;
		this.a = pennHmm.a;
		this.B1 = pennHmm.B1;
		this.B2 = pennHmm.B2;
		this.B3 = pennHmm.B3;
		this.log = pennHmm.log;
	}

	public int getN() {
		return N;
	}

	public double[][] getA() {
		return a;
	}

	public void setA(double[][] a) {
		this.a = a;
	}

	public Logger getLog() {
		return log;
	}

	public BStatus getB1() {
		return B1;
	}

	public double[] getPi() {
		return pi;
	}

	public BStatus getB2() {
		return B2;
	}

	public BStatus getB3() {
		return B3;
	}

	/**
	 * take log of pi values<br>
	 * WARNING : modifies internal pi array;
	 */
	private void logPi() {
		double[] loggedPi = new double[pi.length];
		for (int i = 0; i < loggedPi.length; i++) {
			double pt = pi[i];
			if (pt == 0) {/* eliminate problems with zero probability */
				pt = NOT_ZERO_PI;
			}
			loggedPi[i] = Math.log(pt);
		}
		this.pi = loggedPi;
	}

	/**
	 * @param state
	 *            copynumber state
	 * @param bStatus
	 *            B1, B2...
	 * @param o
	 *            the intensity observation (LRR value)
	 * @return
	 */
	private static double b1iot(int state, BStatus bStatus, double o) {
		double p = 0;
		p = bStatus.getB_uf();
		p += (1 - bStatus.getB_uf()) * bStatus.getGaussians()[state].probability(new ObservationReal(o));

		if (p == 0) {
			p = Float.MIN_VALUE;

		}
		return Math.log(p);
	}

	private static double b2iot(int state, BStatus bStatus, double pfb, double b) {
		double p = 0;

		double uf = bStatus.getB_uf();
		p = uf;
		// double mean0 = bStatus.getB_mean()[1];//note the index starts at 1
		// double mean25 = bStatus.getB_mean()[2];
		// double mean33 = bStatus.getB_mean()[3];
		// double mean50 = bStatus.getB_mean()[4];
		// double mean50_state1 = bStatus.getB_mean()[5];
		// double sd0 = bStatus.getB_sd()[1];
		// double sd25 = bStatus.getB_sd()[2];
		// double sd33 = bStatus.getB_sd()[3];
		// double sd50 = bStatus.getB_sd()[4];
		// double sd50_state1 = bStatus.getB_sd()[5];

		// if (state == 1) {
		// if (b==0) {
		// p+= (1-uf) * cdf_normal (0, mean50_state1, sd50_state1);
		// } else if (b==1) {
		// p+= (1-uf) * cdf_normal (0, mean50_state1, sd50_state1);
		// } else {
		// p+= (1-uf) * pdf_normal (b, mean50_state1, sd50_state1);
		// }
		// }
		ObservationReal o = new ObservationReal(b);
		if (state == 0) {
			OpdfGaussian opdfGaussian = bStatus.getGaussians()[4];

			if (b == 0) {
				p += (1 - uf) * opdfGaussian.cdf(new ObservationReal(0));
			} else if (b == 1) {
				p += (1 - uf) * opdfGaussian.cdf(new ObservationReal(0));
			} else {
				p += (1 - uf) * opdfGaussian.probability(o);
			}

		} else if (state == 1) {
			if (b == 0) {
				p += (1 - uf) * (1 - pfb) / 2;
			} else if (b == 1) {
				p += (1 - uf) * pfb / 2;
			} else {
				OpdfGaussian opdfGaussian = bStatus.getGaussians()[0];
				OpdfGaussian opdfGaussianMinus = new OpdfGaussian(1 - bStatus.getB_mean()[0], Math.pow(bStatus.getB_sd()[0], 2));

				p += (1 - uf) * (1 - pfb) * opdfGaussian.probability(o);
				p += (1 - uf) * pfb * opdfGaussianMinus.probability(o);
			}
		} else if (state == 2) {
			if (b == 0) {
				p += (1 - uf) * (1 - pfb) * (1 - pfb) / 2;
			} else if (b == 1) {
				p += (1 - uf) * pfb * pfb / 2;
			} else {
				OpdfGaussian opdfGaussian = bStatus.getGaussians()[0];
				OpdfGaussian opdfGaussianMinus = new OpdfGaussian(1 - bStatus.getB_mean()[0], Math.pow(bStatus.getB_sd()[0], 2));
				OpdfGaussian opdfGaussian5 = bStatus.getGaussians()[3];

				p += (1 - uf) * (1 - pfb) * (1 - pfb) * opdfGaussian.probability(o);
				p += (1 - uf) * 2 * pfb * (1 - pfb) * opdfGaussian5.probability(o);
				p += (1 - uf) * pfb * pfb * opdfGaussianMinus.probability(o);
			}
		} else if (state == 3) {
			if (b == 0) {
				p += (1 - uf) * (1 - pfb) / 2;
			} else if (b == 1) {
				p += (1 - uf) * pfb / 2;
			} else {
				OpdfGaussian opdfGaussian = bStatus.getGaussians()[1 - 1];

				p += (1 - uf) * (1 - pfb) * opdfGaussian.probability(o);
				p += (1 - uf) * pfb * opdfGaussian.probability(o);
			}
		} else if (state == 4) {
			if (b == 0) {
				p += (1 - uf) * (1 - pfb) * (1 - pfb) * (1 - pfb) / 2;
			} else if (b == 1) {
				p += (1 - uf) * pfb * pfb * pfb / 2;
			} else {
				OpdfGaussian opdfGaussian = bStatus.getGaussians()[0];
				OpdfGaussian opdfGaussianMinus = new OpdfGaussian(1 - bStatus.getB_mean()[0], Math.pow(bStatus.getB_sd()[0], 2));
				OpdfGaussian opdfGaussian33 = bStatus.getGaussians()[2];
				OpdfGaussian opdfGaussianMinus33 = new OpdfGaussian(1 - bStatus.getB_mean()[2], Math.pow(bStatus.getB_sd()[2], 2));

				p += (1 - uf) * (1 - pfb) * (1 - pfb) * (1 - pfb) * opdfGaussian.probability(o);
				p += (1 - uf) * 3 * (1 - pfb) * (1 - pfb) * pfb * opdfGaussian33.probability(o);
				p += (1 - uf) * 3 * (1 - pfb) * pfb * pfb * opdfGaussianMinus33.probability(o);
				p += (1 - uf) * pfb * pfb * pfb * opdfGaussianMinus.probability(o);
			}
		} else if (state == 5) {
			if (b == 0) {
				p += (1 - uf) * (1 - pfb) * (1 - pfb) * (1 - pfb) * (1 - pfb) / 2;
			} else if (b == 1) {
				p += (1 - uf) * pfb * pfb * pfb * pfb / 2;
			} else {
				OpdfGaussian opdfGaussian = bStatus.getGaussians()[0];
				OpdfGaussian opdfGaussianMinus = new OpdfGaussian(1 - bStatus.getB_mean()[0], Math.pow(bStatus.getB_sd()[0], 2));
				OpdfGaussian opdfGaussian25 = bStatus.getGaussians()[1];
				OpdfGaussian opdfGaussianMinus25 = new OpdfGaussian(1 - bStatus.getB_mean()[1], Math.pow(bStatus.getB_sd()[1], 2));
				OpdfGaussian opdfGaussian5 = bStatus.getGaussians()[3];

				p += (1 - uf) * (1 - pfb) * (1 - pfb) * (1 - pfb) * (1 - pfb) * opdfGaussian.probability(o);
				p += (1 - uf) * 4 * (1 - pfb) * (1 - pfb) * (1 - pfb) * pfb * opdfGaussian25.probability(o);
				p += (1 - uf) * 6 * (1 - pfb) * (1 - pfb) * pfb * pfb * opdfGaussian5.probability(o);
				p += (1 - uf) * 4 * (1 - pfb) * pfb * pfb * pfb * opdfGaussianMinus25.probability(o);
				p += (1 - uf) * pfb * pfb * pfb * pfb * opdfGaussianMinus.probability(o);
			}
		}
		if (p == 0)
			p = Float.MIN_VALUE;
		return Math.log(p);
	}

	public static PennHmm loadPennHmm(String hmmFile, Logger log) {
		if (!Files.exists(hmmFile)) {
			log.reportFileNotFound(hmmFile);
			return null;
		} else {
			try {
				BufferedReader reader = Files.getAppropriateReader(hmmFile);
				String M = reader.readLine().trim();
				int m = -1;
				int n = -1;
				double[][] a = null;
				double[][] b = null;

				if (!M.startsWith("M=")) {
					log.reportTimeError("cannot read M annotation from HMM file");
					return null;
				} else {
					m = Integer.parseInt(M.replaceAll("M=", ""));
				}
				String N = reader.readLine().trim();
				if (!N.startsWith("N=")) {
					log.reportTimeError("cannot read N annotation from HMM file");
					return null;
				} else {
					n = Integer.parseInt(N.replaceAll("N=", ""));
				}
				String aTag = reader.readLine().trim();
				if (!aTag.startsWith("A:")) {
					log.reportTimeError("cannot read A: annotation from HMM file");
					return null;
				} else {
					a = loadMatrix(reader, m, n);
				}
				String bTag = reader.readLine().trim();
				if (!bTag.startsWith("B:")) {
					log.reportTimeError("cannot read B: annotation from HMM file");
					return null;
				} else {
					b = loadMatrix(reader, m, n);
				}
				double[] pi = loadTwoLineDoubleArray("pi", reader, log);
				BStatus B1 = loadBstatus("B1", reader, log);
				BStatus B2 = loadBstatus("B2", reader, log);
				BStatus B3 = loadBstatus("B3", reader, log);
				reader.close();

				return new PennHmm(m, n, pi, a, B1, B2, B3, log);
			} catch (FileNotFoundException fnfe) {
				log.reportError("Error: file \"" + hmmFile + "\" not found in current directory");
				return null;
			} catch (IOException ioe) {
				log.reportError("Error reading file \"" + hmmFile + "\"");
				return null;
			}
		}
	}

	public static ViterbiResult ViterbiLogNP_CHMM(PennHmm pennHmm, double[] o1, double[] o2, double[] pfb, int[] snpdist, boolean[] copyNumberOnlyDef) {
		if (o1.length != o2.length || o1.length != pfb.length || o1.length != snpdist.length || o1.length != copyNumberOnlyDef.length) {
			String error = "BUG: mismatched array lengths";
			pennHmm.getLog().reportTimeError(error);
			throw new IllegalArgumentException(error);
		}

		double[][] biot = new double[pennHmm.getN()][o1.length];
		int[] q = new int[o1.length];
		double[][] delta = new double[o1.length][pennHmm.getN()];
		int[][] psi = new int[o1.length][pennHmm.getN()];
		PennHmm pennHmmLog = new PennHmm(pennHmm);
		/* 0. Preprocessing */
		pennHmmLog.logPi();
		int cnCount = 0;
		int snpCount = 0;
		for (int i = 0; i < pennHmmLog.getN(); i++) {
			for (int t = 0; t < o1.length; t++) {
				if (copyNumberOnlyDef[t]) {
					biot[i][t] = b1iot(i, pennHmmLog.getB3(), o1[t]);
					cnCount++;
				} else {
					double bioTmp = b1iot(i, pennHmmLog.getB1(), o1[t]);
					bioTmp += b2iot(i, pennHmmLog.getB2(), pfb[t], o2[t]);
					biot[i][t] = bioTmp;
					snpCount++;
				}
			}
		}

		pennHmm.getLog().reportTimeInfo("Encountered " + snpCount + " snp probes and " + cnCount + " cn probes");

		for (int i = 0; i < pennHmmLog.getN(); i++) {
			delta[0][i] = pennHmmLog.getPi()[i] + biot[i][0];
		}
		for (int t = 1; t < o1.length; t++) {
			PennHmm converted = convertHMMTransition(pennHmmLog, snpdist[t - 1]);
			for (int j = 0; j < converted.getN(); j++) {
				double maxval = -1 * VITHUGE;
				int maxvalind = 0;
				for (int i = 0; i < converted.getN(); i++) {
					double val = delta[t - 1][i] + Math.log(converted.getA()[i][j]);
					if (val > maxval) {
						maxval = val;
						maxvalind = i;
					}
				}
				delta[t][j] = maxval + biot[j][t];
				psi[t][j] = maxvalind;
			}
		}

		double pprob = -1 * VITHUGE;
		q[o1.length - 1] = 1;
		for (int i = 0; i < pennHmmLog.getN(); i++) {
			if (delta[o1.length - 1][i] > pprob) {
				pprob = delta[o1.length - 1][i];
				q[o1.length - 1] = i;
			}
		}
		for (int t = o1.length - 1; t >= 0; t--) {
			q[t] = psi[t][q[t]];
		}
		return new ViterbiResult(q);
	}

	private static class ViterbiResult{
		private int[] q;

		public ViterbiResult(int[] q) {
			super();
			this.q = q;
		}
	
	}


	private static int[][] getSNPDist(Project proj) {
		MarkerSet markerSet = proj.getMarkerSet();
		int[][] chrPos = markerSet.getPositionsByChr();
		int[][] snpDists = new int[chrPos.length][];
		for (int i = 0; i < chrPos.length; i++) {
			int[] distsTmp = new int[chrPos[i].length];
			if (distsTmp.length > 0) {
				distsTmp[0] = chrPos[i][0];
				for (int j = 1; j < distsTmp.length - 1; j++) {
					int dist = chrPos[i][j + 1] - chrPos[i][j];
					distsTmp[j] = dist > 0 ? dist : 1;
				}
			}
			snpDists[i] = distsTmp;
		}
		return snpDists;

	}

	private static BStatus loadBstatus(String bTag, BufferedReader reader, Logger log) throws IOException {
		double[] bmean = loadTwoLineDoubleArray(bTag + "_mean:", reader, log);

		if (bmean != null) {
			double[] bSd = loadTwoLineDoubleArray(bTag + "_sd:", reader, log);
			if (bSd != null) {
				double[] b_uf = loadTwoLineDoubleArray(bTag + "_uf:", reader, log);
				if (b_uf != null) {
					BStatus bStatus = new BStatus(bmean, bSd, b_uf[0]);
					return bStatus;
				}
			}
		}

		return null;
	}

	private static double[] loadTwoLineDoubleArray(String tag, BufferedReader reader, Logger log) throws IOException {
		String atag = reader.readLine();
		if (!atag.startsWith(tag)) {
			log.reportTimeError("cannot read " + tag + " annotation from HMM file");
			return null;
		} else {
			String[] line = reader.readLine().trim().split("[\\s]+");
			double[] data = Array.toDoubleArray(line);
			return data;
		}
	}

	private static double[][] loadMatrix(BufferedReader reader, int m, int n) throws IOException {
		double[][] a;
		a = new double[n][m];
		for (int i = 0; i < n; i++) {
			String[] line = reader.readLine().trim().split("[\\s]+");
			double[] tmp = Array.toDoubleArray(line);
			a[i] = tmp;
		}
		return a;
	}

	/**
	 * Stores the B1* etc entries for the data distributions of each state
	 *
	 */
	private static class BStatus {
		private double[] b_mean;
		private double[] b_sd;
		private OpdfGaussian[] gaussians;
		private double b_uf;

		public BStatus(double[] b_mean, double[] b_sd, double b_uf) {
			super();
			this.b_mean = b_mean;
			this.b_sd = b_sd;
			this.b_uf = b_uf;
			this.gaussians = new OpdfGaussian[b_mean.length];
			for (int i = 0; i < b_mean.length; i++) {
				gaussians[i] = new OpdfGaussian(b_mean[i], Math.pow(b_sd[i], 2));
			}
		}

		public OpdfGaussian[] getGaussians() {
			return gaussians;
		}

		public double[] getB_mean() {
			return b_mean;
		}

		public double[] getB_sd() {
			return b_sd;
		}

		public double getB_uf() {
			return b_uf;
		}

	}

	public static PennHmm convertHMMTransition(PennHmm pennHmm, int dist)
	/* this subroutine convert HMM transition probabilities using the P=Pref*(1-exp(-d/D)) formula for off-diagonal cells in the matrix */
	{
		PennHmm converted = new PennHmm(pennHmm);
		int i, j;
		double D = STATE_CHANGE;
		double offdiagonal_sum = 0;
		double[][] tmpA = new double[pennHmm.getA().length][pennHmm.getA()[0].length];
		for (i = 0; i < converted.getN(); i++) {
			offdiagonal_sum = 0;
			for (j = 0; j < converted.getN(); j++) {
				if (i != j) {
					if (i == 3) {
						tmpA[i][j] = pennHmm.getA()[i][j] * (1 - Math.exp(-dist / D / 1000)) / (1 - Math.exp(-5000 / D / 1000));
					} else {
						tmpA[i][j] = pennHmm.getA()[i][j] * (1 - Math.exp(-dist / D)) / (1 - Math.exp(-5000 / D));
					}
					if (tmpA[i][j] > 1) {
						pennHmm.getLog().reportTimeWarning("WARNING: Off-diagonal cell A[%i][%i] (%f to %f by %i) in transition matrix is over boundary of 1 (HMM model is not optimized). Assign 0.999 as the value instead.\n" + i + "\t" + j + "\t" + pennHmm.getA()[i][j] + "\t" + tmpA[i][j] + "\t" + dist);
						tmpA[i][j] = 0.999; /* maximum possible off-diagonal value (since state3 frequency is 0.999) */
					}
					offdiagonal_sum += tmpA[i][j];
				}
			}
			if (offdiagonal_sum >= 1) {
				for (j = 0; j < converted.getN(); j++) {
					tmpA[i][j] /= (offdiagonal_sum / 0.999);
				}
				offdiagonal_sum = 0.999;
			}
			tmpA[i][i] = 1 - offdiagonal_sum;
		}
		converted.setA(tmpA);
		return converted;
	}

	private static boolean[] determineCNOnly(Project proj, String[] markers) {
		boolean[] copyNumberOnlyDef = new boolean[markers.length];
		ARRAY array = proj.getArrayType();
		for (int i = 0; i < markers.length; i++) {
			copyNumberOnlyDef[i] = array.isCNOnly(markers[i]);
		}
		return copyNumberOnlyDef;
	}

	public static void test(Project proj, String hmmFile) {
		PennHmm pennHmm = loadPennHmm(hmmFile, new Logger());
		if (!Files.exists(proj.CUSTOM_PFB_FILENAME.getValue())) {
			PennCNV.populationBAF(proj);
		}
		PFB pfb = PFB.loadPFB(proj, proj.CUSTOM_PFB_FILENAME.getValue());
		GcModel gcModel = GcAdjustor.GcModel.populateFromFile(proj.GC_MODEL_FILENAME.getValue(false, false), false, proj.getLog());
		Sample sample = proj.getFullSampleFromRandomAccessFile(proj.getSamples()[0]);
		GcAdjustor gcAdjustor = GcAdjustor.getComputedAdjustor(proj, sample, gcModel, true, true, true, false);
		double[] lrrs = gcAdjustor.getCorrectedIntensities();
		double[] bafs = Array.toDoubleArray(sample.getBAFs());

		MarkerSet markerSet = proj.getMarkerSet();
		int[][] chrIndices = markerSet.getIndicesByChr();
		int[][] snpdists = getSNPDist(proj);
		for (int i = 0; i < chrIndices.length; i++) {
			if (chrIndices[i].length > 10) {
				double[] lrrChr = Array.subArray(lrrs, chrIndices[i]);
				double[] bafsChr = Array.subArray(bafs, chrIndices[i]);
				double[] pfbsChr = Array.subArray(pfb.getPfbs(), chrIndices[i]);
				String[] markersChr = Array.subArray(markerSet.getMarkerNames(), chrIndices[i]);
				boolean[] cnDef = determineCNOnly(proj, markersChr);
				System.out.println(Array.booleanArraySum(cnDef));
				System.out.println(pfbsChr[0]);
				System.out.println(bafsChr[0]);
				System.out.println(lrrChr[0]);
				System.out.println(snpdists[i][0]);

				ViterbiLogNP_CHMM(pennHmm, lrrChr, bafsChr, pfbsChr, snpdists[i], cnDef);
			}
		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "C:/bin/pennCNV/penncnv/penncnv/lib/hhall.hmm";
		Project proj = new Project("C:/workspace/Genvisis/projects/OSv2_hg19.properties", false);

		test(proj, filename);
	}

}
