package cnv.qc;

import java.io.FileWriter;
import java.io.PrintWriter;

import cnv.filesys.MarkerData;
import cnv.filesys.Pedigree;
import cnv.filesys.Project;
import cnv.manage.MDL;

/**
 * @author Kitty Check for mendelian errors in a pedigree , taken from http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#mendel
 */

public class MendelErrors {

	private byte chr;
	private int offSex;

	private byte offGenotype;
	private byte moGenotype;
	private byte faGenotype;

	/**
	 * @param chr
	 *            Only used for chr23 and male offspring
	 * @param offSex
	 *            Only used for chr23 and male offspring
	 * @param offGenotype
	 *            -1,0,1,2 genotype
	 * @param moGenotype
	 *            -1,0,1,2 genotype
	 * @param faGenotype
	 *            -1,0,1,2 genotype
	 */
	public MendelErrors(byte chr, int offSex, byte offGenotype, byte faGenotype, byte moGenotype) {
		super();
		this.chr = chr;
		this.offSex = offSex;
		this.offGenotype = offGenotype;
		this.faGenotype = faGenotype;
		this.moGenotype = moGenotype;
	}

	public MendelErrorCheck checkMendelError() {
		if (offGenotype == -1) {
			return new MendelErrorCheck(-1, false, false);
		}
		if (moGenotype == -1 && faGenotype == -1) {
			return new MendelErrorCheck(-1, false, false);
		}
		if (chr < 23) {
			// note these are not in error code order
			if (faGenotype == 0 && moGenotype == 0 && offGenotype == 1) {
				return new MendelErrorCheck(1, true, true);
			}
			if (faGenotype == 2 && moGenotype == 2 && offGenotype == 1) {
				return new MendelErrorCheck(2, true, true);
			}
			if (faGenotype == 2 && moGenotype == 2 && offGenotype == 0) {
				return new MendelErrorCheck(5, true, true);
			}
			if (faGenotype == 0 && moGenotype == 0 && offGenotype == 2) {
				return new MendelErrorCheck(8, true, true);
			}
			if (faGenotype == 2 && offGenotype == 0) {
				return new MendelErrorCheck(3, true, false);
			}
			if (moGenotype == 2 && offGenotype == 0) {
				return new MendelErrorCheck(4, false, true);
			}
			if (faGenotype == 0 && offGenotype == 2) {
				return new MendelErrorCheck(6, true, false);
			}
			if (moGenotype == 0 && offGenotype == 2) {
				return new MendelErrorCheck(7, false, true);
			}

		} else if (chr == 23) {
			if (offSex == 1 && moGenotype == 0 && offGenotype == 2) {
				return new MendelErrorCheck(9, false, true);
			}
			if (offSex == 1 && moGenotype == 2 && offGenotype == 0) {
				return new MendelErrorCheck(10, false, true);
			}
		} else if (chr == 26) {

			if (moGenotype == 2 && offGenotype != 2) {
				return new MendelErrorCheck(11, false, true);
			}
			if (moGenotype == 0 && offGenotype != 0) {
				return new MendelErrorCheck(12, false, true);
			}
		}
		return new MendelErrorCheck(-1, false, false);
	}

	/**
	 * @author Kitty Stores the mendel error, if any, and where it occurred
	 */
	public static class MendelErrorCheck {

		/**
		 * 
		 * -1 is no error<br>
		 * Code Pat , Mat -> Offspring<br>
		 * 
		 * 1 AA , AA -> AB <br>
		 * 2 BB , BB -> AB<br>
		 * 
		 * 3 BB , ** -> AA<br>
		 * 4 ** , BB -> AA<br>
		 * 5 BB , BB -> AA<br>
		 * 
		 * 6 AA , ** -> BB<br>
		 * 7 ** , AA -> BB<br>
		 * 8 AA , AA -> BB<br>
		 * 
		 * 9 ** , AA -> BB (X chromosome male offspring)<br>
		 * 10 ** , BB -> AA (X chromosome male offspring)<br>
		 * 
		 * 11 **,AA -> !AA (Mitochondrial, maternal) <br>
		 * 12 **,BB -> !BB (Mitochondrial, maternal)
		 * 
		 */
		private int errorCode;
		private boolean faMendelError;
		private boolean moMendelError;

		public MendelErrorCheck(int errorCode, boolean faMendelError, boolean moMendelError) {
			super();
			this.errorCode = errorCode;
			this.faMendelError = faMendelError;
			this.moMendelError = moMendelError;
		}

		public boolean hasError() {
			return faMendelError || moMendelError;
		}

		public int getErrorCode() {
			return errorCode;
		}

		public boolean hasFaMendelError() {
			return faMendelError;
		}

		public boolean hasMoMendelError() {
			return moMendelError;
		}
	}

	public static void detectMendelMarkers(Project proj) {
		Pedigree pedigree = proj.loadPedigree();
		if (pedigree == null) {
			return;
		} else {
			String output = proj.PROJECT_DIRECTORY.getValue() + "mendelErrorMarkers.txt";
			try {
				PrintWriter writer = new PrintWriter(new FileWriter(output));
				writer.println("MarkerName\tNumMendelErrors");
				boolean[] samplesToCheck = proj.getSamplesToInclude(null);
				MDL mdl = new MDL(proj, proj.getMarkerNames(), 2, 100);
				while (mdl.hasNext()) {
					MarkerData markerData = mdl.next();
					MendelErrorCheck[] mendelErrorChecks = Pedigree.PedigreeUtils.checkMendelErrors(pedigree, markerData, samplesToCheck, null, null, 0);
					int num = 0;
					for (int i = 0; i < mendelErrorChecks.length; i++) {
						if (mendelErrorChecks[i].getErrorCode() > 0) {
							num++;
						}
					}
					if (num > 0) {
						writer.println(markerData.getMarkerName() + "\t" + num);
						proj.getLog().reportTimeInfo(markerData.getMarkerName() + "\t" + num);
					}

				}

				writer.close();
			} catch (Exception e) {
				proj.getLog().reportError("Error writing to " + output);
				proj.getLog().reportException(e);
			}

		}

	}

	public static void main(String[] args) {
		int numArgs = args.length;
		String filename = "C:/workspace/Genvisis/projects/Poynter_PCs.properties";
		String logfile = null;

		String usage = "\n" + "cnv.qc.MendelErrors requires 0-1 arguments\n" + "   (1) filename (i.e. proj=" + filename + " (default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("proj=")) {
				filename = args[i].split("=")[1];
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			Project proj = new Project(filename, false);
			detectMendelMarkers(proj);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
