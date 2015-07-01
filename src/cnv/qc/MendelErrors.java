package cnv.qc;

/**
 * @author Kitty Check for mendelian errors in a trio , taken from http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#mendel
 */
public class MendelErrors {

	private byte chr;
	private int offSex;

	private byte offGenotype;
	private byte moGenotype;
	private byte faGenotype;

	public MendelErrors(byte chr, int offSex, byte offGenotype,
			byte moGenotype, byte faGenotype) {
		super();
		this.chr = chr;
		this.offSex = offSex;
		this.offGenotype = offGenotype;
		this.moGenotype = moGenotype;
		this.faGenotype = faGenotype;
	}

	public MendelErrorCheck checkMendelError() {
		if (offGenotype == -1) {
			return new MendelErrorCheck(-1, false, false);

		}
		if (moGenotype == -1 && faGenotype == -1) {
			return new MendelErrorCheck(-1, false, false);
		}
		if (faGenotype == 0 && moGenotype == 0 && offGenotype == 1) {
			return new MendelErrorCheck(1, true, true);
		}
		if (faGenotype == 2 && moGenotype == 2 && offGenotype == 1) {
			return new MendelErrorCheck(2, true, true);
		}
		if (faGenotype == 2 && offGenotype == 0) {
			return new MendelErrorCheck(3, true, false);
		}
		if (moGenotype == 2 && offGenotype == 0) {
			return new MendelErrorCheck(4, false, true);
		}
		if (faGenotype == 2 && moGenotype == 2 && offGenotype == 0) {
			return new MendelErrorCheck(5, true, true);
		}
		if (faGenotype == 0 && offGenotype == 2) {
			return new MendelErrorCheck(6, true, false);
		}
		if (moGenotype == 0 && offGenotype == 2) {
			return new MendelErrorCheck(7, false, true);
		}
		if (faGenotype == 0 && moGenotype == 0 && offGenotype == 2) {
			return new MendelErrorCheck(8, true, true);
		}
		if (chr == 23 && offSex == 1 && moGenotype == 0 && offGenotype == 2) {
			return new MendelErrorCheck(9, false, true);
		}
		if (chr == 23 && offSex == 1 && moGenotype == 2 && offGenotype == 1) {
			return new MendelErrorCheck(10, false, true);
		}
		return new MendelErrorCheck(-1, false, false);

	}

	/**
	 * @author Kitty Stores the mendel error, if any, and where it occurred
	 */
	public static class MendelErrorCheck {

		// The error codes are as follows:
		// -1 = no error
		// Code Pat , Mat -> Offspring
		//
		// 1 AA , AA -> AB
		// 2 BB , BB -> AB
		//
		// 3 BB , ** -> AA
		// 4 ** , BB -> AA
		// 5 BB , BB -> AA
		//
		// 6 AA , ** -> BB
		// 7 ** , AA -> BB
		// 8 AA , AA -> BB
		//
		// 9 ** , AA -> BB (X chromosome male offspring)
		// 10 ** , BB -> AA (X chromosome male offspring)

		private int errorCode;
		private boolean faMendelError;
		private boolean maMendelError;

		public MendelErrorCheck(int errorCode, boolean faMendelError,
				boolean maMendelError) {
			super();
			this.errorCode = errorCode;
			this.faMendelError = faMendelError;
			this.maMendelError = maMendelError;
		}

		public boolean hasError() {
			return faMendelError || maMendelError;
		}

		public int getErrorCode() {
			return errorCode;
		}

		public boolean isFaMendelError() {
			return faMendelError;
		}

		public boolean isMaMendelError() {
			return maMendelError;
		}
	}
}
