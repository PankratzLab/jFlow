package org.genvisis.seq.manage;

import java.util.ArrayList;
import java.util.List;

import org.genvisis.common.Array;
import org.genvisis.common.Logger;
import org.genvisis.common.Sort;
import org.genvisis.seq.analysis.Blast.BlastResults;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

public class CigarOps {

	/**
	 * @param cigars
	 * @return an array of {@link Cigar} ordered by their length matched to the reference
	 */
	public static int[] sortByRefMatchLength(Cigar[] cigars) {
		int[] refLengthMatch = new int[cigars.length];
		for (int i = 0; i < cigars.length; i++) {
			int reftmp = getRefLength(cigars[i]);
			refLengthMatch[i] = reftmp;
		}
		int[] order = Sort.quicksort(refLengthMatch, Sort.DESCENDING);

		return order;
	}

	public static int getRefLength(Cigar cigar) {
		int reftmp = 0;
		for (CigarElement cigarElement : cigar.getCigarElements()) {
			if (cigarElement.getOperator() == CigarOperator.EQ) {
				reftmp += cigarElement.getLength();
			}
		}
		return reftmp;
	}

	/**
	 * @param cigars
	 * @return String[] of the cigar representations
	 */
	public static String[] toStringArray(Cigar[] cigars) {
		String[] cigarStrings = new String[cigars.length];
		for (int i = 0; i < cigarStrings.length; i++) {
			cigarStrings[i] = toString(cigars[i]);
		}
		return cigarStrings;
	}

	public static String toString(Cigar cigar) {
		return cigar.toString();
	}

	public static Cigar getConstantCigar(int length, CigarOperator ci) {
		ArrayList<CigarElement> cigarElements = new ArrayList<CigarElement>();
		cigarElements.add(new CigarElement(length, ci));

		return new Cigar(cigarElements);
	}

	/**
	 * @param blastResults
	 *            must have the btop entry present or will fail
	 * @param initialSequencLength
	 *            used to check the resulting length of the cigar string
	 * @param log
	 * @return the {@link Cigar } representation of the btop in {@link BlastResults}
	 */
	public static Cigar convertBtopToCigar(BlastResults blastResults, int initialSequencLength, Logger log) {
		String btop = blastResults.getBtop();
		ArrayList<CigarElement> cigarElements = new ArrayList<CigarElement>();
		Cigar cigar = null;

		if (btop == null) {
			String error = "Blast results must contain a \"btop\" entry in order to use this function...";
			log.reportTimeError(error);
			throw new IllegalArgumentException(error);
		} else if (btop.equals("NA")) {
			cigarElements.add(new CigarElement(initialSequencLength, CigarOperator.X));
			cigar = new Cigar(cigarElements);
		} else {
			if (isAllMatched(blastResults, initialSequencLength)) {// perfect alignment, completely equals the reference...
				cigarElements.add(new CigarElement(initialSequencLength, CigarOperator.EQ));
				cigar = new Cigar(cigarElements);
			} else if (isInt(btop)) {// partial alignment of all matching bps
				int alignmentLength = Integer.parseInt(btop);
				if (blastResults.getQstart() - 1 > 0) {//
					cigarElements.add(new CigarElement(blastResults.getQstart() - 1, CigarOperator.X));
				}
				cigarElements.add(new CigarElement(alignmentLength, CigarOperator.EQ));
				if (blastResults.getQstop() != initialSequencLength) {
					cigarElements.add(new CigarElement(initialSequencLength - blastResults.getQstop(), CigarOperator.X));
				}
				cigar = new Cigar(cigarElements);
				if (cigar.getReadLength() != initialSequencLength) {
					String error = "INT ONLY REP: Cigar length representation of " + cigar.getReadLength() + " did not equal the query length of " + initialSequencLength;
					error += "\n BLAST:  " + Array.toStr(blastResults.getResults());
					error += "\n CIGAR:  " + cigar.toString();
					log.reportTimeError(error);
					throw new IllegalArgumentException(error);
				}
			} else {// has gaps and or mismatches
				if (blastResults.getQstart() - 1 > 0) {//
					cigarElements.add(new CigarElement(blastResults.getQstart() - 1, CigarOperator.X));
				}
				String[] btopBroken = breakUpBtop(btop, log);
				for (int i = 0; i < btopBroken.length; i++) {
					if (isInt(btopBroken[i])) {
						cigarElements.add(new CigarElement(Integer.parseInt(btopBroken[i]), CigarOperator.EQ));
					} else {
						if (btopBroken[i].length() == 2) {

							if (btopBroken[i].startsWith("-")) {// query has deletion compared to ref
								cigarElements.add(new CigarElement(1, CigarOperator.D));
								// Insertion consumes read bases, subject = ref query= probe

							} else if (btopBroken[i].endsWith("-")) {// query has insertion compared to ref
								cigarElements.add(new CigarElement(1, CigarOperator.I));
							} else {
								cigarElements.add(new CigarElement(1, CigarOperator.X));
							}
						} else {
							String error = "Non int btop strings were supposed to have length two..";
							log.reportTimeError(error);
							throw new IllegalArgumentException(error);
						}
					}
				}
				if (blastResults.getQstop() != initialSequencLength) {
					cigarElements.add(new CigarElement(initialSequencLength - blastResults.getQstop(), CigarOperator.X));
				}
				cigar = new Cigar(cigarElements);
				if (cigar.getReadLength() != initialSequencLength) {
					String error = "STRING INT REP: Cigar length representation of " + cigar.getReadLength() + " did not equal the query length of " + initialSequencLength;
					error += "\n BLAST:  " + Array.toStr(blastResults.getResults());
					error += "\n CIGAR:  " + cigar.toString();

					log.reportTimeError(error);
					throw new IllegalArgumentException(error);
				}
			}
		}

		if (cigar != null && blastResults.getSstart() > blastResults.getSstop()) {

			// flip the strand of the cigar
			List<CigarElement> cigarElementsStrandCurrent = cigar.getCigarElements();
			ArrayList<CigarElement> cigarElementsStrandFlip = new ArrayList<CigarElement>();
			for (int i = cigarElementsStrandCurrent.size() - 1; i >= 0; i--) {
				cigarElementsStrandFlip.add(cigarElementsStrandCurrent.get(i));
			}
			// System.out.println(cigar.toString());
			cigar = new Cigar(cigarElementsStrandFlip);
			// System.out.println(cigar.toString());
		}

		int currentReadLength = cigar.getReadLength();
		int currentRefLength= cigar.getPaddedReferenceLength();
		//System.out.println(btop + "\t" + "HI1" + cigar);

		Cigar original = cigar;
		cigar = uniquify(cigar);

		if (currentReadLength != cigar.getReadLength() || cigar.getPaddedReferenceLength() != currentRefLength) {
			String error = "Could not properly uniqify " + original.toString() + ", came out as " + cigar.toString();
			System.out.println(currentReadLength + "\t" + cigar.getCigarElements().size() + "\t" + cigar.toString());
			System.out.println(currentReadLength + "\t" + cigar.getCigarElements().size() + "\t" + original.toString());
			log.reportTimeError(error);
			throw new IllegalStateException(error);
		}

		if (cigar != null && cigar.getReadLength() != initialSequencLength) {
			String error = "Cigar length representation of " + cigar.getReadLength() + " did not equal the query length of " + initialSequencLength;
			log.reportTimeError(error);
			throw new IllegalArgumentException(error);
		}

		// log.reportTimeInfo("BTOP -> " + blastResults.getBtop());
		// log.reportTimeInfo(" : CIGAR " + cigar.toString());
		return cigar;
	}

	/**
	 * @param cigar
	 * @return a new cigar with any identical {@link CigarElement} in a row combined
	 */
	private static Cigar uniquify(Cigar cigar) {
		ArrayList<CigarElement> unique = new ArrayList<CigarElement>();
		if (cigar.getCigarElements().size() == 1) {
			// System.out.println("HI2\t" + cigar);
			return cigar;
		} else {
			// System.out.println("HI3\t" + cigar);

			int i = 1;
			CigarElement tmp = cigar.getCigarElements().get(0);
			while (i < cigar.getCigarElements().size()) {
				CigarElement current = cigar.getCigarElements().get(i);
				//System.out.println(current.getLength() + "\t" + current.getOperator());
				i++;
				if (tmp.getOperator() == current.getOperator()) {
					tmp = new CigarElement(tmp.getLength() + current.getLength(), tmp.getOperator());
				} else {
					unique.add(tmp);
					tmp = current;
				}
			}
			unique.add(tmp);
			return new Cigar(unique);
		}
	}

	/**
	 * Break up the btop entry into Integer strings and String strings for later processing
	 * 
	 * @param btop
	 * @param log
	 * @return
	 */
	private static String[] breakUpBtop(String btop, Logger log) {
		ArrayList<String> btopBroken = new ArrayList<String>();
		String currentInt = null;
		String currentString = null;
		for (int i = 0; i < btop.length(); i++) {
			if (currentString != null && currentString.length() == 2) {
				btopBroken.add(currentString);
				currentString = null;
			}

			if (isInt(btop.charAt(i) + "")) {// integer representing reference match
				if (currentString != null) {
					if (currentString.length() != 2) {
						String error = "Internal error, unaccounted for string breakup length";
						log.reportTimeError(error);
						throw new IllegalStateException(error);
					}
					btopBroken.add(currentString);
					currentString = null;
				}
				if (currentInt == null) {
					currentInt = "";
				}
				currentInt += btop.charAt(i);
			}

			else {
				if (currentInt != null) {
					btopBroken.add(currentInt);
					currentInt = null;
				}
				if (currentString == null) {
					currentString = "";
				}

				currentString += btop.charAt(i);
			}
		}
		if (currentInt != null && currentString != null) {
			String error = "Internal error, unaccounted for breakup";
			log.reportTimeError(error);
			throw new IllegalStateException(error);
		}
		if (currentInt != null) {
			btopBroken.add(currentInt);
		}
		if (currentString != null) {
			btopBroken.add(currentString);
		}
		return Array.toStringArray(btopBroken);
	}

	private static boolean isAllMatched(BlastResults blastResults, int initialSequencLength) {
		return blastResults.getAlignmentLength() == initialSequencLength && blastResults.getGapOpens() == 0 && blastResults.getMismatches() == 0;
	}

	private static boolean isInt(String potentialInt) {
		try {
			Integer.parseInt(potentialInt);
			return true;
		} catch (NumberFormatException nfe) {
			return false;
		}
	}

}
