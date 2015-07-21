package seq.manage;

import java.util.ArrayList;

import common.Sort;
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
			int reftmp = 0;
			for (CigarElement cigarElement : cigars[i].getCigarElements()) {
				if (cigarElement.getOperator() == CigarOperator.EQ) {
					reftmp += cigarElement.getLength();
				}
			}
			refLengthMatch[i] = reftmp;
		}
		int[] order = Sort.quicksort(refLengthMatch, Sort.DESCENDING);

		return order;
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

}
