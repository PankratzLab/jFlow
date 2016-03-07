package seq.manage;

import htsjdk.tribble.annotation.Strand;

public class StrandOps {

	public static String[] flipIfNeeded(String[] b, Strand strand, boolean ignoreInvalidAlleles) {
		String[] flipped = new String[b.length];
		for (int i = 0; i < flipped.length; i++) {
			flipped[i] = flipIfNeeded(b[i], strand, ignoreInvalidAlleles);
		}
		return flipped;
	}

	public static String flipsIfNeeded(String b, Strand strand, boolean ignoreInvalidAlleles) {
		return flipsIfNeeded(b, strand, ignoreInvalidAlleles, false);
	}

	public static String flipsIfNeeded(String b, Strand strand, boolean ignoreInvalidAlleles, boolean reverse) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < b.length(); i++) {
			sb.append(flipIfNeeded(b.charAt(i) + "", strand, ignoreInvalidAlleles));
		}
		return reverse ? sb.reverse().toString() : sb.toString();
	}

	public static String flipIfNeeded(String b, Strand strand, boolean ignoreInvalidAlleles) {
		if (strand == Strand.NEGATIVE) {
			if (b.equals("A")) {
				return "T";
			} else if (b.equals("G")) {
				return "C";
			} else if (b.equals("C")) {
				return "G";
			} else if (b.equals("T")) {
				return "A";
			} else {
				if (!ignoreInvalidAlleles) {
					throw new IllegalArgumentException("Invalid base for strand flip " + b);
				}
				return b;
			}
		} else {
			return b;
		}
	}
}
