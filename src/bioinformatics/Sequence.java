package bioinformatics;

import common.ext;

public class Sequence {
//	public static final String[] SENSE = {"A", "C", "G", "T", "I", "D"};
//	public static final String[] ANTISENSE = {"T", "G", "C", "A", "I", "D"};
	public static final char[] ALLELES = {'A', 'C', 'G', 'T'};
	public static final char[] SENSE = {'A', 'C', 'G', 'T', 'I', 'D', 'a', 'c', 'g', 't', 'i', 'd'};
	public static final char[] ANTISENSE = {'T', 'G', 'C', 'A', 'I', 'D', 't', 'g', 'c', 'a', 'i', 'd'};

	public static char flip(char allele) {
		return ANTISENSE[ext.indexOfChar(allele, SENSE)];		
	}

	public static String flip(String allele) {
		return flip(allele.charAt(0))+"";		
	}

	public static boolean validAllele(String allele) {
		return allele.length()==1 && ext.indexOfChar(allele.charAt(0), SENSE)>=0;		
	}
}
