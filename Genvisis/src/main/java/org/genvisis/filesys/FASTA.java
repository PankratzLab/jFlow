package org.genvisis.filesys;

import org.genvisis.common.ext;

public final class FASTA {
	
	private static final String INDEX_EXT = ".fai";
	private static final String DICTIONARY_EXT = ".dict";

	private FASTA() {
		// prevent instantiation of utility class
	}
	
	public static final String getIndex(String fasta) {
		return fasta + INDEX_EXT;
	}
	
	public static final String getDictionary(String fasta) {
		return ext.addToRoot(fasta, DICTIONARY_EXT);
	}
}
