package org.genvisis.one.JL.mica;

import java.io.File;

import org.genvisis.common.Logger;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;

public class HLANucFast {

	public static void main(String[] args) {
		String fa = "/Users/Kitty/git/Analysis/HLA/hla_nuc.fasta";
		Logger log = new Logger();
		File f = new File(fa);
		FastaSequenceFile hlaFa = new FastaSequenceFile(f, false);

		ReferenceSequence ref = hlaFa.nextSequence();
		while (ref != null) {
			if (ref.getName().contains("MICA")) {
				log.report(ref.getName() + "\t" + ref.getBases().length);
			}
			ref = hlaFa.nextSequence();
		}
		hlaFa.close();

	}

}
