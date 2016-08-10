package org.genvisis.one.JL.aging;

import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.seq.analysis.MitoSeqCN;
import org.genvisis.sra.SRARunTable;

/**
 * 
 *
 */
public class TestPipe {

	private TestPipe() {

	}

	public static void main(String[] args) {
		SRARunTable sraRunTable = SRARunTable.load("/Volumes/Beta/data/aric_sra/prep/SraRunTable.txt", new Logger());
		String targetBam = "/Volumes/Beta/data/aric_sra/bams/SRR1654226.bam";
		
		
		String sraSamp = ext.rootOf(targetBam);
		System.out.println(sraRunTable.get(sraSamp).toString());
		
	//	MitoSeqCN.run(fileOfBams, outDir, captureBed, referenceGenomeFasta, params, numthreads);
		
		
		
		
	}

}
