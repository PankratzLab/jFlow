package org.genvisis.one.jl.aging;

import java.io.File;

import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.seq.SeqVariables.ASSEMBLY_NAME;
import org.genvisis.seq.SeqVariables.PLATFORM;
import org.genvisis.seq.analysis.genage.Pipeline;
import org.genvisis.seq.manage.BamOps;
import org.genvisis.sra.SRARunTable;
import org.genvisis.sra.SRASample;

/**
 * 
 * Done
 */
public class TestPipe {

	private TestPipe() {

	}

	
	
	public static void main(String[] args) {

		SRARunTable sraRunTable = SRARunTable.load("/Volumes/Beta/data/aric_sra/prep/SraRunTable.txt", new Logger());
		String outDir = "/Volumes/Beta/data/aric_sra/test/";
		String refGenome = "/Volumes/Beta/ref/GRCh37_canon.fa";
		new File(outDir).mkdirs();

		Logger log = new Logger(outDir + ".log");
		String captureBed = "/Volumes/Beta/ref/VCRome_2_1_hg19_capture_targets.bed";
		String targetBam = "/Volumes/Beta/data/aric_sra/bams/SRR1654226.bam";
		BamOps.verifyIndex(targetBam, log);
		String sraSamp = ext.rootOf(targetBam);

		SRASample current = sraRunTable.get(sraSamp);
		log.reportTimeInfo(sraRunTable.get(sraSamp).toString());
		String bamList = outDir + sraSamp + ".bamList";
		Files.write(targetBam, bamList);
		if (current.getaName() == ASSEMBLY_NAME.GRCH37 && current.getPlatform() == PLATFORM.ILLUMINA) {
//			Pipeline.pipeline(targetBam, outDir, refGenome, captureBed, current.getaType(), current.getaName(), 1, log);
		} else {
			log.reportTimeWarning("Skipping sample " + current.toString());
		}

	}

}
