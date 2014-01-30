package affy;

//import java.io.*;
//import java.util.Vector;
//import common.Array;
//import common.Files;
//import common.CmdLine;
//import common.ext;

// target sketch? 
public class AffyPowerTools {

	/**
	 * apt-probeset-genotype -c lib/GenomeWideSNP_6.cdf -a birdseed --read-models-birdseed lib/GenomeWideSNP_6.birdseed.models --special-snps lib/GenomeWideSNP_6.specialSNPs --out-dir apt --cel-files listfile call rate >96% SD of normalized intensity (LRR) <0.35 the GC wave factor of LRR was between −0.02 and 0.02 were accepted CNV call count <80
	 * 
	 */

	public static void main(String[] args) {
		System.out.println("affy");
		String base = "/lustre/pankrat2/normDat/";
		// String intitialCelList = base+"lists/bigListFinal.txt";
		// String pbsOutDir = "/home/pankrat2/lanej/PBS/";
		String affyChunk = "cels";
		String batchChunk = "SNP_";
		String genvisisSource = "/lustre/pankrat2/normDat/output/ARICGenvisis/00src/";
		String javaClass = "affy.AffySNP6Tables";
		String project = "/home/pankrat2/lanej/projects/dbGaP_ARIC_11908.properties";

		int numBatches = 1;
		int numJobs = 8;
		int lineBuffer = 20;
		int memory = 14999;
		// double wallTime = 96.00;

		String affyQCfolderOut = "QCOut";
		String affyGenofolderOut = "genoTypeOut";
		String affySumfolderOut = "summarOut";
		String affyGenoClusterFolderOut = "clusterOut";
		String affyDetectFolder = "indDetect";
		String affyChunkProbe = "probes";
		String quantNorm = "quant-norm.pm-only.med-polish.expr.summary.txt";

		String finalCelList = base + "lists/bigListFinal.txt";
		String lists = base + "lists/";

		// String snpProbesetList = "AllGenoTypingSnps.txt";
		// String cnProbesetList = "AllCopyNumberProbes.txt";
		String park = "/home/pankrat2/lanej/park.jar -Xmx" + memory / 1024 + "G ";
		String sexFile = lists + "file_sex.txt";
		String affyLib = base + "affy/lib/";
		String outDir = base + "output/";
		String pennCNVbin = base + "pennCNV/bin/";
		String pennCNVlib = base + "pennCNV/lib/";
		// int totalMemory = memory * numJobs;
		String affyCDF = affyLib + "GenomeWideSNP_6.cdf";
		String affyBirdseedModel = affyLib + "GenomeWideSNP_6.birdseed-v2.models";
		String affySpecialSnps = affyLib + "GenomeWideSNP_6.specialSNPs";
		String affyHapMapQuant = affyLib + "hapmap.quant-norm.normalization-target.txt";
		// String affyRefInput = affyLib + "GenomeWideSNP_6.hapmap270.na33.r1.a5.ref";
		String affyChrX = affyLib + "GenomeWideSNP_6.chrXprobes";
		String affyChrY = affyLib + "GenomeWideSNP_6.chrYprobes";
		String affyQCA = affyLib + "GenomeWideSNP_6.r2.qca";
		String affyQCC = affyLib + "GenomeWideSNP_6.r2.qcc";
		// String annoFile = affyLib + "GenomeWideSNP_6.na33.annot.db";
		String affyGenoQC = affyLib + "apt-geno-qc";
		String affyGenotype = affyLib + "apt-probeset-genotype";
		String affySummarize = affyLib + "apt-probeset-summarize";
		String affyChpToText = affyLib + "apt-chp-to-txt";
		// String qcOut = outDir + "qcOut.txt";
		// String birdseedReport = "birdseed-v2.report.txt";
		String birdseedCalls = "birdseed-v2.calls.txt";
		String birdseedConf = "birdseed-v2.confidences.txt";
		String locFile = pennCNVlib + "affygw6.hg18.pfb";
		String genoCluster = "gw6.genocluster";
		String genoLrrBaf = "gw6.lrr_baf.txt";
		String pennSplitOut = "signalListFile.txt";
		String pennHmm = pennCNVlib + "affygw6.hmm";
		String detect_cnv = pennCNVbin + "penncnv/";

		// QC params
		// double callRateCut = .95;

		String[] affyQCJobs = new String[numJobs];
		String[] affyGenoJobs = new String[numJobs];
		String[] affySumJobs = new String[numJobs];
		String[] affyGenoClusterJobs = new String[numJobs];
		// String[] affyCNClusterJobs = new String[numJobs];
		String[] affySNP6Tables = new String[numJobs];
		String[] LRRBAFJobs = new String[numJobs];
		String[] startKColumn = new String[numJobs];
		String[] stopKColumn = new String[numJobs];
		String[] detectCNVs = new String[numJobs];
		String[] chpToTxt = new String[numJobs];

		for (int j = 0; j < numBatches; j++) {
			String batchID = "";
			if (j < 10) {
				batchID = batchChunk + "B0" + j + "_";
			}
			if (j >= 10) {
				batchID = batchChunk + "B" + j + "_";
			}

			// String affyGenoPBS = pbsOutDir + "GT" + batchID + ".pbs";
			// String affySumPBS = pbsOutDir + "summarOut" + batchID + ".pbs";
			// String affyPennGenoClustPBS = pbsOutDir + "PenGenoClust" + batchID + ".pbs";
			// String affyChpToTxtPBS = pbsOutDir + "ChpToTxt" + batchID + ".pbs";
			// String affyPennCNClustPBS = pbsOutDir + "PenCNClust" + batchID + ".pbs";
			// String affyPennGenoLRRBAFPBS = pbsOutDir + "PennLRR_BAF" + batchID + ".pbs";
			// String affySNP6TablesPBS = pbsOutDir + "AS6T" + batchID + ".pbs";
			//
			// String affyQCPBS = pbsOutDir + "QC" + batchID + ".pbs";
			// String affyKColStartPBS = pbsOutDir + "KColStart" + batchID + ".pbs";
			// String affyKColStopPBS = pbsOutDir + "KColStop" + batchID + ".pbs";
			// String affyDetectPBS = pbsOutDir + "JLDetectB0" + j + ".pbs";

			for (int i = 0; i < numJobs; i++) {
				String jobID = "";
				String detectBatch = "";
				if (i < 10) {
					jobID = batchID + "0" + i;
					detectBatch = "0" + i;
				}
				if (i >= 10) {
					jobID = batchID + i;
					detectBatch = "" + i;
				}

				String qc = affyGenoQC + " -c " + affyCDF + " --qca-file " + affyQCA + " -qcc-file " + affyQCC + " --chrX-probes " + affyChrX + " --chrY-probes " + affyChrY + " --out-file " + outDir + affyQCfolderOut + "/" + jobID + ".qc --cel-files " + lists + affyChunk + jobID;
				affyQCJobs[i] = qc;

				String genotypeCommand = affyGenotype + " -c " + affyCDF + " --cc-chp-output --table-output true -a birdseed-v2 --set-gender-method cn-probe-chrXY-ratio --read-models-birdseed " + affyBirdseedModel + " --special-snps " + affySpecialSnps + " -out-dir " + outDir + affyGenofolderOut + jobID + "/ --cel-files " + finalCelList + " --chrX-probes " + affyChrX + " --chrY-probes " + affyChrY + " --probeset-ids " + lists + affyChunkProbe + jobID;
				affyGenoJobs[i] = genotypeCommand;

				String aptChptoTxt = affyChpToText + " " + outDir + affyGenofolderOut + jobID + "/cc-chp*/*.chp -o " + outDir + "ARICGenvisis/00src/" + affyGenofolderOut + jobID + "/cc-chp";
				chpToTxt[i] = aptChptoTxt;

				String summarizeCommand = affySummarize + " --cdf-file " + affyCDF + " --analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true --target-sketch " + affyHapMapQuant + " --out-dir " + outDir + affySumfolderOut + jobID + "/ --cel-files " + finalCelList + " --probeset-ids " + lists + affyChunkProbe + jobID;

				affySumJobs[i] = summarizeCommand;

				String affySNP6TablesCommand = "/soft/java/jdk1.7.0_45/bin/java -cp " + park + " " + javaClass + " proj=" + project + " calls=" + outDir + affyGenofolderOut + jobID + "/" + birdseedCalls + " conf=" + outDir + affyGenofolderOut + jobID + "/" + birdseedConf + " sig=" + outDir + affySumfolderOut + jobID + "/" + quantNorm + " out=" + genvisisSource + jobID + "/ -" + batchChunk + " numLinesBuffer=" + lineBuffer;

				affySNP6Tables[i] = affySNP6TablesCommand;
				// Include sex file later
				String generate_affy_geno_cluster = "perl " + pennCNVbin + "generate_affy_geno_cluster.pl " + outDir + affyGenofolderOut + jobID + "/" + birdseedCalls + " " + outDir + affyGenofolderOut + jobID + "/" + birdseedConf + " " + outDir + affySumfolderOut + jobID + "/" + quantNorm + " " + batchChunk + " -locfile " + locFile + " -sexfile " + sexFile + " -out " + outDir + affyGenoClusterFolderOut + jobID + "/" + genoCluster;

				affyGenoClusterJobs[i] = generate_affy_geno_cluster;

				// String generate_affy_CN_cluster = "perl " + pennCNVbin+"generate_affy_CN_cluster.pl "+
				// outDir + affyGenofolderOut + jobID +"/" + birdseedCalls
				// + " " + outDir + affyGenofolderOut + jobID +"/" + birdseedConf + " "
				// + outDir + affySumfolderOut + jobID +"/" +quantNorm + " -locfile " +
				// locFile + " -sexfile "+ sexFile + " -out " + outDir + affyGenoClusterFolderOut + jobID + "/" + genoCluster;
				//
				// affyCNClusterJobs[i] = generate_affy_CN_cluster;

				String lrrBafCalc = "perl " + pennCNVbin + "normalize_affy_geno_cluster.pl " + outDir + affyGenoClusterFolderOut + jobID + "/" + genoCluster + " " + outDir + affySumfolderOut + jobID + "/" + quantNorm + " -locfile " + locFile + " -out " + outDir + affyGenoClusterFolderOut + jobID + "/" + genoLrrBaf;

				LRRBAFJobs[i] = lrrBafCalc;

				String splitSignalFileStart = "perl " + pennCNVbin + "kcolumn.pl " + outDir + affyGenoClusterFolderOut + jobID + "/" + genoLrrBaf + " split 2 -tab -head 3 --name --start_split 5001 -out " + outDir + "Kcol/Kcol" + jobID + "/gw6_split --filenameout " + outDir + "Kcol/Kcol" + jobID + "/" + pennSplitOut;

				startKColumn[i] = splitSignalFileStart;

				String splitSignalFileStop = "perl " + pennCNVbin + "kcolumn.pl " + outDir + affyGenoClusterFolderOut + jobID + "/" + genoLrrBaf + " split 2 -tab -head 3 --name --end_split 5000 -out " + outDir + "Kcol/Kcol" + jobID + "/gw6_split --filenameout " + outDir + "Kcol/Kcol" + jobID + "/" + pennSplitOut;
				stopKColumn[i] = splitSignalFileStop;

				String detectCNV = "perl " + detect_cnv + "detect_cnv.pl --test -hmm " + pennHmm + " -pfb " + locFile + " -conf -gcmodel " + pennCNVlib + "affygw6.hg18.gcmodel -log " + outDir + affyDetectFolder + detectBatch + "/gw6.log -out " + outDir + affyDetectFolder + detectBatch + "/gw6.rawcnv" + " -list " + lists + "indDetect" + detectBatch;

				detectCNVs[i] = detectCNV;
				// System.out.println(qc);
				// System.out.println(generate_affy_geno_cluster);
				System.out.println(detectCNV);
			}

			// Files.qsubMultiple(affyQCPBS , affyQCJobs, numJobs, memory ,totalMemory, wallTime);
			// Files.qsubMultiple(affyGenoPBS , affyGenoJobs, numJobs, memory ,totalMemory, wallTime);
			// Files.qsubMultiple(affyChpToTxtPBS , chpToTxt, numJobs, memory ,totalMemory, wallTime);

			// Files.qsubMultiple(affySumPBS , affySumJobs, numJobs, memory ,totalMemory, wallTime);
			// Files.qsubMultiple(affySNP6TablesPBS , affySNP6Tables, numJobs, memory ,totalMemory, wallTime);
			// Files.qsubMultiple(affyPennGenoClustPBS , affyGenoClusterJobs, numJobs, memory ,totalMemory, wallTime);
			// Files.qsubMultiple(affyPennCNClustPBS , affyGenoClusterJobs, numJobs, memory ,totalMemory, wallTime);

			// Files.qsubMultiple(affyPennGenoLRRBAFPBS , LRRBAFJobs, numJobs, memory ,totalMemory, wallTime);

			// Files.qsubMultiple(affyKColStartPBS , startKColumn, numJobs, memory ,totalMemory, wallTime);
			// Files.qsubMultiple(affyKColStopPBS , stopKColumn, numJobs, memory ,totalMemory, wallTime);
			// Files.qsubMultiple(affyDetectPBS, detectCNVs, numJobs, memory, totalMemory, wallTime);
		}
	}

}

// works
// String qcCommand = affyGenoQC + " -c " + affyCDF + " --qca-file " + affyQCA + " -qcc-file " +
// affyQCC + " --chrX-probes " + affyChrX +" --chrY-probes " + affyChrY + " --cel-files "
// + intitialCelList + " --out-file " + qcOut;
// String qcCommand = affyGenoQC + " -c " + affyCDF + " --qca-file " + affyQCA + " -qcc-file " +
// affyQCC + " --chrX-probes " + affyChrX +" --chrY-probes " + affyChrY + " --out-file $wdir/qcOut.qc" + " *.CEL";
// //works
// String genotypeCommand = affyGenotype + " --summaries --write-models -c " + affyCDF + " -a birdseed-v2 --read-models-birdseed " +
// affyBirdseedModel + " --special-snps " + affySpecialSnps + " -out-dir " +
// outDir + " --cel-files " + intitialCelList + " --chrX-probes " + affyChrX +" --chrY-probes " +
// affyChrY ;
// works
// String summarizeCommand = "apt-probeset-summarize --cdf-file " + affyCDF +
// " --analysis quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true --target-sketch "+
// affyHapMapQuant + " --out-dir " + outDir + " --cel-files " + finalCelList;
// works
// String copyNumberWorkflow = "apt-copynumber-workflow -v 1 --cdf-file " + affyCDF +
// " --chrX-probes " + affyChrX +" --chrY-probes " + affyChrY + " --special-snps "
// + affySpecialSnps + " --annotation-file " + annoFile + " --reference-input " +
// affyRefInput + " -out-dir " + outDir + " --cel-files " + finalCelList + " --text-output true";
// works
// String generate_affy_geno_cluster = "perl " + affyPennCnv+"generate_affy_geno_cluster.pl "+
// birdseedCalls + " " + birdseedConf + " " + quantNorm + " -locfile " +
// locFile + " -sexfile  " + sexFile + " -out " + genoCluster;
// works
// String lrrBafCalc = "perl " + affyPennCnv+"normalize_affy_geno_cluster.pl "+
// genoCluster + " " + quantNorm + " -locfile " + locFile +
// " -out " + genoLrrBaf;
// //works
// String splitSignalFileStart = "perl " + pennCNV+"kcolumn.pl "+
// genoLrrBaf + " split 2 -tab -head 3 --name_by_header -out gw6 --filenameout " +pennSplitOut ;
// //works
// String detectCNV = "perl " + pennCNV+"detect_cnv.pl --test -hmm " + pennHmm +
// " -pfb " + locFile + " -list " + pennSplitOut + " -log gw6.log -out gw6.rawcnv";
//

// System.out.println(qcCommand);
// CmdLine.run(qcCommand , outDir );
//
// genFinalCelList(qcOut , finalCelList , callRateCut);
// System.out.println(genotypeCommand);
// CmdLine.run(genotypeCommand , outDir);
//
// System.out.println(summarizeCommand);
// CmdLine.run(summarizeCommand , outDir);
//
// genSexFile(birdseedReport, sexFile);
// //System.out.println(copyNumberWorkflow);
// //CmdLine.run(copyNumberWorkflow , outDir);
// System.out.println(generate_affy_geno_cluster);
// CmdLine.run(generate_affy_geno_cluster , outDir);

// System.out.println(lrrBafCalc);
// CmdLine.run(lrrBafCalc , outDir);
// run command from output directory

// System.out.println(splitSignalFileStart);
// CmdLine.run(splitSignalFileStart , outDir);

// System.out.println(detectCNV);
// CmdLine.run(detectCNV , outDir);