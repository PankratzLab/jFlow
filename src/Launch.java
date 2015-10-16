import java.io.*;

import one.SuperNovo;
import cnv.analysis.FilterCalls;
import cnv.analysis.MeanLRR;
import cnv.analysis.cnvTrio;
import cnv.qc.CNVConcordance;
import cnv.qc.CNVFilter;
import cnv.qc.CNVTrioFilter;
import link.Heritability;
import link.TrimFam;
import mining.Transformations;
import bioinformatics.*;
import common.*;
import parse.*;
import seq.Vcf;
import seq.manage.VCF;
import seq.manage.VCFOps;
import gwas.*;
import db.*;

public class Launch {
	public static final String[] LAUNCH_TYPES = { "lookup - using a list of keys, pull data from multiple files", 
		"dummy - create a dummy dataset with N instances of each pattern", 
		"counts", 
		"miss", 
		"indep", 
		"genes", 
		"filterSNPs - filters SNP positions based on a set of regions with start and end positions", 
		"filterByLists - filter unique IDs via a keeps file and a removes file", 
		"plink - convert a PLINK text data set (plink.ped/plink.map) into a tab-delimited spreadsheet", 
		"simpleM - generate a simpleM script to estimate the number of effective tests", 
		"score", 
		"parse - use GenParser to edit a text file", 
		"ucsc", 
		"split", 
		"cat - concatenate the specified files", 
		"rename - rename the specified files", 
		"db", 
		"merge", 
		"mergeSNPs", 
		"trimFam - use the TrimFam algorithm to reduce the complexity of a set of pedigrees down to critical individuals and those needed to connect them", 
		"freq - computes weighted allele frequency", 
		"uniform - creates a hits control file where each file listed has the same column names, only with a different prefix", 
		"metal - perform a meta-analysis using METAL", 
		"transform", 
		"forest", 
		"unique", 
		"dir", 
		"copy", 
		"meta", 
		"gwaf", 
		"sas - merge results from a series of dumped sas.xln files in different folders", 
		"results - merge map and frequency information into a final results file", 
		"vcf - lookup chr pos ref alt and return allele counts and frequencies", 
		"FilterDB - filter based on column names, thresholds and error messages", 
		"filterCNVs - calls FilterCalls to apply size/score/span limits", 
		"MeanLRR - compute mean LRRs for specific regions, then analyze or export to a text file", 
		"gwas.Qc - runs full QC protocol using PLINK", 
		"descriptive - summarize a phenotype file", 
		"phenoPrep - transform trait, reorder ids, and deal with outliers", 
		"bestTransformation - determine transformations that improve normality (i.e., minimzes skewness and kurtosis)", 
		"peakat - takes the first or last N lines of a file, or counts the lines", 
		"grep - filters a file line by line depending on the presence/absence of inclusion/exclusion criteria", 
		CNVTrioFilter.COMMAND_CNV_TRIO_CRF + CNVTrioFilter.COMMAND_CNV_TRIO_CRF_DESCRIPTION, 
		VCF.VCF_INIT, 
		VCF.VCF_COMMAND, 
		CNVFilter.COMMAND_CNV_FILTER_CRF, 
		CNVFilter.COMMAND_CNV_FILTER_DESCRIPTION, 
		CNVConcordance.COMMAND_CNV_CONCORDANCE, 
		CNVConcordance.COMMAND_CNV_CONCORDANCE_DESCRIPTION, 
		VCFOps.COMMAND_VCF_OPS_EXTRACT,
		VCFOps.COMMAND_VCF_EXTRACT_DESCRIPTION, 
		"replaceAll - replace Strings in a file using a list of replacements",
		};

	public static void run(String filename, Logger log) throws Elision {
		String temp;

		try {
			if (new File(filename).length() == 0) {
				temp = ext.rootOf(filename);
			} else {
				temp = Files.getReader(filename, false, true, log, true).readLine();
			}
			
			temp = temp.trim();
			
			log.report("Launching option '"+temp+"'");

			if (temp == null || temp.equals("")) {
				log.reportError("Below is a list of valid launch types:");
				log.reportError(Array.toStr(LAUNCH_TYPES, "\n"));
			} else if (temp.equals("lookup")) {
				LookupTable.fromParameters(filename, log);
			} else if (temp.equals("dummy")) {
				DummyDataset.createFromParameters(filename, log);
			} else if (temp.equals("counts")) {
				DummyDataset.createReverseFromParameters(filename, log);
			} else if (temp.equals("miss")) {
				MarkerQC.parseParameters(filename, log, true);
			} else if (temp.equals("indep")) {
				IndependentSNPs.selectFromParameters(filename, log);
			} else if (temp.equals("genes")) {
				MapGenesToSNPs.filter(filename, log);
			} else if (temp.equals("filterSNPs")) {
				FilterSNPsByRegion.fromParameters(filename, log);
			} else if (temp.equals("filterByLists")) {
				FilterByLists.fromParameters(filename, log);
			} else if (temp.equals("plink")) {
				CreateDatabaseFromPlink.createDatabase(filename, log);
			} else if (temp.equals("simpleM")) {
				CreateDatabaseFromPlink.createCountsMatrix(filename, log);
			} else if (temp.equals("score")) {
				CreateDatabaseFromPlink.createScore(filename, log);
			} else if (temp.equals("parse")) {
				GenParser.parseFile(filename, log);
			} else if (temp.equals("ucsc")) {
				UCSCtrack.describeFromParameters(filename, log);
			} else if (temp.equals("split")) {
				Files.splitFileFromParameters(filename, log);
			} else if (temp.equals("dir")) {
				Files.summarizeDirectoryFromParameters(filename, log);
			} else if (temp.equals("copy")) {
				Files.copySpecificFiles(filename, log);
			} else if (temp.equals("cat")) {
				Files.catFilesFromParameters(filename, log);
			} else if (temp.equals("rename")) {
				Files.renameFilesFromParameters(filename, log);
			} else if (temp.equals("db")) {
				CreateMarkerDatabase.createFromParameters(filename, log);
			} else if (temp.equals("merge")) {
				Files.mergeFromParameters(filename, log);
			} else if (temp.equals("mergeSNPs")) {
				Files.mergeSNPlistsFromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("trimFam")) {
				TrimFam.fromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("freq")) {
				Metal.calculateWeightedAlleleFrequencyFromParamters(filename, log);
			} else if (temp.equalsIgnoreCase("uniform")) {
				Metal.generateUniformsFromParamters(filename, log);
			} else if (temp.equalsIgnoreCase("metal")) {
				Metal.fromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("transform")) {
				Transformations.fromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("forest")) {
				ForestPlot.fromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("unique")) {
				Unique.fromParamters(filename, log);
			} else if (temp.equalsIgnoreCase("match")) {
				MatchSamples.matchFromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("meta")) {
				Conditional.metaAllRegions(filename, log);
			} else if (temp.equalsIgnoreCase("gwaf")) {
				CreateDatabaseFromPlink.gwafFromParamters(filename, log);
			} else if (temp.equalsIgnoreCase("sas")) {
				DumpSAS.mergeFromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("results")) {
				ResultsPackager.createFromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("vcf")) {
				Vcf.createFromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("heritability")) {
				Heritability.fromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("FilterDB")) {
				FilterDB.fromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("filterCNVs")) {
				FilterCalls.fromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("MeanLRR")) {
				MeanLRR.fromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("phenotype")) {
				SummarizePhenotype.filesFromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("descriptive")) {
				SummarizePhenotype.summarizeFromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("hitWindows")) {
				HitWindows.fromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("phenoPrep")) {
				PhenoPrep.fromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("bestTransformation")) {
				PhenoPrep.summarizeFromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("gzip")) {
				Zip.fromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("SuperNovo")) {
				SuperNovo.fromParameters(filename, log);
			} else if (temp.equalsIgnoreCase(CNVTrioFilter.COMMAND_CNV_TRIO_CRF)) {
				cnvTrio.fromParameters(filename, log);
			} else if (temp.equalsIgnoreCase(VCF.VCF_INIT)) {
				VCF.fromParameters(filename, log);
			} else if (temp.equalsIgnoreCase(CNVFilter.COMMAND_CNV_FILTER_CRF)) {
				CNVFilter.fromParameters(filename, log);
			}else if (temp.equalsIgnoreCase(CNVConcordance.COMMAND_CNV_CONCORDANCE)) {
				CNVConcordance.fromParameters(filename, log);
			} else if (temp.equalsIgnoreCase(VCFOps.COMMAND_VCF_OPS_EXTRACT)) {
				VCFOps.fromParameters(filename, VCFOps.UTILITY_TYPE.EXTRACT_SEGMENTS_ANNOTATION, log);
			}else if (temp.equalsIgnoreCase("transpose")) {
				Files.transposeFromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("gwas.Qc")) {
				Qc.fromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("peakat")) {
				widgets.peakAt.fromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("grep")) {
				widgets.grep.fromParameters(filename, log);
			} else if (temp.equalsIgnoreCase("replaceAll")) {
				Files.replaceAllFromParameters(filename, log);
			} else {
				log.reportError("Error - '"+temp+"' is an invalid launch type, options include:");
				log.reportError(Array.toStr(LAUNCH_TYPES, "\n"));
			}
		} catch (Exception e) {
			log.reportError(e.getMessage());
			log.reportException(e);
			ext.waitForResponse();
		}
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String filename = "transform.crf";
		boolean suppress = false;
		boolean create = false;

		String usage = "\n"+
		"park.crfDB requires 0-1 arguments\n"+
		"   (1) filename (i.e. "+filename+" (default)\n"+
		"   (2) suppress Windows stdin (i.e. -suppress (not the default)\n"+
		"   (3) create a blank file before running if it does not exist (i.e. -create (not the default)\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("-suppress")) {
				suppress = true;
				numArgs--;
			} else if (args[i].startsWith("-create")) {
				create = true;
				numArgs--;
			} else {
				filename = args[i];
				numArgs--;
			}
		}
		if (numArgs!=0) {
			System.err.println(usage);
			System.exit(1);
		}
		
		if (!Files.exists(filename)) {
			if (create) {
				new File(filename).createNewFile();
			} else {
				System.err.println("Error - file '"+filename+"' does not exist; add -create to the commandline to autogenerate");
				ext.waitForResponse();
				return;
			}
		}
				
		try {
			run(filename, new Logger(ext.rootOf(filename)+".log"));
		} catch (Exception e) {
			e.printStackTrace();
		}
		if (!suppress) {
			ext.waitForResponse("Press any key to close this window");
		}
	}
}
