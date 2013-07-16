import java.io.*;

import cnv.analysis.FilterCalls;
import cnv.analysis.MeanLRR;

import link.Heritability;
import link.TrimFam;
import mining.Transformations;
import bioinformatics.*;
import common.*;
import parse.*;
import seq.Vcf;
import gwas.*;
import db.*;


public class Launch {
	public static final String[] LAUNCH_TYPES = {"hits", "dummy", "counts", "miss", "indep", "genes", "filterSNPs - filters SNP positions based on a set of regions with start and end positions", "filterByLists - filter unique IDs via a keeps file and a removes file", "plink", "simpleM", "score", "parse", "ucsc", "split", "cat", "db", "merge", "mergeSNPs", "trimFam", "freq - computes weighted allele frequency", "uniform - creates a hits control file where each file listed has the same column names, only with a different prefix", "metal", "transform", "forest", "unique", "dir", "copy", "meta", "gwaf", "sas - merge results from a series of dumped sas.xln files in different folders", "results - merge map and frequency information into a final results file", "vcf - lookup chr pos ref alt and return allele counts and frequencies", "FilterDB - filter based on column names, thresholds and error messages", "filterCNVs - calls FilterCalls to apply size/score/span limits", "MeanLRR - compute mean LRRs for specific regions, then analyze or export to a text file"};

	public static void run(String filename, Logger log) throws Elision {
		BufferedReader reader;
		String temp;

		try {
			reader = Files.getReader(filename, false, true, log, true);
			temp = reader.readLine();
			
			log.report("Launching option '"+temp+"'");

			if (temp == null || temp.equals("")) {
				log.reportError("Below is a list of valid launch types:");
				log.reportError(Array.toStr(LAUNCH_TYPES, "\n"));
			} else if (temp.equals("hits")) {
				LookupTable.fromParameters(filename, log);
			} else if (temp.equals("dummy")) {
				DummyDataset.createFromParameters(filename, log);
			} else if (temp.equals("counts")) {
				DummyDataset.createReverseFromParameters(filename, log);
			} else if (temp.equals("miss")) {
				MarkerQC.parseParameters(filename, log);
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
				Files.splitFileFromParamters(filename, log);
			} else if (temp.equals("dir")) {
				Files.summarizeDirectoryFromParameters(filename, log);
			} else if (temp.equals("copy")) {
				Files.copySpecificFiles(filename, log);
			} else if (temp.equals("cat")) {
				Files.catFilesFromParamters(filename, log);
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
			} else {
				log.reportError("Error - '"+temp+"' is an invalid launch type, options include:");
				log.reportError(Array.toStr(LAUNCH_TYPES, "\n"));
			}
		} catch (Exception e) {
			log.reportError(e.getMessage());
			log.reportException(e);
		}
	}

	public static void main(String[] args) throws IOException {
		int numArgs = args.length;
		String filename = "transform.crf";
		boolean suppress = false;

		String usage = "\n"+
		"park.crfDB requires 0-1 arguments\n"+
		"   (1) filename (i.e. "+filename+" (default)\n"+
		"   (2) suppress Windows stdin (i.e. -suppress (not the default)\n"+
		"";

		for (int i = 0; i<args.length; i++) {
			if (args[i].equals("-h")||args[i].equals("-help")||args[i].equals("/h")||args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("-suppress")) {
				suppress = true;
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
		try {
			run(filename, new Logger(ext.rootOf(filename)+".log"));
		} catch (Exception e) {
			e.printStackTrace();
		}
		if (!suppress && System.getProperty("os.name").startsWith("Windows")) {
			System.err.println();
			System.err.println("Press any key to close this window");
			new BufferedReader(new InputStreamReader(System.in)).readLine();
		}
	}
}
