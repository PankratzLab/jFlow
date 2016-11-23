package org.genvisis.one.JL.mica;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;

import org.genvisis.CLI;
import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.ExcelConverter;
import org.genvisis.common.HashVec;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.common.ExcelConverter.ExcelConversionParams;
import org.genvisis.seq.manage.GenotypeOps;
import org.genvisis.seq.manage.VCFOps;

import com.google.common.io.Files;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

public class MicaPlink {
	private static final String[] MPERM = new String[] { "EMP1", "EMP2" };
	private static final String[] ASSOC = new String[] { "C_A", "C_U" };
	private static final String[] LOC = new String[] { "BP", "A1", "A2" };

	private static void run(String vcf, String fam, String covars) {
		String outDir = ext.parseDirectoryOfFile(vcf) + "micaPlink/";
		new File(outDir).mkdirs();
		Logger log = new Logger(outDir + "log.log");
		int[] gqs = new int[] { -5, 20, 50 };
		String[] covarsH = Array.subArray(org.genvisis.common.Files.getHeaderOfFile(covars, log), 2);
		ArrayList<Result> results = new ArrayList<>();
		// , 20, 50
		for (int gq : gqs) {
			VCFFileReader reader = new VCFFileReader(new File(vcf), false);
			String root = outDir + VCFOps.getAppropriateRoot(vcf, true) + "_GQ_" + gq;
			String filtVcf = root + ".vcf";
			VariantContextWriter writer = VCFOps.initWriterWithHeader(reader, filtVcf, VCFOps.DEFUALT_WRITER_OPTIONS,
					new Logger());

			for (VariantContext vc : reader) {

				VariantContextBuilder builder = new VariantContextBuilder(vc);
				ArrayList<Genotype> genotypes = new ArrayList<>();
				for (Genotype g : vc.getGenotypes()) {
					GenotypeBuilder gb = new GenotypeBuilder(g);
					if (g.getGQ() < gq) {
						gb.alleles(GenotypeOps.getNoCall());

					}
					genotypes.add(gb.make());
				}

				GenotypesContext bc = GenotypesContext.create(genotypes);
				builder.genotypes(bc);
				writer.add(builder.make());
			}
			reader.close();
			writer.close();
			CmdLine.run("plink2 --vcf " + filtVcf + " --out " + root + " --make-bed --keep " + covars, outDir);

			try {
				Files.copy(new File(fam), new File(root + ".fam"));
			} catch (IOException e) {
				log.reportException(e);
			}
			String run = "plink2 --assoc counts mperm=100000 --freqx --bfile " + root + " --out " + root + " --covar "
					+ covars + " --covar-name " + Array.toStr(covarsH, ",");
			CmdLine.run(run, outDir);
			String t1 = "cat " + root + ".assoc.mperm | tr -s ' ' '\t'|cut -f 2- > " + root + ".assoc.mperm.tabs";
			String t2 = "cat " + root + ".assoc | tr -s ' ' '\t' |cut -f 2-> " + root + ".assoc.tabs";
			org.genvisis.common.Files.write(t1 + "\n" + t2, root + ".sh");
			org.genvisis.common.Files.chmod(root + ".sh");
			CmdLine.run(root + ".sh", outDir);

			Hashtable<String, String> mperm = HashVec.loadFileToHashString(root + ".assoc.mperm.tabs", "SNP", MPERM,
					"\t");

			Hashtable<String, String> assoc = HashVec.loadFileToHashString(root + ".assoc.tabs", "SNP", ASSOC, "\t");

			Hashtable<String, String> loc = HashVec.loadFileToHashString(root + ".assoc.tabs", "SNP", LOC, "\t");
			results.add(new Result(mperm, assoc, loc, gq));
		}
		String rsFile = ext.parseDirectoryOfFile(vcf) + "rsIdsOfInterest.txt";
		Hashtable<String, String> rsOfInterest = HashVec.loadFileToHashString(rsFile, 0, new int[] { 1, 2 }, "\t",
				true);
	
		HashSet<String> rsHaves = new HashSet<>();

		StringBuilder builder = new StringBuilder(
				"SNP\t" + Array.toStr(LOC) + "\tRSID\tRSID_OF_INTEREST\tPopFreqMax\tFunction");
		for (Result result : results) {
			builder.append("\t" + Array.toStr(Array.tagOn(ASSOC, "GQ_" + result.GQ, null)) + "\t"
					+ Array.toStr(Array.tagOn(MPERM, "GQ_" + result.GQ, null)));
		}
		builder.append("\n");
		for (String snp : results.get(0).assoc.keySet()) {
			String rs = parseRSID(snp);
			if (rsOfInterest.containsKey(rs)) {
				rsHaves.add(rs);
			}
			builder.append(snp + "\t" + results.get(0).loc.get(snp) + "\t" + rs + "\t" + rsOfInterest.containsKey(rs)
					+ "\t" + Array.toStr(getPopFunc(snp)));
			for (Result result : results) {
				builder.append("\t" + result.assoc.get(snp) + "\t" + result.mperm.get(snp));
			}
			builder.append("\n");
		}
		org.genvisis.common.Files.write(builder.toString(), outDir + "parsedResults.txt");
		ArrayList<ExcelConversionParams> excelFile = new ArrayList<ExcelConversionParams>();
		ExcelConversionParams parsedResults = new ExcelConversionParams(outDir + "parsedResults.txt", "\t",
				"ParsedResults");
		excelFile.add(parsedResults);
		ExcelConversionParams famTab = new ExcelConversionParams(fam, "\t", "FamFile");
		excelFile.add(famTab);

		ExcelConversionParams covarTab = new ExcelConversionParams(covars, "\t", "Covars");
		excelFile.add(covarTab);

		StringBuilder builderRS = new StringBuilder("RS_ID\tGenotyped\tCHR\tPOS\n");
		for (String rs : rsOfInterest.keySet()) {
			builderRS.append(rs + "\t" + rsHaves.contains(rs) + "\t" + rsOfInterest.get(rs) + "\n");
		}
		String rsSummary = ext.addToRoot(rsFile, ".have");
		org.genvisis.common.Files.write(builderRS.toString(), rsSummary);
		ExcelConversionParams rsTab = new ExcelConversionParams(rsSummary, "\t", "RsIdsOfInterests");
		excelFile.add(rsTab);
		
		ExcelConversionParams regions = new ExcelConversionParams(ext.parseDirectoryOfFile(vcf) + "targetedRegions.txt", "\t", "TargetRegions");
		excelFile.add(regions);

		String output = outDir + VCFOps.getAppropriateRoot(vcf, true) + "_summary.xlsx";
		new ExcelConverter(excelFile, output, log).convert(true);
	}

	private static String[] getPopFunc(String snpID) {
		String[] tmp = snpID.split("--");
		return new String[] { tmp[tmp.length - 1], tmp[tmp.length - 3] };
	}

	private static String parseRSID(String snpID) {
		if (snpID.startsWith("rs")) {
			return snpID.split("--")[0];
		}
		return "NA";
	}

	private static class Result {
		private Hashtable<String, String> mperm;
		private Hashtable<String, String> assoc;
		private Hashtable<String, String> loc;
		private int GQ;

		public Result(Hashtable<String, String> mperm, Hashtable<String, String> assoc, Hashtable<String, String> loc,
				int gQ) {
			super();
			this.mperm = mperm;
			this.assoc = assoc;
			this.loc = loc;
			GQ = gQ;
		}

	}

	public static void main(String[] args) {
		CLI c = new CLI(MicaPlink.class);

		c.addArgWithDefault("vcf", "vcf to annotate with default methods", "a.vcf");
		c.addArgWithDefault("fam", "fam file", "a.fam");
		c.addArgWithDefault("covars", "covar file", "a.covar");

		c.parseWithExit(args);

		run(c.get("vcf"), c.get("fam"), c.get("covars"));

	}

}
