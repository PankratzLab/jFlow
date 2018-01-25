package org.genvisis.gwas.results.impl;

import java.io.IOException;
import java.util.Map;

import org.genvisis.CLI;
import org.genvisis.cnv.plots.ManhattanPlot;
import org.genvisis.cnv.plots.QQPlot;
import org.genvisis.common.Files;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.gwas.HitWindows;
import org.genvisis.gwas.results.files.AbstractColumnFilter;
import org.genvisis.gwas.results.files.AbstractFileParserFactory;
import org.genvisis.gwas.results.files.AliasedFileColumn;
import org.genvisis.gwas.results.files.Aliases;
import org.genvisis.gwas.results.files.Aliases.MultipleAliasStrategy;
import org.genvisis.gwas.results.files.ColumnFilter;
import org.genvisis.gwas.results.files.FileColumn;
import org.genvisis.gwas.results.files.FileLink;
import org.genvisis.gwas.results.files.FileParser;
import org.genvisis.gwas.results.files.FileParserFactory;
import org.genvisis.gwas.results.files.StandardFileColumns;

/**
 * Parse Plink association test results, either linear or logistic. <br />
 * <br />
 * Can also pull in information from a plink.frq file (A1, A2, and MAF). <br />
 * <br />
 * Also runs the results through {@link HitWindows} and generates a {@link ManhattanPlot} and
 * {@link QQPlot}.
 */
public class PlinkAssocParser {

	public void run(String resultsFile, String freqFile, String outFile, boolean hits, boolean man,
									boolean qq) {
		FileColumn<Integer> chr = StandardFileColumns.chr("CHR");
		FileColumn<String> snp = StandardFileColumns.snp("SNP");
		FileColumn<Integer> pos = StandardFileColumns.pos("BP");
		FileColumn<Double> beta = StandardFileColumns.beta("BETA");
		FileColumn<Double> se = StandardFileColumns.stdErr("SE");
		FileColumn<Double> p = StandardFileColumns.pVal("P");

		final FileColumn<String> test = new AliasedFileColumn("TEST",
																													new Aliases(new String[] {"TEST"},
																																			MultipleAliasStrategy.FAIL,
																																			false));

		ColumnFilter addFilter = new AbstractColumnFilter(test) {
			@Override
			public boolean filter(Map<FileColumn<?>, String> values) {
				return values.get(test).equals("ADD");
			}
		};

		FileColumn<String> snp1 = StandardFileColumns.snp("SNP");
		FileColumn<String> a1 = StandardFileColumns.a1("A1");
		FileColumn<String> a2 = StandardFileColumns.a2("A2");
		FileColumn<Double> maf = StandardFileColumns.alleleFreq("MAF");
		FileLink link = freqFile != null && Files.exists(freqFile)
																															 ? FileLink.setup(freqFile).keys(snp1)
																																				 .values(a1, a2, maf)
																															 : null;
		AbstractFileParserFactory factory = FileParserFactory.setup(resultsFile, snp, chr, pos, beta,
																																se, p)
																												 .filter(addFilter);
		if (link != null) {
			factory.link(link);
		}
		FileParser parser = factory.build();

		try {
			parser.parseToFile(outFile, "\t");

			if (hits) {
				HitWindows.main(new String[] {"file=" + outFile,
																			"out=" + ext.rootOf(outFile, false) + ".hits"});
			}

			if (man) {
				// manhattan plot
				ManhattanPlot mp = new ManhattanPlot(null);
				mp.loadFileAuto(outFile);
				mp.waitForData();
				mp.getManPan().setSize(800, 600);
				mp.screenshot(ext.rootOf(outFile, false) + "_manPlot.png");
			}

			if (qq) {
				// qqplot
				QQPlot qqPlot = QQPlot.loadPvals(new String[] {outFile}, "Q-Q Plot", false, true,
																				 false, -1, -1, false,
																				 Float.MAX_VALUE,
																				 new Logger());
				qqPlot.screenCap(ext.rootOf(outFile, false) + "_qqPlot.png");
			}

		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	public static void main(String[] args) {
		String fileArg = "file";
		String fileDesc = "Association results file (assoc.linear or assoc.logistic)";
		String freqArg = "freq";
		String freqDesc = "Frequency file";
		String outArg = "out";
		String outDesc = "Output filename";
		String hitsArg = "hits";
		String hitsDesc = "Run HitWindows";
		String manArg = "man";
		String manDesc = "Create a ManhattanPlot";
		String qqArg = "qq";
		String qqDesc = "Create a QQPlot";

		CLI cli = new CLI(PlinkAssocParser.class);
		cli.addArg(fileArg, fileDesc, true);
		cli.addArg(freqArg, freqDesc, false);
		cli.addArg(outArg, outDesc, false);
		cli.addArg(hitsArg, hitsDesc, false);
		cli.addArg(manArg, manDesc, false);
		cli.addArg(qqArg, qqDesc, false);

		cli.parseWithExit(args);

		String file = cli.get(fileArg);
		String freq = cli.has(freqArg) ? cli.get(freqArg) : null;
		String out = cli.has(outArg) ? cli.get(outArg)
																 : ext.rootOf(ext.rootOf(fileArg, false)) + ".results";
		boolean hits = cli.has(hitsArg) ? Boolean.parseBoolean(cli.get(hitsArg)) : true;
		boolean man = cli.has(manArg) ? Boolean.parseBoolean(cli.get(manArg)) : true;
		boolean qq = cli.has(qqArg) ? Boolean.parseBoolean(cli.get(qqArg)) : true;

		new PlinkAssocParser().run(file, freq, out, hits, man, qq);
	}
}
