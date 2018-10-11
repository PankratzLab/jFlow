package org.pankratzlab.gwas.parse;

import java.io.IOException;
import org.genvisis.cnv.plots.ManhattanPlot;
import org.genvisis.cnv.plots.PlotUtilities;
import org.genvisis.cnv.plots.QQPlot;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.parsing.AbstractColumnFilter;
import org.pankratzlab.common.parsing.AbstractFileParserFactory;
import org.pankratzlab.common.parsing.AliasedFileColumn;
import org.pankratzlab.common.parsing.Aliases;
import org.pankratzlab.common.parsing.ColumnFilter;
import org.pankratzlab.common.parsing.DataLine;
import org.pankratzlab.common.parsing.FileColumn;
import org.pankratzlab.common.parsing.FileLink;
import org.pankratzlab.common.parsing.FileParser;
import org.pankratzlab.common.parsing.FileParserFactory;
import org.pankratzlab.common.parsing.StandardFileColumns;
import org.pankratzlab.common.parsing.Aliases.MultipleAliasStrategy;
import org.pankratzlab.utils.gwas.windows.HitWindows;
import org.pankratzlab.common.CLI;

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
                  boolean qq, boolean af) {
    FileColumn<Byte> chr = StandardFileColumns.chr("CHR");
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
      public boolean filter(DataLine values) {
        return values.hasValid(test) && values.getString(test).equals("ADD");
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
    try (FileParser parser = factory.build()) {
      parser.parseToFile(outFile, "\t");

      if (hits) {
        PlotUtilities.runHitWindows(outFile);
      }

      if (man) {
        PlotUtilities.createManPlotScreenshot(outFile);
      }

      if (qq) {
        PlotUtilities.createQQPlotScreenshot(outFile);
      }

      if (af) {
        PlotUtilities.createAFPlotScreenshot(outFile);
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
    String afArg = "af";
    String afDesc = "Create an AFPlot";

    CLI cli = new CLI(PlinkAssocParser.class);
    cli.addArg(fileArg, fileDesc, true);
    cli.addArg(freqArg, freqDesc, false);
    cli.addArg(outArg, outDesc, false);
    cli.addArg(hitsArg, hitsDesc, false);
    cli.addArg(manArg, manDesc, false);
    cli.addArg(qqArg, qqDesc, false);
    cli.addArg(afArg, afDesc, false);

    cli.parseWithExit(args);

    String file = cli.get(fileArg);
    String freq = cli.has(freqArg) ? cli.get(freqArg) : null;
    String out = cli.has(outArg) ? cli.get(outArg)
                                 : ext.rootOf(ext.rootOf(fileArg, false)) + ".results";
    boolean hits = cli.has(hitsArg) ? Boolean.parseBoolean(cli.get(hitsArg)) : true;
    boolean man = cli.has(manArg) ? Boolean.parseBoolean(cli.get(manArg)) : true;
    boolean qq = cli.has(qqArg) ? Boolean.parseBoolean(cli.get(qqArg)) : true;
    boolean af = cli.has(afArg) ? Boolean.parseBoolean(cli.get(afArg)) : true;

    new PlinkAssocParser().run(file, freq, out, hits, man, qq, af);
  }
}
