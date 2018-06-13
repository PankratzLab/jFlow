package org.genvisis.common.parsing.impl;

import java.io.IOException;
import org.genvisis.CLI;
import org.genvisis.cnv.plots.PlotUtilities;
import org.genvisis.common.ext;
import org.genvisis.common.parsing.AliasedFileColumn;
import org.genvisis.common.parsing.Aliases;
import org.genvisis.common.parsing.Aliases.MultipleAliasStrategy;
import org.genvisis.common.parsing.ColumnFilter;
import org.genvisis.common.parsing.DoubleWrapperColumn;
import org.genvisis.common.parsing.FileColumn;
import org.genvisis.common.parsing.FileLink;
import org.genvisis.common.parsing.FileParser;
import org.genvisis.common.parsing.FileParserFactory;
import org.genvisis.common.parsing.FixedValueColumn;
import org.genvisis.common.parsing.IntegerFilter;
import org.genvisis.common.parsing.IntegerWrapperColumn;
import org.genvisis.common.parsing.StandardFileColumns;
import org.genvisis.stats.Maths.COMPARISON;

public class PlinkCNVParser {

  private String inputFileRoot = "";
  private String outFile = "";

  private int ncnvFloor = 5;

  private boolean hits = true;
  private boolean man = true;
  private boolean qq = true;

  public PlinkCNVParser(String inputFileRoot, String outputFile) {
    this.inputFileRoot = inputFileRoot;
    this.outFile = outputFile;
  }

  public PlinkCNVParser setNCNVFloor(int floor) {
    this.ncnvFloor = floor;
    return this;
  }

  public PlinkCNVParser runHitWindows(boolean run) {
    hits = run;
    return this;
  }

  public PlinkCNVParser runManhattanPlot(boolean run) {
    man = run;
    return this;
  }

  public PlinkCNVParser runQQPlot(boolean run) {
    qq = run;
    return this;
  }

  public void run() {
    FileColumn<Byte> chr = StandardFileColumns.chr("CHR");
    FileColumn<String> snp = StandardFileColumns.snp("SNP");
    FileColumn<Integer> pos = StandardFileColumns.pos("BP");

    FileColumn<Integer> NCNV = new IntegerWrapperColumn(new AliasedFileColumn("NCNV",
                                                                              new Aliases(new String[] {"NCNV"},
                                                                                          MultipleAliasStrategy.FAIL,
                                                                                          false)));

    FileColumn<Double> M0 = new DoubleWrapperColumn(new AliasedFileColumn("M0",
                                                                          new Aliases(new String[] {"M0"},
                                                                                      MultipleAliasStrategy.FAIL,
                                                                                      false)));
    FileColumn<Double> M1 = new DoubleWrapperColumn(new AliasedFileColumn("M1",
                                                                          new Aliases(new String[] {"M1"},
                                                                                      MultipleAliasStrategy.FAIL,
                                                                                      false)));
    ColumnFilter ncnvFilter = new IntegerFilter(NCNV, COMPARISON.GT, ncnvFloor);

    FileColumn<String> snp1 = StandardFileColumns.snp("SNP");

    FileColumn<Double> EMP1 = new DoubleWrapperColumn(new AliasedFileColumn("P",
                                                                            new Aliases(new String[] {"EMP1"},
                                                                                        MultipleAliasStrategy.FAIL,
                                                                                        false)));
    FileColumn<Double> EMP2 = new DoubleWrapperColumn(new AliasedFileColumn("EMP2",
                                                                            new Aliases(new String[] {"EMP2"},
                                                                                        MultipleAliasStrategy.FAIL,
                                                                                        false)));

    FileLink link = FileLink.setup(inputFileRoot + ".qt.summary.mperm").keys(snp1).values(EMP1,
                                                                                          EMP2);

    try (FileParser parser = FileParserFactory.setup(inputFileRoot + ".qt.summary", snp, chr, pos,
                                                     NCNV, EMP1, EMP2, M0, M1)
                                              .filter(ncnvFilter).link(link).build()) {

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

    } catch (IOException e) {
      e.printStackTrace();
    }

    // write new map file:
    chr = StandardFileColumns.chr("CHR");
    snp = StandardFileColumns.snp("SNP");
    FileColumn<String> zeroCol = new FixedValueColumn("", "0");
    pos = StandardFileColumns.pos("BP");

    try (FileParser parser = FileParserFactory.setup(outFile, chr, snp, zeroCol, pos).build()) {
      parser.parseToFile(outFile + ".map", "\t");
    } catch (IOException e) {
      e.printStackTrace();
    }

  }

  public static void main(String[] args) {
    String fileArg = "file";
    String fileDesc = "PLINK Permutation results file root (/path/to/directory/with/plinkRoot)";
    String outArg = "out";
    String outDesc = "Output filename";
    String hitsArg = "hits";
    String hitsDesc = "Run HitWindows";
    String manArg = "man";
    String manDesc = "Create a ManhattanPlot";
    String qqArg = "qq";
    String qqDesc = "Create a QQPlot";
    String ncnvArg = "ncnv";
    String ncnvDesc = "Minimum acceptable value in NCNV column";

    CLI cli = new CLI(PlinkCNVParser.class);
    cli.addArg(fileArg, fileDesc, true);
    cli.addArg(outArg, outDesc, false);
    cli.addArg(hitsArg, hitsDesc, false);
    cli.addArg(manArg, manDesc, false);
    cli.addArg(qqArg, qqDesc, false);
    cli.addArgWithDefault(ncnvArg, ncnvDesc, 5);

    cli.parseWithExit(args);

    String file = cli.get(fileArg);
    String out = cli.has(outArg) ? cli.get(outArg)
                                 : ext.rootOf(ext.rootOf(fileArg, false)) + ".results";
    boolean hits = cli.has(hitsArg) ? Boolean.parseBoolean(cli.get(hitsArg)) : true;
    boolean man = cli.has(manArg) ? Boolean.parseBoolean(cli.get(manArg)) : true;
    boolean qq = cli.has(qqArg) ? Boolean.parseBoolean(cli.get(qqArg)) : true;
    int floor = cli.getI(ncnvArg);

    new PlinkCNVParser(file, out).runHitWindows(hits).runManhattanPlot(man).runQQPlot(qq)
                                 .setNCNVFloor(floor).run();
  }

}
