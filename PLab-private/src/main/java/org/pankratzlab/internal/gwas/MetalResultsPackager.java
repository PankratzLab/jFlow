package org.pankratzlab.internal.gwas;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.pankratzlab.common.CLI;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.parsing.AbstractColumnFilter;
import org.pankratzlab.common.parsing.AbstractFileParserFactory;
import org.pankratzlab.common.parsing.AliasedFileColumn;
import org.pankratzlab.common.parsing.ColumnFilter;
import org.pankratzlab.common.parsing.DataLine;
import org.pankratzlab.common.parsing.FileColumn;
import org.pankratzlab.common.parsing.FileLink;
import org.pankratzlab.common.parsing.FileParser;
import org.pankratzlab.common.parsing.FileParserFactory;
import org.pankratzlab.common.parsing.IntegerWrapperColumn;
import org.pankratzlab.common.parsing.NumberWrapperColumn;
import org.pankratzlab.common.parsing.RoundedDoubleWrapperColumn;
import org.pankratzlab.common.parsing.StandardFileColumns;

public class MetalResultsPackager {

  String mapFile;
  List<String> addlCols;

  String inputFile;
  int nonMissFilter;
  double mafFilter;
  String outputFile;

  Logger logger;

  public MetalResultsPackager(String inputFile, String outputFile, int nonMiss, double mafThresh,
                              String mapFile, String... addlMapCols) {
    this.inputFile = inputFile;
    this.outputFile = outputFile;
    this.nonMissFilter = nonMiss;
    this.mafFilter = mafThresh;

    this.mapFile = mapFile;
    this.addlCols = Arrays.asList(addlMapCols);

    this.logger = new Logger();
  }

  private void run() throws IOException {
    AliasedFileColumn colMarkerName = new AliasedFileColumn("MarkerName");
    AliasedFileColumn colDirection = new AliasedFileColumn("Direction");
    AliasedFileColumn colFreq1 = new AliasedFileColumn("Freq1");
    NumberWrapperColumn<Double> colMAF = new NumberWrapperColumn<Double>(new AliasedFileColumn("MAF",
                                                                                               "Freq1"),
                                                                         (s) -> {
                                                                           double d = Double.parseDouble(s);
                                                                           if (d > 0.5) {
                                                                             return 1 - d;
                                                                           }
                                                                           return d;
                                                                         });
    RoundedDoubleWrapperColumn colMAFRounded = new RoundedDoubleWrapperColumn(colMAF, 4);
    FileColumn<String> colElse = StandardFileColumns.allExcept("\t", colMarkerName, colDirection,
                                                               colFreq1);

    AliasedFileColumn colRsID = new AliasedFileColumn("MarkerName", "rs_id");
    IntegerWrapperColumn colChr = new IntegerWrapperColumn(new AliasedFileColumn("chr"));
    IntegerWrapperColumn colPos = new IntegerWrapperColumn(new AliasedFileColumn("pos",
                                                                                 "position"));
    FileColumn<?>[] valueColumns = new FileColumn[addlCols.size() + 2];
    valueColumns[0] = colChr;
    valueColumns[1] = colPos;
    for (int i = 0; i < addlCols.size(); i++) {
      valueColumns[2 + i] = new AliasedFileColumn(addlCols.get(i));
    }

    ColumnFilter filterDirection = new AbstractColumnFilter(colDirection) {

      @Override
      public boolean filter(DataLine values) {
        if (values.hasValid(colDirection)) {
          String direction = values.getUnsafe(colDirection);
          int count = 0;
          for (int i = 0; i < direction.length(); i++) {
            if (direction.charAt(i) != '?') {
              count++;
            }
          }
          return count >= nonMissFilter;
        }
        return false;
      }
    };

    ColumnFilter filterMAF = new AbstractColumnFilter(colMAF) {

      @Override
      public boolean filter(DataLine values) {
        if (values.has(colMAF)) {
          if (values.hasValid(colMAF)) {
            return values.getUnsafe(colMAF).doubleValue() > mafFilter;
          }
          return false;
        }
        return true;
      }
    };

    FileLink mapLink = FileLink.setup(mapFile).keys(colRsID).values(valueColumns);
    FileParser parser;
    AbstractFileParserFactory factory = FileParserFactory.setup(inputFile, colMarkerName,
                                                                colDirection, colElse)
                                                         .optionalColumns(colFreq1, colMAF,
                                                                          colMAFRounded)
                                                         .link(mapLink);
    if (nonMissFilter >= 0) {
      factory.filter(filterDirection);
    }
    if (mafFilter >= 0) {
      factory.filter(filterMAF);
    }
    logger.reportTime("Beginning parsing...");
    parser = factory.build();
    logger.reportTime("Loaded Map Info...");
    List<FileColumn<?>> outputColumnOrder = new ArrayList<>();
    outputColumnOrder.add(colMarkerName);
    outputColumnOrder.add(colChr);
    outputColumnOrder.add(colPos);
    outputColumnOrder.add(colFreq1);
    outputColumnOrder.add(colMAFRounded);
    outputColumnOrder.add(colDirection);
    outputColumnOrder.add(colElse);
    for (int i = 2; i < valueColumns.length; i++) {
      outputColumnOrder.add(valueColumns[i]);
    }

    logger.reportTime("Parsing results...");
    parser.parseToFile(outputFile, "\t", outputColumnOrder);
    logger.reportTime("Done!");
  }

  public static void main(String[] args) throws IOException {
    CLI cli = new CLI(MetalResultsPackager.class);

    cli.addArg("input", "Input File Path", true);
    cli.addArg("output", "Output File Path", true);
    cli.addArg("numStudies", "Threshold number of non-missing studies (-1 to disable)", false);
    cli.addArg("mafThreshold", "MAF threshold (if present; -1 to disable)", false);
    cli.addArg("mapFile", "Map file with chr/pos information", true);
    cli.addArg("addlMapCols", "Additional map file columns to add", false);

    cli.parseWithExit(args);

    String in = cli.get("input");
    String out = cli.get("output");
    int numS = cli.has("numStudies") ? cli.getI("numStudies") : -1;
    double maf = cli.has("mafThreshold") ? cli.getD("mafThreshold") : -1d;
    String map = cli.get("mapFile");
    String[] addlCols = cli.has("addlMapCols") ? cli.get("addlMapCols").split(",") : new String[0];

    MetalResultsPackager mrp = new MetalResultsPackager(in, out, numS, maf, map, addlCols);
    mrp.run();
  }

}
