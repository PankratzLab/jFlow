package org.genvisis.one.JL;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;

import org.genvisis.cnv.annotation.AnnotationFileLoader.QUERY_ORDER;
import org.genvisis.cnv.annotation.AnnotationParser;
import org.genvisis.cnv.annotation.MarkerAnnotationLoader;
import org.genvisis.cnv.annotation.MarkerBlastAnnotation;
import org.genvisis.cnv.annotation.MarkerGCAnnotation;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.ExtProjectDataParser;
import org.genvisis.cnv.manage.ExtProjectDataParser.ProjectDataParserBuilder;
import org.genvisis.cnv.qc.MarkerMetrics;
import org.genvisis.common.Array;
import org.genvisis.common.Files;
import org.genvisis.common.ext;
import org.genvisis.filesys.Segment;
import org.genvisis.seq.manage.ReferenceGenome;
import org.genvisis.stats.Histogram.DynamicAveragingHistogram;
import org.genvisis.stats.Rscript.COLUMNS_MULTIPLOT;
import org.genvisis.stats.Rscript.PLOT_DEVICE;
import org.genvisis.stats.Rscript.RScatter;
import org.genvisis.stats.Rscript.RScatters;
import org.genvisis.stats.Rscript.SCATTER_TYPE;

public class GcLook {
  public enum CROSS_HYBE_FILTER {
                                 ALL(-1), ALIGN_25(25), ALIGN_30(30), ALIGN_35(35), ALIGN_45(45);

    private final int minTally;

    private CROSS_HYBE_FILTER(int buffer) {
      minTally = buffer;

    }

    public int getMinTally() {
      return minTally;
    }

  }
  public enum KB_BUFFER {
                         DESIGN(-1);
    // , BP_50(50), BP_100(100), BP_250(250), BP_500(500), BP_1000(1000), BP_10000(10000);

    private final int buffer;

    private KB_BUFFER(int buffer) {
      this.buffer = buffer;

    }

    public int getBuffer() {
      return buffer;
    }

  }

  private static final int[] CHRS = new int[] {-1, 26};
  // private static final int[] CHRS = new int[] { -1, 26, 23, 24 };
  private static final int[][] QC_GROUPINGS =
      new int[][] {{1}, {2, 3, 4}, {5, 6}, {7, 8, 9}, {10, 11, 12}, {23}};

  private static final String[] QC_TITLES =
      new String[] {"CallRate", "MeanClusterTheta", "DiffTheta", "SDClusterTheta", "MeanClusterR",
                    "LRR_SD"};

  private static final String GC_CONTENT = "GC_Content";

  public static void gcQCSummary(Project proj) {

    if (!Files.exists(proj.MARKER_METRICS_FILENAME.getValue())) {
      MarkerMetrics.fullQC(proj, null, null, false, 12);
    }

    MarkerSet markerSet = proj.getMarkerSet();
    String[] markerNames = proj.getMarkerNames();
    byte[] chrs = markerSet.getChrs();
    int[] pos = markerSet.getPositions();
    String dir = proj.PROJECT_DIRECTORY.getValue() + "GC_LOOK/";
    new File(dir).mkdirs();
    String out = dir + "gcLook.txt";

    ReferenceGenome referenceGenome =
        new ReferenceGenome(proj.REFERENCE_GENOME_FASTA_FILENAME.getValue(), proj.getLog());

    ProjectDataParserBuilder builder = new ExtProjectDataParser.ProjectDataParserBuilder();
    builder.separator("\t");
    builder.sampleBased(false);
    builder.requireAll(true);
    builder.dataKeyColumnName(MarkerMetrics.FULL_QC_HEADER[0]);

    ExtProjectDataParser parser;
    try {
      parser = builder.build(proj, proj.MARKER_METRICS_FILENAME.getValue());
      parser.loadData();
      ArrayList<String> titles = new ArrayList<String>();

      DynamicAveragingHistogram[][][] dHistograms =
          new DynamicAveragingHistogram[parser.getNumericData().length][CHRS.length][KB_BUFFER.values().length];
      for (int qcMetric = 0; qcMetric < parser.getNumericData().length; qcMetric++) {
        for (int chrIndex = 0; chrIndex < CHRS.length; chrIndex++) {
          for (int bufferIndex = 0; bufferIndex < KB_BUFFER.values().length; bufferIndex++) {
            String title = parser.getNumericDataTitles()[qcMetric] + "_chr"
                           + (CHRS[chrIndex] >= 0 ? CHRS[chrIndex] : "All") + "_"
                           + KB_BUFFER.values()[bufferIndex] + "_bp";
            titles.add(title);
            dHistograms[qcMetric][chrIndex][bufferIndex] = new DynamicAveragingHistogram(0, 1, 2);
            dHistograms[qcMetric][chrIndex][bufferIndex].setTitle(title);
          }
        }
      }
      if (!Files.exists(out)) {
        MarkerAnnotationLoader markerAnnotationLoader =
            new MarkerAnnotationLoader(proj, null, proj.BLAST_ANNOTATION_FILENAME.getValue(),
                                       proj.getMarkerSet(), true);
        markerAnnotationLoader.setReportEvery(500000);
        MarkerGCAnnotation[] gcAnnotations =
            MarkerGCAnnotation.initForMarkers(proj, markerNames,
                                              markerAnnotationLoader.getMarkerSet(),
                                              markerAnnotationLoader.getIndices());
        MarkerBlastAnnotation[] blastResults = MarkerBlastAnnotation.initForMarkers(markerNames);

        ArrayList<AnnotationParser[]> parsers = new ArrayList<AnnotationParser[]>();
        parsers.add(gcAnnotations);
        parsers.add(blastResults);

        markerAnnotationLoader.fillAnnotations(null, parsers, QUERY_ORDER.ONE_PER_IN_ORDER);

        try {

          PrintWriter writer = new PrintWriter(new FileWriter(out));
          PrintWriter writerSeparate = new PrintWriter(new FileWriter(out));

          writer.print(GC_CONTENT);
          for (int qcMetric = 0; qcMetric < parser.getNumericData().length; qcMetric++) {
            for (int chrIndex = 0; chrIndex < CHRS.length; chrIndex++) {
              for (int bufferIndex = 0; bufferIndex < KB_BUFFER.values().length; bufferIndex++) {
                writer.print("\t" + dHistograms[qcMetric][chrIndex][bufferIndex].getTitle());
              }
            }
          }

          writer.println();
          for (int i = 0; i < markerNames.length; i++) {
            if (i % 10000 == 0) {
              proj.getLog().reportTimeInfo("parsed " + i);
            }
            for (int bufferIndex = 0; bufferIndex < KB_BUFFER.values().length; bufferIndex++) {
              double gc = Double.NaN;
              if (KB_BUFFER.values()[bufferIndex] != KB_BUFFER.DESIGN) {
                Segment markerSegment = new Segment(chrs[i], pos[i], pos[i]);
                markerSegment =
                    markerSegment.getBufferedSegment(KB_BUFFER.values()[bufferIndex].getBuffer());
                gc = referenceGenome.getGCContentFor(markerSegment);
              } else {
                try {
                  gc = Double.parseDouble(gcAnnotations[i].getAnnotations()[0].getData());

                } catch (NumberFormatException nfe) {

                }
              }
              for (int qcMetric = 0; qcMetric < parser.getNumericData().length; qcMetric++) {
                for (int chrIndex = 0; chrIndex < CHRS.length; chrIndex++) {
                  if (CHRS[chrIndex] < 0 || CHRS[chrIndex] == chrs[i]) {
                    dHistograms[qcMetric][chrIndex][bufferIndex].addDataPair(gc,
                                                                             parser.getNumericData()[qcMetric][i]);
                  }
                }
              }
            }
          }
          for (int qcMetric = 0; qcMetric < parser.getNumericData().length; qcMetric++) {
            for (int chrIndex = 0; chrIndex < CHRS.length; chrIndex++) {
              for (int bufferIndex = 0; bufferIndex < KB_BUFFER.values().length; bufferIndex++) {
                dHistograms[qcMetric][chrIndex][bufferIndex].average();
              }
            }
          }

          for (int i = 0; i < dHistograms[0][0][0].getCounts().length; i++) {
            writer.print(dHistograms[0][0][0].getBins()[i]);
            for (int qcMetric = 0; qcMetric < parser.getNumericData().length; qcMetric++) {
              for (int chrIndex = 0; chrIndex < CHRS.length; chrIndex++) {
                for (int bufferIndex = 0; bufferIndex < KB_BUFFER.values().length; bufferIndex++) {
                  writer.print("\t"
                               + dHistograms[qcMetric][chrIndex][bufferIndex].getAverages()[i]);

                }
              }
            }
            writer.println();
          }

          writer.close();
          writerSeparate.close();
        } catch (Exception e) {
          proj.getLog().reportError("Error writing to " + out);
          proj.getLog().reportException(e);
        }
      }
      ArrayList<RScatter> rScatters = new ArrayList<RScatter>();
      for (int l = 0; l < QC_GROUPINGS.length; l++) {
        for (int chrIndex = 0; chrIndex < CHRS.length; chrIndex++) {
          for (int bufferIndex = 0; bufferIndex < KB_BUFFER.values().length; bufferIndex++) {
            String groupPlot = ext.rootOf(out, false) + "_" + QC_TITLES[l] + "_chr"
                               + (CHRS[chrIndex] >= 0 ? CHRS[chrIndex] : "All") + "_"
                               + KB_BUFFER.values()[bufferIndex] + "bp";
            String title = "n=" + Array.sum(dHistograms[0][chrIndex][bufferIndex].getCounts());
            ArrayList<String> ys = new ArrayList<String>();
            for (int k = 0; k < QC_GROUPINGS[l].length; k++) {
              ys.add(dHistograms[QC_GROUPINGS[l][k]][chrIndex][bufferIndex].getTitle());
            }
            String[] yColumns = ys.toArray(new String[ys.size()]);
            RScatter rScatterGroupAvg =
                new RScatter(out, groupPlot + ".rscript", ext.removeDirectoryInfo(groupPlot),
                             groupPlot + ".pdf", GC_CONTENT, yColumns, SCATTER_TYPE.POINT,
                             proj.getLog());
            rScatterGroupAvg.setxLabel(GC_CONTENT);
            rScatterGroupAvg.setyLabel(QC_TITLES[l]);
            rScatterGroupAvg.setOverWriteExisting(false);
            rScatterGroupAvg.setFontsize(6);
            rScatterGroupAvg.setTitle(title);
            rScatterGroupAvg.setyMin(0);
            rScatterGroupAvg.execute();
            rScatters.add(rScatterGroupAvg);
          }
        }
        RScatters rScattersAll =
            new RScatters(rScatters.toArray(new RScatter[rScatters.size()]), out + ".rscript",
                          out + ".pdf", COLUMNS_MULTIPLOT.COLUMNS_MULTIPLOT_1, PLOT_DEVICE.PDF,
                          proj.getLog());

        rScattersAll.execute();
      }
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    }
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = null;
    String usage = "\n" + "one.JL.GcLook requires 0-1 arguments\n"
                   + "   (1) project filename (i.e. proj=" + filename + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("proj=")) {
        filename = arg.split("=")[1];
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      Project proj = new Project(filename, false);
      gcQCSummary(proj);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
