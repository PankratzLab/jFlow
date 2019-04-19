package org.pankratzlab.common.bioinformatics;

import java.io.PrintWriter;
import java.util.Vector;

import org.pankratzlab.common.Aliases;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.common.filesys.GeneData;
import org.pankratzlab.common.filesys.GeneSet;
import org.pankratzlab.common.filesys.Positions;
import org.pankratzlab.common.filesys.Segment;

public class GeneDensityInRegion {

  public static void computeDensity(String region, int window) {
    PrintWriter writer;
    GeneSet geneSet;
    int[] loc;
    Vector<GeneData> inRegion;
    GeneData[] genes;
    GeneData regionAsGene;
    Vector<Segment> segments;
    int sum;
    int[][] exonBoundaries;

    // TODO not clear if this should be pulling from resources?
    geneSet = GeneSet.load(Aliases.getPathToFileInReferenceDirectory(GeneSet.REFSEQ_DB, true,
                                                                     new Logger()));
    genes = geneSet.getSet();
    loc = Positions.parseUCSClocation(region);
    regionAsGene = new GeneData("", new String[0], (byte) loc[0], true, (byte) 0, loc[1], loc[2],
                                new int[][] {}, (byte) 0, false);
    inRegion = new Vector<>();
    for (GeneData gene : genes) {
      if (gene.overlaps(regionAsGene)) {
        inRegion.add(gene);
      }
    }
    genes = new GeneSet(inRegion).getSet();
    try {
      writer = Files.openAppropriateWriter(ext.replaceAllWith(region, ":", "_") + ".xln");
      writer.println("Gene\tAssession #'s\tChr\tStart\tStop\tNumExons");
      for (GeneData gene : genes) {
        writer.println(gene.getGeneName() + "\t"
                       + ArrayUtils.toStr(gene.getNcbiAssessionNumbers(), "|") + "\t"
                       + gene.getChr() + "\t" + gene.getStart() + "\t" + gene.getStop() + "\t"
                       + gene.getExonBoundaries().length);
      }
      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing results");
      e.printStackTrace();
    }

    try {
      writer = Files.openAppropriateWriter(ext.replaceAllWith(region, ":", "_") + ".out");
      writer.println("Total region: " + ext.prettyUpDistance(regionAsGene.getSize(), 1));
      writer.println("Number of non-redundant RefSeq genes: " + genes.length);

      segments = new Vector<>();
      for (GeneData gene : genes) {
        segments.add(new Segment(gene.getStart(), gene.getStop()));
      }
      Segment.mergeOverlapsAndSort(segments);
      sum = 0;
      for (int i = 0; i < segments.size(); i++) {
        sum += segments.elementAt(i).getSize();
      }
      writer.println("Genes span: " + ext.prettyUpDistance(sum, 1) + " ("
                     + (int) ((double) sum / (double) regionAsGene.getSize() * 100) + "%)");

      segments = new Vector<>();
      for (GeneData gene : genes) {
        segments.add(new Segment(gene.getStart() - window, gene.getStop() + window));
      }
      Segment.mergeOverlapsAndSort(segments);
      sum = 0;
      for (int i = 0; i < segments.size(); i++) {
        sum += segments.elementAt(i).getSize();
      }
      writer.println("Including " + ext.prettyUpDistance(window, 0)
                     + " up- and down-stream, the genes span: " + ext.prettyUpDistance(sum, 1)
                     + " (" + (int) ((double) sum / (double) regionAsGene.getSize() * 100) + "%)");

      segments = new Vector<>();
      for (GeneData gene : genes) {
        exonBoundaries = gene.getExonBoundaries();
        for (int[] exonBoundarie : exonBoundaries) {
          segments.add(new Segment(exonBoundarie[0], exonBoundarie[1]));
        }
      }
      Segment.mergeOverlapsAndSort(segments);
      sum = 0;
      for (int i = 0; i < segments.size(); i++) {
        sum += segments.elementAt(i).getSize();
      }
      writer.println("Number of exons: " + segments.size());
      writer.println("Exons span: " + ext.prettyUpDistance(sum, 1) + " ("
                     + (int) ((double) sum / (double) regionAsGene.getSize() * 100) + "%)");

      for (GeneData gene : genes) {
        segments.add(new Segment(gene.getStart() - window, gene.getStart()));
        segments.add(new Segment(gene.getStop(), gene.getStop() + window));
      }
      Segment.mergeOverlapsAndSort(segments);
      sum = 0;
      for (int i = 0; i < segments.size(); i++) {
        sum += segments.elementAt(i).getSize();
      }
      writer.println("Including " + ext.prettyUpDistance(window, 0)
                     + " both 5' and 3', the total non-overlapping distance is: "
                     + ext.prettyUpDistance(sum, 1) + " ("
                     + (int) ((double) sum / (double) regionAsGene.getSize() * 100) + "%)");

      writer.close();
    } catch (Exception e) {
      System.err.println("Error writing results");
      e.printStackTrace();
    }

    Files.more(ext.replaceAllWith(region, ":", "_") + ".out");

  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String region = "chr2:221735000-240816800"; // Original 19.1 Mb
    // String region = "chr2:227241128-235934686"; // 8.7 Mb
    int window = 2000;

    String usage = "\n" + "bioinformatics.GeneDensityInRegion requires 0-1 arguments\n"
                   + "   (1) region (i.e. region=" + region + " (default))\n"
                   + "   (2) window in bases (i.e. win=" + window + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("region=")) {
        region = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("window=")) {
        window = Integer.parseInt(arg.split("=")[1]);
        numArgs--;
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }
    try {
      computeDensity(region, window);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }
}
