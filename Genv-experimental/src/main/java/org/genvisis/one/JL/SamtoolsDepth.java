package org.genvisis.one.JL;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import org.genvisis.cnv.LocusSet;
import org.genvisis.seq.manage.BEDFileReader;
import org.genvisis.seq.manage.BEDFileReader.BEDFeatureSeg;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.CmdLine;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.HashVec;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;
import org.pankratzlab.shared.filesys.Positions;
import org.pankratzlab.shared.filesys.Segment;
import org.genvisis.seq.manage.BamOps;
import org.genvisis.seq.manage.BedOps;

/**
 * Use samtools to get bp resolution of DOC
 */
public class SamtoolsDepth {

  private static boolean runDepth(String bam, String out, Logger log) {
    List<String> command = new ArrayList<>();
    command.add("samtools");
    command.add("depth");
    command.add(bam);
    command.add(">");
    command.add(out);
    String c = out + ".sh";
    CmdLine.prepareBatchForCommandLine(out + ".sh", true, log, ArrayUtils.toStringArray(command));
    return CmdLine.runCommandWithFileChecks(new String[] {c}, "", new String[] {bam},
                                            new String[] {out}, true, false, false, log);

  }

  public static void main(String[] args) {
    String outDir = "/Volumes/Beta/data/Cushings/ATM/depth/";
    String inputDir = "/Volumes/Beta/data/Cushings/ATM/out/";
    String captureTargets = "/Volumes/Beta/data/Cushings/ATM/geneATM.txt";
    String output = outDir + "fullDepth.txt";

    new File(outDir).mkdirs();
    Logger log = new Logger(outDir + "depth.log");

    BedOps.verifyBedIndex(captureTargets, log);
    BEDFileReader readerCapture = new BEDFileReader(captureTargets, false);
    LocusSet<BEDFeatureSeg> set = readerCapture.loadAll(log);

    String[] bams = Files.listFullPaths(inputDir, ".bam");
    log.reportTimeInfo(bams.length + " bams");

    PrintWriter writer = Files.getAppropriateWriter(output);
    writer.println("SAMPLE\tEXON_NUMBER\tCHR\tPOS\tCOV");
    for (String bam : bams) {
      String out = outDir + ext.rootOf(bam) + ".depth";
      runDepth(bam, out, log);
      String[] data = HashVec.loadFileToStringArray(out, false, null, false);
      String sample = BamOps.getSampleName(bam, log);
      for (String d : data) {
        String[] parse = d.split("\t");
        Segment seg = new Segment(Positions.chromosomeNumber(parse[0]), Integer.parseInt(parse[1]),
                                  Integer.parseInt(parse[1]));
        int[] olap = set.getOverlappingIndices(seg);
        if (olap != null && olap.length > 0) {
          for (int en : olap) {
            writer.println(sample + "\t" + olap[0] + "\t" + seg.getChr() + "\t" + seg.getStart()
                           + "\t" + parse[2]);
          }

        }

      }

    }
  }

}
