package org.genvisis.one.JL;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import org.genvisis.cnv.LocusSet;
import org.genvisis.cnv.LocusSet.TO_STRING_TYPE;
import org.genvisis.cnv.filesys.CNVariant;
import org.genvisis.cnv.filesys.MarkerDetailSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.filesys.CNVariant.CNVBuilder;
import org.pankratzlab.common.ArrayUtils;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;

/**
 * Creating output to prototype WBS in R
 */
public class WBS {

  public static void main(String[] args) {
    Project proj = new Project("/Users/Kitty/.genvisis/projects/LLFS.properties");

    System.out.println(proj.getSamples().length);
    Sample samp = proj.getFullSampleFromRandomAccessFile("10007060");
    MarkerDetailSet md = proj.getMarkerSet();
    System.out.println(ArrayUtils.toStr(ArrayUtils.subArray(samp.getLRRs(),
                                                            new int[] {100, 200, 5000})));
    // System.out.println(ArrayUtils.mean(md.getPositions()));
    // chr4:34,769,681-34,834,776

    String outDir = "/Volumes/Beta/data/LLFS/cnvsWBS/";
    String wbsCalls = "/Volumes/Beta/data/LLFS/cnvsWBS/10007060.wbs.cps.txt";
    new File(outDir).mkdirs();
    List<MarkerDetailSet.Marker> markers = md.markersAsList();
    Map<MarkerDetailSet.Marker, Integer> map = md.getMarkerIndexMap();
    PrintWriter writer = Files.getAppropriateWriter(outDir + samp.getSampleName() + ".wbs.gz");
    writer.println("MARKER\tCHR\tPOS\tLRR");
    for (MarkerDetailSet.Marker mark : markers) {
      if (mark.getChr() > 0 && mark.getChr() == 4) {
        writer.println(mark.getName() + "\t" + mark.getChr() + "\t" + mark.getPosition() + "\t"
                       + samp.getLRRs()[map.get(mark)]);
      }
    }

    try {
      BufferedReader reader = Files.getAppropriateReader(wbsCalls);
      reader.readLine();
      List<CNVariant> wbscnvs = new ArrayList<>();
      int num = 0;
      int alt = 3;
      while (reader.ready()) {
        String[] line = reader.readLine().trim().split("\t");
        CNVBuilder builder = new CNVBuilder();
        if (num == 0) {
          num++;
          builder.start(0);
        } else {
          builder.start(wbscnvs.get(wbscnvs.size() - 1).getStop());
        }
        builder.stop(Integer.parseInt(line[2]));
        builder.cn(alt);
        if (alt == 3) {
          alt = 1;
        } else {
          alt = 3;
        }
        builder.chr(Byte.parseByte(line[1]));
        builder.familyID(samp.getSampleName());
        builder.individualID(samp.getSampleName());
        wbscnvs.add(builder.build());
      }

      new LocusSet<>(wbscnvs, true, new Logger()).writeRegions(ext.addToRoot(wbsCalls, "wbsCNVs"),
                                                               TO_STRING_TYPE.REGULAR, true,
                                                               new Logger());
    } catch (FileNotFoundException e) {
      e.printStackTrace();
    } catch (IOException e) {
      e.printStackTrace();
    }
  }

}
