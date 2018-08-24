package org.genvisis.one.spencer.files;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.stream.Collectors;
import org.genvisis.common.Files;
import org.genvisis.common.GenomicPosition;
import org.genvisis.common.PSF;
import com.google.common.collect.BiMap;
import com.google.common.collect.HashBiMap;

public class SnpChecking {

  public static void main(String[] args) {
    try (BufferedReader heightReader = Files.getAppropriateReader("F:\\OSv2_ImputationResults\\Mirabello Height\\OSv2_Mirabello_Height_MR.txt")) {
      heightReader.readLine();
      heightReader.readLine();
      BiMap<String, GenomicPosition> rsids = heightReader.lines()
                                                         .map(l -> l.split(PSF.Regex.GREEDY_WHITESPACE))
                                                         .collect(Collectors.toMap(c -> c[0],
                                                                                   c -> new GenomicPosition(Byte.parseByte(c[1]),
                                                                                                            Integer.parseInt(c[2])),
                                                                                   (a, b) -> {
                                                                                     throw new IllegalStateException();
                                                                                   },
                                                                                   HashBiMap::create));
      try (BufferedReader reader1000g = Files.getAppropriateReader("F:\\OSv2_ImputationResults\\Mirabello Height\\1000GP_Phase3_combined.legend.gz")) {
        reader1000g.readLine();
        while (reader1000g.ready()) {
          String[] colds = reader1000g.readLine().split(PSF.Regex.GREEDY_WHITESPACE);
          String rsid = colds[0].split(":")[0];
          GenomicPosition pos = new GenomicPosition(Byte.parseByte(colds[1]),
                                                    Integer.parseInt(colds[2]));
          if (rsids.containsKey(rsid)) {
            if (rsids.get(rsid).equals(pos)) {
              rsids.remove(rsid);
            } else {
              System.out.println("Mismatched Position for " + rsid + ": " + pos.toString() + ", "
                                 + rsids.get(rsid));
            }
          } else if (rsids.containsValue(pos)) {
            System.out.println("No rsid in 1000g for " + rsid + Arrays.toString(colds));
          }
        }

        System.out.println();
        System.out.println(rsids.size() + " unmatched rsids:");
        rsids.keySet().forEach(System.out::println);
        System.out.println();
      }
    } catch (FileNotFoundException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }
}
