package org.genvisis.one.JL.beta;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import org.genvisis.common.Files;
import org.genvisis.gwas.windows.BasicHit;
import org.genvisis.gwas.windows.GeneralHitWindowDetector;
import org.genvisis.gwas.windows.HitWindow;

public class HitWindowTest {

  private static List<BasicHit> loadHits(String betaFile) {

    ArrayList<BasicHit> hits = new ArrayList<>();
    try {
      BufferedReader reader = Files.getAppropriateReader(betaFile);
      reader.readLine();
      while (reader.ready()) {
        String[] line = reader.readLine().trim().split("\t");
        try {
          BasicHit test = new BasicHit(line[0], Byte.parseByte(line[1]), Integer.parseInt(line[2]),
                                       Double.parseDouble(line[3]));
          hits.add(test);
        } catch (NumberFormatException nfe) {

        }

      }

      reader.close();

    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }

    return hits;
  }

  public static void main(String[] args) {
    String betaFile = "/Users/Kitty/Downloads/input.txt";
    List<BasicHit> hits = loadHits(betaFile);
    System.out.println(hits.size());
    // 0.00001 0.0001 500000

    GeneralHitWindowDetector<BasicHit> generalHitWindowDetector = new GeneralHitWindowDetector<BasicHit>(hits,
                                                                                                         500000,
                                                                                                         0.0001,
                                                                                                         0.00001);
    while (generalHitWindowDetector.hasNext()) {
      HitWindow<BasicHit> hit = generalHitWindowDetector.next();
      System.out.println(hit.toString());
    }

  }
}
