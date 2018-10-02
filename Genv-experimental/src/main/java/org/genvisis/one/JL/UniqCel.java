package org.genvisis.one.JL;

import java.util.ArrayList;
import org.pankratzlab.common.Files;
import org.pankratzlab.common.Logger;
import org.pankratzlab.common.ext;

public class UniqCel {

  public static void main(String[] args) {
    String dir = "/scratch.global/lanej/Affy6_1000g/affy6/";
    String p1 = dir + "phase1_2/";
    String[] celsP1 = Files.list(p1, ".CEL");

    String p2 = dir + "phase3/";
    String[] celsP2 = Files.list(p2, ".CEL");

    ArrayList<String> finalCels = new ArrayList<>();
    for (String element : celsP1) {
      finalCels.add(p1 + element);
    }
    for (String element : celsP2) {
      if (ext.indexOfStr(element, celsP1) >= 0) {
        System.out.println(element + " is duplicate");
        String dupRename = ext.addToRoot(element, "_2");
        if (!Files.exists(dupRename)) {
          Files.copyFileUsingFileChannels(p2 + element, p2 + dupRename, new Logger());
          finalCels.add(p2 + dupRename);
          System.out.println("renamed to " + dupRename);
        } else {
          throw new IllegalArgumentException();
        }
      } else {
        finalCels.add(p1 + element);

      }
    }
    Files.writeIterable(finalCels, dir + "g1000Cels.txt");

  }
}
