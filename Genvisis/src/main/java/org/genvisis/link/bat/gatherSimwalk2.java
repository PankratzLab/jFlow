package org.genvisis.link.bat;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.StringTokenizer;
import java.util.Vector;

public class gatherSimwalk2 {
  public static double LOD_THRESHOLD = 2.0;

  public static void main(String[] args) throws IOException {
    try {
      if (args.length == 1) {
        new gatherSimwalk2(Integer.valueOf(args[0]).intValue(),
                           Integer.valueOf(args[0]).intValue());
      } else if (args.length == 2) {
        new gatherSimwalk2(Integer.valueOf(args[0]).intValue(),
                           Integer.valueOf(args[1]).intValue());
      } else {
        System.out.println("Expecting 1-2 arguments: how many models, or a range of models");
      }
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public gatherSimwalk2(int start, int stop) throws IOException {
    BufferedReader reader = null;
    PrintWriter writer = null;
    PrintWriter peaks = null;
    StringTokenizer st;
    String chrome, temp;
    Vector<String> vLines = new Vector<String>();
    Vector<String> mrkrNames = new Vector<String>();
    String vLine;
    int lineNumber;
    double pos;
    // double max_pos;
    double lod, max_lod;

    peaks = new PrintWriter(new FileWriter("linkage peaks.xls"));
    for (int i = start; i <= stop; i++) {
      // peaks.print("\tModel"+(i+1)+"\tcM");
      peaks.print("\tModel" + i);
    }
    peaks.println();
    for (int chromosome = 1; chromosome <= 23; chromosome++) {
      chrome =
          (Integer.valueOf(chromosome + "").intValue() < 10) ? "0" + chromosome : "" + chromosome;
      peaks.print(chromosome);
      writer = new PrintWriter(new FileWriter("summary" + chrome + ".xls"));
      mrkrNames.removeAllElements();
      vLines.removeAllElements();
      vLines.add("");
      for (int i = start; i <= stop; i++) {
        try {
          max_lod = -999;
          // max_pos = -1;
          if (chromosome < 23) {
            reader = new BufferedReader(new FileReader("chr" + chromosome + "/Model" + i + "/SCORE-"
                                                       + chrome + ".ALL"));
            do {
              temp = reader.readLine();
            } while (!temp.startsWith("Marker"));
            reader.readLine();
            reader.readLine();
            lineNumber = 0;
            temp = "";
            do {
              for (int j = 0; j < 11; j++) {
                st = new StringTokenizer(reader.readLine());
                pos = Double.valueOf(st.nextToken()).doubleValue();
                temp = st.nextToken();
                if (temp.equals("NaN")) {
                  lod = -998;
                } else {
                  lod = Double.valueOf(temp).doubleValue();
                }
                if (lod > max_lod) {
                  max_lod = lod;
                  // max_pos = pos;
                }
                vLine = vLines.elementAt(lineNumber);
                if (i == start) {
                  vLine = pos + "\t" + lod;
                  vLines.add("");
                } else {
                  vLine += "\t" + lod;
                }
                vLines.removeElementAt(lineNumber);
                vLines.insertElementAt(vLine, lineNumber);
                lineNumber++;
              }
              temp = reader.readLine();
              mrkrNames.add(temp);
            } while (!temp.startsWith(" "));
          }

          if (max_lod > LOD_THRESHOLD) {
            // peaks.print("\t"+max_lod+"\t"+max_pos);
            peaks.print("\t" + max_lod);
          } else {
            // peaks.print("\t"+max_lod+"\t");
            peaks.print("\t" + max_lod);
          }
          reader.close();
        } catch (Exception e) {
          System.err.println("Missed: chr" + chromosome + "/Model" + i + "/SCORE-" + chrome
                             + ".ALL");
        }
      }
      peaks.println();
      for (int i = start; i <= stop; i++) {
        writer.print("\tModel" + i);
      }
      writer.println();
      while (vLines.size() > 1) {
        for (int i = 0; i < 11; i++) {
          writer.println(vLines.elementAt(0));
          vLines.removeElementAt(0);
        }
        writer.println(mrkrNames.elementAt(0));
        mrkrNames.removeElementAt(0);
      }

      // for (int i=0; i<vLines.size()-1; i++) {
      // temp = vLines.elementAt(i);
      // writer.println(temp);
      // }
      writer.close();
    }
    peaks.close();
  }
}
