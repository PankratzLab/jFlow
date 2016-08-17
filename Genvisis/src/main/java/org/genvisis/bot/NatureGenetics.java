package org.genvisis.bot;

import java.util.HashSet;

import org.genvisis.common.Logger;
import org.genvisis.common.ext;

public class NatureGenetics {
  private final Logger log;
  private final AhkBot bot;

  private final String[] lines;
  private int index;

  private final HashSet<String> hash;

  public NatureGenetics(AhkBot bot, Logger log) {
    this.bot = bot;
    this.log = log;

    lines = ext.getClipboard().split("\n");
    index = 0;
    hash = new HashSet<String>();
  }

  public void next() {
    String[] line, location;

    try {
      if (index < lines.length) {
        line = lines[index].split("\t", -1);
        if (!hash.contains(line[0] + "\t" + line[2])) {
          bot.send(line[0], true);
          bot.send("\t\t");
          bot.send(line[2], true);
          bot.send("\t");
          bot.send(line[6], true);
          bot.send("\t");
          bot.send("\t");
          bot.send(line[4], true);
          bot.send("\t");
          bot.send(line[3], true);
          bot.send("\t");
          bot.send("\t");
          bot.send("\t");
          location = line[5].split(", ");
          bot.send(location[0], true);
          bot.send("\t");
          if (location.length == 3) {
            bot.send(location[1], true);
          }
          bot.send("\t");
          if (location.length == 3 && location[2].equals("USA")) {
            bot.send(" united s\n", false);
            bot.send("\t");

            // bot.sleep(200);
            // bot.send("n", false);
            // bot.sleep(200);
            // bot.send("[DownArrow]", false);
            // bot.sleep(200);
            // bot.send("[DownArrow]", false);
          } else {
            log.report(line[0] + " " + line[2] + ": " + location[1]);
          }
          bot.send("\t");
          bot.send("00000");
          bot.send("\t");
          bot.send("\t");
          bot.send("\t");
          bot.send("\t");
          bot.send("\t");
          bot.send("\t");
          if (line[7].toUpperCase().equals("X")) {
            bot.send(" ");
          }
          bot.send("\t");
          if (line[10].toUpperCase().equals("X")) {
            bot.send(" ");
          }
          bot.send("\t");
          if (line[12].toUpperCase().equals("X")) {
            bot.send(" ");
          }
          bot.send("\t");
          if (line[11].toUpperCase().equals("X")) {
            bot.send(" ");
          }
          bot.send("\t");
          bot.send("\t");
          if (line[13].toUpperCase().equals("X")) {
            bot.send(" ");
          }
          bot.send("\t");
          if ((line[8] + "" + line[9]).length() > 0) {
            bot.send(" ");
          }
          bot.send("\t");
          if (line[8].toUpperCase().equals("X")) {
            bot.send("recruited and assessed participants\n", true);
          }
          if (line[9].toUpperCase().equals("X")) {
            bot.send("generated genotyping data\n", true);
          }

          for (int i = 0; i < 14; i++) {
            bot.send("\t");
          }
          bot.send("\n");

          hash.add(line[0] + "\t" + line[2]);
        } else {
          log.report(line[0] + " " + line[2] + " is already inserted");
        }
        index++;
      } else {
        log.report("Error - reached the end of the author list");
      }
    } catch (Exception e) {
      log.report(e.getMessage());
    }
  }
}
