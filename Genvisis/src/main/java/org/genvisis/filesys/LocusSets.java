package org.genvisis.filesys;

import java.util.ArrayList;

import org.genvisis.common.Logger;
import org.genvisis.filesys.LocusSet.TO_STRING_TYPE;

public class LocusSets<T extends Segment> {

  public enum UTILITY_TYPE {
    UNIONIZE;
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String[] filenames = null;
    String output = null;
    // String logfile = null;
    UTILITY_TYPE uType = UTILITY_TYPE.UNIONIZE;

    String usage = "\n" + "cnv.var.LocusSets requires 0-1 arguments\n";
    usage +=
        "   (1) comma delimited list of filenames with chr\tstart\tstop (i.e. files= (no default))\n"
            + "";
    usage += "   (2) output file (i.e. out= ( no default))\n" + "";

    usage += "   (3) utility type (i.e. utility=" + uType + " (default))\n" + "";
    usage += "    	 utility options are:";
    for (int i = 0; i < UTILITY_TYPE.values().length; i++) {
      usage += UTILITY_TYPE.values()[i];
    }

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("files=")) {
        filenames = arg.split("=")[1].split(",");
        numArgs--;
      } else if (arg.startsWith("out=")) {
        output = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("utility=")) {
        uType = UTILITY_TYPE.valueOf(arg.split("=")[1]);
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }
    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }

    Logger log = new Logger(null);
    switch (uType) {
      case UNIONIZE:
        unionize(filenames, output, log);
        break;
      default:
        System.err.println("invalid utility type " + uType);
        break;

    }

  }

  public static void unionize(String[] files, String output, Logger log) {

    ArrayList<LocusSet<Segment>> sets = new ArrayList<LocusSet<Segment>>();
    for (String file : files) {
      LocusSet<Segment> tmp = LocusSet.loadSegmentSetFromFile(file, 0, 1, 2, 0, true, true, 0, log);
      log.reportTimeInfo("Loaded " + tmp.getLoci().length + " segments from " + file);
      sets.add(tmp.mergeOverlapping());
    }

    LocusSets<Segment> lSets = new LocusSets<Segment>(sets, log);
    LocusSet<Segment> union = lSets.getUnion();
    long unionSize = union.getBpCovered();
    log.reportTimeInfo("Union contains " + unionSize + " bp total");
    for (int i = 0; i < sets.size(); i++) {
      long curCov = sets.get(i).getBpCovered();
      double corPerecent = (double) unionSize / curCov;
      log.reportTimeInfo("Regions in " + files[i] + " had " + curCov + " bp covered and "
          + corPerecent + " percent was in the union");
    }
    union.writeRegions(output, TO_STRING_TYPE.REGULAR, false, log);

  }

  private final ArrayList<LocusSet<T>> locusSets;

  private final Logger log;

  public LocusSets(ArrayList<LocusSet<T>> locusSets, Logger log) {
    super();
    this.locusSets = locusSets;
    this.log = log;
  }

  public LocusSet<Segment> getUnion() {
    int maxNum = -1;
    int maxIndex = 1;
    ArrayList<Segment> union = new ArrayList<Segment>();
    for (int i = 0; i < locusSets.size(); i++) {
      if (locusSets.get(i).getLoci().length > maxNum) {
        maxNum = locusSets.get(i).getLoci().length;
        maxIndex = i;
      }
    }
    for (int i = 0; i < locusSets.get(maxIndex).getLoci().length; i++) {
      T tmp = locusSets.get(maxIndex).getLoci()[i];
      ArrayList<T> finalOverLap = new ArrayList<T>();
      for (int j = 0; j < locusSets.size(); j++) {
        if (j != maxIndex) {
          T[] curOverlap = locusSets.get(j).getOverLappingLoci(tmp);
          if (curOverlap != null) {

            for (T element : curOverlap) {

              if (finalOverLap.size() == 0) {
                finalOverLap.add(element);

              } else {
                boolean found = false;
                for (int k2 = 0; k2 < finalOverLap.size(); k2++) {
                  if (finalOverLap.get(k2).equals(element)) {
                    found = true;
                  }
                }
                if (!found) {
                  finalOverLap.add(element);
                }
              }
            }
          }
        }
      }
      for (int j = 0; j < finalOverLap.size(); j++) {
        union.add(tmp.getIntersection(finalOverLap.get(j), log));
      }
    }
    LocusSet<Segment> finalUnion =
        new LocusSet<Segment>(union.toArray(new Segment[union.size()]), true, log) {

          /**
           * 
           */
          private static final long serialVersionUID = 1L;

        };

    finalUnion = finalUnion.mergeOverlapping();
    for (int i = 0; i < finalUnion.getLoci().length; i++) {
      for (int j = 0; j < locusSets.size(); j++) {

        if (locusSets.get(j).getOverlappingIndices(finalUnion.getLoci()[i]) == null) {

          String error = "Internal error developing union of the dataset";
          log.reportTimeError(error);
          throw new IllegalStateException(error);
        } else {

        }
      }
    }

    return finalUnion;
  }
}
