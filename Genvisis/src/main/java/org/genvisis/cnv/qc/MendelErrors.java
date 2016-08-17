package org.genvisis.cnv.qc;

import java.io.FileWriter;
import java.io.PrintWriter;

import org.genvisis.cnv.filesys.MarkerData;
import org.genvisis.cnv.filesys.Pedigree;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.MDL;
import org.genvisis.common.Logger;

/**
 * @author Kitty Check for Mendelian errors in a pedigree, taken from
 *         http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#mendel
 */

public class MendelErrors {

  /**
   * @author Kitty Stores the Mendel error, if any, and where it occurred
   */
  public static class MendelErrorCheck {

    /**
     * 
     * -1 is no error<br>
     * Code Pat , Mat -> Offspring<br>
     * 
     * 1 AA , AA -> AB <br>
     * 2 BB , BB -> AB<br>
     * 
     * 3 BB , ** -> AA<br>
     * 4 ** , BB -> AA<br>
     * 5 BB , BB -> AA<br>
     * 
     * 6 AA , ** -> BB<br>
     * 7 ** , AA -> BB<br>
     * 8 AA , AA -> BB<br>
     * 
     * 9 ** , AA -> BB (X chromosome male offspring)<br>
     * 10 ** , BB -> AA (X chromosome male offspring)<br>
     * 
     * 11 **,AA -> !AA (Mitochondrial, maternal) <br>
     * 12 **,BB -> !BB (Mitochondrial, maternal)
     * 
     */
    private final int errorCode;
    private String error;
    private boolean faMendelError;
    private boolean moMendelError;

    public MendelErrorCheck(int errorCode) {
      super();
      this.errorCode = errorCode;
      switch (errorCode) {
        case -1:
          faMendelError = false;
          moMendelError = false;
          error = null;
          break;
        case 1:
          faMendelError = true;
          moMendelError = true;
          error = "AA , AA -> AB";
          break;
        case 2:
          faMendelError = true;
          moMendelError = true;
          error = "BB , BB -> AB";
          break;
        case 3:
          faMendelError = true;
          moMendelError = false;
          error = "BB , ** -> AA";
          break;
        case 4:
          faMendelError = false;
          moMendelError = true;
          error = "** , BB -> AA";
          break;
        case 5:
          faMendelError = true;
          moMendelError = true;
          error = "BB , BB -> AA";
          break;
        case 6:
          faMendelError = true;
          moMendelError = false;
          error = "AA , ** -> BB";
          break;
        case 7:
          faMendelError = false;
          moMendelError = true;
          error = "** , AA -> BB";
          break;
        case 8:
          faMendelError = true;
          moMendelError = true;
          error = "AA , AA -> BB";
          break;
        case 9:
          faMendelError = false;
          moMendelError = true;
          error = "** , AA -> BB (X chromosome male offspring)";
          break;
        case 10:
          faMendelError = false;
          moMendelError = true;
          error = "** , BB -> AA (X chromosome male offspring)";
          break;
        case 11:
          faMendelError = false;
          moMendelError = true;
          error = "**,AA -> !AA (Mitochondrial, maternal)";
          break;
        case 12:
          faMendelError = false;
          moMendelError = true;
          error = "**,BB -> !BB (Mitochondrial, maternal)";
          break;
        default:
          System.err.println("Warning - Unrecognized Mendel Error Code (" + errorCode + ")");
          faMendelError = false;
          moMendelError = false;
          error = null;
          break;
      }
    }

    public String getError() {
      return error;
    }

    public int getErrorCode() {
      return errorCode;
    }

    public boolean hasError() {
      return faMendelError || moMendelError;
    }

    public boolean hasFaMendelError() {
      return faMendelError;
    }

    public boolean hasMoMendelError() {
      return moMendelError;
    }
  }

  public static void detectMendelMarkers(Project proj) {
    Pedigree pedigree = proj.loadPedigree();
    Logger log;

    log = proj.getLog();
    if (pedigree == null) {
      return;
    } else {
      String output = proj.PROJECT_DIRECTORY.getValue() + "mendelErrorMarkers.txt";
      try {
        PrintWriter writer = new PrintWriter(new FileWriter(output));
        writer.println("MarkerName\tNumMendelErrors");
        boolean[] samplesToCheck = proj.getSamplesToInclude(null);
        MDL mdl = new MDL(proj, proj.getMarkerSet(), proj.getMarkerNames(), 2, 100);
        while (mdl.hasNext()) {
          MarkerData markerData = mdl.next();
          MendelErrorCheck[] mendelErrorChecks = Pedigree.PedigreeUtils.checkMendelErrors(pedigree,
              markerData, samplesToCheck, null, null, 0, log);
          int num = 0;
          for (MendelErrorCheck mendelErrorCheck : mendelErrorChecks) {
            if (mendelErrorCheck.getErrorCode() > 0) {
              num++;
            }
          }
          if (num > 0) {
            writer.println(markerData.getMarkerName() + "\t" + num);
            proj.getLog().reportTimeInfo(markerData.getMarkerName() + "\t" + num);
          }

        }

        writer.close();
      } catch (Exception e) {
        proj.getLog().reportError("Error writing to " + output);
        proj.getLog().reportException(e);
      }

    }

  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String filename = "C:/workspace/Genvisis/projects/Poynter_PCs.properties";

    String usage = "\n" + "cnv.qc.MendelErrors requires 0-1 arguments\n"
        + "   (1) filename (i.e. proj=" + filename + " (default))\n" + "";

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
      detectMendelMarkers(proj);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  private final byte chr;
  private final int offSex;

  private final byte offGenotype;

  private final byte moGenotype;

  private final byte faGenotype;

  /**
   * @param chr Only used for chr23 and male offspring
   * @param offSex Only used for chr23 and male offspring
   * @param offGenotype -1,0,1,2 genotype
   * @param moGenotype -1,0,1,2 genotype
   * @param faGenotype -1,0,1,2 genotype
   */
  public MendelErrors(byte chr, int offSex, byte offGenotype, byte faGenotype, byte moGenotype) {
    super();
    this.chr = chr;
    this.offSex = offSex;
    this.offGenotype = offGenotype;
    this.faGenotype = faGenotype;
    this.moGenotype = moGenotype;
  }

  public MendelErrorCheck checkMendelError() {
    if (offGenotype == -1) {
      return new MendelErrorCheck(-1);
    }
    if (moGenotype == -1 && faGenotype == -1) {
      return new MendelErrorCheck(-1);
    }
    if (chr < 23) {
      // note these are not in error code order
      if (faGenotype == 0 && moGenotype == 0 && offGenotype == 1) {
        return new MendelErrorCheck(1);
      }
      if (faGenotype == 2 && moGenotype == 2 && offGenotype == 1) {
        return new MendelErrorCheck(2);
      }
      if (faGenotype == 2 && moGenotype == 2 && offGenotype == 0) {
        return new MendelErrorCheck(5);
      }
      if (faGenotype == 0 && moGenotype == 0 && offGenotype == 2) {
        return new MendelErrorCheck(8);
      }
      if (faGenotype == 2 && offGenotype == 0) {
        return new MendelErrorCheck(3);
      }
      if (moGenotype == 2 && offGenotype == 0) {
        return new MendelErrorCheck(4);
      }
      if (faGenotype == 0 && offGenotype == 2) {
        return new MendelErrorCheck(6);
      }
      if (moGenotype == 0 && offGenotype == 2) {
        return new MendelErrorCheck(7);
      }

    } else if (chr == 23) {
      if (offSex == 1 && moGenotype == 0 && offGenotype == 2) {
        return new MendelErrorCheck(9);
      }
      if (offSex == 1 && moGenotype == 2 && offGenotype == 0) {
        return new MendelErrorCheck(10);
      }
    } else if (chr == 26) {

      if (moGenotype == 2 && offGenotype != 2) {
        return new MendelErrorCheck(11);
      }
      if (moGenotype == 0 && offGenotype != 0) {
        return new MendelErrorCheck(12);
      }
    }
    return new MendelErrorCheck(-1);
  }
}
