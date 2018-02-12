package org.genvisis.cnv.hmm;

import java.io.FileNotFoundException;
import org.genvisis.cnv.filesys.MarkerSetInfo;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.manage.ExtProjectDataParser;
import org.genvisis.cnv.manage.ExtProjectDataParser.ProjectDataParserBuilder;

/**
 * @author lane0212 Handles the pfb data for the hmm in reduced format
 */
public class PFB {

  private final Project proj;
  private final double[] pfbs;

  public PFB(PFB pfb) {
    proj = pfb.proj;
    pfbs = pfb.pfbs;
  }

  private PFB(Project proj, double[] pfbst) {
    super();
    this.proj = proj;
    pfbs = pfbst;
    MarkerSetInfo markerSet = proj.getMarkerSet();
    if (pfbs.length != markerSet.getMarkerNames().length) {
      String error = "Found " + pfbs.length + " pfb entries, but the project has "
                     + markerSet.getMarkerNames().length + " markers";
      proj.getLog().reportError(error);
      throw new IllegalArgumentException(error);
    } else {
      this.proj.getLog().reportTimeInfo("Loaded " + pfbst.length + " pfb entries");
    }
    int problems = 0;
    for (int i = 0; i < pfbs.length; i++) { // what PennCNV does
      if (!Double.isNaN(pfbs[i]) && pfbs[i] >= 0 && pfbs[i] <= 1) {
        if (pfbs[i] < 0.01) {
          pfbs[i] = 0.01;
        }
        if (pfbs[i] > .99) {
          pfbs[i] = .99;
        }
      } else if (!Double.isNaN(pfbs[i]) && pfbs[i] < 1) {
        if (!proj.ARRAY_TYPE.getValue().isCNOnly(markerSet.getMarkerNames()[i])) {
          problems++;
        }
      }
    }
    if (problems > 0) {
      proj.getLog().reportTimeWarning(problems + " markers " + " had a pfb value less than 1 ");
    }
  }

  public double[] getPfbs() {
    return pfbs;
  }

  public static PFB loadPFB(Project proj) {
    return loadPFB(proj, proj.CUSTOM_PFB_FILENAME.getValue());
  }

  public static PFB loadPFB(Project proj, String fullPathToPfb) {

    ProjectDataParserBuilder builder = new ProjectDataParserBuilder();
    builder.dataKeyColumnName("Name");
    builder.numericDataTitles(new String[] {"PFB"});
    builder.sampleBased(false);
    builder.requireAll(true);
    builder.treatAllNumeric(false);
    builder.verbose(false);

    try {
      ExtProjectDataParser extProjectDataParser = builder.build(proj, fullPathToPfb);
      extProjectDataParser.determineIndicesFromTitles();
      extProjectDataParser.loadData();
      double[] pfbs = extProjectDataParser.getNumericDataForTitle("PFB");
      return new PFB(proj, pfbs);
    } catch (FileNotFoundException e) {
      e.printStackTrace();
      proj.getLog().reportFileNotFound(fullPathToPfb);
      return null;
    }

  }

}
