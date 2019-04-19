package org.genvisis.cnv;

import java.io.File;
import java.util.Calendar;
import java.util.Date;

import org.pankratzlab.common.LauncherManifest;
import org.pankratzlab.common.PSF;

/**
 * Specialization of {@link LauncherManifest} for Genvisis.
 */
public final class GenvisisManifest {

  {
    LauncherManifest.setLaunchClass(Launch.class);
  }

  private static File validateFile(File file) {
    if (!file.exists() || !file.getAbsolutePath().endsWith(".jar")) {
      file = new File("../" + PSF.Java.GENVISIS);
    }
    return file;
  }

  /**
   * @return The manifest for the jar file used to launch this application
   */
  public static LauncherManifest loadGenvisisManifest() {
    return LauncherManifest.validateAndLoadManifest(GenvisisManifest::validateFile);
  }

  /**
   * @return An informative string about the current genvisis version
   */
  public static String getGenvisisInfo() {
    try {
      LauncherManifest manifest = loadGenvisisManifest();// until it always works
      return "Genvisis, " + manifest.getVersion().getNormalVersion() + "\n"
             + manifest.getCopyright() + "\n\n" + (new Date());
    } catch (Exception e) {
      return "Genvisis, v0.0.0\n(c)2009-" + Calendar.getInstance().get(Calendar.YEAR)
             + " Nathan Pankratz, GNU General Public License, v2\n\n" + (new Date());
    }
  }

  public static void main(String[] args) {
    GenvisisManifest.loadGenvisisManifest();
  }
}
