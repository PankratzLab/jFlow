package org.genvisis.common;

import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Scanner;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JTextArea;

/**
 * @author lane0212, fairly specific class that hopefully can check and update genvisis
 *
 */
public class HttpUpdate {
  private static final String CHANGELOG =
                                        "https://github.com/npankrat/Genvisis/blob/master/CHANGELOG.md";
  public static final String REMOTE_JAR = "http://genvisis.org/genvisis.jar";
  public static final String UNDETERMINED_VERSION = "v.-1.-1.-1";

  public static RemoteJarStatus getRemoteJarVersion(String remoteJar, Logger log) {
    CHECK_STATUS status = CHECK_STATUS.OTHER_ERROR;
    String version = UNDETERMINED_VERSION;
    if (HttpDownloadUtility.canDownload(remoteJar, log)) {
      URL url;
      try {
        url = new URL(remoteJar);
        ZipInputStream zip = new ZipInputStream(url.openStream());
        ZipEntry entry;
        do {
          entry = zip.getNextEntry();
        } while (entry != null && !entry.getName().startsWith("META-INF/MANIFEST.MF"));
        if (entry != null) {
          Scanner sc = new Scanner(zip);
          while (sc.hasNextLine()) {
            String line = sc.nextLine();
            if (line.startsWith(CurrentManifest.IMPLEMENTATION_VERSION)) {
              version = line.replaceAll(CurrentManifest.IMPLEMENTATION_VERSION + ":", "")
                            .replaceAll(" ", "");
              status = CHECK_STATUS.OK;
              break;
            }
          }
          if (version.equals("undetermined")) {
            status = CHECK_STATUS.LINE_NOT_FOUND;
          }
          sc.close();
        } else {
          status = CHECK_STATUS.FILE_NOT_FOUND;
        }
        zip.close();
      } catch (MalformedURLException e) {
        status = CHECK_STATUS.HTTP_ERROR;
        // TODO Auto-generated catch block
        e.printStackTrace();
      } catch (IOException e) {
        status = CHECK_STATUS.FILE_NOT_FOUND;
        // TODO Auto-generated catch block
        e.printStackTrace();
      }

    } else {
      status = CHECK_STATUS.HTTP_ERROR;
    }

    return new RemoteJarStatus(new Version(version), remoteJar, status);
  }

  public static String checkGenvisisVersion(Logger log) {
    try {
      RemoteJarStatus remoteJarStatus = getRemoteJarVersion(REMOTE_JAR, log);
      CurrentManifest currentManifest = CurrentManifest.loadGenvisisManifest();
      if (currentManifest.getVersion().isLessThan(remoteJarStatus.getVersion())) {
        return "Your version of Genvisis (" + currentManifest.getVersion().getVersion()
               + ") is not up to date. Latest version is "
               + remoteJarStatus.getVersion().getVersion();
      } else if (currentManifest.getVersion().isGreaterThan(remoteJarStatus.getVersion())) {
        if (!remoteJarStatus.getVersion().getVersion().equals(UNDETERMINED_VERSION)) {
          String fun = "Looks like you are using a bleeding edge version of genvisis ("
                       + currentManifest.getVersion().getVersion() + ")\n";
          fun += "The latest released version at " + REMOTE_JAR + " is "
                 + remoteJarStatus.getVersion().getVersion() + ", good luck to ya";
          return fun;
        } else {
          return "Unable to determine genvisis version";

        }
      } else {
        return "Genvisis (" + currentManifest.getVersion().getVersion() + ") is up to date";
      }
    } catch (Exception e) {
      return "Unable to determine genvisis version";
    }
  }

  public enum CHECK_STATUS {
                            OK, HTTP_ERROR, FILE_NOT_FOUND, LINE_NOT_FOUND, OTHER_ERROR;
  }

  public static class Version {
    private final String version;
    private final int major;
    private final int minor;
    private final int patch;

    public Version(String version) {
      super();
      this.version = version;
      int[] tmp = parse();
      major = tmp[0];
      minor = tmp[1];
      patch = tmp[2];
    }

    public String getVersion() {
      return version;
    }

    public int getMajor() {
      return major;
    }

    public int getMinor() {
      return minor;
    }

    public int getPatch() {
      return patch;
    }

    public boolean isLessThan(Version other) {
      if (major < other.getMajor()) {
        return true;
      } else if (minor < other.getMinor()) {
        return true;
      } else if (patch < other.getPatch()) {
        return true;
      }
      return false;
    }

    public boolean isGreaterThan(Version other) {
      if (major > other.getMajor()) {
        return true;
      } else if (minor > other.getMinor()) {
        return true;
      } else if (patch > other.getPatch()) {
        return true;
      }
      return false;
    }

    private int[] parse() {
      int[] v = Array.toIntArray(version.replaceAll("v", "").split("\\."));
      return v;

    }

  }

  public static class RemoteJarStatus {
    private final Version version;
    private final CHECK_STATUS status;
    private final String jarChecked;

    public RemoteJarStatus(Version version, String jarChecked, CHECK_STATUS status) {
      super();
      this.version = version;
      this.status = status;
      this.jarChecked = jarChecked;
    }

    public Version getVersion() {
      return version;
    }

    public CHECK_STATUS getStatus() {
      return status;
    }

    public String getJarChecked() {
      return jarChecked;
    }

  }

  public static class UpdateInfo extends JFrame {

    /**
     *
     */
    private static final long serialVersionUID = 1L;
    private JTextArea infoPane;
    private JButton ok;
    private JButton cancel;
    private JPanel pan1;
    private JPanel pan2;
    private final RemoteJarStatus remoteJarStatus;
    private final CurrentManifest manifest;
    private final String newFileDir;
    private final String newJarFile;
    private boolean upToDate;
    private final Logger log;

    public UpdateInfo(String directoryToSave, RemoteJarStatus remoteJarStatus, Logger log) {
      newFileDir = directoryToSave;
      manifest = CurrentManifest.loadGenvisisManifest();
      new File(newFileDir).mkdirs();
      newJarFile = newFileDir
                   + ext.addToRoot(PSF.Java.GENVISIS, remoteJarStatus.getVersion().getVersion());
      this.remoteJarStatus = remoteJarStatus;
      initComponents();
      this.log = log;
    }

    private String generateMessage() {
      StringBuilder builder = new StringBuilder();
      if (manifest.getVersion().isLessThan(remoteJarStatus.getVersion())) {
        builder.append("The current version of Genvisis is " + manifest.getVersion().getVersion()
                       + "\n");
        builder.append("You can check out the latest change log at " + CHANGELOG + "\n");
        builder.append("Would you like to update to " + remoteJarStatus.getVersion().getVersion()
                       + "\n");
        builder.append("The latest version will be saved to " + newJarFile);
        upToDate = false;
      } else {
        builder.append("Genvisis is up to date " + remoteJarStatus.getVersion().getVersion()
                       + "\n");
        builder.append("If you think this is an error, you can just download "
                       + remoteJarStatus.getJarChecked() + "\n");
        upToDate = true;
      }
      return builder.toString();
    }

    private void initComponents() {

      setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);
      setTitle("New Update Found");
      pan1 = new JPanel();
      pan1.setLayout(new BorderLayout());

      pan2 = new JPanel();
      pan2.setLayout(new FlowLayout());

      infoPane = new JTextArea();
      infoPane.setEditable(false);
      infoPane.setText(generateMessage());

      if (!upToDate) {
        ok = new JButton("Update");

        cancel = new JButton("Cancel");
      } else {
        ok = new JButton("OK");

        cancel = new JButton("");
        cancel.setVisible(false);
      }
      ok.addActionListener(new ActionListener() {

        @Override
        public void actionPerformed(ActionEvent e) {
          update();
        }
      });
      cancel.addActionListener(new ActionListener() {

        @Override
        public void actionPerformed(ActionEvent e) {
          UpdateInfo.this.dispose();
        }
      });
      pan2.add(ok);
      pan2.add(cancel);
      pan1.add(pan2, BorderLayout.SOUTH);
      pan1.add(infoPane, BorderLayout.CENTER);

      this.add(pan1);
      pack();
      setVisible(true);
    }

    private void update()

    {
      boolean success = false;
      if (!upToDate) {
        if (HttpDownloadUtility.canDownload(remoteJarStatus.getJarChecked(), log)) {
          try {
            HttpDownloadUtility.downloadFile(remoteJarStatus.getJarChecked(), newJarFile, true,
                                             log);
            success = Files.exists(newJarFile);
            if (success) {
              log.reportTimeInfo("New version of Genvisis can be found at " + newJarFile);
            }
          } catch (IOException e) {
            log.reportTimeError("Could not download " + remoteJarStatus.getJarChecked());
            e.printStackTrace();
          }
        } else {
          log.reportTimeError("Could not download " + remoteJarStatus.getJarChecked());
        }
      }
      dispose();
    }
  }

  /**
   * @param remoteJar the full url
   * @param dir directory to save any updates
   * @param log
   */
  public static void update(String remoteJar, String dir, Logger log) {
    RemoteJarStatus remoteJarStatus = getRemoteJarVersion(remoteJar, new Logger());
    if (remoteJarStatus.getStatus() == CHECK_STATUS.OK) {
      new UpdateInfo(dir, remoteJarStatus, new Logger());
    } else {
      log.reportTimeError("Unable to fetch updates");
    }
  }

  public static void main(String[] args) {
    long time = System.currentTimeMillis();
    RemoteJarStatus remoteJarStatus = getRemoteJarVersion("http://genvisis.org/genvisis_dev.jar",
                                                          new Logger());
    new Logger().reportTimeElapsed(time);
    if (remoteJarStatus.getStatus() == CHECK_STATUS.OK) {

    }
  }
}
