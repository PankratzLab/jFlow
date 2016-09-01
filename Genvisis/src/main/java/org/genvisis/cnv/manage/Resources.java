package org.genvisis.cnv.manage;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.compress.compressors.gzip.GzipCompressorInputStream;
import org.genvisis.common.Files;
import org.genvisis.common.HttpDownloadUtility;
import org.genvisis.common.Logger;
import org.genvisis.seq.manage.VCFOps;

/**
 * All resources are ultimately a string path. Need to be able to ask for a resource via API, or by
 * string Some resources are nested within other resources
 */
public class Resources {

  public static final String DEFAULT_URL = "http://genvisis.org/rsrc/";
  public static final String DEFAULT_LOCAL_DIR = "resources" + File.separator;
  public static final String BIN_DIR = "bin";
  public static final String GENOME_DIR = "Genome";

  // TODO, a big TODO
  // need to add web-based download, and local file structure
  // could probably do this like project properties...

  public static MiniMac miniMac(Logger log) {
    return new MiniMac(log);
  }

  public static class MiniMac extends AbstractResourceFactory {
    public MiniMac(Logger log) {
      super(BIN_DIR + "/Minimac3", log);
    }

    public Resource getMiniMac3() {
      return getTarGzResource("Minimac3.v1.0.14.tar.gz");
    }

    @Override
    public List<Resource> getResources() {
      List<Resource> resources = new ArrayList<Resource>();
      resources.add(getMiniMac3());
      return resources;
    }
  }

  public static Shapeit shapeit(Logger log) {
    return new Shapeit(log);
  }

  public static class Shapeit extends AbstractResourceFactory {
    public Shapeit(Logger log) {
      super(DEFAULT_LOCAL_DIR + BIN_DIR + File.separator + "shapeit" + File.separator
            + "shapeit.tar.gz",
            "https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz",
            log);
    }

    public Resource getShapeit() {
      return getResource("");
    }

    @Override
    public List<Resource> getResources() {
      List<Resource> resources = new ArrayList<Resource>();
      resources.add(getShapeit());
      return resources;
    }
  }

  public static class Chr extends AbstractResourceFactory {
    private String build;

    public Chr(GENOME_BUILD build, Logger log) {
      super(GENOME_DIR + "/" + build.getBuild() + "/chr", log);
      this.build = build.getBuild();
    }

    public Resource getGeneticMap(CHROMASOME c) {
      if (CHROMASOME.CX_PAR.equals(c)) {
        getResource(getPath("genetic_map_", c.toString() + "2", ".txt.gz"));
        return getResource(getPath("genetic_map_", c.getLabel() + "1", ".txt.gz"));
      }
      return getResource(getPath("genetic_map_", c.getLabel(), ".txt.gz"));
    }

    public Resource getG1Kphase3v5RefPanel(CHROMASOME c) {
      return getResource(getPath("1000genomes_ref_panel_Phase3v5_", c.getLabel(), ".m3vcf.gz"));
    }

    private String getPath(String prefix, String chromasome, String suffix) {
      return prefix + build + "_chr" + chromasome + suffix;
    }

    @Override
    public List<Resource> getResources() {
      List<Resource> resources = new ArrayList<Resource>();
      for (CHROMASOME c : CHROMASOME.values()) {
        resources.add(getGeneticMap(c));
        resources.add(getG1Kphase3v5RefPanel(c));
      }
      return resources;
    }
  }

  public static Genome genome(GENOME_BUILD build, Logger log) {
    return new Genome(build, log);
  }

  public static class Genome extends AbstractResourceFactory {
    private GENOME_BUILD build;

    public Genome(GENOME_BUILD build, Logger log) {
      super(GENOME_DIR, log);
      this.build = build;
    }

    public Resource getModelBase() {
      return getResource(getPath() + "_gc5Base.txt");
    }

    public Resource getDBSNP() {
      return getVCFResource(getPath() + "_dbSnp147.vcf.gz");
    }

    private String getPath() {
      String b = build.getBuild();
      return new StringBuilder().append(b).append("/").append(b).toString();
    }

    public Chr chr() {
      return new Chr(build, log());
    }

    @Override
    public List<Resource> getResources() {
      List<Resource> resources = new ArrayList<Resource>();
      resources.add(getModelBase());
      resources.add(getDBSNP());
      return resources;
    }
  }

  public static MitoCN mitoCN(Logger log) {
    return new MitoCN(log);
  }

  public static class MitoCN extends AbstractResourceFactory {
    public MitoCN(Logger log) {
      super("MitoCN", log);
    }

    public Resource getWhiteWBC() {
      return getResource("Whites_WBC_TOTAL_SingleSNPmatched.final.beta");
    }

    public Resource getBlackWBC() {
      return getResource("Blacks_WBC_TOTAL_SingleSNPmatched.final.beta");
    }

    public Resource getTotalWBC() {
      return getResource("WBC_TOTAL_SingleSNPmatched.final.beta");
    }

    @Override
    public List<Resource> getResources() {
      List<Resource> resources = new ArrayList<Resource>();
      resources.add(getWhiteWBC());
      resources.add(getBlackWBC());
      resources.add(getTotalWBC());
      return resources;
    }
  }

  /**
   * Illumina Bundle TODO
   */

  public static AffySnp6 affy(Logger log) {
    return new AffySnp6(log);
  }

  public static class AffySnp6 extends AbstractResourceFactory {
    public AffySnp6(Logger log) {
      super("Arrays/AffySnp6", log);
    }

    public Resource getMarkerPositions(GENOME_BUILD build) {
      return getResource(build.getBuild() + "_markerPositions.txt");
    }

    public Resource getHMM() {
      return getResource("affygw6.hmm");
    }

    public Resource getABLookup() {
      return getResource("AB_lookup.dat");
    }

    @Override
    public List<Resource> getResources() {
      List<Resource> resources = new ArrayList<Resource>();
      resources.add(getHMM());
      resources.add(getABLookup());

      for (GENOME_BUILD build : GENOME_BUILD.values()) {
        resources.add(getMarkerPositions(build));
      }
      return resources;
    }
  }

  private abstract static class AbstractResourceFactory implements ResourceFactory {
    private final String localPath;
    private final String remotePath;
    private final Logger log;

    public AbstractResourceFactory(String subPath, Logger log) {
      this(DEFAULT_LOCAL_DIR + subPath + File.separator, DEFAULT_URL + subPath + "/", log);
    }

    public AbstractResourceFactory(String localPath, String url, Logger log) {
      this.localPath = localPath;
      this.remotePath = url;
      this.log = log;
    }

    protected Logger log() {
      return log;
    }

    protected Resource getTarGzResource(String rsrc) {
      return new TarGzResource(rsrc, localPath, remotePath, log);
    }

    protected Resource getVCFResource(String rsrc) {
      return new VCFResource(rsrc, localPath, remotePath, log);
    }

    protected Resource getResource(String rsrc) {
      return new DefaultResource(rsrc, localPath, remotePath, log);
    }
  }

  private static interface ResourceFactory {
    List<Resource> getResources();
  }

  public static class TarGzResource extends AbstractResource {
    private String unzippedDir;
    private String unzippedPath;

    public TarGzResource(String resourceName, String path, String url, Logger log) {
      super(resourceName, path, url, log);
      unzippedDir = resourceName.substring(0, resourceName.lastIndexOf(".tar.gz",
                                                                       resourceName.length()));
      unzippedDir += File.separator;
      unzippedPath = path + unzippedDir;
    }

    @Override
    public String get() {
      if (isLocallyAvailable(unzippedPath)) {
        return unzippedPath;
      }

      String path = super.get();
      if (path != null) {
        return extractTarGz(path, unzippedPath);
      }

      return null;
    }

    private String extractTarGz(String targzPath, String destination) {
      final int BUFFER = 2048;

      try {
        FileInputStream fin = new FileInputStream(targzPath);
        BufferedInputStream in = new BufferedInputStream(fin);
        GzipCompressorInputStream gzIn = new GzipCompressorInputStream(in);
        TarArchiveInputStream tarIn = new TarArchiveInputStream(gzIn);

        TarArchiveEntry entry = null;

        /** Read the tar entries using the getNextEntry method **/

        while ((entry = (TarArchiveEntry) tarIn.getNextEntry()) != null) {

          /** If the entry is a directory, create the directory. **/

          if (entry.isDirectory()) {

            File f = new File(destination + entry.getName());
            f.mkdirs();
          }
          /**
           * If the entry is a file,write the decompressed file to the disk and close destination
           * stream.
           **/
          else {
            int count;
            byte[] data = new byte[BUFFER];

            FileOutputStream fos = new FileOutputStream(destination + entry.getName());
            BufferedOutputStream dest = new BufferedOutputStream(fos, BUFFER);

            while ((count = tarIn.read(data, 0, BUFFER)) != -1) {
              dest.write(data, 0, count);
            }
            dest.close();
          }
        }

        /** Close the input stream **/
        tarIn.close();
        new File(targzPath).delete();
        return destination;
      } catch (Exception e) {
        log().reportError("Failed to extract: " + targzPath);
        log().reportException(e);
        return null;
      }
    }
  }

  public static class VCFResource extends AbstractResource {
    public VCFResource(String resourceName, String path, String url, Logger log) {
      super(resourceName, path, url, log);
    }

    @Override
    public String get() {
      String path = super.get();
      if (path != null) {
        String indexPath = VCFOps.getIndex(getName());
        if (get(indexPath) == null) {
          log().reportError("Warning: no index found for vcf file: " + path);
        }
      }
      return path;
    }
  }

  public static class DefaultResource extends AbstractResource {
    public DefaultResource(String resourceName, String path, String url, Logger log) {
      super(resourceName, path, url, log);
    }
  }

  public abstract static class AbstractResource implements Resource {

    private final String localPath;
    private final String remotePath;
    private final String rsrc;
    private final Logger log;

    /**
     *
     * @param localPath This can be used to create fully qualified locations i.e
     *        /home/usr/resources, or relative i.e resources/<br>
     *        Thinking this will be set by launch properties
     * @param subPath The path of the resource within local path and url
     * @param url Typically {@link Resources#DEFAULT_URL}
     */
    public AbstractResource(String resourceName, String path, String url, Logger log) {
      localPath = path;
      remotePath = url;
      this.log = log;
      rsrc = resourceName;
    }

    public String getName() {
      return rsrc;
    }

    @Override
    public String getLocalPath() {
      return localPath + rsrc;
    }

    public Logger log() {
      return log;
    }

    protected boolean isLocallyAvailable(String file) {
      return Files.exists(file);
    }

    protected boolean isRemotelyAvailable(String file) {
      return HttpDownloadUtility.canDownload(file, log);
    }

    private boolean downloadResource(String url, String downloadPath) {
      if (isRemotelyAvailable(url)) {
        try {
          HttpDownloadUtility.downloadFile(url, downloadPath, true, log);
          return true;
        } catch (IOException e) {
          log.reportTimeError("Could not retrieve resource from " + downloadPath + " and save it to"
                              + localPath);
          log.reportException(e);
        }
      } else {
        log.reportTimeError("Resource is not available for download: " + url);
      }
      return false;
    }

    protected String get(String file) {
      String path = localPath + file;
      String url = remotePath + file;

      if (isLocallyAvailable(path)) {
        return path;
      }
      log.report("Resource is not available at " + path + ", will attempt to download from " + url);

      if (!downloadResource(url, path)) {
        log.reportError("Download failed for: " + url);
      } else if (!isLocallyAvailable(path)) {
        log.reportError("Downloaded resource cannot be found at " + path);
      } else {
        return path;
      }
      return null;
    }

    @Override
    public String get() {
      return get(rsrc);
    }

    @Override
    public boolean isAvailable() {
      return isAvailable(false);
    }

    @Override
    public boolean isAvailable(boolean showHint) {
      boolean isAvailable = isLocallyAvailable(localPath + rsrc) || isRemotelyAvailable(remotePath + rsrc);

      if (!isAvailable) {
        log.reportTimeError("Could not find local file " + getLocalPath()
                            + " and could not download it from " + remotePath + rsrc
                            + " please manually download and save to " + getLocalPath());
      }

      return isAvailable;
    }
  }

  public static interface Resource {
    boolean isAvailable();

    boolean isAvailable(boolean showHint);

    String get();

    String getLocalPath();
  }

  public enum GENOME_BUILD {
                            HG19("hg19", 37), HG18("hg18", 36);

    private final String build;
    private final int buildInt;

    private GENOME_BUILD(String build, int buildInt) {
      this.build = build;
      this.buildInt = buildInt;
    }

    public String getBuild() {
      return build;
    }

    public int getBuildInt() {
      return buildInt;
    }
  }

  public enum CHROMASOME {
                          C1("1"), C2("2"), C3("3"), C4("4"), C5("5"), C6("6"), C7("7"), C8("8"), C9("9"), C10("10"), C11("11"), C12("12"), C13("13"), C14("14"), C15("15"), C16("16"), C17("17"), C18("18"), C19("19"), C20("20"), C21("21"), C22("22"), CX_PAR("X_PAR"), CX_nonPAR("X_nonPAR");

    private String label;

    private CHROMASOME(String c) {
      label = c;
    }

    public String getLabel() {
      return label;
    }
  }
}
