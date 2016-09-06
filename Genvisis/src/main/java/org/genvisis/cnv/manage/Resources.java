package org.genvisis.cnv.manage;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
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
 * Static utility class for accessing {@link Resource} instances.
 */
public final class Resources {

  public static final String DEFAULT_URL = "http://genvisis.org/rsrc/";
  public static final String DEFAULT_LOCAL_DIR = "resources" + File.separator;
  public static final String BIN_DIR = "bin";
  public static final String GENOME_DIR = "Genome";

  private Resources() {
    // prevent instantiation of utility class
  }

  /**
   * Helper method for chaining resource calls
   */
  public static MiniMac miniMac(Logger log) {
    return new MiniMac(log);
  }

  /**
   * MiniMac resource container
   */
  public static class MiniMac extends AbstractResourceFactory {
    public MiniMac(Logger log) {
      super(BIN_DIR + "/Minimac3", log, MiniMac.class);
    }

    /**
     * @return A resource for the MiniMac3 app
     */
    public Resource getMiniMac3() {
      return getTarGzResource("Minimac3.v1.0.14.tar.gz");
    }
  }

  /**
   * Helper method for chaining resource calls
   */
  public static Shapeit shapeit(Logger log) {
    return new Shapeit(log);
  }

  /**
   * Shapeit resource container
   */
  public static class Shapeit extends AbstractResourceFactory {
    public Shapeit(Logger log) {
      super(DEFAULT_LOCAL_DIR + BIN_DIR + File.separator + "shapeit" + File.separator
            + "shapeit.tar.gz",
            "https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz",
            log, Shapeit.class);
    }

    /**
     * @return A resource for the shapeit app
     */
    public Resource getShapeit() {
      return getResource("");
    }
  }

  /**
   * Chromasome-related resource container. Always requires a related {@link GENOME_BUILD}.
   */
  public static class Chr extends AbstractResourceFactory {
    private String build;
    private CHROMASOME c;

    public Chr(GENOME_BUILD build, CHROMASOME c, Logger log) {
      super(GENOME_DIR + "/" + build.getBuild() + "/chr", log, Chr.class);
      this.build = build.getBuild();
      this.c = c;
    }

    /**
     * @return The genetic map for the requested {@link CHROMASOME}
     */
    public Resource getGeneticMap() {
      String prefix = "genetic_map_";
      String extension = ".txt.gz";

      if (CHROMASOME.CX_PAR.equals(c)) {
        getResource(getPath(prefix, c.toString() + "2", extension));
        return getResource(getPath(prefix, c.getLabel() + "1", extension));
      }
      return getResource(getPath(prefix, c.getLabel(), extension));
    }

    /**
     * @return The G1K ref for the requested {@link CHROMASOME}
     */
    public Resource getG1Kphase3v5RefPanel() {
      return getResource(getPath("1000genomes_ref_panel_Phase3v5_", c.getLabel(), ".m3vcf.gz"));
    }

    /**
     * Helper method to format the resource path
     */
    private String getPath(String prefix, String chromasome, String suffix) {
      return prefix + build + "_chr" + chromasome + suffix;
    }
  }

  /**
   * Helper method for chaining resource calls
   */
  public static Genome genome(GENOME_BUILD build, Logger log) {
    return new Genome(build, log);
  }

  /**
   * Container for {@link GENOME_BUILD}-related resources.
   */
  public static class Genome extends AbstractResourceFactory {
    private GENOME_BUILD build;

    public Genome(GENOME_BUILD build, Logger log) {
      super(GENOME_DIR, log, Genome.class);
      this.build = build;
    }

    /**
     * @return The RefSeq.gtrack for this {@link GENOME_BUILD}
     */
    public Resource getGTrack() {
      return getResource(build.getBuild() + "/RefSeq_" + build.getBuild() + ".gtrack");
    }

    /**
     * @return The GC5 base for this {@link GENOME_BUILD}
     */
    public Resource getModelBase() {
      return getResource(getPath() + "_gc5Base.txt");
    }

    /**
     * @return The DB Snp for this {@link GENOME_BUILD}
     */
    public Resource getDBSNP() {
      return getVCFResource(getPath() + "_dbSnp147.vcf.gz");
    }

    /**
     * @return The "genes##.xln" for this {@link GENOME_BUILD}
     */
    public Resource getGenes() {
      return getResource(build.getBuild() + "/genes" + build.getBuildInt() + ".xln");
    }

    /**
     * Helper method for formatting resource path
     */
    private String getPath() {
      String b = build.getBuild();
      return new StringBuilder().append(b).append("/").append(b).toString();
    }

    /**
     * Helper method for chaining resource calls
     */
    public Chr chr(CHROMASOME c) {
      return new Chr(build, c, log());
    }
  }

  /**
   * Helper method for chaining resource calls
   */
  public static MitoCN mitoCN(Logger log) {
    return new MitoCN(log);
  }

  /**
   * Container for MitoCN resources
   */
  public static class MitoCN extends AbstractResourceFactory {
    public MitoCN(Logger log) {
      super("MitoCN", log, MitoCN.class);
    }

    /**
     * @return White WBC data
     */
    public Resource getWhiteWBC() {
      return getResource("Whites_WBC_TOTAL_SingleSNPmatched.final.beta");
    }

    /**
     * @return Black WBC data
     */
    public Resource getBlackWBC() {
      return getResource("Blacks_WBC_TOTAL_SingleSNPmatched.final.beta");
    }

    /**
     * @return Total WBC data
     */
    public Resource getTotalWBC() {
      return getResource("WBC_TOTAL_SingleSNPmatched.final.beta");
    }
  }

  /**
   * Illumina Bundle TODO
   */

  /**
   * Helper method for chaining resource calls
   */
  public static AffySnp6 affy(Logger log) {
    return new AffySnp6(log);
  }

  /**
   * Container for Affy resources.
   */
  public static class AffySnp6 extends AbstractResourceFactory {
    public AffySnp6(Logger log) {
      super("Arrays/AffySnp6", log, AffySnp6.class);
    }

    /**
     * @return Marker positions for the specified {@link GENOME_BUILD}
     */
    public Resource getMarkerPositions(GENOME_BUILD build) {
      return getResource(build.getBuild() + "_markerPositions.txt");
    }

    /**
     * @return HMM file
     */
    public Resource getHMM() {
      return getResource("affygw6.hmm");
    }

    /**
     * @return ABLookup file
     */
    public Resource getABLookup() {
      return getResource("AB_lookup.dat");
    }
  }

  /**
   * Abstract superclass for containers that create {@link Resource} instances.
   * <p>
   * Note: the {@link #getResources()} implementation uses reflection to find all zero-parameter
   * methods that return a {@link Resource} in the class list this factory is constructed with. A
   * list of classes is used to facilitate inheritance between {@link ResourceFactory} classes. But,
   * if you do not use this standard method format for accessing resources you should override this
   * method.
   * </p>
   */
  private abstract static class AbstractResourceFactory implements ResourceFactory {
    private final String localPath;
    private final String remotePath;
    private final Logger log;
    private final Class<?>[] classes;

    public AbstractResourceFactory(String subPath, Logger log, Class<?>... classes) {
      this(DEFAULT_LOCAL_DIR + subPath + File.separator, DEFAULT_URL + subPath + "/", log, classes);
    }

    public AbstractResourceFactory(String localPath, String url, Logger log, Class<?>... classes) {
      this.localPath = localPath;
      this.remotePath = url;
      this.log = log;
      this.classes = classes;
    }

    @Override
    public List<Resource> getResources() {
      List<Resource> resources = new ArrayList<Resource>();

      for (Class<?> c : classes) {
        for (Method m : c.getDeclaredMethods()) {
          if (m.getReturnType().equals(Resources.Resource.class) && m.getParameterCount() == 0) {
            try {
              resources.add((Resource) m.invoke(this));
            } catch (Exception e) {
              log.reportError("Failed to add resource: " + m.getName());
            }
          }
        }
      }
      return resources;
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

  /**
   * A ResourceFactory holds all the configuration for building {@link Resource} instances.
   * Typically this is just the local and remote path for a resource folder.
   */
  private static interface ResourceFactory {
    List<Resource> getResources();
  }

  /**
   * Resource that may be .tar.gz compressed
   */
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

  /**
   * VCF resource with an accompanying index file
   */
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

  /**
   * Basic resource class
   */
  public static class DefaultResource extends AbstractResource {
    public DefaultResource(String resourceName, String path, String url, Logger log) {
      super(resourceName, path, url, log);
    }
  }

  /**
   * Abstract {@link Resource} superclass
   */
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
      boolean isAvailable = isLocallyAvailable(localPath + rsrc)
                            || isRemotelyAvailable(remotePath + rsrc);

      if (!isAvailable) {
        log.reportTimeError("Could not find local file " + getLocalPath()
                            + " and could not download it from " + remotePath + rsrc
                            + " please manually download and save to " + getLocalPath());
      }

      return isAvailable;
    }
  }

  /**
   * A resource is a general-purpose file used by Genvisis but not shipped with the Genvisis core.
   * Resources are typically available remotely and thus can be downloaded automatically if missing.
   * Use the {@link #isAvailable()} method to check if this resource can be obtained, and
   * {@link #get()} to return the local path to the resource.
   */
  public static interface Resource {
    /**
     * @return {@code true} if this resource can be found locally or remotely.
     */
    boolean isAvailable();

    /**
     * As {@link #isAvailable()} but can display a hint to the user on how to manually download, if
     * the resource is not available.
     *
     * @return {@code true} if this resource can be found locally or remotely.
     */
    boolean isAvailable(boolean showHint);

    /**
     * Ensure this resource is available locally, downloading it if necessary.
     *
     * @return The local path to this resource.
     */
    String get();

    /**
     * Unlike {@link #get()}, this method will not download a remote resource.
     *
     * @return The local path to this resource.
     */
    String getLocalPath();
  }

  /**
   * Supported genome reference builds
   */
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

  /**
   * Supported chromasomes
   */
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
