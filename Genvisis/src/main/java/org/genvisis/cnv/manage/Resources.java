package org.genvisis.cnv.manage;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.compress.compressors.gzip.GzipCompressorInputStream;
import org.genvisis.cnv.LaunchProperties;
import org.genvisis.cnv.LaunchProperties.LaunchKey;
import org.genvisis.common.Files;
import org.genvisis.common.HttpDownloadUtility;
import org.genvisis.common.Logger;
import org.genvisis.filesys.FASTA;
import org.genvisis.seq.manage.BedOps;
import org.genvisis.seq.manage.VCFOps;

/**
 * Static utility class for accessing {@link Resource} instances.
 */
public final class Resources {

	public static final String DEFAULT_URL = "http://genvisis.org/rsrc/";
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
			return getTarGzResource(localPath() + "Minimac3", remotePath() + "Minimac3.v1.0.14.tar.gz");
		}
	}

	/**
	 * Helper method for chaining resource calls
	 */
	public static Eigenstrat eigenstrat(Logger log) {
		return new Eigenstrat(log);
	}

	/**
	 * Eigenstrat resource container
	 */
	public static class Eigenstrat extends AbstractResourceFactory {
		public Eigenstrat(Logger log) {
			super(BIN_DIR + "/Eigenstrat", log, Eigenstrat.class);
		}

		/**
		 * @return A resource for the Eigenstrat program
		 */
		public Resource getEigenstrat() {
			return getTarGzResource("eigenstrat");
		}

		/**
		 * @return A resource for the EigenstratQTL program
		 */
		public Resource getEigenstratQTL() {
			return getTarGzResource("eigenstratQTL");
		}
	}

	/**
	 * Helper method for chaining resource calls
	 */
	public static Convertf convertf(Logger log) {
		return new Convertf(log);
	}

	/**
	 * Convertf resource container
	 */
	public static class Convertf extends AbstractResourceFactory {
		public Convertf(Logger log) {
			super(BIN_DIR + "/Convertf", log, Convertf.class);
		}

		/**
		 * @return A resource for the convertf app
		 */
		public Resource getConvertf() {
			return getTarGzResource("convertf");
		}
	}

	/**
	 * Helper method for chaining resource calls
	 */
	public static Smartpca smartpca(Logger log) {
		return new Smartpca(log);
	}

	/**
	 * Smartpca resource container
	 */
	public static class Smartpca extends AbstractResourceFactory {
		public Smartpca(Logger log) {
			super(BIN_DIR + "/Smartpca", log, Smartpca.class);
		}

		/**
		 * @return A resource for the MiniMac3 app
		 */
		public Resource getSmartpca() {
			return getTarGzResource("smartpca");
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
			super(LaunchProperties.get(LaunchKey.RESOURCES_DIR) + BIN_DIR + File.separator + "shapeit", "", log, Shapeit.class);
		}

		/**
		 * @return A resource for the shapeit app
		 */
		public Resource getShapeit() {
			return getTarGzResource(localPath(),
															"https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz");
		}
	}

	/**
	 * Chromosome-related resource container. Always requires a related {@link GENOME_BUILD}.
	 */
	public static class Chr extends AbstractResourceFactory {
		private final String build;
		private final CHROMOSOME c;

		public Chr(GENOME_BUILD build, CHROMOSOME c, Logger log) {
			super(GENOME_DIR + "/" + build.getBuild() + "/chr", log, Chr.class);
			this.build = build.getBuild();
			this.c = c;
		}

		/**
		 * @return The genetic map for the requested {@link CHROMOSOME}
		 */
		public Resource getGeneticMap() {
			String prefix = "genetic_map_";
			String extension = ".txt.gz";

			if (CHROMOSOME.CX_PAR.equals(c)) {
				getResource(getPath(prefix, c.toString() + "2", extension));
				return getResource(getPath(prefix, c.getLabel() + "1", extension));
			}
			return getResource(getPath(prefix, c.getLabel(), extension));
		}

		/**
		 * @return The G1K ref for the requested {@link CHROMOSOME}
		 */
		public Resource getG1Kphase3v5RefPanel() {
			return getResource(getPath("1000genomes_ref_panel_Phase3v5_", c.getLabel(), ".m3vcf.gz"));
		}

		/**
		 * Helper method to format the resource path
		 */
		private String getPath(String prefix, String chromosome, String suffix) {
			return prefix + build + "_chr" + chromosome + suffix;
		}
	}

	/**
	 * Helper {@link Pathway} accessor for method chaining
	 */
	public static Pathway path(Logger log) {
		return new Pathway(log);
	}

	/**
	 * Container for pathway database resources
	 */
	public static class Pathway extends AbstractResourceFactory {
		public Pathway(Logger log) {
			super("Pathways", log, Pathway.class);
		}

		public Resource getKegg() {
			return getResource("kegg.ser");
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
		private final GENOME_BUILD build;

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
		 * @return The 100-mer mappability track for this {@link GENOME_BUILD}, note we do not download
		 *         the index - generate if needed with {@link BedOps#verifyBedIndex(String, Logger)}
		 */
		public Resource get100MerMappabilityTrack() {
			return getResource(getPath() + "_wgEncodeCrgMapabilityAlign100mer.bedGraph");
		}

		/**
		 * @return b138 DB .ser resource
		 */
		public Resource getB138() {
			return getResource(build.getBuild() + "/b138_" + build.getBuildInt() + "_3.ser");
		}

		/**
		 * @return The "genes##.xln" for this {@link GENOME_BUILD}
		 */
		public Resource getGenes() {
			return getResource(build.getBuild() + "/genes" + build.getBuildInt() + ".xln");
		}

		/**
		 * @return The reference genome FASTA for this {@link GENOME_BUILD}
		 */
		public Resource getFASTA() {
			return getFASTAResource(getPath() + ".fa");
		}

		/**
		 * @return .dat list of known problematic regions (e.g. centromere, chromosome ends)
		 */
		public Resource getProblematicRegions() {
			return getResource(build.getBuild() + "/problematicRegions_" + build.getBuild() + ".dat");
		}

		/**
		 * Helper method for formatting resource path: formatted "{build}/{build}", e.g. for "hg19/hg19.fa"
		 */
		private String getPath() {
			String b = build.getBuild();
			return new StringBuilder().append(b).append("/").append(b).toString();
		}

		/**
		 * Helper method for chaining resource calls
		 */
		public Chr chr(CHROMOSOME c) {
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
	 * Container for build-specific Affy resources
	 */
	public static class AffyGenomes extends AbstractResourceFactory {
		private final GENOME_BUILD build;

		public AffyGenomes(GENOME_BUILD build, Logger log) {
			super(AffySnp6.DIR, log, AffyGenomes.class);
			this.build = build;
		}

		/**
		 * @return Marker positions for the specified {@link GENOME_BUILD}
		 */
		public Resource getMarkerPositions() {
			return getResource(build.getBuild() + "_markerPositions.txt");
		}

		/**
		 * @return gcmodel file
		 */
		public Resource getGcmodel() {
			return getResource("affygw6." + build.getBuild() + ".gcmodel");
		}

		/**
		 * @return pfb file
		 */
		public Resource getHmm() {
			return getResource("affygw6." + build.getBuild() + ".pfb");
		}
	}

	/**
	 * Container for Affy resources.
	 */
	public static class AffySnp6 extends AbstractResourceFactory {
		private static final String DIR = "Arrays/AffySnp6";

		public AffySnp6(Logger log) {
			super(DIR, log, AffySnp6.class);
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

		public AffyGenomes genome(GENOME_BUILD build) {
			return new AffyGenomes(build, log());
		}
	}

	/**
	 * Helper method for chaining resource calls
	 */
	public static CNV cnv(Logger log) {
		return new CNV(log);
	}

	/**
	 * Get genome-specific CNV resources
	 */
	public static class CNVGenomes extends AbstractResourceFactory {
		private final GENOME_BUILD build;

		public CNVGenomes(GENOME_BUILD build, Logger log) {
			super(CNV.DIR, log, CNVGenomes.class);
			this.build = build;
		}

		public Resource getAllPfb() {
			return getPfb("all");
		}

		public Resource getAllGcmodel() {
			return getGc("all");
		}

		public Resource get550Pfb() {
			return getPfb("550");
		}

		public Resource get550Gcmodel() {
			return getGc("550");
		}

		private Resource getPfb(String cnv) {
			return get(cnv, "pfb");
		}

		private Resource getGc(String cnv) {
			return get(cnv, "gcmodel");
		}

		private Resource get(String cnv, String extension) {
			StringBuilder sb = new StringBuilder();
			sb.append("hh").append(cnv).append(".").append(build.build).append(".").append(extension);
			return getResource(sb.toString());
		}
	}


	/**
	 * Get general annotation resources
	 */
	public static class Annotation extends AbstractResourceFactory {
		private static final String DIR = "Annotation";

		public Annotation(Logger log) {
			super(DIR, log, Annotation.class);
		}

		public Resource getGDI() {
			return getResource("GDI.txt");
		}


		/**
		 * Helper method for chaining resource calls.
		 */
		public Annotation annotation() {
			return new Annotation(log());
		}
	}

	/**
	 * Helper method for chaining resource calls
	 */
	public static Annotation annotation(Logger log) {
		return new Annotation(log);
	}



	/**
	 * Get general CNV resources
	 */
	public static class CNV extends AbstractResourceFactory {
		private static final String DIR = "CNV";

		public CNV(Logger log) {
			super(DIR, log, CNV.class);
		}

		public Resource getAllHmm() {
			return getResource("hhall.hmm");
		}

		public Resource get550Hmm() {
			return getResource("hh550.hmm");
		}

		/**
		 * Helper method for chaining resource calls.
		 */
		public CNVGenomes genome(GENOME_BUILD build) {
			return new CNVGenomes(build, log());
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
			this(LaunchProperties.get(LaunchKey.RESOURCES_DIR) + subPath + File.separator, DEFAULT_URL + subPath + "/", log, classes);
		}

		public AbstractResourceFactory(String localPath, String url, Logger log, Class<?>... classes) {
			this.localPath = localPath;
			remotePath = url;
			this.log = log;
			this.classes = classes;
		}

		@Override
		public List<Resource> getResources() {
			List<Resource> resources = new ArrayList<Resource>();

			for (Class<?> c : classes) {
				for (Method m : c.getDeclaredMethods()) {
					if (m.getReturnType().equals(Resources.Resource.class)
							&& m.getParameterTypes().length == 0) {
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

		protected String localPath() {
			return localPath;
		}

		protected String remotePath() {
			return remotePath;
		}

		protected Logger log() {
			return log;
		}

		protected Resource getTarGzResource(String rsrc) {
			return getTarGzResource(localPath() + rsrc, remotePath() + rsrc + ".tar.gz");
		}

		protected Resource getTarGzResource(String unzippedPath, String remotePath) {
			return new TarGzResource(unzippedPath, localPath(), remotePath, log);
		}

		protected Resource getVCFResource(String rsrc) {
			return getVCFResource(localPath + rsrc, remotePath + rsrc);
		}

		protected Resource getVCFResource(String localPath, String remotePath) {
			return new VCFResource(localPath, remotePath, log);
		}

		protected Resource getFASTAResource(String rsrc) {
			return getFASTAResource(localPath + rsrc, remotePath + rsrc);
		}

		protected Resource getFASTAResource(String localPath, String remotePath) {
			return new FASTAResource(localPath, remotePath, log);
		}

		protected Resource getResource(String rsrc) {
			return getResource(localPath + rsrc, remotePath + rsrc);
		}

		protected Resource getResource(String localPath, String remotePath) {
			return new DefaultResource(localPath, remotePath, log);
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
	 * Resource that is .tar.gz compressed
	 */
	public static class TarGzResource extends AbstractResource {
		private final String unzippedPath;

		/**
		 * @param path Local path to unzipped indicator file
		 * @param url Remote tar.gz path
		 */
		public TarGzResource(String unzippedPath, String extractionDir, String url, Logger log) {
			super(extractionDir + new File(url).getName(), url, log);
			this.unzippedPath = unzippedPath;
		}

		@Override
		public String get() {
			if (isLocallyAvailable(unzippedPath)) {
				return unzippedPath;
			}

			String path = super.get();
			if (path != null) {
				return extractTarGz(path);
			}

			return null;
		}

		private String extractTarGz(String targzPath) {
			final int BUFFER = 2048;

			try {
				FileInputStream fin = new FileInputStream(targzPath);
				BufferedInputStream in = new BufferedInputStream(fin);
				GzipCompressorInputStream gzIn = new GzipCompressorInputStream(in);
				TarArchiveInputStream tarIn = new TarArchiveInputStream(gzIn);

				TarArchiveEntry entry = null;
				String destination = new File(targzPath).getParent() + File.separator;

				/** Read the tar entries using the getNextEntry method **/
				Set<TarArchiveEntry> dirs = new HashSet<TarArchiveEntry>();
				Set<TarArchiveEntry> files = new HashSet<TarArchiveEntry>();

				while ((entry = (TarArchiveEntry) tarIn.getNextEntry()) != null) {
					// Sort out the directories and files. Ensure the directories are all created first.
					if (entry.isDirectory()) {
						dirs.add(entry);
					} else {
						files.add(entry);
					}
				}

				/** If the entry is a directory, create the directory. **/
				for (TarArchiveEntry e : dirs) {
					File f = new File(destination + e.getName());
					f.mkdirs();
				}

				/**
				 * If the entry is a file,write the decompressed file to the disk and close destination
				 * stream.
				 **/
				for (TarArchiveEntry e : files) {
					int count;
					byte[] data = new byte[BUFFER];
					String fileName = destination + e.getName();
					File f = new File(fileName).getParentFile();
					f.mkdirs();

					FileOutputStream fos = new FileOutputStream(fileName);
					BufferedOutputStream dest = new BufferedOutputStream(fos, BUFFER);

					while ((count = tarIn.read(data, 0, BUFFER)) != -1) {
						dest.write(data, 0, count);
					}
					dest.close();
				}

				/** Close the input stream **/
				tarIn.close();
				new File(targzPath).delete();
				return unzippedPath;
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
		private Resource index;

		public VCFResource(String path, String url, Logger log) {
			super(path, url, log);
			index = new DefaultResource(VCFOps.getIndex(path), VCFOps.getIndex(url), log);
		}

		@Override
		public String get() {
			String path = super.get();
			if (path != null && index.get() == null) {
				log().reportTimeWarning("no index found for vcf file: " + path);
			}
			return path;
		}
	}

	/**
	 * FASTA Resource with accompanying index and dictionary files
	 */
	public static class FASTAResource extends AbstractResource {
		private Resource index;
		private Resource dictionary;

		public FASTAResource(String path, String url, Logger log) {
			super(path, url, log);
			index = new DefaultResource(FASTA.getIndex(path), FASTA.getIndex(url), log);
			dictionary = new DefaultResource(FASTA.getDictionary(path), FASTA.getDictionary(url), log);
		}

		public String get() {
			String path = super.get();
			if (path != null && index.get() == null) {
				log().reportTimeWarning("no index found for FASTA file: " + path);
			}
			if (path != null && dictionary.get() == null) {
				log().reportTimeWarning("no dictionary found for FASTA file: " + path);
			}
			return path;
		}
	}

	/**
	 * Basic resource class
	 */
	public static class DefaultResource extends AbstractResource {
		public DefaultResource(String path, String url, Logger log) {
			super(path, url, log);
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
		public AbstractResource(String path, String url, Logger log) {
			localPath = path;
			remotePath = url;
			this.log = log;
			rsrc = new File(path).getName();
		}

		public String getName() {
			return rsrc;
		}

		@Override
		public String getLocalPath() {
			return localPath;
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
					log.reportError("Could not retrieve resource from "	+ url + " and save it to"
													+ downloadPath);
					log.reportException(e);
				}
			} else {
				log.reportError("Resource is not available for download: " + url);
			}
			return false;
		}

		@Override
		public String get() {
			if (isLocallyAvailable(localPath)) {
				return localPath;
			}
			log.report("Resource is not available at "	+ localPath + ", will attempt to download from "
									+ remotePath);

			if (!downloadResource(remotePath, localPath)) {
				log.reportError("Download failed for: " + remotePath);
			} else if (!isLocallyAvailable(localPath)) {
				log.reportError("Downloaded resource cannot be found at " + localPath);
			} else {
				return localPath;
			}
			return null;
		}

		@Override
		public String getAbsolute() {
			return new File(get()).getAbsolutePath();
		}

		@Override
		public boolean isAvailable() {
			return isAvailable(false);
		}

		@Override
		public boolean isAvailable(boolean showHint) {
			boolean isAvailable = isLocallyAvailable(localPath) || isRemotelyAvailable(remotePath);

			if (!isAvailable) {
				log.reportError("Could not find local file "	+ localPath
												+ " and could not download it from " + remotePath
												+ " please manually download and save to " + localPath);
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
		 * As {@link #get()} but always returns a fully qualified path.
		 *
		 * @see #get()
		 * @return The absolute path to this resource.
		 */
		String getAbsolute();

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
	 * Supported chromosomes
	 */
	public enum CHROMOSOME {
													C1("1"), C2("2"), C3("3"), C4("4"), C5("5"), C6("6"), C7("7"), C8("8"), C9("9"), C10("10"), C11("11"), C12("12"), C13("13"), C14("14"), C15("15"), C16("16"), C17("17"), C18("18"), C19("19"), C20("20"), C21("21"), C22("22"), CX_PAR("X_PAR"), CX_nonPAR("X_nonPAR");

		private String label;

		private CHROMOSOME(String c) {
			label = c;
		}

		public String getLabel() {
			return label;
		}
	}
}
