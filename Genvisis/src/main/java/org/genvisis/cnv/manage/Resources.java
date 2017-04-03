package org.genvisis.cnv.manage;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import org.apache.commons.codec.digest.DigestUtils;
import org.apache.commons.compress.archivers.ArchiveStreamFactory;
import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.compress.utils.IOUtils;
import org.genvisis.CLI;
import org.genvisis.cnv.LaunchProperties;
import org.genvisis.cnv.LaunchProperties.DefaultLaunchKeys;
import org.genvisis.cnv.LaunchProperties.LaunchKey;
import org.genvisis.common.AbstractStartupCheck;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HttpDownloadUtility;
import org.genvisis.common.Logger;
import org.genvisis.common.PSF;
import org.genvisis.common.StartupCheck;
import org.genvisis.common.StartupValidation;
import org.genvisis.filesys.FASTA;
import org.genvisis.seq.manage.BedOps;
import org.genvisis.seq.manage.VCFOps;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;

/**
 * Static utility class for accessing {@link Resource} instances.
 */
public final class Resources {

	private static final double MS_IN_DAY = 8.64e+7;
	public static final String DEFAULT_URL = "http://genvisis.org/rsrc/";
	public static final String BIN_DIR = "bin";
	public static final String GENOME_DIR = "Genome";
	private static Set<Resource> allResources = null;

	public static LaunchKey LAST_RESOURCE_CHECK = new LaunchKey("0", false, "LAST_RESOURCE_CHECK");

	private Resources() {
		// prevent instantiation of utility class
	}

	/**
	 * {@link StartupCheck} to verify and update the versions of local resources.
	 */
	public static class ResourceVersionCheck extends AbstractResourceCheck {

		@Override
		protected String warningHeader() {
			return "WARNING: the following resources may be out of date. Recommend deletion of local copies:";
		}

		@Override
		protected void doCheck() {
			String lastCheck = LaunchProperties.get(LAST_RESOURCE_CHECK);
			// Only check once a day
			if (lastCheck == null || System.currentTimeMillis() - Long.parseLong(lastCheck) > MS_IN_DAY) {
				System.out.println("Validating local resources... ");
				for (Resource rsrc : listAll()) {
					String localMD5 = rsrc.getMD5Local();
					if (localMD5 == null) {
						continue;
					}
					String remoteMD5 = rsrc.getMD5Remote();

					if (remoteMD5 != null && !remoteMD5.equals(localMD5)) {
						addMessage(rsrc.getLocalPath());
					}
				}
				System.out.println("Local resource check complete.");
				LaunchProperties.put(LAST_RESOURCE_CHECK, String.valueOf(System.currentTimeMillis()));
			}
		}

		@Override
		public boolean requiresRemote() {
			return true;
		}
	}

	/**
	 * {@link StartupCheck} for ensuring things in the resources directory are actually resources.
	 */
	public static class LocalResourceCheck extends AbstractResourceCheck {
		/**
		 * Check all local resources against a list of known resources. Report any resources present
		 * locally not in the known list.
		 */
		@Override
		protected void doCheck() {
			Set<String> knownPaths = new HashSet<String>();
			for (Resource r : listAll()) {
				knownPaths.add(r.getLocalPath());
			}
			for (String localFile : localResources()) {
				if (!knownPaths.contains(localFile)) {
					addMessage(localFile);
				}
			}
		}

		@Override
		protected String warningHeader() {
			return "WARNING: the following local file(s) in "
						 + LaunchProperties.get(DefaultLaunchKeys.RESOURCES_DIR) + " are not tracked:";
		}

		@Override
		public boolean requiresRemote() {
			return false;
		}
	}

	/**
	 * Abstract superclass for resource-based {@link StartupCheck}s that want to check
	 * {@link #localResources());
	 */
	private abstract static class AbstractResourceCheck extends AbstractStartupCheck {
		private Set<String> localResources;

		public AbstractResourceCheck() {
			// Build the list of files in the local resource dir
			localResources = new HashSet<String>();
			String resourceDir = LaunchProperties.get(DefaultLaunchKeys.RESOURCES_DIR);
			if (Files.exists(resourceDir)) {
				for (String resource : Files.listAllFilesInTree(resourceDir, false)) {
					localResources.add(new File(resourceDir + resource).getAbsolutePath());
				}
			}
		}

		protected Set<String> localResources() {
			return localResources;
		}
	}

	/**
	 * @return A list of all available resource instances.
	 */
	public static Set<Resource> listAll() {
		if (allResources == null) {
			initResourceList();
		}
		return allResources;
	}

	/**
	 * Synchronized initialization method to ensure resource list is only generated once.
	 */
	private static synchronized void initResourceList() {
		if (allResources == null) {
			Set<Resource> resources = new HashSet<Resource>();
			resources.addAll(path(null).getResources());
			resources.addAll(hapMap(null).getResources());
			resources.addAll(mitoCN(null).getResources());
			resources.addAll(cnv(null).getResources());
			resources.addAll(affy(null).getResources());

			resources.addAll(miniMac(null).getResources());
			resources.addAll(shapeit(null).getResources());
			resources.addAll(eigenstrat(null).getResources());
			resources.addAll(convertf(null).getResources());
			resources.addAll(smartpca(null).getResources());

			for (GENOME_BUILD build : GENOME_BUILD.values()) {
				resources.addAll(genome(build, null).getResources());
				resources.addAll(cnv(null).genome(build).getResources());
				resources.addAll(affy(null).genome(build).getResources());
				for (CHROMOSOME c : CHROMOSOME.values()) {
					resources.addAll(genome(build, null).chr(c).getResources());
				}
			}
			allResources = Collections.unmodifiableSet(resources);
		}
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
			return getTarGzResource("Minimac3", remotePath() + "Minimac3.v2.0.1.tar.gz", "Minimac3");
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
			super(LaunchProperties.get(DefaultLaunchKeys.RESOURCES_DIR) + BIN_DIR + "/shapeit/", "", log,
						Shapeit.class);
		}

		/**
		 * @return A resource for the shapeit app
		 */
		public Resource getShapeit() {
			return getTarGzResource("bin/shapeit",
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
		 * Helper method for formatting resource path: formatted "{build}/{build}", e.g. for
		 * "hg19/hg19.fa"
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
	public static HapMap hapMap(Logger log) {
		return new HapMap(log);
	}

	public static class HapMap extends AbstractResourceFactory {

		public HapMap(Logger log) {
			super("HapMap", log, HapMap.class);
		}

		/**
		 * @return PLINK dataset of Unambiguous HapMap founders
		 */
		public Resource getUnambiguousHapMapFounders() {
			String plinkroot = "unambiguousHapMapFounders";
			Collection<String> indicators = PSF.Plink.getPlinkBedBimFamSet(plinkroot);
			indicators.add("CEUFounders.txt");
			return getTarGzResource(indicators, plinkroot);
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
			return getResource("GDI_percentile.txt");
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

		/**
		 * Use this constructor when resource strucutre is consistent between remote and local paths
		 *
		 * @param subPath Relative path to this resource's top-level for both local (to the RESOURCES
		 *        directory) and remote (where resources are hosted) locations.
		 */
		public AbstractResourceFactory(String subPath, Logger log, Class<?>... classes) {
			this(LaunchProperties.get(DefaultLaunchKeys.RESOURCES_DIR) + subPath + File.separator,
					 DEFAULT_URL + subPath + "/", log, classes);
		}

		/**
		 * Use this constructor when URL and local dir must be significantly different.
		 *
		 * @param localPath Fully qualified local path to a base directory for this resource
		 * @param url Fully qualified URL
		 */
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

		/**
		 * Use this method when the resource a) follows consistent structure on local and remote dirs
		 * and b) does not require Make to be called
		 *
		 * @param rsrc Relative (from local dir) path to the "essential" archived file
		 */
		protected Resource getTarGzResource(String rsrc) {
			return getTarGzResource(rsrc, remotePath() + rsrc + ".tar.gz");
		}

		/**
		 * Use this method when the remote and local path are not consistent, but Make does not need to
		 * be called.
		 * 
		 * @see #getTarGzResource(String)
		 * @param remotePath Absolute path to remote download location
		 */
		protected Resource getTarGzResource(String rsrc, String remotePath) {
			return getTarGzResource(rsrc, remotePath, null);
		}

		/**
		 * Use this when Make needs to be called in the extracted directory
		 *
		 * @see #getTarGzResource(String, String)
		 * @param makeFile Relative (from local) path to make file, or null if this resource does not
		 *        need to be built.
		 */
		protected Resource getTarGzResource(String rsrc, String remotePath, String makePath) {
			String makeDir = makePath == null ? null : localPath() + makePath;
			return new TarGzResource(localPath() + rsrc, localPath(), remotePath, log, makeDir);
		}

		/**
		 * Use this method when the resource a) has an indicator that does not match the return resource
		 * and b) the return resource has consistent structure with the remote resource and c) does not
		 * require Make to be called
		 *
		 *
		 * @param rsrc Relative (from local dir) path to the resource path to return
		 * @param indicatorFiles Relative (from local dir) path to extracted files that indicate success
		 */
		protected Resource getTarGzResource(Collection<String> indicatorFiles, String rsrc) {
			return getTarGzResource(indicatorFiles, rsrc, remotePath() + rsrc + ".tar.gz", null);
		}

		/**
		 * Use this method when the resource [has an indicator that does not match the return resource
		 * and the remote and local path are not consistent, but Make does not need to be called.
		 *
		 * @see #getTarGzResource(Collection, String)
		 * @param remotePath Absolute path to remote download location
		 */
		protected Resource getTarGzResource(Collection<String> indicatorFiles, String rsrc,
																				String remotePath) {
			return getTarGzResource(indicatorFiles, rsrc, remotePath, null);
		}

		/**
		 * Use this when Make needs to be called in the extracted directory and indicator files do not
		 * match the returned resource
		 *
		 * @see #getTarGzResource(Collection, String, String)
		 * @param makeFile Relative (from local) path to make file, or null if this resource does not
		 *        need to be built.
		 * 
		 */
		protected Resource getTarGzResource(Collection<String> indicatorFiles, String rsrc,
																				String remotePath, String makePath) {

			String makeDir = makePath == null ? null : localPath() + makePath;
			Collection<String> unzippedPaths = Lists.newArrayListWithCapacity(indicatorFiles.size());
			for (String indicatorFile : indicatorFiles) {
				unzippedPaths.add(localPath() + indicatorFile);
			}
			return new TarGzResource(unzippedPaths, localPath() + rsrc, localPath(), remotePath, log,
															 makeDir);
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
	// FIXME: this class is currently doing too much.
	// It would be nice to split out things that need to be built and things that need to be chmodded
	// into some sort of merge-able stack,
	// e.g. by constructing multiple Resources wrapping each other.
	public static class TarGzResource extends AbstractResource {
		private final Collection<String> unzippedPaths;
		private final String returnPath;
		private final String makeDir;

		/**
		 * 
		 * @param unzippedPath single local path to indicate success and return as Resource
		 * @param extractionDir directory to extract to
		 * @param url Remote tar.gz path
		 * @param log
		 * @param makeDir directory to call make, null to not make
		 */
		public TarGzResource(String unzippedPath, String extractionDir, String url, Logger log,
												 String makeDir) {
			this(Collections.singleton(unzippedPath), unzippedPath, extractionDir, url, log, makeDir);
		}

		/**
		 * @param path
		 * @param url Remote tar.gz path
		 */
		/**
		 * 
		 * @param unzippedPaths local paths to unzipped indicator files
		 * @param returnPath local path to return as Resource
		 * @param extractionDir directory to extract to
		 * @param url Remote tar.gz path
		 * @param log
		 * @param makeDir directory to call make, null to not make
		 */
		public TarGzResource(Iterable<String> unzippedPaths, String returnPath, String extractionDir,
												 String url, Logger log, String makeDir) {
			super(extractionDir + new File(url).getName(), url, log);
			this.unzippedPaths = ImmutableList.copyOf(unzippedPaths);
			this.returnPath = returnPath;
			this.makeDir = makeDir;
		}



		@Override
		public String get() {
			boolean available = true;
			for (String unzippedPath : unzippedPaths) {
				if (!isLocallyAvailable(unzippedPath)) {
					available = false;
					break;
				}
			}
			if (available)
				return returnPath;

			String path = super.get();
			if (path != null) {
				String extracted = extractTarGz(path);
				if (extracted != null) {
					// Run make if the zipped archive requires compilation
					if (makeDir != null && !CmdLine.runDefaults("make -s", makeDir)) {
						return null;
					}
					if (Files.exists(extracted)) {
						// If return path is a file, ensure extracted file can be used
						Files.chmod(extracted);
					}
				}
				return extracted;
			}

			return null;
		}

		private String extractTarGz(String targzPath) {
			final File inputFile = new File(targzPath);
			final String destination = inputFile.getParent() + File.separator;
			inputFile.deleteOnExit();

			// from: http://stackoverflow.com/a/7556307/1027800
			try {
				// Unzip the downloaded targz
				// Not all files end with .tar.gz.. some are .tgz.. so handle both cases
				String unzippedName = inputFile.getName().substring(0,
																														inputFile.getName().lastIndexOf('.'));
				if (!unzippedName.endsWith(".tar")) {
					unzippedName += ".tar";
				}
				final File unzippedTar = new File(destination, unzippedName);

				final GZIPInputStream in = new GZIPInputStream(new FileInputStream(inputFile));
				final FileOutputStream out = new FileOutputStream(unzippedTar);
				IOUtils.copy(in, out);
				in.close();
				out.close();
				unzippedTar.deleteOnExit();

				// Extract all entries in the tar
				final InputStream is = new FileInputStream(unzippedTar);
				final TarArchiveInputStream tarStream = (TarArchiveInputStream) new ArchiveStreamFactory().createArchiveInputStream("tar",
																																																														is);
				TarArchiveEntry entry = null;
				while ((entry = (TarArchiveEntry) tarStream.getNextEntry()) != null) {
					final File outputFile = new File(destination, entry.getName());
					if (entry.isDirectory()) {
						if (!outputFile.exists() && !outputFile.mkdirs()) {
							throw new IllegalStateException(String.format("Couldn't create directory %s.",
																														outputFile.getAbsolutePath()));
						}
					} else {
						final File parentDir = outputFile.getParentFile();
						if (!parentDir.exists() && !parentDir.mkdirs()) {
							throw new IllegalStateException(String.format("Couldn't create directory %s.",
																														parentDir));
						} else {
							final OutputStream outputFileStream = new FileOutputStream(outputFile);
							IOUtils.copy(tarStream, outputFileStream);
							outputFileStream.close();
						}
					}
				}
				tarStream.close();
				for (String unzippedPath : unzippedPaths) {
					if (!Files.exists(unzippedPath)) {
						log().reportError("Resource extracted but did not find " + unzippedPath);
						return null;
					}
				}
				return returnPath;
			} catch (Exception e) {
				log().reportError("Failed to extract: " + targzPath);
				log().reportException(e);
				// Remove the unzippedPaths as these are the markers of success
				for (String unzippedPath : unzippedPaths) {
					new File(unzippedPath).deleteOnExit();
				}
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

		@Override
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
			if (log == null) {
				this.log = new Logger();
			} else {
				this.log = log;
			}
			rsrc = new File(path).getName();
		}

		public String getName() {
			return rsrc;
		}

		@Override
		public String getLocalPath() {
			return localPath;
		}

		@Override
		public String getMD5Local() {
			String md5 = null;
			if (isLocallyAvailable(localPath)) {
				try {
					String md5Path = localPath + ".md5";
					if (isLocallyAvailable(md5Path)) {
						// If we already have a local .md5, we just need to read it and return.
						BufferedReader reader = Files.getAppropriateReader(md5Path);
						md5 = reader.readLine().trim();
						reader.close();
					} else {
						// Compute the md5 for this resource
						log.report("Generating md5 checksum for resource: " + absolutePath());
						FileInputStream fis = new FileInputStream(new File(localPath));
						md5 = DigestUtils.md5Hex(fis);
						fis.close();
						// Write out the computed md5
						PrintWriter writer = Files.getAppropriateWriter(md5Path);
						writer.println(md5);
						writer.flush();
						writer.close();
					}
				} catch (FileNotFoundException e) {
					log.reportException(e);
				} catch (IOException e) {
					log.reportException(e);
				}
			}
			return md5;
		}

		@Override
		public String getMD5Remote() {
			String md5 = null;
			try {
				md5 = HttpDownloadUtility.readFileAsHexString(remotePath + ".md5");
			} catch (IOException e) {
				// If remote MD5 isn't available, that can be OK.
			}
			return md5;
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

					// verify the download
					String md5Remote = getMD5Remote();
					if (md5Remote == null) {
						log.reportTimeWarning("Remote MD5 not available.");
						return true;
					}
					if (getMD5Local().equals(getMD5Remote())) {
						return true;
					}

					// md5 did not match
					log.reportError("Local md5 checksum does not match remote for resource: " + url);
				} catch (IOException e) {
					log.reportError("Could not retrieve resource from " + url + " and save it to"
													+ downloadPath);
					log.reportException(e);
				}
			} else {
				log.reportError("Resource is not available for download: " + url);
			}
			return false;
		}

		private String absolutePath() {
			return new File(localPath).getAbsolutePath();
		}

		@Override
		public String get() {
			// Ensure resource validation is complete
			StartupValidation.passed();

			if (isLocallyAvailable(localPath)) {
				return localPath;
			}
			log.report("Resource is not available at " + absolutePath()
								 + ", will attempt to download from " + remotePath);

			if (!downloadResource(remotePath, localPath)) {
				log.reportError("Download failed for: " + remotePath);
			} else if (!isLocallyAvailable(localPath)) {
				log.reportError("Downloaded resource cannot be found at " + absolutePath());
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
				log.reportError("Could not find local file " + absolutePath()
												+ " and could not download it from " + remotePath
												+ " please manually download and save to " + absolutePath());
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

		/**
		 * @return The computed md5 checksum for this resource, or null if the resource is not locally
		 *         available.
		 */
		String getMD5Local();

		/**
		 * @return The "official" md5 checksum for this resource, or null if lacking internet
		 *         connectivity.
		 */
		String getMD5Remote();
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
		C1("1"),
		C2("2"),
		C3("3"),
		C4("4"),
		C5("5"),
		C6("6"),
		C7("7"),
		C8("8"),
		C9("9"),
		C10("10"),
		C11("11"),
		C12("12"),
		C13("13"),
		C14("14"),
		C15("15"),
		C16("16"),
		C17("17"),
		C18("18"),
		C19("19"),
		C20("20"),
		C21("21"),
		C22("22"),
		CX_PAR("X_PAR"),
		CX_nonPAR("X_nonPAR");

		private String label;

		private CHROMOSOME(String c) {
			label = c;
		}

		public String getLabel() {
			return label;
		}
	}

	public static void main(String... args) {
		CLI c = new CLI(Resources.class);
		final String localCheck = "local";

		c.addFlag(localCheck, "Check local resources for untracked files.");
		c.parseWithExit(args);

		if (c.has(localCheck)) {
			List<String> checkResults = new LocalResourceCheck().check();
			for (String s : checkResults) {
				System.err.println(s);
			}
		}
	}
}
