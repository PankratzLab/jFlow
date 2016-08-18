package org.genvisis.cnv.manage;

import java.io.IOException;

import org.genvisis.common.CmdLine;
import org.genvisis.common.Files;
import org.genvisis.common.HttpDownloadUtility;
import org.genvisis.common.Logger;
import org.genvisis.common.ext;
import org.genvisis.seq.manage.VCFOps;

public class Resources {

	public static final String DEFAULT_URL = "http://genvisis.org/rsrc/";
	public static final String DEFUALT_LOCAL_DIR_BASE = "resources/";

	public static final String BIN_SUB_DIR = "bin/";
	public static final String GENOME_SUB_DIR = "Genome/";
	public static final String MITO_SUB_DIR = "MitoCN/";

	public static final String ARRAY_SUB_DIR = "Arrays/";

	public static final String CHR_SUB_DIR = "chr/";

	// TODO, a big TODO
	// need to add web-based download, and local file structure
	// could probably do this like project properties...

	public static String getLocalDirBase() {
		return DEFUALT_LOCAL_DIR_BASE;
	}

	public enum BIN_RESOURCE_TYPE {

		SHAPEIT("shapeit/bin/shapeit",
				"https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz",
				false,
				"shapeit/shapeit.tar.gz"),

		MINIMAC3("Minimac3/bin/Minimac3-omp",
				DEFAULT_URL + BIN_SUB_DIR + "Minimac3/Minimac3.v1.0.14.tar.gz",
				false,
				"Minimac3.tar.gz",
				"Minimac3/");

		private String localSubPath;
		private String url;
		private boolean windows;
		private boolean tarGz;
		private String tarGzSubPath;
		private boolean make;
		private String makeSubDir;

		private BIN_RESOURCE_TYPE(String localSubPath, String url, boolean windows) {
			this(localSubPath, url, windows, null, null);
		}

		private BIN_RESOURCE_TYPE(String localSubPath, String url, boolean windows, String tarGzSubPath) {
			this(localSubPath, url, windows, tarGzSubPath, null);
		}

		/**
		 * @param localSubPath
		 *            local subpath of binary
		 * @param url
		 *            full url to retrieve
		 * @param windows
		 *            true if binary is supported on windows
		 * @param tarGzSubPath
		 *            if resource is archived, path to download tar.gz
		 * @param makeSubDir
		 *            if resource needs to be built, directory to call make from after extracting
		 */
		private BIN_RESOURCE_TYPE(String localSubPath, String url, boolean windows, String tarGzSubPath, String makeSubDir) {
			this.localSubPath = localSubPath;
			this.url = url;
			this.windows = windows;
			this.tarGz = tarGzSubPath != null;
			this.tarGzSubPath = tarGzSubPath;
			this.make = makeSubDir != null;
			this.makeSubDir = makeSubDir;
		}

		public Resource getResource() {
			final String localBinDir = getLocalDirBase() + BIN_SUB_DIR;
			return new Resource(localBinDir + localSubPath, url) {
				public boolean downloadResource(Logger log) {
					if (!windows && Files.isWindows()) {
						log.reportTimeError("Requested binary resource is not supported on Windows");
						return false;
					}
					if (!tarGz)
						return super.downloadResource(log);
					if (!downloadResource(localBinDir + tarGzSubPath, log))
						return false;
					if (tarGz && !extractTarGz(log)) {
						log.reportTimeError(tarGzSubPath + " could not be extracted");
						return false;
					}
					if (make && !makeBinary(log))
						return false;
					return true;
				}

				private boolean extractTarGz(Logger log) {
					// TODO use Apache Commons or other Java utility to allow on Windows and not use command line
					String file = ext.removeDirectoryInfo(tarGzSubPath);
					String dir = ext.parseDirectoryOfFile(localBinDir + tarGzSubPath);
					log.report("Extracting " + file);
					if (!CmdLine.runDefaults("tar -xzf " + file, dir))
						return false;
					log.report("Extracted to " + dir);
					log.report("Removing " + file);
					CmdLine.runDefaults("rm " + file, dir);
					return true;
				}

				private boolean makeBinary(Logger log) {
					String dir = localBinDir + makeSubDir;
					log.report("Building " + dir);
					if (!CmdLine.runDefaults("make -s", dir))
						return false;
					log.report("Built to " + getFullLocalPath());
					return true;
				}
			};
		}

	}

	public enum GENOME_CHROMOSOME_RESOURCE_TYPE {
		GENETIC_MAP("genetic_map_", ".txt.gz", DEFAULT_URL),
		G1K_PHASE3v5_REF_PANEL("1000genomes_ref_panel_Phase3v5_", ".m3vcf.gz", DEFAULT_URL);

		private String namePrefix;
		private String nameSuffix;
		private String url;

		/**
		 * @param namePrefix
		 * @param nameSuffix
		 * @param url
		 */
		private GENOME_CHROMOSOME_RESOURCE_TYPE(String namePrefix, String nameSuffix, String url) {
			this.namePrefix = namePrefix;
			this.nameSuffix = nameSuffix;
			this.url = url;
		}

		public Resource getResource(GENOME_BUILD build, String chr) {
			String resourceSubPath = GENOME_RESOURCE_TYPE.getGenomeBuildSubDir(build) + CHR_SUB_DIR  + namePrefix + build.getBuild() + "_" + "chr" + chr + nameSuffix;
			return new Resource(getLocalDirBase(), resourceSubPath, url) {
			};
		}

	}

	public enum GENOME_RESOURCE_TYPE {
		/**
		 * A gc5base file, for constructing gc-models
		 */
		GC5_BASE("", "_gc5Base.txt", DEFAULT_URL),
		DB_SNP("", "_dbSnp147.vcf.gz", DEFAULT_URL), ;

		private String namePrefix;
		private String nameSuffix;
		private String url;

		/**
		 * @param namePrefix
		 * @param nameSuffix
		 * @param url
		 */
		private GENOME_RESOURCE_TYPE(String namePrefix, String nameSuffix, String url) {
			this.namePrefix = namePrefix;
			this.nameSuffix = nameSuffix;
			this.url = url;
		}

		public Resource getResource(GENOME_BUILD build) {
			String resourceSubPath = getGenomeBuildSubDir(build) + namePrefix + build.getBuild() + nameSuffix;
			switch (this) {
			case DB_SNP:
				return new VCFResource(getLocalDirBase(), resourceSubPath, url);
			default:
				return new Resource(getLocalDirBase(), resourceSubPath, url) {
				};
			}
		}
		
		protected static String getGenomeBuildSubDir (GENOME_BUILD build) {
			return GENOME_SUB_DIR + build.getBuild() + "/";
		}
	}



	public enum MITO_RESOURCE_TYPE {
		/**
		 * A gc5base file, for constructing gc-models
		 */
		WHITE_WBC_BETA("", "Whites_WBC_TOTAL_SingleSNPmatched.final.beta", DEFAULT_URL),
		BLACK_WBC_BETA("", "WBC_TOTAL_SingleSNPmatched.final.beta", DEFAULT_URL),
		ALL_WBC_BETA("", "Blacks_WBC_TOTAL_SingleSNPmatched.final.beta", DEFAULT_URL);

		private String namePrefix;
		private String nameSuffix;
		private String url;

		/**
		 * @param namePrefix
		 * @param nameSuffix
		 * @param url
		 */
		private MITO_RESOURCE_TYPE(String namePrefix, String nameSuffix, String url) {
			this.namePrefix = namePrefix;
			this.nameSuffix = nameSuffix;
			this.url = url;
		}

		public Resource getResource() {
			String resourceSubPath = MITO_SUB_DIR + namePrefix + nameSuffix;
			return new Resource(getLocalDirBase(), resourceSubPath, url) {

			};
		}

	}

	public enum ARRAY_RESOURCE_TYPE {

		/**
		 * Affy Bundle
		 */
		AFFY_SNP6_MARKER_POSITIONS("AffySnp6/", "_markerPositions.txt", DEFAULT_URL, true),
		AFFY_SNP6_HMM("AffySnp6/", "affygw6.hmm", DEFAULT_URL, false),
		AFFY_SNP6_ABLOOKUP("AffySnp6/", "AB_lookup.dat", DEFAULT_URL, false);
		
		/**
		 * Illumina Bundle TODO
		 */
		
		
		private String namePrefix;
		private boolean genomeBuildSpecific;// some Array resources are, some aren't
		private String nameSuffix;
		private String url;

		/**
		 * @param namePrefix
		 * @param nameSuffix
		 * @param url
		 */
		private ARRAY_RESOURCE_TYPE(String namePrefix, String nameSuffix, String url, boolean genomeBuildSpecific) {
			this.namePrefix = namePrefix;
			this.nameSuffix = nameSuffix;
			this.url = url;
			this.genomeBuildSpecific = genomeBuildSpecific;
		}

		public Resource getResource(GENOME_BUILD build) {
			String resourceSubPath = ARRAY_SUB_DIR + namePrefix + (genomeBuildSpecific ? build.getBuild() : "") + nameSuffix;
			return new Resource(getLocalDirBase(), resourceSubPath, url) {
			};
		}
	}

	public enum GENOME_BUILD {

		HG19("hg19", 37),
		HG18("hg18", 36);

		private String build;
		private int buildInt;

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

	public static abstract class Resource {

		private String fullLocalPath;
		private String fullUrl;

		/**
		 * 
		 * @param localPath
		 *            This can be used to create fully qualified locations i.e /home/usr/resources, or relative i.e resources/<br>
		 *            Thinking this will be set by launch properties
		 * @param resourceSubPath
		 *            The path of the resource within local path and url
		 * @param url
		 *            Typically {@link Resources#DEFAULT_URL}
		 */
		private Resource(String localPath, String resourceSubPath, String url) {
			this(localPath + resourceSubPath, url + resourceSubPath);
		}

		/**
		 * 
		 * @param fullLocalPath
		 *            full path to resource on local system
		 * @param fullUrl
		 *            full url path to resource on internet
		 */
		private Resource(String fullLocalPath, String fullUrl) {
			super();
			this.fullLocalPath = fullLocalPath;
			this.fullUrl = fullUrl;
		}

		private boolean isLocallyAvailable() {
			return Files.exists(fullLocalPath);
		}

		public boolean validateWithHint(Logger log) {
			boolean available = isLocallyAvailable() || isRemotelyAvailable(log);
			if (!available) {
				log.reportTimeError("Could not find local file " + fullLocalPath + " and could not download it from " + fullUrl + " please manually download and save to " + fullLocalPath);

			}

			return available;
		}

		private boolean isRemotelyAvailable(Logger log) {
			return HttpDownloadUtility.canDownload(fullUrl, log);
		}

		protected String getFullLocalPath() {
			return fullLocalPath;
		}

		protected boolean downloadResource(String downloadPath, Logger log) {
			if (isRemotelyAvailable(log)) {
				try {
					HttpDownloadUtility.downloadFile(fullUrl, downloadPath, true, log);
					return true;
				} catch (IOException e) {
					log.reportTimeError("Could not retrieve resource from " + downloadPath + " and save it to" + fullLocalPath);
					e.printStackTrace();
				}
			} else
				log.reportTimeError("Resource is not available for download");
			return false;
		}

		public boolean downloadResource(Logger log) {
			return downloadResource(fullLocalPath, log);
		}

		public boolean isAvailable(Logger log) {
			return isLocallyAvailable() || isRemotelyAvailable(log);
		}

		/**
		 * @param log
		 * @return the local path (immediately if available, or after downloading to the local path)
		 */
		public String getResource(Logger log) {
			if (isLocallyAvailable())
				return fullLocalPath;
			log.report("Resource is not available at " + fullLocalPath + ", will attempt to download from " + fullUrl);
			if (downloadResource(log)) {
				if (isLocallyAvailable())
					return fullLocalPath;
				log.reportError("Downloaded resource cannot be found at " + fullLocalPath);
			}
			return null;
		}
	}

	private static class VCFResource extends Resource {
		private Resource index;

		private VCFResource(String fullLocalPath, String fullUrl) {
			super(fullLocalPath, fullUrl);
			if (!fullLocalPath.endsWith(".vcf") && !fullLocalPath.endsWith(".vcf.gz")) {
				throw new IllegalArgumentException("This should only be used for vcf files");
			}
			this.index = new Resource(VCFOps.getIndex(fullLocalPath), VCFOps.getIndex(fullUrl)) {

			};
		}

		private VCFResource(String localPath, String resourceSubPath, String url) {
			this(localPath + resourceSubPath, url + resourceSubPath);
		}

		@Override
		public String getResource(Logger log) {
			index.getResource(log);
			return super.getResource(log);
		}

	}

	private static void test() {
		Logger log = new Logger();

		for (GENOME_RESOURCE_TYPE gType : GENOME_RESOURCE_TYPE.values()) {
			for (GENOME_BUILD gb : GENOME_BUILD.values()) {
				log.reportTimeInfo("Testing " + gType + "\t" + gb);
				Resource gResource = gType.getResource(gb);
				gResource.downloadResource(log);
				log.report("Available? " + gResource.isAvailable(log));
				gResource.getResource(log);
			}
		}

		for (BIN_RESOURCE_TYPE bType : BIN_RESOURCE_TYPE.values()) {
			log.reportTimeInfo("Testing " + bType.toString());
			Resource bResource = bType.getResource();
			bResource.downloadResource(log);
			log.report("Available? " + bResource.isAvailable(log));
			bResource.getResource(log);
		}
		for (ARRAY_RESOURCE_TYPE atype : ARRAY_RESOURCE_TYPE.values()) {
			for (GENOME_BUILD gb : GENOME_BUILD.values()) {
				log.reportTimeInfo("Testing " + atype + "\t" + gb);
				Resource bResource = atype.getResource(gb);
				bResource.downloadResource(log);
				log.report("Available? " + bResource.isAvailable(log));
				bResource.getResource(log);
			}
		}
	}

	public static void main(String[] args) {
		test();
	}
}
