package cnv.manage;

import java.io.File;
import java.io.IOException;

import common.CmdLine;
import common.Files;
import common.HttpDownloadUtility;
import common.Logger;
import common.ext;

public class Resources {

	public static final String DEFAULT_URL = "http://genvisis.org/rsrc/resources/";
	public static final String DEFUALT_LOCAL_DIR_BASE = "resources/";
	
    public static final String GENOME_SUB_DIR = "Genome/";
    public static final String BIN_SUB_DIR = "bin/";

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
				"shapeit.tar.gz"),
				
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
		 * @param localSubPath local subpath of binary
		 * @param url full url to retrieve
		 * @param windows true if binary is supported on windows
		 * @param tarGzSubPath if resource is archived, path to download tar.gz
		 * @param makeSubDir if resource needs to be built, directory to call make from after extracting
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
					if (!windows && System.getProperty("os.name").startsWith("Windows")) {
						log.reportTimeError("Requested binary resource is not supported on Windows");
						return false;
					}
					if (!tarGz) return super.downloadResource(log);
					if (!downloadResource(localBinDir + tarGzSubPath, log)) return false;
					if (tarGz && !extractTarGz(log)) {
						log.reportTimeError(tarGzSubPath + " could not be extracted");
						return false;
					}
					if (make && !makeBinary(log)) return false;
					return true;
				}
				
				private boolean extractTarGz(Logger log) {
					//TODO use Apache Commons or other Java utility to allow on Windows and not use command line
					String file = ext.removeDirectoryInfo(tarGzSubPath);
					String dir = ext.parseDirectoryOfFile(localBinDir + tarGzSubPath);
					log.report("Extracting " + file);
					if (!CmdLine.runDefaults("tar -xzf " + file, dir)) return false;
					log.report("Extracted to " + dir);
					log.report("Removing " + file);
					CmdLine.runDefaults("rm " + file, dir);
					return true;
				}
				
				private boolean makeBinary(Logger log) {
					String dir = localBinDir + makeSubDir;
					log.report("Building " + dir);
					if (!CmdLine.runDefaults("make -s", dir)) return false;
					log.report("Built to " + getFullLocalPath());
					return true;
				}
			};
		}

	}

	public enum GENOME_RESOURCE_TYPE {
		/**
		 * A gc5base file, for constructing gc-models
		 */
		GC5_BASE("", "_gc5Base.txt" ,DEFAULT_URL);
		
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
			String resourceSubPath = GENOME_SUB_DIR + build.getBuild() + "/" + namePrefix + build.getBuild() + nameSuffix;
			return new Resource(getLocalDirBase(), resourceSubPath, url) { };
		}
		
	}

	public enum ARRAY_RESOURCE_TYPE {
		// Project proj = new Project()
		// SNP6, specific Illumnina, etc
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
		 * @param localPath This can be used to create fully qualified locations i.e /home/usr/resources, or relative i.e resources/<br>
		 * 					Thinking this will be set by launch properties
		 * @param resourceSubPath The path of the resource within local path and url
		 * @param url Typically {@link Resources#DEFAULT_URL}
		 */
		private Resource(String localPath, String resourceSubPath, String url) {
			this(localPath + resourceSubPath, url + resourceSubPath);
		}
		
		/**
		 * 
		 * @param fullLocalPath full path to resource on local system
		 * @param fullUrl full url path to resource on internet
		 */
		private Resource(String fullLocalPath, String fullUrl) {
			super();
			this.fullLocalPath = fullLocalPath;
			this.fullUrl = fullUrl;
		}

		private boolean isLocallyAvailable() {
			return Files.exists(fullLocalPath);
		}

		private boolean isRemotelyAvailable(Logger log) {
			return HttpDownloadUtility.canDownload(fullUrl, log);
		}

		protected String getFullLocalPath() {
			return fullLocalPath;
		}
		
		protected boolean downloadResource(String downloadPath, Logger log) {
			if (isRemotelyAvailable(log)){
				try {
					HttpDownloadUtility.downloadFile(fullUrl, downloadPath, true, log);
					return true;
				} catch (IOException e) {
					log.reportTimeError("Could not retrieve resource from " + downloadPath + " and save it to" + fullLocalPath);
					e.printStackTrace();
				}
			} else log.reportTimeError("Resource is not available for download");
			return false;
		}
		
		public boolean downloadResource(Logger log) {
			return downloadResource(fullUrl, log);
		}

		public boolean isAvailable(Logger log) {
			return isLocallyAvailable() || isRemotelyAvailable(log);
		}

		/**
		 * @param log
		 * @return the local path (immediately if available, or after downloading to the local path)
		 */
		public String getResource(Logger log) {
			if (isLocallyAvailable()) return fullLocalPath;
			log.report("Resource is not available at " + fullLocalPath + ", will attempt to download from " + fullUrl);
			if (downloadResource(log)) {
				if (isLocallyAvailable()) return fullLocalPath;
				log.reportError("Downloaded resource cannot be found at " + fullLocalPath);
			}
			return null;
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
	}

	public static void main(String[] args) {
		test();
	}
}
