package cnv.manage;

import java.io.File;
import java.io.IOException;

import common.Files;
import common.HttpDownloadUtility;
import common.Logger;

public class Resources {

	public static final String DEFAULT_URL = "http://genvisis.org/rsrc/resources/";
	public static final String DEFUALT_LOCAL_DIR_BASE = "resources/";

	// TODO, a big TODO
	// need to add web-based download, and local file structure
	// could probably do this like project properties...

	public static String getLocalDirBase() {
		return DEFUALT_LOCAL_DIR_BASE;
	}
	
	
	public enum BIN_RESOURCE_TYPE {
		 
		SHAPEIT,
		MINIMAC3;
		
		//TODO figure out how bin fits in
		
//		private static final String BIN_SUB_DIR = "bin/";
//		
//		private String localSubPath;
//		private String webSubPath;
//		private String url;

	}

	public enum GENOME_RESOURCE_TYPE {
		/**
		 * A gc5base file, for constructing gc-models
		 */
		GC5_BASE("", "_gc5Base.txt" ,DEFAULT_URL);
		
		private static final String GENOME_SUB_DIR = "Genome/";
		
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

		public boolean isAvailable(Logger log) {
			return isLocallyAvailable() || isRemotelyAvailable(log);
		}

		/**
		 * @param log
		 * @return the local path (immediately if available, or after downloading to the local path)
		 */
		public String getResource(Logger log) {
			if (Files.exists(fullLocalPath)) {
				return fullLocalPath;
			} else if (HttpDownloadUtility.canDownload(fullUrl, log)) {
				try {
					HttpDownloadUtility.downloadFile(fullUrl, fullLocalPath, true, log);
					return fullLocalPath;
				} catch (IOException e) {
					log.reportTimeError("Could not retrieve file from " + fullUrl + " and save it to" + fullLocalPath);
					e.printStackTrace();
				}
			}
			log.reportTimeError("Could not find "+fullLocalPath+" and could not download from "+fullUrl);

			return null;

		}

	}

	private static void test(String tmpDir) {
		new File(tmpDir).delete();
		new File(tmpDir).mkdirs();
		Logger log = new Logger();

		for (GENOME_RESOURCE_TYPE gType : GENOME_RESOURCE_TYPE.values()) {
			for (GENOME_BUILD gb : GENOME_BUILD.values()) {
				log.reportTimeInfo(gType + "\t" + gb);
				Resource gResource = gType.getResource(gb);
				System.out.println(gResource.isAvailable(log));
				gResource.getResource(log);
			}
		}
	}

	public static void main(String[] args) {
		String outDir = "C:/tmpresourceExtraStuff/";
		test(outDir);
	}
}
