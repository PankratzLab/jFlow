package cnv.manage;

import java.io.IOException;
import java.util.ResourceBundle;

import common.Files;
import common.HttpDownloadUtility;
import common.Logger;

public class Resources {

	public static final String DEFAULT_URL = "http://genvisis.org/rsrc/resources/";
	public static final String DEFUALT_LOCAL_DIR_BASE = "resources/";
	public static final String GENOME_SUB_DIR = "Genome/";

	// TODO, a big TODO
	// need to add web-based download, and local file structure
	// could probably do this like project properties...

	public enum BIN_RESOURCE_TYPE {
		/**
		 * place holder...
		 */
		PLINK;
	}

	public enum GENOME_RESOURCE_TYPE {
		/**
		 * A gc5base file, for constructing gc-models
		 */
		GC5_BASE;
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

	public static Resource getGenomeResource(GENOME_RESOURCE_TYPE type, GENOME_BUILD build) {
		return getGenomeResource(type, DEFUALT_LOCAL_DIR_BASE, DEFAULT_URL, build);
	}

	private static Resource getGenomeResource(GENOME_RESOURCE_TYPE type, String localDirBase, String url, GENOME_BUILD build) {
		switch (type) {
		case GC5_BASE:
			return getGCBaseResource(localDirBase, GENOME_SUB_DIR, url, build);
		default:
			throw new IllegalArgumentException("Invalid resource selection " + type);
		}
	}

	private static Resource getGCBaseResource(String localDirBase, String genomeSubDir, String url, GENOME_BUILD build) {
		Resource gcBaseResource = new Resource(localDirBase, genomeSubDir + build.getBuild() + "/" + build.getBuild() + "_gc5Base.txt", DEFAULT_URL) {
		};
		return gcBaseResource;
	}

	public static abstract class Resource {
		/**
		 * This can be used to create fully qualified locations i.e /home/usr/resources, or relative i.e resources/<br>
		 * Thinking this will be set by launch properties
		 */
		private String localPath;
		/**
		 * The path of the resource within the local path
		 */
		private String resourceSubPath;

		private String fullLocalPath;
		private String fullUrl;

		/**
		 * Typically {@link Resources#DEFAULT_URL}
		 */
		private String url;

		private Resource(String localPath, String resourceSubPath, String url) {
			super();
			this.localPath = localPath;
			this.resourceSubPath = resourceSubPath;
			this.fullLocalPath = localPath + resourceSubPath;
			this.url = url;
			this.fullUrl = url + resourceSubPath;
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
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			return null;
			// else if(url check){
			// download the file to the local path
			// report download in progress
			// }
		}

	}
}
