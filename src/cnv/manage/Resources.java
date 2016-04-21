package cnv.manage;

import cnv.filesys.Project;
import common.Files;

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
		Resource gcBaseResource = new Resource(localDirBase + genomeSubDir + build.getBuild() + "/" + build.getBuild() + "_gc5Base.txt", DEFAULT_URL) {
		};
		return gcBaseResource;
	}

	public static abstract class Resource {
		private String localPath;
		private String url;

		private Resource(String localPath, String url) {
			super();
			this.localPath = localPath;
			this.url = url;
		}

		public boolean isAvailable() {
			return Files.exists(localPath);// TODO, check of URL dl
		}

		public String getResource() {
			if (Files.exists(localPath)) {
				return localPath;
			}
			return null;
			// else if(url check){
			// download the file to the local path
			// report download in progress
			// }
		}

	}
}
