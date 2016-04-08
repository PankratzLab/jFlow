package cnv.manage;

import common.Files;

public class Resources {
	// TODO, a big TODO
	// need to add web-based download, and local file structure
	// could probably do this like project properties...

	public enum GENOME_BUILD {

		HG19("hg19"),
		HG18("hg18");

		private String build;

		private GENOME_BUILD(String build) {
			this.build = build;
		}

		public String getBuild() {
			return build;
		}

	}

	public enum ARRAY_TYPE {
		// SNP6, specific Illumnina, etc
	}

	public static Resource getGCBaseResource(GENOME_BUILD build) {
		Resource gcBaseResource = new Resource("resources/Genome/" + build.getBuild() + "/" + build.getBuild() + "_gc5Base.txt", null) {
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
			//
			// }
		}

	}
}
