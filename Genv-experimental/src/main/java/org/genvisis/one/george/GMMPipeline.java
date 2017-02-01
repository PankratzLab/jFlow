package org.genvisis.one.george;

//import org.apache.spark.ml.clustering.GaussianMixture;
//import org.apache.spark.ml.clustering.GaussianMixtureModel;
//import org.apache.spark.sql.Dataset;
//import org.apache.spark.sql.Row;
//import org.apache.spark.sql.SparkSession;

public class GMMPipeline {
	String fcs_dirpath;
	
	public GMMPipeline(String fcs_dirpath) {
		this.fcs_dirpath = fcs_dirpath;
	}
	// Fits GMM on fcs datasets converted to csvs
	public void run() {
//		SparkSession spark = SparkSession
//				.builder()
//				.appName("cytometry analysis pipeline")
//				.config("spark.some.config.option", "some-value")
//				.getOrCreate();
		
		StringBuilder out = new StringBuilder();
		
		// get list of csv data sets
		Tuple<String[],String[]> tup = FileManipulator.filesWithExt(".csv", fcs_dirpath);
		String[] fileListExt = tup.x; // just filenames
		String[] pathListExt = tup.y; // file paths
		
		// run gmm model
		for (int i = 0; i < pathListExt.length; i++) {
//			Dataset<Row> dataset = spark.read().csv(pathListExt[i]);
			
			// choose model with lowest bic NOT IMPLEMENTED
			for (int n_components = 3; n_components < 8; n_components++) {
//				GaussianMixture mix = new GaussianMixture().setK(n_components);
//				GaussianMixtureModel model = mix.fit(dataset);
				
				// bic = -2 * num_samples * average P(X | Z)  + num free parameters * log(num_samples)
				// where num free parameters = n_components * n_features * (n_features + 1) / 2
			}
			
		}
		// plot
	}

	public static void main(String[] args) {
		GMMPipeline pipeline = new GMMPipeline("./fcs");
		pipeline.run();
	}
}
