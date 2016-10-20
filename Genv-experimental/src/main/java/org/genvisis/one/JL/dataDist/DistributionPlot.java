package org.genvisis.one.JL.dataDist;

import org.genvisis.CLI;
import org.genvisis.cnv.filesys.ABLookup;
import org.genvisis.cnv.filesys.Project;

import javafx.application.Application;

import javafx.scene.Scene;
import javafx.scene.chart.BarChart;
import javafx.scene.chart.CategoryAxis;
import javafx.scene.chart.NumberAxis;

import javafx.stage.Stage;

/**
 * Class aimed at detecting meta-parameters for CNVs via LRR and BAF clusters
 *
 */
public class DistributionPlot {
	private DistributionPlot() {

	}

	/**
	 * @param proj
	 *            Project to look at
	 */
	public static void launch(Project proj) {
		BaCha.proj = proj;
		BaCha.launch(BaCha.class, new String[] { "DF" });
	}

	/**
	 * Widget application
	 *
	 */
	public static class BaCha extends Application {
		private static Project proj;

		/**
		 * Must be a public constructor
		 */
		public BaCha() {
			// called by launch
		}

		@Override
		public void start(Stage primaryStage) throws Exception {
			primaryStage.setTitle(proj.getPropertyFilename());
			final CategoryAxis xAxis = new CategoryAxis();
			final NumberAxis yAxis = new NumberAxis();
			final BarChart<String, Number> bc = new BarChart<String, Number>(xAxis, yAxis);
			primaryStage.setScene(new Scene(bc, 450, 450));
			primaryStage.show();
		}

	}

	public static void main(String[] args) {
		CLI c = new CLI(ABLookup.class);
		c.addArgWithDefault("proj", "project properties file", null);
		c.parseWithExit(args);
		launch(new Project(c.get("proj"), false));
	}

}
