package org.genvisis.one.JL.dataDist;

import org.genvisis.CLI;
import org.genvisis.cnv.filesys.ABLookup;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.common.Array;
import org.genvisis.stats.BasicHistogram;

import javafx.application.Application;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.collections.FXCollections;
import javafx.scene.Scene;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.ScatterChart;
import javafx.scene.chart.XYChart;
import javafx.scene.control.ComboBox;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;
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
		private ScatterChart<Number, Number> lrrScat;
		private ScatterChart<Number, Number> bafScat;

		/**
		 * Must be a public constructor
		 */
		public BaCha() {
			// called by launch
		}

		@Override
		public void start(Stage primaryStage) throws Exception {
			primaryStage.setTitle(proj.PROJECT_NAME.getValue());
			String sampleName = proj.getSamples()[10];
			int[] subset = proj.getAutosomalMarkerIndices();

			lrrScat = initialize();
			lrrScat.setTitle("LRR values");
			lrrScat.getData().add(update(sampleName, subset, false));

			bafScat = initialize();
			bafScat.setTitle("BAF values");

			bafScat.getData().add(update(sampleName, subset, true));

			final ComboBox<String> comboBox = new ComboBox<String>(
					FXCollections.observableArrayList(proj.getSamples()));
			comboBox.getSelectionModel().selectedItemProperty().addListener(new ChangeListener<String>() {
				public void changed(ObservableValue<? extends String> observable, String oldValue, String newValue) {
					XYChart.Series<Number, Number> newDataLrr = update(newValue, subset, false);
					lrrScat.setData(FXCollections.observableArrayList());
					lrrScat.getData().add(newDataLrr);

					XYChart.Series<Number, Number> newDataBaf = update(newValue, subset, true);
					bafScat.setData(FXCollections.observableArrayList());
					bafScat.getData().add(newDataBaf);

				}
			});
			VBox root = new VBox();

			HBox hBox = new HBox();

			hBox.getChildren().addAll(lrrScat, bafScat);
			root.getChildren().addAll(comboBox, hBox);

			primaryStage.setScene(new Scene(root, 1000, 1000));
			primaryStage.show();
		}

		private static ScatterChart<Number, Number> initialize() {
			NumberAxis yAxis = new NumberAxis();
			yAxis.setAutoRanging(true);
			NumberAxis xAxis = new NumberAxis();
			xAxis.setAutoRanging(true);
			ScatterChart<Number, Number> bc = new ScatterChart<>(xAxis, yAxis);
			bc.setHorizontalGridLinesVisible(true);
			bc.setVerticalGridLinesVisible(false);
			bc.setAnimated(false);
			return bc;
		}

		private XYChart.Series<Number, Number> update(String sampleName, int[] subset, boolean baf) {
			Sample samp = proj.getFullSampleFromRandomAccessFile(sampleName);
			float[] dataF = baf ? samp.getBAFs() : samp.getLRRs();
			double[] data = Array.subArray(Array.toDoubleArray(dataF), subset);
			return getHistogram(data, samp.getSampleName());
		}

		private static XYChart.Series<Number, Number> getHistogram(double[] data, String name) {

			XYChart.Series<Number, Number> countsBin = new XYChart.Series<>();
			countsBin.setName(name);
			BasicHistogram histogram = BasicHistogram.getHistogram(100, data);
			for (int j = 0; j < histogram.getCounts().length; j++) {
				double y = Math.log(histogram.getCounts()[j]);
				if (!Double.isFinite(y)) {
					y = 0;
				}
				if (Double.isFinite(histogram.getBinMin()[j])) {
					countsBin.getData().add(new XYChart.Data<Number, Number>(histogram.getBinMin()[j], y));
					countsBin.getData().add(new XYChart.Data<Number, Number>(histogram.getBinMax()[j], y));
				}
			}
			return countsBin;
		}

	}

	public static void main(String[] args) {
		CLI c = new CLI(ABLookup.class);
		c.addArgWithDefault("proj", "project properties file", null);
		c.parseWithExit(args);
		launch(new Project(c.get("proj"), false));
	}

}
