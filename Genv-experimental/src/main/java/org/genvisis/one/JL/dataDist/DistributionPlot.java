package org.genvisis.one.JL.dataDist;

import java.util.ArrayList;
import java.util.Hashtable;

import org.genvisis.CLI;
import org.genvisis.cnv.filesys.ABLookup;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.Array;
import org.genvisis.common.Matrix;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.LocusSet;
import org.genvisis.stats.BasicHistogram;

import javafx.application.Application;
import javafx.beans.binding.Bindings;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.collections.FXCollections;
import javafx.geometry.Side;
import javafx.scene.Scene;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.ScatterChart;
import javafx.scene.chart.XYChart;
import javafx.scene.control.ComboBox;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Pane;
import javafx.scene.layout.VBox;
import javafx.scene.shape.Circle;

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
		BaCha.launch(BaCha.class);
	}

	/**
	 * Widget application
	 *
	 */
	public static class BaCha extends Application {
		private static Project proj;
		private ScatterChart<Number, Number> lrrScat;
		private ScatterChart<Number, Number> bafScat;
		private Plot plot;

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
			Sample samp = proj.getFullSampleFromRandomAccessFile(sampleName);
			LocusSet<CNVariant> cnvs = CNVariant.loadLocSet(proj.CNV_FILENAMES.getValue()[0], proj.getLog());
			Hashtable<String, LocusSet<CNVariant>> inds = CNVariant.breakIntoInds(cnvs, proj.getLog());
			float[] bafs = samp.getBAFs();
			float[] lrrs = samp.getLRRs();
			Axes axes = new DistributionPlot().new Axes(400, 300, 0, 1, 1, -10, 10, .1);

			plot = new DistributionPlot().new Plot(8, 0.1, axes);
			plot.replot(bafs, lrrs, axes);

			lrrScat = initialize();
			lrrScat.setTitle("LRR values");
			SampleData sampleData = proj.getSampleData(0, false);
			lrrScat.getData().add(update(sampleData, sampleName, inds, bafs, lrrs, subset, false));

			bafScat = initialize();
			bafScat.setTitle("BAF values");

			bafScat.getData().add(update(sampleData, sampleName, inds, bafs, lrrs, subset, true));

			final ComboBox<String> comboBox = new ComboBox<String>(
					FXCollections.observableArrayList(proj.getSamples()));
			comboBox.getSelectionModel().selectedItemProperty().addListener(new ChangeListener<String>() {
				public void changed(ObservableValue<? extends String> observable, String oldValue, String newValue) {
					Sample samp = proj.getFullSampleFromRandomAccessFile(newValue);
					float[] bafs = samp.getBAFs();
					float[] lrrs = samp.getLRRs();
					plot.replot(bafs, lrrs, axes);
					XYChart.Series<Number, Number> newDataLrr = update(sampleData, newValue, inds, bafs, lrrs, subset,
							false);
					lrrScat.setData(FXCollections.observableArrayList());
					lrrScat.getData().add(newDataLrr);

					XYChart.Series<Number, Number> newDataBaf = update(sampleData, newValue, inds, bafs, lrrs, subset,
							true);
					bafScat.setData(FXCollections.observableArrayList());
					bafScat.getData().add(newDataBaf);

				}
			});

			VBox root = new VBox();

			HBox hBox = new HBox();

			hBox.getChildren().addAll(lrrScat, bafScat);
			root.getChildren().addAll(comboBox, hBox, plot);

			primaryStage.setScene(new Scene(root, 1000, 1000));
			primaryStage.show();
		}

		private static ScatterChart<Number, Number> initialize() {
			NumberAxis yAxis = new NumberAxis();
			yAxis.setAutoRanging(true);
			yAxis.setLabel("Log count");
			NumberAxis xAxis = new NumberAxis();
			xAxis.setAutoRanging(true);
			xAxis.setLabel("Value");

			ScatterChart<Number, Number> bc = new ScatterChart<>(xAxis, yAxis);
			bc.setHorizontalGridLinesVisible(true);
			bc.setVerticalGridLinesVisible(false);
			bc.setAnimated(false);
			return bc;
		}

		private XYChart.Series<Number, Number> update(SampleData sampleData, String sampleName,
				Hashtable<String, LocusSet<CNVariant>> inds, float[] bafs, float[] lrrs, int[] subset, boolean baf) {
			String look = sampleData.lookup(sampleName)[1];
			System.out.println(look);
			if (inds.containsKey(look)) {
				System.out.println("DHFD");
			}
			double[] data = Array.subArray(Array.toDoubleArray(baf ? bafs : lrrs), subset);
			return getHistogram(data, sampleName);
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

	class Axes extends Pane {
		private NumberAxis xAxis;
		private NumberAxis yAxis;

		public Axes(int width, int height, double xLow, double xHi, double xTickUnit, double yLow, double yHi,
				double yTickUnit) {
			setMinSize(Pane.USE_PREF_SIZE, Pane.USE_PREF_SIZE);
			setPrefSize(width, height);
			setMaxSize(Pane.USE_PREF_SIZE, Pane.USE_PREF_SIZE);

			xAxis = new NumberAxis(xLow, xHi, xTickUnit);
			xAxis.setSide(Side.BOTTOM);
			xAxis.setMinorTickVisible(false);
			xAxis.setPrefWidth(width);
			xAxis.setLayoutY((double) height / 2);

			yAxis = new NumberAxis(yLow, yHi, yTickUnit);
			yAxis.setSide(Side.LEFT);
			yAxis.setMinorTickVisible(false);
			yAxis.setPrefHeight(height);

			yAxis.layoutXProperty().bind(Bindings.subtract((width / 2) + 1, yAxis.widthProperty()));

			getChildren().setAll(xAxis, yAxis);
		}

		public NumberAxis getXAxis() {
			return xAxis;
		}

		public NumberAxis getYAxis() {
			return yAxis;
		}
	}

	class Plot extends Pane {
		public Plot(double xMax, double xInc, Axes axes) {
			setMinSize(Pane.USE_PREF_SIZE, Pane.USE_PREF_SIZE);
			setPrefSize(axes.getPrefWidth(), axes.getPrefHeight());
			setMaxSize(Pane.USE_PREF_SIZE, Pane.USE_PREF_SIZE);

		}

		private void replot(float[] xs, float[] ys, Axes axes) {
			ArrayList<Circle> circles = new ArrayList<>();
			double[][] data = new double[][] { Array.toDoubleArray(xs), Array.toDoubleArray(ys) };
			double[][] dist = Matrix.transpose(data);

			for (int i = 0; i < ys.length; i++) {
				if (Math.random() < .10) {
					Circle circle = new Circle(mapX(xs[i], axes), mapY(ys[i], axes), 4);

					circles.add(circle);
				}

			}

			getChildren().setAll(axes);
			getChildren().addAll(circles);
		}

		private double mapX(double x, Axes axes) {
			double tx = axes.getPrefWidth() / 2;
			double sx = axes.getPrefWidth() / (axes.getXAxis().getUpperBound() - axes.getXAxis().getLowerBound());

			return x * sx + tx;
		}

		private double mapY(double y, Axes axes) {
			double ty = axes.getPrefHeight() / 2;
			double sy = axes.getPrefHeight() / (axes.getYAxis().getUpperBound() - axes.getYAxis().getLowerBound());

			return -y * sy + ty;
		}
	}

	public static void main(String[] args) {
		CLI c = new CLI(ABLookup.class);
		c.addArgWithDefault("proj", "project properties file", null);
		c.parseWithExit(args);
		launch(new Project(c.get("proj"), false));
	}

}
