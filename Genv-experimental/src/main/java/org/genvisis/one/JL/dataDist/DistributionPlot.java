package org.genvisis.one.JL.dataDist;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;

import org.genvisis.CLI;
import org.genvisis.cnv.filesys.ABLookup;
import org.genvisis.cnv.filesys.MarkerSet;
import org.genvisis.cnv.filesys.Project;
import org.genvisis.cnv.filesys.Sample;
import org.genvisis.cnv.var.SampleData;
import org.genvisis.common.ArrayUtils;
import org.genvisis.filesys.CNVariant;
import org.genvisis.filesys.LocusSet;
import org.genvisis.filesys.Segment;
import org.genvisis.stats.BasicHistogram;

import com.google.common.primitives.Ints;

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
		private ScatterChart<Number, Number> bafVLrrScat;

		/**
		 * Must be a public constructor
		 */
		public BaCha() {
			// called by launch
		}

		@Override
		public void start(Stage primaryStage) throws Exception {
			primaryStage.setTitle(proj.PROJECT_NAME.getValue());
			MarkerSet markerSet = proj.getMarkerSet();
			LocusSet<CNVariant> cnvs = CNVariant.loadLocSet(proj.CNV_FILENAMES.getValue()[0], proj.getLog());
			final Hashtable<String, LocusSet<CNVariant>> inds = CNVariant.breakIntoInds(cnvs, proj.getLog());

			lrrScat = initialize("LRR bin", "Log count");
			lrrScat.setTitle("LRR values");
			final SampleData sampleData = proj.getSampleData(0, false);

			bafScat = initialize("BAF bin", "Log count");
			bafScat.setTitle("BAF values");

			bafVLrrScat = initialize("BAF", "LRR");
			bafVLrrScat.setTitle("BAF v LRR");

			final ComboBox<String> comboBox = new ComboBox<>(
					FXCollections.observableArrayList(proj.getSamples()));
			comboBox.getSelectionModel().selectedItemProperty().addListener(new ChangeListener<String>() {
				public void changed(ObservableValue<? extends String> observable, String oldValue, String newValue) {
					Sample samp = proj.getFullSampleFromRandomAccessFile(newValue);
					float[] bafs = samp.getBAFs();
					float[] lrrs = samp.getLRRs();
					LocusSet<CNVariant> ind = inds.get(sampleData.lookup(newValue)[1]);
					ArrayList<Integer> copyNumbers = new ArrayList<>();

					if (ind != null) {
						HashSet<Integer> markersInCnv = new HashSet<>();
						for (int i = 0; i < markerSet.getChrs().length; i++) {
							if (markerSet.getChrs()[i] > 0 && markerSet.getChrs()[i] < 23) {
								Segment seg = new Segment(markerSet.getChrs()[i], markerSet.getPositions()[i]-1,
										markerSet.getPositions()[i] + 1);
								int[] olaps = ind.getOverlappingIndices(seg);
								if (olaps != null && olaps.length > 0) {
									markersInCnv.add(i);
									copyNumbers.add(ind.getLoci()[olaps[0]].getCN());
								}
							}
						}
						int[] finals = Ints.toArray(markersInCnv);
						proj.getLog().report(newValue+"\t"+ind.getLoci().length + " cnvs over " + finals.length + " markers");
						bafs = ArrayUtils.subArray(bafs, finals);
						lrrs = ArrayUtils.subArray(lrrs, finals);
						if (copyNumbers.size() != bafs.length) {
							throw new IllegalStateException();
						}

					}
					lrrScat.setData(FXCollections.observableArrayList());
					bafScat.setData(FXCollections.observableArrayList());
					bafVLrrScat.setData(FXCollections.observableArrayList());

					if (ind != null) {
						XYChart.Series<Number, Number> newDataLrr = update(sampleData, newValue, bafs, lrrs, false);
						XYChart.Series<Number, Number> newDataBaf = update(sampleData, newValue, bafs, lrrs, true);
						bafScat.getData().add(newDataBaf);
						lrrScat.getData().add(newDataLrr);

						ArrayList<XYChart.Series<Number, Number>> cnSeries = new ArrayList<>();
						for (int i = 0; i < 5; i++) {
							cnSeries.add(new XYChart.Series<Number, Number>());
							cnSeries.get(i).setName("CN-" + i);
						}
						if (!copyNumbers.isEmpty()) {
							for (int j = 0; j < lrrs.length; j++) {
								if (Double.isFinite(bafs[j]) && Double.isFinite(lrrs[j])) {
									cnSeries.get(copyNumbers.get(j)).getData()
											.add(new XYChart.Data<Number, Number>(bafs[j], lrrs[j]));
								}
							}
						}
						bafVLrrScat.getData().addAll(cnSeries);
					}
				}
			});

			VBox root = new VBox();
			HBox hBox = new HBox();

			hBox.getChildren().addAll(lrrScat, bafScat);
			root.getChildren().addAll(comboBox, hBox, bafVLrrScat);

			primaryStage.setScene(new Scene(root, 1000, 1000));
			primaryStage.show();
		}

		private static ScatterChart<Number, Number> initialize(String x, String y) {
			NumberAxis yAxis = new NumberAxis();
			yAxis.setAutoRanging(true);
//			yAxis.setUpperBound(1);
//			yAxis.setLowerBound(-1);

			yAxis.setLabel(y);
			NumberAxis xAxis = new NumberAxis();
			xAxis.setAutoRanging(true);
			xAxis.setLabel(x);

			ScatterChart<Number, Number> bc = new ScatterChart<>(xAxis, yAxis);
			bc.setHorizontalGridLinesVisible(true);
			bc.setVerticalGridLinesVisible(false);
			bc.setAnimated(false);
			return bc;
		}

		private XYChart.Series<Number, Number> update(SampleData sampleData, String sampleName, float[] bafs,
				float[] lrrs, boolean baf) {

			double[] data = ArrayUtils.toDoubleArray(baf ? bafs : lrrs);
			return getHistogram(data, sampleName);
		}

		private static XYChart.Series<Number, Number> getHistogram(double[] data, String name) {

			XYChart.Series<Number, Number> countsBin = new XYChart.Series<Number, Number>();
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
		c.addArgWithDefault(CLI.ARG_PROJ, CLI.DESC_PROJ, null);
		c.parseWithExit(args);
		launch(new Project(c.get(CLI.ARG_PROJ), false));
	}

}
