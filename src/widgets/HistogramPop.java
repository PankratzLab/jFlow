package widgets;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;

import javax.imageio.ImageIO;

import org.apache.commons.math3.random.EmpiricalDistribution;

import common.Array;
import common.Files;
import common.ext;
import javafx.application.Application;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.embed.swing.SwingFXUtils;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.scene.Scene;
import javafx.scene.chart.AreaChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
import javafx.scene.control.Button;
import javafx.scene.control.Slider;
import javafx.scene.image.WritableImage;
import javafx.scene.layout.VBox;
import javafx.stage.Stage;
import javafx.stage.StageStyle;

public class HistogramPop extends Application {
	private static final int DEFAULT_BIN = 100;

	private ParseResult parseResult;
	final Slider slider = new Slider(1, 100, 1);
	private Scene scene;
	/**
	 * @return get the data from the clipboard
	 */
	private static String[] getData() {
		return ext.getClipboard().trim().split("\\n");
	}

	private ParseResult parseData() {
		String[] data = getData();
		ArrayList<ArrayList<Double>> ds = new ArrayList<ArrayList<Double>>();
		ArrayList<String> titles = new ArrayList<String>();
		for (int i = 0; i < data.length; i++) {
			String[] all = data[i].split("\t");
			for (int j = 0; j < all.length; j++) {
				try {
					double d = Double.parseDouble(all[j]);
					if (Double.isFinite(d)) {
						if (ds.size() <= j) {
							for (int k = ds.size() - 1; k < j; k++) {
								ds.add(new ArrayList<Double>());
							}
						}
						ds.get(j).add(d);
					}
				} catch (NumberFormatException nfe) {
					if (i == 0) {
						titles.add(all[j]);
					} else {
						titles.add("Histogram " + (j + 1));

					}
				}
			}
		}
		double[][] dataParsed = new double[ds.size()][];
		for (int i = 0; i < dataParsed.length; i++) {
			dataParsed[i] = Array.toDoubleArray(ds.get(i));
		}
		return new ParseResult(Array.toStringArray(titles), dataParsed);
	}

	private final class ChangeListenerImplementation implements ChangeListener {
		private final Stage stage;

		private ChangeListenerImplementation(Stage stage) {
			this.stage = stage;
		}

		@Override
		public void changed(@SuppressWarnings("rawtypes") ObservableValue arg0, Object arg1, Object arg2) {
			show(stage, (int) slider.getValue());

		}
	}

	private static class ParseResult {
		private double[][] data;
		private String[] titles;

		public ParseResult(String[] titles, double[][] data) {
			super();
			this.data = data;
			this.titles = titles;
		}
	}

	private static BasicHistogram[] getHistogram(int binCount, double[][] data) {
		BasicHistogram[] basicHistograms = new BasicHistogram[data.length];
		for (int i = 0; i < basicHistograms.length; i++) {
			long[] counts = new long[binCount];
			double[] binMax = new double[binCount];
			double[] binMin = new double[binCount];

			EmpiricalDistribution distribution = new EmpiricalDistribution(binCount);
			distribution.load(data[i]);
			int k = 0;
			for (org.apache.commons.math3.stat.descriptive.SummaryStatistics stats : distribution.getBinStats()) {
				counts[k] = stats.getN();
				binMax[k] = stats.getMax();
				binMin[k] = stats.getMin();
				k++;
			}
			basicHistograms[i] = new BasicHistogram(counts, binMax, binMin);
		}
		return basicHistograms;
	}

	private static class BasicHistogram {
		private long[] counts;
		private double[] binMax;
		private double[] binMin;

		public BasicHistogram(long[] counts, double[] binMax, double[] binMin) {
			super();
			this.counts = counts;
			this.binMax = binMax;
			this.binMin = binMin;
		}

	}

	@SuppressWarnings("unchecked")
	@Override
	public void start(final Stage stage) {
		this.parseResult = parseData();
		stage.initStyle(StageStyle.UNIFIED);
		slider.setShowTickMarks(true);
		slider.setValue(DEFAULT_BIN);
		slider.setShowTickLabels(true);
		slider.setMajorTickUnit(10f);
		slider.setMinorTickCount(10);
		slider.setSnapToTicks(true);
		slider.setBlockIncrement(1);
		slider.valueProperty().addListener(new ChangeListenerImplementation(stage));
		
		
		show(stage, DEFAULT_BIN);

	}

	public void show(Stage stage, int numBins) {
		BasicHistogram[] histograms = getHistogram((int) slider.getValue(), parseResult.data);
		final NumberAxis xAxis = new NumberAxis();
		final NumberAxis yAxis = new NumberAxis();
		xAxis.setLabel("Bins");
		yAxis.setLabel("Counts");
		DecimalFormat twoDForm = new DecimalFormat("#.##");
		final AreaChart<Number, Number> areaChart =
				new AreaChart<Number, Number>(xAxis, yAxis);
		areaChart.setCreateSymbols(true);
		areaChart.setHorizontalGridLinesVisible(false);
		areaChart.setVerticalGridLinesVisible(false);
		for (int i = 0; i < histograms.length; i++) {
			BasicHistogram histogram = histograms[i];
			XYChart.Series<Number, Number> countsBin = new XYChart.Series<Number, Number>();
			countsBin.setName(parseResult.titles[i]);
			for (int j = 0; j < histogram.binMax.length; j++) {
				if (Double.isFinite(histogram.binMin[j]) && Double.isFinite(histogram.binMax[j])) {
					countsBin.getData().add(new XYChart.Data<Number, Number>(Double.valueOf(twoDForm.format(histogram.binMin[j])), histogram.counts[j]));
					countsBin.getData().add(new XYChart.Data<Number, Number>(Double.valueOf(twoDForm.format(histogram.binMax[j])), histogram.counts[j]));
				}
			}
			areaChart.getData().add(countsBin);
		}
		final Button button1 = new Button();
		button1.setText("Screen Shot");
		button1.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent e) {
				button1.setVisible(false);
				WritableImage snapShot = scene.snapshot(null);
				String file = "./histograms/clipboard_histogram.png";
				int cnt = 1;
				while (Files.exists(file)) {
					file = "./histograms/clipboard_histogram_" + cnt++ + ".png";
				}
				try {

					ImageIO.write(SwingFXUtils.fromFXImage(snapShot, null), "png", new File(file));
				} catch (IOException io) {
					// TODO Auto-generated catch block
					io.printStackTrace();
				}
				button1.setVisible(true);
			}
		});
		areaChart.setAnimated(true);
		VBox root = new VBox();
		root.setStyle("-fx-background-color: white");
		root.getChildren().add(areaChart);
		root.getChildren().add(slider);
		root.getChildren().add(button1);
		this.scene = new Scene(root, 1600, 400);

		stage.setScene(scene);
		stage.setScene(scene);
		stage.sizeToScene();
		stage.show();
		
		
	}

	public static void main(String[] args) {

		launch(args);
	}

}
