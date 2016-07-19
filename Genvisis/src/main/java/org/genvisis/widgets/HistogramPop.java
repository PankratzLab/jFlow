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
import javafx.scene.control.ToolBar;
import javafx.scene.image.WritableImage;
import javafx.scene.layout.VBox;
import javafx.stage.Stage;
import javafx.stage.StageStyle;

/**
 * Javafx histogram widget
 *
 */
public class HistogramPop extends Application {
	private static final int DEFAULT_BIN = 30;
	private ParseResult parseResult;
	final Slider slider = new Slider(0, 100, 1);
	private ToolBar toolBar;
	private Scene scene;
	private AreaChart<Number, Number> areaChart;
	private VBox root;
	private NumberAxis xAxis;
	private NumberAxis yAxis;

	/**
	 * @return get the data from the clipboard
	 */
	private static String[] getData() {
		return ext.getClipboard().trim().split("\\n");
	}

	/**
	 * Parse the clipboard data
	 */
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

	/**
	 * Little data container
	 *
	 */
	private static class ParseResult {
		private double[][] data;
		private String[] titles;

		private ParseResult(String[] titles, double[][] data) {
			super();
			this.data = data;
			this.titles = titles;
		}
	}

	/**
	 * Using apache commons for speed, can change later
	 */
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

		private BasicHistogram(long[] counts, double[] binMax, double[] binMin) {
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
		slider.valueProperty().addListener(new ChangeListenerImplementation());
		show(stage);

	}

	private void show(Stage stage) {
		this.xAxis = new NumberAxis();
		this.yAxis = new NumberAxis();
		yAxis.setAutoRanging(true);
		xAxis.setAutoRanging(true);

		areaChart = new AreaChart<Number, Number>(xAxis, yAxis);
		updateChart();

		final Button button1 = new Button();
		paramaterizeButton(button1);
		this.toolBar = new ToolBar(button1);
		this.root = new VBox();
		root.setStyle("-fx-background-color: white");
		root.getChildren().addAll(toolBar, areaChart, slider, button1);
		this.scene = new Scene(root, 1600, 600);
		stage.setScene(scene);
		stage.show();

	}

	@SuppressWarnings("rawtypes")
	private final class ChangeListenerImplementation implements ChangeListener {
		private ChangeListenerImplementation() {
		}

		@Override
		public void changed(ObservableValue arg0, Object arg1, Object arg2) {
			updateChart();
		}
	}

	private void updateChart() {
		areaChart.setCreateSymbols(false);
		areaChart.setHorizontalGridLinesVisible(true);
		areaChart.setVerticalGridLinesVisible(false);
		areaChart.setAnimated(false);
		BasicHistogram[] histograms = getHistogram(Math.max(1, (int) slider.getValue()), parseResult.data);
		addDataToChart(histograms, new DecimalFormat("#.##"));
		areaChart.getXAxis().setLabel("Bins " + (int) slider.getValue());

	}

	/**
	 * @param button1
	 *            give this button a screen shot
	 */
	private void paramaterizeButton(final Button button1) {
		button1.setText("Screen Shot");
		button1.setOnAction(new EventHandler<ActionEvent>() {
			@Override
			public void handle(ActionEvent e) {
				button1.setVisible(false);
				slider.setVisible(false);
				toolBar.setVisible(false);
				WritableImage snapShot = scene.snapshot(null);
				String dir = "./histograms/";
				String file = dir + "clipboard_histogram.png";
				int cnt = 1;

				while (Files.exists(file)) {
					file = "./histograms/clipboard_histogram_" + cnt++ + ".png";
				}
				try {
					ImageIO.write(SwingFXUtils.fromFXImage(snapShot, null), "png", new File(file));
				} catch (IOException io) {
					io.printStackTrace();
				}
				toolBar.setVisible(true);
				button1.setVisible(true);
				slider.setVisible(true);

			}
		});
	}

	private void addDataToChart(BasicHistogram[] histograms, DecimalFormat twoDForm) {
		boolean replace = areaChart.getData().size() > 0;
		for (int i = 0; i < histograms.length; i++) {
			BasicHistogram histogram = histograms[i];
			XYChart.Series<Number, Number> countsBin = new XYChart.Series<Number, Number>();
			countsBin.setName(parseResult.titles[i]);
			for (int j = 0; j < histogram.binMax.length; j++) {
				if (Double.isFinite(histogram.binMin[j]) && Double.isFinite(histogram.binMax[j])) {
					countsBin.getData().add(new XYChart.Data<Number, Number>(
							Double.valueOf(twoDForm.format(histogram.binMin[j])), histogram.counts[j]));
					countsBin.getData().add(new XYChart.Data<Number, Number>(
							Double.valueOf(twoDForm.format(histogram.binMax[j])), histogram.counts[j]));
				}
			}
			if (replace && i < areaChart.getData().size()) {
				 areaChart.getData().set(i, countsBin);
				// areaChart.getData().add(countsBin);
			} else {
				areaChart.getData().add(countsBin);
			}
		}
	}

	public static void main(String[] args) {
		launch(args);
	}

}
