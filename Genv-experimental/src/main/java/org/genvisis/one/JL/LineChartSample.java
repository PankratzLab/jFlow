package org.genvisis.one.JL;

import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

import org.genvisis.common.Array;
import org.genvisis.common.HashVec;
import org.genvisis.common.ext;

import javafx.application.Application;
import javafx.application.Platform;
import javafx.embed.swing.SwingFXUtils;
import javafx.scene.Scene;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
import javafx.scene.chart.XYChart.Series;
import javafx.scene.image.WritableImage;
import javafx.stage.Stage;

public class LineChartSample extends Application {

  // http://docs.oracle.com/javafx/2/charts/line-chart.htm

  public static void main(String[] args) {
    System.setProperty("java.awt.headless", "true");
    // new LineChartSample(new double[] { 0, 4, 2 });
    launch(args);
  }

  // private double[] data;
  //
  //
  // public LineChartSample(double[] data) {
  // super();
  // this.data = data;
  // }
  @SuppressWarnings({"rawtypes", "unchecked"})
  @Override
  public void start(Stage stage) {
    String file = "C:/data/misc/betaPlots/Whites_WBC_TOTAL_SingleSNPmatched.final_beta_summary.txt";
    double[] pcs =
        Array.toDoubleArray(HashVec.loadFileToStringArray(file, true, new int[] {0}, false));
    double[] pvals =
        Array.toDoubleArray(HashVec.loadFileToStringArray(file, true, new int[] {22}, false));

    double[] correl =
        Array.toDoubleArray(HashVec.loadFileToStringArray(file, true, new int[] {15}, false));
    String[] types = HashVec.loadFileToStringArray(file, true, new int[] {21}, false);
    String[] uniqTypes = Array.unique(types);
    XYChart.Series[] all = new XYChart.Series[uniqTypes.length];
    for (int i = 0; i < all.length; i++) {
      if (uniqTypes[i].contains("NATURAL")) {
        all[i] = new XYChart.Series();
        all[i].setName(uniqTypes[i]);
      }
    }

    stage.setTitle("Beta plots");
    final NumberAxis xAxis = new NumberAxis();
    final NumberAxis yAxis = new NumberAxis();
    xAxis.setLabel("PC");
    yAxis.setLabel("Beta Correlation of inverse mito estimates (p<.05)");

    final XYChart<Number, Number> lineChart = new LineChart<Number, Number>(xAxis, yAxis);

    lineChart.setTitle(ext.rootOf(file));

    for (int i = 0; i < pcs.length; i++) {
      if (pvals[i] == 0.05) {
        if (types[i].contains("NATURAL")) {
          all[ext.indexOfStr(types[i], uniqTypes)].getData()
              .add(new XYChart.Data(pcs[i], correl[i]));
        }
      }
    }

    for (Series element : all) {
      if (element != null) {
        lineChart.getData().add(element);
      }
    }

    lineChart.setAnimated(false);
    Scene scene = new Scene(lineChart, 800, 600);

    // stage.setScene(scene);
    // stage.show();
    WritableImage snapShot = scene.snapshot(null);

    try {
      ImageIO.write(SwingFXUtils.fromFXImage(snapShot, null), "png",
          new File("C:/data/misc/betaPlots/test.png"));
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    Platform.exit();
  }
}
