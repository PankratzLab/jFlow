package org.genvisis.one.JL;

import java.util.List;
import java.util.Map;

import org.genvisis.cnv.manage.temp;

import javafx.application.Application;
import javafx.event.EventHandler;
import javafx.scene.Cursor;
import javafx.scene.Node;
import javafx.scene.Scene;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.ScatterChart;
import javafx.scene.chart.XYChart;
import javafx.scene.control.ContextMenu;
import javafx.scene.effect.DropShadow;
import javafx.scene.input.MouseButton;
import javafx.scene.input.MouseEvent;
import javafx.scene.layout.VBox;
import javafx.stage.Stage;
import javafx.stage.StageStyle;

public class TestPlot extends Application {
	private DropShadow ds = new DropShadow();
	private ContextMenu contextMenu;

	private XYChart.Series selectedSeries;

	@Override
	public void start(Stage stage) {
		stage.setTitle("Scatter Chart Sample");
		stage.initStyle(StageStyle.UNIFIED);
		final NumberAxis xAxis = new NumberAxis(0, 10, 1);
		final NumberAxis yAxis = new NumberAxis(-100, 500, 100);
		final ScatterChart<Number, Number> sc = new
				ScatterChart<Number, Number>(xAxis, yAxis);
		xAxis.setLabel("Age (years)");
		yAxis.setLabel("Returns to date");
		sc.setTitle("Investment Overview");
		Parameters parameters = getParameters();
		Map<String, String> namedParameters = parameters.getNamed();
		List<String> rawArguments = parameters.getRaw();
		List<String> unnamedParameters = parameters.getUnnamed();

		System.out.println("\nnamedParameters -");
		for (Map.Entry<String, String> entry : namedParameters.entrySet())
			System.out.println(entry.getKey() + " : " + entry.getValue());

		System.out.println("\nrawArguments -");
		for (String raw : rawArguments)
			System.out.println(raw);

		System.out.println("\nunnamedParameters -");
		for (String unnamed : unnamedParameters)
			System.out.println(unnamed);
		
		XYChart.Series<Number,Number> series1 = new XYChart.Series<Number,Number>();
		series1.setName("Equities");
		series1.getData().add(new XYChart.Data<Number,Number>(4.2, 193.2));
		series1.getData().add(new XYChart.Data<Number,Number>(2.8, 33.6));
		series1.getData().add(new XYChart.Data<Number,Number>(6.2, 24.8));
		series1.getData().add(new XYChart.Data<Number,Number>(1, 14));
		series1.getData().add(new XYChart.Data<Number,Number>(1.2, 26.4));
		series1.getData().add(new XYChart.Data(4.4, 114.4));
		series1.getData().add(new XYChart.Data(8.5, 323));
		series1.getData().add(new XYChart.Data(6.9, 289.8));
		series1.getData().add(new XYChart.Data(9.9, 287.1));
		series1.getData().add(new XYChart.Data(0.9, -9));
		series1.getData().add(new XYChart.Data(3.2, 150.8));
		series1.getData().add(new XYChart.Data(4.8, 20.8));
		series1.getData().add(new XYChart.Data(7.3, -42.3));
		series1.getData().add(new XYChart.Data(1.8, 81.4));
		series1.getData().add(new XYChart.Data(7.3, 110.3));
		series1.getData().add(new XYChart.Data(2.7, 41.2));

		XYChart.Series<Number,Number> series2 = new XYChart.Series<Number,Number>();
		series2.setName("Mutual funds");
		series2.getData().add(new XYChart.Data(5.2, 229.2));
		series2.getData().add(new XYChart.Data(2.4, 37.6));
		series2.getData().add(new XYChart.Data(3.2, 49.8));
		series2.getData().add(new XYChart.Data(1.8, 134));
		series2.getData().add(new XYChart.Data(3.2, 236.2));
		series2.getData().add(new XYChart.Data(7.4, 114.1));
		series2.getData().add(new XYChart.Data(3.5, 323));
		series2.getData().add(new XYChart.Data(9.3, 29.9));
		series2.getData().add(new XYChart.Data(8.1, 287.4));

		sc.getData().addAll(series1, series2);
		sc.setStyle("-fx-background-color: white");
		Scene scene = new Scene(sc, 500, 400);
		VBox root = new VBox();
		root.setStyle("-fx-background-color: white");

		stage.setScene(scene);
		stage.show();
	}

	private void applyMouseEvents(final XYChart.Series series) {

		final Node node = series.getNode();

		node.setOnMouseEntered(new EventHandler<MouseEvent>() {

			@Override
			public void handle(MouseEvent arg0) {
				node.setEffect(ds);
				node.setCursor(Cursor.HAND);
			}
		});

		node.setOnMouseExited(new EventHandler<MouseEvent>() {

			@Override
			public void handle(MouseEvent arg0) {
				node.setEffect(null);
				node.setCursor(Cursor.DEFAULT);
			}
		});

		node.setOnMouseReleased(new EventHandler<MouseEvent>() {

			@Override
			public void handle(MouseEvent mouseEvent) {
				if (mouseEvent.getButton().equals(MouseButton.SECONDARY)) {
					contextMenu.show(node, mouseEvent.getScreenX() + 1, mouseEvent.getScreenY() + 1);
					// Set as selected
					selectedSeries = series;
					System.out.println("Selected Series data " + selectedSeries.getData());
				}
			}
		});
	}

	public static void test(String[] args) {
		launch(args);
	}

	public static void main(String[] args) {
		test(args);
	}
}