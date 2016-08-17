package org.genvisis.one.JL;

import java.io.BufferedReader;

import org.controlsfx.control.table.TableFilter;
import org.genvisis.common.Files;

import javafx.application.Application;
import javafx.beans.property.DoubleProperty;
import javafx.beans.property.SimpleDoubleProperty;
import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.StringProperty;
import javafx.beans.value.ObservableValue;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.Scene;
import javafx.scene.control.Label;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableColumn.CellDataFeatures;
import javafx.scene.control.TableView;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.HBox;
import javafx.stage.Stage;
import javafx.util.Callback;

public class testTableView extends Application {

  public class Item {
    private final StringProperty name = new SimpleStringProperty(this, "name");
    private final DoubleProperty value1 = new SimpleDoubleProperty(this, "value1");
    private final DoubleProperty value2 = new SimpleDoubleProperty(this, "value2");

    public Item(String name, double value1, double value2) {
      setName(name);
      setValue1(value1);
      setValue2(value2);
    }


    public final java.lang.String getName() {
      return nameProperty().get();
    }

    public final double getValue1() {
      return value1Property().get();
    }

    public final double getValue2() {
      return value2Property().get();
    }

    public final StringProperty nameProperty() {
      return name;
    }

    public final void setName(final java.lang.String name) {
      nameProperty().set(name);
    }

    public final void setValue1(final double value1) {
      value1Property().set(value1);
    }

    public final void setValue2(final double value2) {
      value2Property().set(value2);
    }

    public final DoubleProperty value1Property() {
      return value1;
    }

    public final DoubleProperty value2Property() {
      return value2;
    }


  }

  public static void main(String[] args) {
    launch(args);
  }

  private TableColumn<ObservableList<StringProperty>, String> createColumn(final int columnIndex,
                                                                           String columnTitle) {
    TableColumn<ObservableList<StringProperty>, String> column = new TableColumn<>();
    String title;
    if (columnTitle == null || columnTitle.trim().length() == 0) {
      title = "Column " + (columnIndex + 1);
    } else {
      title = columnTitle;
    }
    column.setText(title);
    column.setCellValueFactory(new Callback<TableColumn.CellDataFeatures<ObservableList<StringProperty>, String>, ObservableValue<String>>() {
      @Override
      public ObservableValue<String> call(CellDataFeatures<ObservableList<StringProperty>, String> cellDataFeatures) {
        ObservableList<StringProperty> values = cellDataFeatures.getValue();
        if (columnIndex >= values.size()) {
          return new SimpleStringProperty("");
        } else {
          return cellDataFeatures.getValue().get(columnIndex);
        }
      }
    });
    return column;
  }

  private BufferedReader getReaderFromUrl(String urlSpec) throws Exception {

    return Files.getAppropriateReader(urlSpec);
  }

  private void populateTable(final TableView<ObservableList<StringProperty>> table,
                             final String urlSpec, final boolean hasHeader) throws Exception {
    table.getItems().clear();
    table.getColumns().clear();
    table.setPlaceholder(new Label("Loading..."));
    // Task<Void> task = new Task<Void>() {
    // @Override
    // protected Void call() throws Exception {
    BufferedReader in = getReaderFromUrl(urlSpec);
    // Header line
    if (hasHeader) {
      final String headerLine = in.readLine();
      final String[] headerValues = headerLine.split("\t");
      // Platform.runLater(new Runnable() {
      // @Override
      // public void run() {
      for (int column = 0; column < headerValues.length; column++) {
        table.getColumns().add(createColumn(column, headerValues[column]));
      }

      // }
      // });
    }

    // Data:

    String dataLine;
    while ((dataLine = in.readLine()) != null) {
      final String[] dataValues = dataLine.split("\t");
      // Platform.runLater(new Runnable() {
      // @Override
      // public void run() {
      // Add additional columns if necessary:
      for (int columnIndex =
          table.getColumns().size(); columnIndex < dataValues.length; columnIndex++) {
        table.getColumns().add(createColumn(columnIndex, ""));
      }
      // Add data to table:
      ObservableList<StringProperty> data = FXCollections.observableArrayList();
      for (String value : dataValues) {
        data.add(new SimpleStringProperty(value));
      }
      table.getItems().add(data);
    }
    // });
    // }

    // return null;
    // }
    // };
    // Thread thread = new Thread(task);
    // thread.setDaemon(true);
    // thread.start();

  }

  @Override
  public void start(Stage primaryStage) {
    final BorderPane root = new BorderPane();
    final TableView<ObservableList<StringProperty>> table = new TableView<>();
    // final TextField urlTextEntry = new TextField();
    // urlTextEntry.setPromptText("Enter path of tab delimited file");
    // final CheckBox headerCheckBox = new CheckBox("Data has header line");
    // urlTextEntry.setOnAction(new EventHandler<ActionEvent>() {
    // @Override
    // public void handle(ActionEvent event) {
    // true);
    // }
    // });

    try {
      populateTable(table, "D:/data/LLFS_GWAS/sexCheck.xln", true);
      TableFilter.forTableView(table).apply();

    } catch (Exception e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }

    HBox controls = new HBox();
    // controls.getChildren().addAll(urlTextEntry, headerCheckBox);
    // HBox.setHgrow(urlTextEntry, Priority.ALWAYS);
    // HBox.setHgrow(headerCheckBox, Priority.NEVER);
    root.setTop(controls);
    root.setCenter(table);
    Scene scene = new Scene(root, 600, 400);
    primaryStage.setScene(scene);
    primaryStage.show();

  }
}
