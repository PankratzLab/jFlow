package org.genvisis.one.ben.imagetag;

import javafx.application.Application;
import javafx.application.Platform;
import javafx.scene.control.Alert;
import javafx.scene.control.Alert.AlertType;
import javafx.scene.control.ButtonType;
import javafx.stage.Stage;

/**
 * JavaFX entry point
 */
public class ImageAnnotatorJavaFXEntry extends Application {

  @Override
  public void start(Stage primaryStage) throws Exception {
    ImageAnnotator.launch();
    // Exit JavaFX
    Platform.exit();
  }

  public static void main(String... args) {
    try {
      Application.launch(ImageAnnotatorJavaFXEntry.class, args);
    } catch (Exception e) {
      e.printStackTrace();
      Platform.runLater(() -> {
        Alert alert = new Alert(AlertType.ERROR, e.getMessage(), ButtonType.CLOSE);
        alert.showAndWait();
      });
    }
  }
}
