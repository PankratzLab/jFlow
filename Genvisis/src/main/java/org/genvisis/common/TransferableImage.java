package org.genvisis.common;

import java.awt.Image;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.io.IOException;

public class TransferableImage implements Transferable {

  Image i;

  public TransferableImage(Image i) {
    this.i = i;
  }

  @Override
  public Object getTransferData(DataFlavor flavor) throws UnsupportedFlavorException, IOException {
    if (flavor.equals(DataFlavor.imageFlavor) && i != null) {
      return i;
    } else {
      throw new UnsupportedFlavorException(flavor);
    }
  }

  @Override
  public DataFlavor[] getTransferDataFlavors() {
    DataFlavor[] flavors = new DataFlavor[1];
    flavors[0] = DataFlavor.imageFlavor;
    return flavors;
  }

  @Override
  public boolean isDataFlavorSupported(DataFlavor flavor) {
    DataFlavor[] flavors = getTransferDataFlavors();
    for (DataFlavor flavor2 : flavors) {
      if (flavor.equals(flavor2)) {
        return true;
      }
    }

    return false;
  }
}
