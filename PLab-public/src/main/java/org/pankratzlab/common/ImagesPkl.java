package org.pankratzlab.common;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;

import javax.imageio.ImageIO;

public class ImagesPkl {

  public static void stitchImages(String dir, String listFile, String outFile, int numImagesPerRow,
                                  Color bgColor, boolean drawInnerBorder,
                                  boolean drawOuterBorder/* , boolean pack */) {
    System.out.print("Loading input file...");
    String list = Files.isRelativePath(listFile) ? dir + listFile : listFile;
    String[] imageFiles = HashVec.loadFileToStringArray(list, false, new int[] {0}, false);
    System.out.println("Complete!");
    stitchImages(dir, imageFiles, outFile, numImagesPerRow, bgColor, drawInnerBorder,
                 drawOuterBorder);
  }

  public static BufferedImage stitchImages(String[] imageFilesWithPaths, int numImagesPerRow,
                                           Color bgColor, boolean drawInnerBorder,
                                           boolean drawOuterBorder) {
    int maxWid = 0;
    int maxHgt = 0;

    final HashMap<String, BufferedImage> images = new HashMap<>();

    System.out.print("Reading image files: <");
    for (String imageFile : imageFilesWithPaths) {
      System.out.print("-");
      try {
        BufferedImage img = ImageIO.read(new File(imageFile));
        maxHgt = Math.max(maxHgt, img.getHeight());
        maxWid = Math.max(maxWid, img.getWidth());
        images.put(imageFile, img);
      } catch (IOException e) {
        e.printStackTrace();
      }
    }
    System.out.println(">");

    BufferedImage finalImage;

    int arrSzCols = numImagesPerRow > 0 ? numImagesPerRow
                                        : (int) Math.ceil(Math.sqrt(imageFilesWithPaths.length));
    int arrSzRows = (int) Math.ceil(imageFilesWithPaths.length
                                    / (numImagesPerRow > 0 ? numImagesPerRow
                                                           : Math.ceil(Math.sqrt(imageFilesWithPaths.length))));
    int BUFFER_SZ = 3;
    int bufferCols = (arrSzCols + 1) * BUFFER_SZ;
    int bufferRows = (arrSzRows + 1) * BUFFER_SZ;

    finalImage = new BufferedImage((maxWid * arrSzCols) + bufferCols,
                                   (maxHgt * arrSzRows) + bufferRows, BufferedImage.TYPE_INT_ARGB);
    if (bgColor != null) {
      finalImage.createGraphics().setColor(bgColor);
      finalImage.createGraphics().fillRect(0, 0, finalImage.getWidth(), finalImage.getHeight());
    }
    System.out.print("Drawing image files: <");
    for (int i = 0; i < imageFilesWithPaths.length; i++) {
      System.out.print("-");
      int y = (i / arrSzCols) * maxHgt;
      int x = (i % arrSzCols) * maxWid;
      Graphics2D graphics = finalImage.createGraphics();
      BufferedImage image = images.get(imageFilesWithPaths[i]);
      if (drawInnerBorder) {
        graphics.setColor(Color.GRAY);
        graphics.drawRect(x, y, image.getWidth(), image.getHeight());
      }
      if (drawOuterBorder) {
        graphics.setColor(Color.BLACK);
        graphics.drawRect(x, y, maxWid, maxHgt);
      }
      graphics.drawImage(image, x, y, null);
    }
    return finalImage;
  }

  public static void stitchImages(String dir, String[] imageFiles, String outFile,
                                  int numImagesPerRow, Color bgColor, boolean drawInnerBorder,
                                  boolean drawOuterBorder/* , boolean pack */) {
    String[] imageFilesWithPaths = new String[imageFiles.length];
    for (int i = 0; i < imageFiles.length; i++) {
      imageFilesWithPaths[i] = Files.isRelativePath(imageFiles[i]) ? dir + imageFiles[i]
                                                                   : imageFiles[i];
    }
    BufferedImage finalImage = stitchImages(imageFilesWithPaths, numImagesPerRow, bgColor,
                                            drawInnerBorder, drawOuterBorder);
    System.out.print(">\nWriting final file...");
    File outputFile = new File(Files.isRelativePath(outFile) ? dir + outFile : outFile);
    // TODO check file is writeable
    try {
      ImageIO.write(finalImage, "png", outputFile);
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }

    System.out.println("Complete!");
  }

}
