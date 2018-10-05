package org.pankratzlab.common;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import javax.imageio.ImageIO;

public class Images {

  public static void stitchImages(String dir, String listFile, String outFile, Color bgColor,
                                  boolean drawInnerBorder,
                                  boolean drawOuterBorder/* , boolean pack */) {
    System.out.print("Loading input file...");
    String[] imageFiles = HashVec.loadFileToStringArray(dir + listFile, false, new int[] {0},
                                                        false);
    System.out.println("Complete!");
    stitchImages(dir, imageFiles, outFile, bgColor, drawInnerBorder, drawOuterBorder);
  }

  public static BufferedImage stitchImages(String[] imageFilesWithPaths, Color bgColor,
                                           boolean drawInnerBorder, boolean drawOuterBorder) {
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
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
    }
    System.out.println(">");

    BufferedImage finalImage;

    int arrSzCols = (int) Math.ceil(Math.sqrt(imageFilesWithPaths.length));
    int arrSzRows = (int) Math.ceil(imageFilesWithPaths.length
                                    / Math.ceil(Math.sqrt(imageFilesWithPaths.length)));
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

  public static void stitchImages(String dir, String[] imageFiles, String outFile, Color bgColor,
                                  boolean drawInnerBorder,
                                  boolean drawOuterBorder/* , boolean pack */) {
    String[] imageFilesWithPaths = new String[imageFiles.length];
    for (int i = 0; i < imageFiles.length; i++) {
      imageFilesWithPaths[i] = dir + imageFiles[i];
    }
    BufferedImage finalImage = stitchImages(imageFilesWithPaths, bgColor, drawInnerBorder,
                                            drawOuterBorder);
    System.out.print(">\nWriting final file...");
    File outputFile = new File(dir + outFile);
    // TODO check file is writeable
    try {
      ImageIO.write(finalImage, "png", outputFile);
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }

    System.out.println("Complete!");
  }

  public static void main(String[] args) {
    int numArgs = args.length;
    String file = "imageList.txt";
    String dir = "D:/data/gedi_gwas/data/corrected/";
    String out = "stitched.png";
    Color bgColor = Color.WHITE;
    boolean outerBorder = true;
    boolean innerBorder = true;

    String usage = "\n" + "common.Images requires 3+ arguments\n" + "   (1) directory (i.e. dir="
                   + dir + " (default))\n" + "   (2) filename (i.e. file=" + file + " (default))\n"
                   + "   (3) output filename (i.e. out=" + out + " (default))\n" + "";

    for (String arg : args) {
      if (arg.equals("-h") || arg.equals("-help") || arg.equals("/h") || arg.equals("/help")) {
        System.err.println(usage);
        System.exit(1);
      } else if (arg.startsWith("file=")) {
        file = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("dir=")) {
        dir = arg.split("=")[1];
        numArgs--;
      } else if (arg.startsWith("out=")) {
        out = arg.split("=")[1];
        numArgs--;
      } else {
        System.err.println("Error - invalid argument: " + arg);
      }
    }

    if (numArgs != 0) {
      System.err.println(usage);
      System.exit(1);
    }

    try {
      stitchImages(dir, file, out, bgColor, innerBorder, outerBorder/* , false */);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
