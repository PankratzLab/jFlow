package org.genvisis.common;

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
    String[] imageFiles =
                        HashVec.loadFileToStringArray(dir + listFile, false, new int[] {0}, false);
    System.out.println("Complete!");
    int maxWid = 0;
    int maxHgt = 0;

    final HashMap<String, BufferedImage> images = new HashMap<String, BufferedImage>();

    System.out.print("Reading image files: <");
    for (String imageFile : imageFiles) {
      System.out.print("-");
      try {
        BufferedImage img = ImageIO.read(new File(dir + imageFile));
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

    // if (pack) {
    // finalImage = packedImage(imageFiles, images, bgColor, drawBorder);
    // } else {
    int arrSzCols = (int) Math.ceil(Math.sqrt(imageFiles.length));
    int arrSzRows = (int) Math.ceil(imageFiles.length / Math.ceil(Math.sqrt(imageFiles.length)));
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
    for (int i = 0; i < imageFiles.length; i++) {
      System.out.print("-");
      int y = (i / arrSzRows) * maxHgt;
      int x = (i % arrSzCols) * maxWid;
      Graphics2D graphics = finalImage.createGraphics();
      BufferedImage image = images.get(imageFiles[i]);
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
    // }

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
  //
  // private static BufferedImage packedImage(final String[] imageFiles, final HashMap<String,
  // BufferedImage> images, Color bgColor, boolean[] drawBorder) {
  // String[] sorted = Array.clone(imageFiles);
  // String[] sortedW = Array.clone(imageFiles);
  // String[] sortedH = Array.clone(imageFiles);
  // Arrays.sort(sorted, new Comparator<String>() {
  // @Override
  // public int compare(String o1, String o2) {
  // // TODO include boilerplate for comparison?
  // BufferedImage i1 = images.get(o1);
  // BufferedImage i2 = images.get(o2);
  // return new Integer(i1.getHeight() + i1.getWidth()).compareTo(i2.getHeight() + i2.getWidth());
  // }
  // });
  // Arrays.sort(sortedW, new Comparator<String>() {
  // @Override
  // public int compare(String o1, String o2) {
  // // TODO include boilerplate for comparison?
  // BufferedImage i1 = images.get(o1);
  // BufferedImage i2 = images.get(o2);
  // return new Integer(i1.getWidth()).compareTo(i2.getWidth());
  // }
  // });
  // Arrays.sort(sortedH, new Comparator<String>() {
  // @Override
  // public int compare(String o1, String o2) {
  // // TODO include boilerplate for comparison?
  // BufferedImage i1 = images.get(o1);
  // BufferedImage i2 = images.get(o2);
  // return new Integer(i1.getHeight()).compareTo(i2.getHeight());
  // }
  // });
  //
  //
  // int BUFFER_SZ = 3;
  //
  // int cols = (int) Math.ceil(Math.sqrt(imageFiles.length));
  // int bufferX = (cols + 1) * BUFFER_SZ;
  // int rows = (int) Math.ceil(imageFiles.length / Math.ceil(Math.sqrt(imageFiles.length)));
  // int bufferY = (rows + 1) * BUFFER_SZ;
  // int maxWdt = images.get(sortedW[0]).getWidth();
  // int maxHgt = images.get(sortedH[0]).getHeight();
  //
  // BufferedImage finalImage = new BufferedImage((maxWdt * cols) + bufferX, (maxHgt * rows) +
  // bufferY, BufferedImage.TYPE_INT_ARGB);
  //
  //
  // if (bgColor != null) {
  // finalImage.createGraphics().setColor(bgColor);
  // finalImage.createGraphics().fillRect(0, 0, finalImage.getWidth(), finalImage.getHeight());
  // }
  // System.out.print("Drawing image files <");
  // for (int i = 0; i < imageFiles.length; i++) {
  // System.out.print("-");
  // int x = (i / rows) * maxWdt;
  // int y = (i % cols) * maxHgt;
  // Graphics2D graphics = finalImage.createGraphics();
  // BufferedImage image = images.get(imageFiles[i]);
  // if (drawBorder != null && drawBorder.length > 1 && drawBorder[1]) {
  // graphics.setColor(Color.GRAY);
  // graphics.drawRect(x, y, image.getWidth(), image.getHeight());
  // }
  // if (drawBorder != null && drawBorder.length > 0 && drawBorder[0]) {
  // graphics.setColor(Color.BLACK);
  // graphics.drawRect(x, y, maxWdt, maxHgt);
  // }
  // graphics.drawImage(image, x, y, null);
  // }
  //
  // return finalImage;
  // }

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
