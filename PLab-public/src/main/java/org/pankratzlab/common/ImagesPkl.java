package org.pankratzlab.common;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;

import javax.imageio.ImageIO;
import javax.imageio.stream.ImageInputStream;

import org.apache.commons.compress.archivers.ArchiveException;
import org.apache.commons.compress.archivers.ArchiveStreamFactory;
import org.apache.commons.compress.archivers.tar.TarArchiveEntry;
import org.apache.commons.compress.archivers.tar.TarArchiveInputStream;
import org.apache.commons.compress.compressors.gzip.GzipCompressorInputStream;

public class ImagesPkl {

  public static ImageInputStream getImageInputStream(String file) throws IOException {
    String[] finds = {".tar", ".tar.gz", ".tgz"};
    String tarFile = null;
    for (String f : finds) {
      int index;
      if ((index = file.indexOf(f + File.separatorChar)) > 0) {
        tarFile = file.substring(0, index + f.length());
        break;
      }
    }
    if (tarFile == null) return ImageIO.createImageInputStream(new File(file));
    return getTarImageInputStream(file);
  }

  public static InputStream getTarInputStream(String file) {
    InputStream is;
    try {
      is = new BufferedInputStream(new FileInputStream(file));
    } catch (FileNotFoundException e1) {
      e1.printStackTrace();
      return null;
    }
    if (file.endsWith("gz")) {
      // either .tar.gz or .tgz
      try {
        is = new GzipCompressorInputStream(is);
      } catch (IOException e) {
        e.printStackTrace();
        return null;
      }
    }
    TarArchiveInputStream tarStream;
    try {
      tarStream = (TarArchiveInputStream) new ArchiveStreamFactory().createArchiveInputStream("tar",
                                                                                              is);
      return tarStream;
    } catch (ArchiveException e) {
      e.printStackTrace();
      return null;
    }
  }

  private static ImageInputStream getTarImageInputStream(String file) throws IOException {
    String[] finds = {".tar", ".tar.gz", ".tgz"};
    String tarFile = null;
    String subFile = null;
    for (String f : finds) {
      int index;
      if ((index = file.indexOf(f + File.separatorChar)) > 0) {
        tarFile = file.substring(0, index + f.length());
        subFile = file.substring(index + f.length() + 1);
        break;
      }
    }
    if (tarFile == null) return null;

    TarArchiveInputStream tarStream = (TarArchiveInputStream) getTarInputStream(tarFile);
    if (tarStream == null) return null;
    TarArchiveEntry entry;
    ImageInputStream iis = null;
    while ((entry = tarStream.getNextTarEntry()) != null) {
      if (new File(entry.getName()).equals(new File(subFile))) {
        byte[] buffer = new byte[1024];
        ByteArrayOutputStream os = new ByteArrayOutputStream();
        int bytesRead;
        while ((bytesRead = tarStream.read(buffer, 0, 1024)) > -1) {
          os.write(buffer, 0, bytesRead);
        }
        iis = ImageIO.createImageInputStream(new ByteArrayInputStream(os.toByteArray()));
        os.close();
        break;
      }
    }
    tarStream.close();
    return iis;
  }

  public static void stitchImages(String dir, String listFile, String outFile, int numImagesPerRow,
                                  Color bgColor, boolean drawInnerBorder,
                                  boolean drawOuterBorder/* , boolean pack */) throws IOException {
    System.out.print("Loading input file...");
    String list = Files.isRelativePath(listFile) ? dir + listFile : listFile;
    String[] imageFiles = HashVec.loadFileToStringArray(list, false, new int[] {0}, false);
    System.out.println("Complete!");
    stitchImages(dir, imageFiles, outFile, numImagesPerRow, bgColor, drawInnerBorder,
                 drawOuterBorder);
  }

  public static BufferedImage stitchImages(String[] imageFilesWithPaths, int numImagesPerRow,
                                           Color bgColor, boolean drawInnerBorder,
                                           boolean drawOuterBorder) throws IOException {
    ImageInputStream[] arr = new ImageInputStream[imageFilesWithPaths.length];
    for (int i = 0; i < arr.length; i++) {
      arr[i] = ImageIO.createImageInputStream(new File(imageFilesWithPaths[i]));
    }
    return stitchImages(arr, numImagesPerRow, bgColor, drawInnerBorder, drawOuterBorder);
  }

  public static BufferedImage stitchImages(ImageInputStream[] imageFilesWithPaths,
                                           int numImagesPerRow, Color bgColor,
                                           boolean drawInnerBorder, boolean drawOuterBorder) {
    int maxWid = 0;
    int maxHgt = 0;

    final HashMap<ImageInputStream, BufferedImage> images = new HashMap<>();

    System.out.print("Reading image files: <");
    for (ImageInputStream imageStream : imageFilesWithPaths) {
      System.out.print("-");
      try {
        BufferedImage img = ImageIO.read(imageStream);
        maxHgt = Math.max(maxHgt, img.getHeight());
        maxWid = Math.max(maxWid, img.getWidth());
        images.put(imageStream, img);
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
                                  boolean drawOuterBorder/* , boolean pack */) throws IOException {
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
