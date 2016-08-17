package org.genvisis.bot;

import java.awt.Color;
import java.awt.MouseInfo;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.Robot;
import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferInt;
import java.io.IOException;
import java.util.Vector;

import org.genvisis.common.Logger;
import org.genvisis.common.Matrix;

public class Ahk {
  public static String colorToHex(Color y) {
    // return "0x"+intToHex(y.getRed())+intToHex(y.getGreen())+intToHex(y.getBlue()); // actual, but
    // doesn't match borrowed AHK code
    String blue = intToHex(y.getBlue());
    String green = intToHex(y.getGreen());
    String red = intToHex(y.getRed());

    blue = (blue.length() == 1 ? "0" : "") + blue;
    green = (green.length() == 1 ? "0" : "") + green;
    red = (red.length() == 1 ? "0" : "") + red;

    return "0x" + blue + green + red;
    // return "0x"+intToHex(y.getBlue())+intToHex(y.getGreen())+intToHex(y.getRed())
  }

  public static int[][] convertToIntPixelMatrix(BufferedImage image) {
    int[][] pixelMatrix;
    int[] pixels;
    int width, height, x, y;


    width = image.getWidth();
    height = image.getHeight();
    pixelMatrix = new int[width][height];
    pixels = convertToIntPixels(image);
    for (int i = 0; i < pixelMatrix.length; i++) {
      x = i % width;
      y = (int) Math.floor(i / width);
      pixelMatrix[x][y] = pixels[i];
    }

    return pixelMatrix;
  }

  public static int[] convertToIntPixels(BufferedImage image) {
    return ((DataBufferInt) image.getRaster().getDataBuffer()).getData();
  }

  public static void doType(Robot bot, int... keyCodes) {
    doType(bot, keyCodes, 0, keyCodes.length);
  }

  public static void doType(Robot bot, int[] keyCodes, int offset, int length) {
    if (length == 0) {
      return;
    }

    bot.keyPress(keyCodes[offset]);
    doType(bot, keyCodes, offset + 1, length - 1);
    bot.keyRelease(keyCodes[offset]);
  }

  public static int[][] findAllMatricesFast(BufferedImage screencapture, int[] matrix, Logger log) {
    int width, dimension; // , height
    int[] pixels;
    Vector<int[]> matches;
    int pX, pY, sX, sY, index, count;

    width = screencapture.getWidth();
    // height = screencapture.getHeight();
    dimension = (int) Math.sqrt(matrix.length);
    pixels = convertToIntPixels(screencapture);

    matches = new Vector<int[]>();
    for (int i = 0; i < pixels.length; i++) {
      sX = i % width;
      sY = (int) Math.floor(i / width);
      count = 0;
      for (int j = 0; j < matrix.length; j++) {
        pX = j % dimension;
        pY = (int) Math.floor(j / dimension);
        index = (sY + pY) * width + (sX + pX);
        if (index >= pixels.length || pixels[index] != matrix[j]) {
          break;
        }
        count++;
      }
      if (count > 0) {
        // log.report("Almost at "+sX+","+sY+" but only "+count);
        // bot.mouseMove(sX, sY);
        // bot.delay(1000);
      }
      if (count == matrix.length) {
        matches.add(new int[] {sX, sY});
      }
    }

    return Matrix.toMatrix(matches);
  }

  public static int[][] findAllMatricesFast(Robot bot, int[] matrix, Logger log) {
    return findAllMatricesFast(bot.createScreenCapture(new Rectangle(Toolkit.getDefaultToolkit()
                                                                            .getScreenSize())),
                               matrix, log);
  }

  public static int[] findMatrix(Robot bot, int ScanStartX, int ScanEndX, int ScanStartY,
                                 int ScanEndY, int Dimension, String Matrix) {
    String[] Pix;
    int Length, count, MatchFound;
    int StartX, StartY, EndX, EndY;
    int[] p;
    int pX, pY, sX, sY;
    String CurrPix;

    Pix = Matrix.split("\\|");
    Length = ScanEndY - ScanStartY;

    // loop for per line of screenY
    for (int A_Index = 1; A_Index <= Length; A_Index++) {
      StartX = ScanStartX;
      EndX = ScanEndX;

      StartY = ScanStartY + A_Index - 1;
      EndY = ScanStartY + A_Index - 1;

      // loop for searching inside a line
      while (true) {
        count = 0;
        MatchFound = 1;

        // change variance and 'Fast' if not suited
        // PixelSearch, pX, pY, %StartX%, %StartY%, %EndX%, %EndY%, %Pix1%, 0, Fast
        p = pixelSearch(bot, StartX, StartY, EndX, EndY, Pix[0]);
        if (p == null) {
          break;
        }
        pX = p[0];
        pY = p[1];
        sX = pX;
        sY = pY;
        // bot.mouseMove(sX, sY);

        // loop for matching
        for (int i = 1; i <= Dimension; i++) {
          for (int j = 1; j <= Dimension; j++) {
            CurrPix = colorToHex(bot.getPixelColor(pX, pY));
            if (!Pix[count].equals(CurrPix)) {
              MatchFound = 0;
            }
            if (MatchFound == 0) {
              break;
            }

            pX++;
            count++;
          }
          if (MatchFound == 0) {
            break;
          }
          pY++;
          pX -= Dimension;
        }

        // match found!
        if (MatchFound == 1) {
          return new int[] {sX, sY};
        } else {
          break;
        }

        // // unreachable code
        // StartX++;
        //
        // DiffX = EndX - StartX;
        // if (DiffX < Dimension) {
        // break;
        // }
      }
      StartY++;
      EndY++;
    }

    return null;
  }

  public static String getClipboard(boolean verbose) {
    Clipboard systemClipboard;
    Transferable contents;

    try {
      systemClipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
      contents = systemClipboard.getContents(null);
      if (contents == null) {
        return ("Clipboard is empty");
      } else {
        try {
          if (contents.isDataFlavorSupported(DataFlavor.stringFlavor)) {
            return (String) contents.getTransferData(DataFlavor.stringFlavor);
          }
        } catch (UnsupportedFlavorException ufe) {
          ufe.printStackTrace();
        } catch (IOException ioe) {
          ioe.printStackTrace();
        }
      }
    } catch (Exception e) {
      if (verbose) {
        e.printStackTrace();
      }
      return ("Clipboard failure");
    }

    return null;
  }

  public static String GetMatrix(Robot bot, int StartX, int StartY, int Dimension) {
    String matrix;

    matrix = "";
    for (int j = 0; j < Dimension; j++) {
      for (int i = 0; i < Dimension; i++) {
        matrix +=
            (i == 0 && j == 0 ? "" : "|") + colorToHex(bot.getPixelColor(StartX + i, StartY + j));
      }
    }

    return matrix;
  }

  public static int[] getMatrixInt(Robot bot, int StartX, int StartY, int dimension) {
    int[] pixels;
    int count;

    count = 0;
    pixels = new int[dimension * dimension];
    for (int j = 0; j < dimension; j++) {
      for (int i = 0; i < dimension; i++) {
        pixels[count++] = bot.getPixelColor(StartX + i, StartY + j).getRGB();
      }
    }

    return pixels;
  }

  // public static int[][] findAllMatricesFast(Robot bot, int[] matrix) {
  // BufferedImage screencapture;
  // int width, dimension; // , height
  // int[][] pixelMatrix;
  // Vector<int[]> matches;
  // int pX, pY, sX, sY, index, count;
  //
  // screencapture = bot.createScreenCapture(new
  // Rectangle(Toolkit.getDefaultToolkit().getScreenSize()));
  // width = screencapture.getWidth();
  //// height = screencapture.getHeight();
  // dimension = (int)Math.sqrt(matrix.length);
  // pixels = convertToIntPixels(screencapture);
  //
  // matches = new Vector<int[]>();
  // for (int i = 0; i < pixelMatrix.length-dimension+1; i++) {
  // for (int j = 0; j < pixelMatrix[i].length-dimension+1; j++) {
  // for (int k = 0; k < matrix.length; k++) {
  // pX = k % dimension;
  // pY = (int)Math.floor(k/dimension);
  // index = (sY+pY)*width+(sX+pX);
  // if (index >= pixels.length || pixels[index] != matrix[k]) {
  // break;
  // }
  // count++;
  // }
  // if (count == matrix.length) {
  // matches.add(new int[] {i % width, });
  // }
  //
  // }
  // }
  //
  // return Matrix.toMatrix(matches);
  // }

  public static String intToHex(int y) {
    String numeric = "0123456789ABCDEF";
    String result;
    int r = y % 16;
    if ((y - r) == 0) {
      result = numeric.substring(r, r + 1);
    } else {
      result = intToHex((y - r) / 16) + numeric.substring(r, r + 1);
    }
    return result;
  }

  public static void mouseClick(Robot bot, int x, int y) {
    mouseClick(bot, x, y, false);
  }

  public static void mouseClick(Robot bot, int x, int y, boolean rightClick) {
    bot.mouseMove(x, y);
    if (rightClick) {
      bot.mousePress(InputEvent.BUTTON3_MASK);
      bot.mouseRelease(InputEvent.BUTTON3_MASK);
    } else {
      bot.mousePress(InputEvent.BUTTON1_MASK);
      bot.mouseRelease(InputEvent.BUTTON1_MASK);
    }
  }

  public static void mouseScrollUp(Robot bot) {
    bot.mouseWheel(-100);
  }

  public static void moveMouseSlowly(Robot bot, int x, int y, int speed) {
    try {
      final Point mouseLoc = MouseInfo.getPointerInfo().getLocation();
      final int currentX = mouseLoc.x;
      final int currentY = mouseLoc.y;
      for (int i = 0; i < 100; i++) {
        int newX = ((x * i) / 100) + (currentX * (100 - i) / 100);
        int newY = ((y * i) / 100) + (currentY * (100 - i) / 100);
        bot.mouseMove(newX, newY);
        bot.delay(speed);
      }
      bot.mouseMove(x, y);
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void outline(Robot bot, int ScanStartX, int ScanEndX, int ScanStartY,
                             int ScanEndY) {
    moveMouseSlowly(bot, ScanStartX, ScanStartY, 25);
    moveMouseSlowly(bot, ScanEndX, ScanEndY, 25);
  }

  public static int[] pixelSearch(Robot bot, int StartX, int StartY, int EndX, int EndY,
                                  String pix) {
    // Logger log;
    //
    // log = new Logger("/home/npankrat/auto/search.log", true);
    // log.report("Searching for '"+pix+"' between "+StartX+","+StartY+" and "+EndX+","+EndY);
    for (int x = StartX; x <= EndX; x++) {
      for (int y = StartY; y <= EndY; y++) {
        // log.report("Found '"+colorToHex(bot.getPixelColor(x, y))+"' at "+x+","+y);
        if (colorToHex(bot.getPixelColor(x, y)).equals(pix)) {
          return new int[] {x, y};
        }
      }
    }
    return null;
  }

  public static void send(Robot bot, String text) {
    for (int i = 0; i < text.length(); i++) {
      type(bot, text.charAt(i));
      sleep(100);
    }
  }

  public static void setClipboard(String text) {
    Toolkit.getDefaultToolkit().getSystemClipboard().setContents(new StringSelection(text), null);
  }

  public static void sleep(int ms) {
    try {
      Thread.sleep(ms);
    } catch (Exception e) {
    }
  }

  public static void type(Robot bot, char character) {
    switch (character) {
      case 'a':
        doType(bot, KeyEvent.VK_A);
        break;
      case 'b':
        doType(bot, KeyEvent.VK_B);
        break;
      case 'c':
        doType(bot, KeyEvent.VK_C);
        break;
      case 'd':
        doType(bot, KeyEvent.VK_D);
        break;
      case 'e':
        doType(bot, KeyEvent.VK_E);
        break;
      case 'f':
        doType(bot, KeyEvent.VK_F);
        break;
      case 'g':
        doType(bot, KeyEvent.VK_G);
        break;
      case 'h':
        doType(bot, KeyEvent.VK_H);
        break;
      case 'i':
        doType(bot, KeyEvent.VK_I);
        break;
      case 'j':
        doType(bot, KeyEvent.VK_J);
        break;
      case 'k':
        doType(bot, KeyEvent.VK_K);
        break;
      case 'l':
        doType(bot, KeyEvent.VK_L);
        break;
      case 'm':
        doType(bot, KeyEvent.VK_M);
        break;
      case 'n':
        doType(bot, KeyEvent.VK_N);
        break;
      case 'o':
        doType(bot, KeyEvent.VK_O);
        break;
      case 'p':
        doType(bot, KeyEvent.VK_P);
        break;
      case 'q':
        doType(bot, KeyEvent.VK_Q);
        break;
      case 'r':
        doType(bot, KeyEvent.VK_R);
        break;
      case 's':
        doType(bot, KeyEvent.VK_S);
        break;
      case 't':
        doType(bot, KeyEvent.VK_T);
        break;
      case 'u':
        doType(bot, KeyEvent.VK_U);
        break;
      case 'v':
        doType(bot, KeyEvent.VK_V);
        break;
      case 'w':
        doType(bot, KeyEvent.VK_W);
        break;
      case 'x':
        doType(bot, KeyEvent.VK_X);
        break;
      case 'y':
        doType(bot, KeyEvent.VK_Y);
        break;
      case 'z':
        doType(bot, KeyEvent.VK_Z);
        break;
      case 'A':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_A);
        break;
      case 'B':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_B);
        break;
      case 'C':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_C);
        break;
      case 'D':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_D);
        break;
      case 'E':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_E);
        break;
      case 'F':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_F);
        break;
      case 'G':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_G);
        break;
      case 'H':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_H);
        break;
      case 'I':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_I);
        break;
      case 'J':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_J);
        break;
      case 'K':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_K);
        break;
      case 'L':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_L);
        break;
      case 'M':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_M);
        break;
      case 'N':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_N);
        break;
      case 'O':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_O);
        break;
      case 'P':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_P);
        break;
      case 'Q':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_Q);
        break;
      case 'R':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_R);
        break;
      case 'S':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_S);
        break;
      case 'T':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_T);
        break;
      case 'U':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_U);
        break;
      case 'V':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_V);
        break;
      case 'W':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_W);
        break;
      case 'X':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_X);
        break;
      case 'Y':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_Y);
        break;
      case 'Z':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_Z);
        break;
      case '`':
        doType(bot, KeyEvent.VK_BACK_QUOTE);
        break;
      case '0':
        doType(bot, KeyEvent.VK_0);
        break;
      case '1':
        doType(bot, KeyEvent.VK_1);
        break;
      case '2':
        doType(bot, KeyEvent.VK_2);
        break;
      case '3':
        doType(bot, KeyEvent.VK_3);
        break;
      case '4':
        doType(bot, KeyEvent.VK_4);
        break;
      case '5':
        doType(bot, KeyEvent.VK_5);
        break;
      case '6':
        doType(bot, KeyEvent.VK_6);
        break;
      case '7':
        doType(bot, KeyEvent.VK_7);
        break;
      case '8':
        doType(bot, KeyEvent.VK_8);
        break;
      case '9':
        doType(bot, KeyEvent.VK_9);
        break;
      case '-':
        doType(bot, KeyEvent.VK_MINUS);
        break;
      case '=':
        doType(bot, KeyEvent.VK_EQUALS);
        break;
      case '~':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_BACK_QUOTE);
        break;
      case '!':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_1);
        break;
      case '@':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_2);
        break;
      case '#':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_3);
        break;
      case '$':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_4);
        break;
      case '%':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_5);
        break;
      case '^':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_6);
        break;
      case '&':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_7);
        break;
      case '*':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_8);
        break;
      case '(':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_9);
        break;
      case ')':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_0);
        break;
      case '_':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_MINUS);
        break;
      case '+':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_EQUALS);
        break;
      case '\t':
        doType(bot, KeyEvent.VK_TAB);
        break;
      case '\n':
        doType(bot, KeyEvent.VK_ENTER);
        break;
      case '[':
        doType(bot, KeyEvent.VK_OPEN_BRACKET);
        break;
      case ']':
        doType(bot, KeyEvent.VK_CLOSE_BRACKET);
        break;
      case '\\':
        doType(bot, KeyEvent.VK_BACK_SLASH);
        break;
      case '{':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_OPEN_BRACKET);
        break;
      case '}':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_CLOSE_BRACKET);
        break;
      case '|':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_BACK_SLASH);
        break;
      case ';':
        doType(bot, KeyEvent.VK_SEMICOLON);
        break;
      case ':':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_SEMICOLON);
        break;
      case '\'':
        doType(bot, KeyEvent.VK_QUOTE);
        break;
      case '"':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_QUOTE);
        break;
      case ',':
        doType(bot, KeyEvent.VK_COMMA);
        break;
      case '<':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_COMMA);
        break;
      case '.':
        doType(bot, KeyEvent.VK_PERIOD);
        break;
      case '>':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_PERIOD);
        break;
      case '/':
        doType(bot, KeyEvent.VK_SLASH);
        break;
      case '?':
        doType(bot, KeyEvent.VK_SHIFT, KeyEvent.VK_SLASH);
        break;
      case ' ':
        doType(bot, KeyEvent.VK_SPACE);
        break;
      default:
        throw new IllegalArgumentException("Cannot type character " + character);
    }
  }
}
