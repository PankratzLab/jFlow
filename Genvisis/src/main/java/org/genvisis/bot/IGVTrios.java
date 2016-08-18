package org.genvisis.bot;

import java.awt.event.KeyEvent;

import org.genvisis.common.Logger;

public class IGVTrios {
  private final AhkBot bot;
  private boolean first;

  public IGVTrios(AhkBot bot, Logger log) {
    this.bot = bot;
  }

  public void next() {
    int x, y;

    // x = bot.getX();
    // y = bot.getY();
    //
    // System.out.println(x+", "+y);
    //
    // bot.mouseClick(x, y+30, false, true);
    // bot.delay(200);
    // bot.mouseClick(x, y+30, true, false);

    // bot.mouseClick(60, 230, false, true);
    // bot.delay(200);

    // bot.mouseClick(60, 260, false, true);

    x = 140;
    if (first) {

      for (y = 200; y < 800; y++) {
        if (bot.getPixelColor(x, y).getRGB() == -2960686) {
          break;
        }

        bot.mouseMove(x, y);
        // System.out.println(bot.getPixelColor(140, y).getRGB());
      }

      bot.mouseClick(x, y + 30, false, true);
      bot.delay(200);
      bot.mouseClick(x, y + 30, true, false);
      bot.delay(500);
      bot.keyPress(KeyEvent.VK_UP);
      bot.delay(100);
      bot.keyPress(KeyEvent.VK_UP);
      bot.delay(100);
      bot.keyPress(KeyEvent.VK_UP);
      bot.delay(100);
      bot.keyPress(KeyEvent.VK_ENTER);
    } else {
      bot.mouseClick(800, 900, true, false);
      bot.delay(500);
      bot.keyPress(KeyEvent.VK_DOWN);
      bot.delay(100);
      bot.keyPress(KeyEvent.VK_DOWN);
      bot.delay(100);
      bot.keyPress(KeyEvent.VK_DOWN);
      bot.delay(100);
      bot.keyPress(KeyEvent.VK_RIGHT);
      bot.delay(100);
      bot.keyPress(KeyEvent.VK_DOWN);
      bot.delay(100);
      bot.keyPress(KeyEvent.VK_DOWN);
      bot.delay(100);
      bot.keyPress(KeyEvent.VK_DOWN);
      bot.delay(100);
      bot.keyPress(KeyEvent.VK_ENTER);
    }

    first = !first;
  }
}
