package cnv.plots;

import java.awt.Color;
import java.util.Random;

public class ColorExt {
	public static Color generateRandomColor(Color mix) {
		Random random = new Random();
		int red = random.nextInt(256);
		int green = random.nextInt(256);
		int blue = random.nextInt(256);

		// mix the color
		if (mix != null) {
			red = (red + mix.getRed()) / 2;
			green = (green + mix.getGreen()) / 2;
			blue = (blue + mix.getBlue()) / 2;
		}

		Color color = new Color(red, green, blue);
		return color;
	}

	public static Color[] generatePastelShades() {
		Color one = generateRandomColor(new Color(255, 255, 255));
		Color two = generateRandomColor(one);
		return new Color[] { one, two };
	}
}
