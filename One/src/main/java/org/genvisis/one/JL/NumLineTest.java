package org.genvisis.one.JL;

import java.io.File;

import org.genvisis.common.Files;


/**
 * @author lane0212
 *Testing line counter
 */
public class NumLineTest {

	public static void main(String[] args) {
		String[] file = new String[] { "1",
				"2",
				"3",
				"4",
				"5"
		};
		Files.writeList(file, "test.txt");
		System.out.println(Files.countLines("test.txt", 0));
		System.out.println(Files.countLines("test.txt", 1));
		System.out.println(Files.countLines("test.txt", 2));

		new File("test.txt").delete();

	}
}
