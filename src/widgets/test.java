package widgets;

import java.io.IOException;

public class test {
	public static void main(String[] args) throws IOException {

		for (int i = 0; i<args.length; i++) {
			System.out.print("\""+args[i]+"\"\t");
		}
		System.out.println();

	}

}
