package one.JL;

import java.util.ArrayList;

import common.Array;
import common.CmdLine;
import common.Logger;

public class bashPG {
	
	
	private static void bashIT() {
		ArrayList<String> b = new ArrayList<String>();
		b.add("/bin/sh");
		b.add("-c");

		b.add("lt");
		b.add("|");
		b.add("grep drw");
		System.out.println(Array.toStr(Array.toStringArray(b)));
		CmdLine.runCommandWithFileChecks(Array.toStringArray(b), "", null, null, false, true, false, new Logger());

	}
	public static void main(String[] args) {
		bashIT();
	}

}
