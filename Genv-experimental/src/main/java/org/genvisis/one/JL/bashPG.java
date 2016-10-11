package org.genvisis.one.JL;

import java.util.ArrayList;

import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Logger;

public class bashPG {


	private static void bashIT() {
		ArrayList<String> b = new ArrayList<String>();
		b.add("/bin/sh");
		b.add("-c");

		b.add("lt");
		b.add("|");
		b.add("grep drw");
		System.out.println(Array.toStr(Array.toStringArray(b)));
		CmdLine.runCommandWithFileChecks(	Array.toStringArray(b), "", null, null, false, true, false,
																			new Logger());

	}

	public static void main(String[] args) {
		bashIT();
	}

}
