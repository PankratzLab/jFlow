package org.genvisis.one.JL;

import java.util.ArrayList;

import org.genvisis.common.ArrayUtils;
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
		System.out.println(ArrayUtils.toStr(ArrayUtils.toStringArray(b)));
		CmdLine.runCommandWithFileChecks(	ArrayUtils.toStringArray(b), "", null, null, false, true, false,
																			new Logger());

	}

	public static void main(String[] args) {
		bashIT();
	}

}
