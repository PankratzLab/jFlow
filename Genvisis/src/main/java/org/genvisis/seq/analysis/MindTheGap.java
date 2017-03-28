/**
 * 
 */
package org.genvisis.seq.analysis;

import java.util.ArrayList;

import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Logger;

/**
 * @author Kitty
 * 
 *         Thing might work https://github.com/GATB/MindTheGap
 * 
 * 
 *         ./bin/MindTheGap find -in
 *         /Users/Kitty/temp/mica/mindthegap/test.r1.fq,/Users/Kitty/temp/mica/mindthegap/test.r2.fq
 *         -ref /Users/Kitty/temp/mica/mindthegap/micaRef.fasta -out example
 * 
 * 
 *         ./bin/MindTheGap fill -graph example2.h5 -bkpt example2.breakpoints -out example2
 *
 */
public class MindTheGap {

	private boolean convertToFasta(	String inputBam, String r1, String r2, String samToFastQLoc,
																	Logger log) {
		String[] inputs = new String[] {inputBam};
		String[] outputs = new String[] {r1, r2};
		ArrayList<String> command = new ArrayList<String>();
		command.add("java");
		command.add("-jar");
		command.add(samToFastQLoc);
		command.add("I=");
		command.add(inputBam);
		command.add("F=");
		command.add(r1);
		command.add("F2=");
		command.add(r2);

		return CmdLine.runCommandWithFileChecks(Array.toStringArray(command), "", inputs, outputs, true,
																						false, false, log);
	}


}
