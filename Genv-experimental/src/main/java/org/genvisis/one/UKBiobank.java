package org.genvisis.one;

import java.io.*;

import org.genvisis.common.*;

public class UKBiobank {
	
	public static void checkChecksums(Logger log) {
		String temp, trav, md5;
		String[] md5s;
		
		md5s = Files.list("./", ".md5");
		log.report("Checking "+md5s.length+" files with the extension .md5");
		
		for (int i = 0; i < md5s.length; i++) {
			try {
				md5 = Files.getHeaderOfFile(md5s[i], log)[0];
				trav = md5s[i].substring(0, md5s[i].length()-4);
				log.report(trav, false, true);
				if (Files.exists(trav)) {
					CmdLine.run("ukb_md5 "+trav, "./", new PrintStream(new File(trav+".log")), null, null, false);
					temp = Files.tail(trav+".log", 1).split("[\\s]+")[1].substring(4);
					log.report("\t"+(md5.equals(temp)?"checks out":"FAILS!!"));
				} else {
					log.report("\tDOES NOT EXIST!!");
				}
			} catch (Exception e) {
				log.reportError("Error with "+md5s[i]);
				log.reportException(e);
			}
		}
	}

	public static void main(String[] args) {
		int numArgs = args.length;
		boolean checkChecksums = false;
		String logfile = null;
		Logger log;

		String usage = "\n" + 
				"org.genvisis.one.UKBiobank requires 0-1 arguments\n" + 
				"   (1) check MD5 checksums (i.e. -checkChecksums (not the default))\n" + "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-h") || args[i].equals("-help") || args[i].equals("/h") || args[i].equals("/help")) {
				System.err.println(usage);
				System.exit(1);
			} else if (args[i].startsWith("-checkChecksums")) {
				checkChecksums = true;
				numArgs--;
			} else if (args[i].startsWith("log=")) {
				logfile = args[i].split("=")[1];
				numArgs--;
			} else {
				System.err.println("Error - invalid argument: " + args[i]);
			}
		}
		if (numArgs != 0) {
			System.err.println(usage);
			System.exit(1);
		}
		try {
			log = new Logger(logfile);
			if (checkChecksums) {
				checkChecksums(log);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
