package org.genvisis.one.JL.ssh;

import java.util.ArrayList;
import java.util.List;

import org.genvisis.common.Array;
import org.genvisis.common.CmdLine;
import org.genvisis.common.Logger;

public class SSH {

	private SSH() {

	}

	/**
	 * JLs remote
	 * 
	 * @param local
	 * @param remoteDir
	 * @param log
	 */
	public static void copyLocalToRemote(String local, String remote, Logger log) {
		copyLocalToRemote(local, "msi", remote, log);
	}

	/**
	 * JL run remote command
	 * 
	 * @param commandToRun
	 * @param log
	 */
	public static void runRemoteCommand(String commandToRun, Logger log) {
		runRemoteCommand("msi", commandToRun, log);
	}

	/**
	 * JLs remote
	 * 
	 * @param local
	 * @param remoteDir
	 * @param log
	 */
	public static void copyRemoteToLocal(String local, String remoteDir, Logger log) {
		copyRemoteToLocal(local, "msi", remoteDir, log);
	}

	private static void copyLocalToRemote(String local, String login, String remote, Logger log) {
		ArrayList<String> command = new ArrayList<String>();
		command.add("scp");
		command.add("-r");
		command.add(local);
		command.add(login + ":" + remote);
		CmdLine.runCommandWithFileChecks(command, "", null, null, true, true, false, false, log);
	}

	private static void copyRemoteToLocal(String local, String login, String remote, Logger log) {
		ArrayList<String> command = new ArrayList<String>();
		command.add("scp");
		command.add("-r");
		command.add(login + ":" + remote);
		command.add(local);
		CmdLine.runCommandWithFileChecks(command, "", null, null, true, true, false, false, log);
	}

	private static void runRemoteCommand(String login, String commandToRun, Logger log) {
		ArrayList<String> command = new ArrayList<String>();
		command.add("ssh");
		command.add(login);
		command.add("'" + commandToRun + "'");
		CmdLine.runCommandWithFileChecks(command, "", null, null, true, true, false, false, log);
	}

}
