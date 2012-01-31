package one;

public class TryThread extends Thread {
	public String str;

	public boolean done;

	public TryThread(String s) {
		str = s;
		setDaemon(true);
	}

	public void run() {
		try {
			Runtime.getRuntime().exec(str).waitFor();
		} catch (Exception e) {
			e.printStackTrace();
		}
		done = true;
	}

	public boolean isDone() {
		return done;
	}
}
