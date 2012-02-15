package common;

/**
 * Signals that an exception of some sort has occurred.
 */
public class Elision extends Exception {
	static final long serialVersionUID = -1;

	/**
	 * Constructs an <code>Elisiono</code> with <code>null</code> as its
	 * error detail message.
	 */
	public Elision() {
		super();
	}

	/**
	 * Constructs an <code>Elision</code> with the specified detail message.
	 * The error message string <code>str</code> can later be retrieved by the
	 * <code>{@link java.lang.Throwable#getMessage}</code> method of class
	 * <code>java.lang.Throwable</code>.
	 * 
	 * @param str
	 *            the detail message.
	 */
	public Elision(String str) {
		super(str);
	}

}
