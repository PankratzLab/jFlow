package org.genvisis;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Option.Builder;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PatternOptionBuilder;
import org.genvisis.common.Logger;

/**
 * Helper class for parsing command line arguments. Typical usage is:
 * <ol>
 * <li>Create a {@link CLI} object with {@link #CLI(String)}. This ensures a standard set of help
 * flags.</li>
 * <li>Build up the parser options by using {@link #addFlag} (e.g. <code>-doThisOption</code>) and
 * {@link #addArg} (e.g., <code>myOpt=someValue</code>)</li>
 * <li>(Optional) If your program includes mutually exclusive commands and you want to ensure a user
 * selects one (and only one), group these options with {@link #addGroup}</li>
 * <li>Parse your the command line with <b>either</b> {@link #parseWithExit}, if you want your
 * application to exit if the command-line usage is printed, or {@link #parse} if you want to handle
 * any exception cases yourself.</li>
 * <li>After parsing, you can interrogate the {@link CLI} with {@link #get(String)} or
 * {@link #has(String)} methods.</li>
 * </ol>
 */
public class CLI {
	/**
	 * Helper enum for supported parameter types. When using Commons-cli, must be a subset of
	 * {@link PatternOptionBuilder} class constants.
	 */
	public enum Arg {
										STRING, NUMBER, FILE
	}

	private static final String HELP_CHAR = "h"; // short string to print usage
	private static final String HELP_STRING = "help"; // long string to print usage
	private static final String ARG_VALUE = "value"; // default value to print for arguments
	private static final boolean REQUIRED = false; // default required state
	private static final Arg TYPE = Arg.STRING; // default validation class
	private static final int OUT_WIDTH = 150; // help output width
	private static final String INDENT = "   "; // indentation for non-required, non-grouped arguments
	private static final String ARG_PREFIX = "";
	private static final char ARG_SEPARATOR = '=';
	private static final String FLAG_PREFIX = "-";

	private final Options options;
	private final String commandName;
	private final Map<String, String> defaults = new HashMap<String, String>();

	private Map<String, String> parsed;

	/**
	 * @see #CLI(String)
	 * @param program Program class for command invocation
	 */
	public CLI(Class<?> program) {
		this(program.getName());
	}

	/**
	 * Create a {@link CLI} object with a default set of command line options (e.g. "-h, -help" to
	 * print usage).
	 * 
	 * @param program Program class for command invocation
	 */
	public CLI(String program) {
		this(program, new Options());
		addFlag(HELP_CHAR, "Print application usage", HELP_STRING);
	}

	/**
	 * @see #CLI(String, Options)
	 * @param program Program class for command invocation
	 * @param o Baseline options to parse in this CLI instance.
	 */
	public CLI(Class<?> program, Options o) {
		this(program.getName(), o);
	}

	/**
	 * Create a {@link CLI} object with the specified command line {@link Option}s as a starting
	 * point.
	 *
	 * @param program Program name for command invocation
	 * @param o Baseline options to parse in this CLI instance.
	 */
	public CLI(String program, Options o) {
		commandName = program;
		options = o;
	}

	/**
	 * @param key Parsed keyword to look up
	 * @return {@code true} if the given key was parsed
	 * @throws IllegalStateException if a {@link #parse} method has not been called yet.
	 */
	public boolean has(String key) {
		checkParse();
		return parsed.containsKey(key);
	}

	/**
	 * @param key Parsed keyword to look up
	 * @return The parsed value for the given key
	 * @throws IllegalStateException if a {@link #parse} method has not been called yet.
	 */
	public String get(String key) {
		checkParse();
		return parsed.get(key);
	}

	/**
	 * @param key Parsed keyword to look up
	 * @return The parsed value for the given key, converted to a {@link File}
	 * @throws IllegalStateException if a {@link #parse} method has not been called yet.
	 */
	public File getF(String key) {
		return new File(get(key));
	}

	/**
	 * @param key Parsed keyword to look up
	 * @return The parsed value for the given key, converted to an {@code int}
	 * @throws IllegalStateException if a {@link #parse} method has not been called yet.
	 */
	public int getI(String key) {
		return Integer.parseInt(get(key));
	}

	/**
	 * @param key Parsed keyword to look up
	 * @return The parsed value for the given key, converted to a {@code double}
	 * @throws IllegalStateException if a {@link #parse} method has not been called yet.
	 */
	public double getD(String key) {
		return Double.parseDouble(get(key));
	}

	/**
	 * Helper method to ensure it is safe to access the parser output.
	 */
	private void checkParse() {
		if (parsed == null) {
			throw new IllegalStateException("Please call parse first");
		}
	}

	/**
	 * Helper method to finish building an {@link Option}, whether from {@link #addFlag} or
	 * {@link #addArg}.
	 */
	private void add(Builder builder, String longName, String desc, boolean required, Arg type) {
		if (longName != null) {
			builder = builder.longOpt(longName);
		}

		// Mark the description so required flags will be clearly indicated in usage text
		String description = (required ? "(required) " : INDENT) + desc;
		Class<?> argType;
		switch (type) {
			case NUMBER:
				argType = PatternOptionBuilder.NUMBER_VALUE;
				break;
			case FILE:
				argType = PatternOptionBuilder.FILE_VALUE;
				break;
			case STRING:
			default:
				argType = PatternOptionBuilder.STRING_VALUE;
		}
		builder = builder.type(argType);
		builder = builder.desc(description).required(required);

		options.addOption(builder.build());
	}

	/**
	 * @see #addArg(String, String, String, boolean, Arg)
	 *
	 * @param name Argument name
	 * @param description Text to use when printing command-line usage
	 */
	public void addArg(String name, String description) {
		addArg(name, description, REQUIRED);
	}

	/**
	 * @see #addArg(String, String, String, boolean, Arg)
	 *
	 * @param name Argument name
	 * @param description Text to use when printing command-line usage
	 * @param required If {@code true}, {@link #parse} will fail if this argument is not explicitly
	 *        set or {@link #addArgWithDefault} was not used.
	 */
	public void addArg(String name, String description, boolean required) {
		addArg(name, description, required, TYPE);
	}

	/**
	 * @see #addArg(String, String, String, boolean, Arg)
	 *
	 * @param name Argument name
	 * @param description Text to use when printing command-line usage
	 * @param type Enforced parameter type from the {@link Arg} enum set.
	 */
	public void addArg(String name, String description, Arg type) {
		addArg(name, description, REQUIRED, type);
	}

	/**
	 * @see #addArg(String, String, String, boolean, Arg)
	 *
	 * @param name Argument name
	 * @param description Text to use when printing command-line usage
	 * @param required If {@code true}, {@link #parse} will fail if this argument is not explicitly
	 *        set or {@link #addArgWithDefault} was not used.
	 * @param type Enforced parameter type from the {@link Arg} enum set.
	 */
	public void addArg(String name, String description, boolean required, Arg type) {
		addArg(name, description, null, required, type);
	}

	/**
	 * @see #addArg(String, String, String, boolean, Arg)
	 *
	 * @param name Argument name
	 * @param description Text to use when printing command-line usage
	 * @param example Example text for this arg.
	 */
	public void addArg(String name, String description, String example) {
		addArg(name, description, example, REQUIRED);
	}


	/**
	 * @see #addArg(String, String, String, boolean, Arg)
	 *
	 * @param name Argument name
	 * @param description Text to use when printing command-line usage
	 * @param example Example text for this arg.
	 * @param required If {@code true}, {@link #parse} will fail if this argument is not explicitly
	 *        set or {@link #addArgWithDefault} was not used.
	 */
	public void addArg(String name, String description, String example, boolean required) {
		addArg(name, description, example, required, TYPE);
	}

	/**
	 * @see #addArg(String, String, String, boolean, Arg)
	 *
	 * @param name Argument name
	 * @param description Text to use when printing command-line usage
	 * @param example Example text for this arg.
	 * @param type Enforced parameter type from the {@link Arg} enum set.
	 */
	public void addArg(String name, String description, String example, Arg type) {
		addArg(name, description, example, REQUIRED, type);
	}

	/**
	 * @see #addArgWithDefault(String, String, String, Arg)
	 *
	 * @param name Argument name
	 * @param description Text to use when printing command-line usage
	 * @param example Default value of this argument. If specified, {@code required} is set to
	 *        {@code false}. If this argument is not present on the command-line, this default value
	 *        will be used.
	 */
	public void addArgWithDefault(String name, String description, String example) {
		addArgWithDefault(name, description, example, TYPE);
	}

	/**
	 * As {@link #addArg(String, String, String, boolean, Arg)} but the {@code example} value will
	 * also be used as the default value. This effectively sets {@code required} to {@code false}.
	 *
	 * @param name Argument name
	 * @param description Text to use when printing command-line usage
	 * @param example Default value of this argument. <i>Default: {@code \<value\>}</i>
	 * @param type Enforced parameter type from the {@link Arg} enum set.
	 */
	public void addArgWithDefault(String name, String description, String example, Arg type) {
		defaults.put(name, example);
		addArg(name, description, example + " (default)", false, type);
	}

	/**
	 * Add an <b>argument</b> (e.g. "key=value") command-line option to this {@link CLI} instance.
	 *
	 * @param name Argument name
	 * @param description Text to use when printing command-line usage
	 * @param example Example text for this arg. <i>Default: {@code \<value\>}</i>
	 * @param required If {@code true}, {@link #parse} will fail if this argument is not explicitly
	 *        set or {@link #addArgWithDefault} was not used. <i>Default: {@code false}</i>
	 * @param type Enforced parameter type from the {@link Arg} enum set. <i>Default:
	 *        {@link Arg#STRING}</i>
	 */
	public void addArg(String name, String description, String example, boolean required, Arg type) {
		String argExample = example == null ? ARG_VALUE : example;

		Builder builder = Option.builder().hasArg().argName(argExample).valueSeparator(ARG_SEPARATOR);

		add(builder, name, description, required, type);
	}

	/**
	 * @see #addFlag(String, String, String, boolean)
	 *
	 * @param name Text to enable this flag on command line
	 * @param description Text to use when printing command-line usage
	 */
	public void addFlag(String name, String description) {
		addFlag(name, description, null);
	}

	/**
	 * @see #addFlag(String, String, String, boolean)
	 *
	 * @param name Text to enable this flag on command line
	 * @param description Text to use when printing command-line usage
	 * @param required If {@code true}, {@link #parse} will fail if this argument is not set.
	 */
	public void addFlag(String name, String description, boolean required) {
		addFlag(name, description, null, required);
	}

	/**
	 * @see #addFlag(String, String, String, boolean)
	 *
	 * @param name Text to enable this flag on command line
	 * @param description Text to use when printing command-line usage
	 * @param longName Alternate long-parse name. If present, this and {@code name} can be used
	 *        interchangeably
	 */
	public void addFlag(String name, String description, String longName) {
		addFlag(name, description, longName, REQUIRED);
	}

	/**
	 * Add a <b>flag</b> (e.g. "-flag") command-line option to this {@link CLI}'s {@link Options} set.
	 *
	 * @param name Text to enable this flag on command line
	 * @param description Text to use when printing command-line usage
	 * @param longName Alternate long-parse name. If present, this and {@code name} can be used
	 *        interchangeably
	 * @param required If {@code true}, {@link #parse} will fail if this argument is not set.
	 *        <i>Default: {@code false}</i>
	 */
	public void addFlag(String name, String description, String longName, boolean required) {
		Builder builder = Option.builder(name);

		add(builder, longName, description, required, TYPE);
	}

	/**
	 * Create and add an {@link OptionGroup} to this {@link CLI}'s {@link Options} instance. An
	 * options group is a mutually-exclusive set of required options: <b>one and only one</b> of the
	 * options in this group will need to be set at runtime.
	 * <p>
	 * Note that adding an option to a group removes any associated default value
	 * </p>
	 *
	 * @param toGroup The name of one or more {@link Option} instances to group together.
	 */
	public void addGroup(String... toGroup) {
		OptionGroup g = new OptionGroup();
		g.setRequired(true);

		if (toGroup.length < 2) {
			throw new IllegalArgumentException("Group creation failed. Can not create an options group with one option.");
		}

		for (String opt : toGroup) {
			Option o = options.getOption(opt);

			if (o == null) {
				throw new IllegalArgumentException("Group creation failed. Option not found: " + opt);
			}

			// Pop the indentation off non-required args
			o.setDescription("(select one) "
												+ o.getDescription().substring(o.isRequired() ? 0 : INDENT.length()));
			g.addOption(o);
			defaults.remove(getName(o));
		}

		g.addOption(options.getOption(HELP_STRING));

		options.addOptionGroup(g);
	}

	/**
	 * @see #parseWithExit(Logger, String...)
	 *
	 * @param args Provided command-line arguments to parse
	 */
	public void parseWithExit(String... args) {
		parseWithExit(new Logger(), args);
	}

	/**
	 * Note: unlike {@link #parse} methods, if the act of parsing does not return normally (e.g. due
	 * to a parsing exception or help flag passed) this method <b>will stop program execution</b>. Use
	 * at your own risk.
	 *
	 * @see #parse(Logger, String...)
	 *
	 * @param log {@link Logger} instance to report any information resulting from parse
	 * @param args Provided command-line arguments to parse
	 */
	public void parseWithExit(Logger log, String... args) {
		try {
			parse(log, args);
		} catch (ParseException e) {
			// If there's no error message then this exception was thrown because a help flag was passed.
			int errorCode = e.getMessage().isEmpty() ? 0 : 1;
			System.exit(errorCode);
		}
	}

	/**
	 * @see #parse(Logger, String...)
	 *
	 * @param args Provided command-line arguments to parse
	 *
	 * @throws ParseException If a problem is encountered during parsing or a help/usage flag is
	 *         parsed.
	 */
	public void parse(String... args) throws ParseException {
		parse(new Logger(), args);
	}

	/**
	 * Performs the actual parsing of command-line options and return parsed values. This is the final
	 * step of command-line argument processing.
	 *
	 * @param log {@link Logger} instance to report any information resulting from parse. <i>Default:
	 *        new Logger()</i>
	 * @param args Provided command-line arguments to parse
	 *
	 * @throws ParseException If a problem is encountered during parsing or a help/usage flag is
	 *         parsed.
	 */
	public void parse(Logger log, String... args) throws ParseException {
		formatArgs(args);

		try {
			CommandLineParser parser = new DefaultParser();
			CommandLine cl = parser.parse(options, args);

			// First see if the help flags are present
			if (cl.hasOption(HELP_CHAR) || cl.hasOption(HELP_STRING)) {
				throw new ParseException("");
			}

			parsed = new HashMap<String, String>();
			for (Option o : options.getOptions()) {
				String opt = getName(o);

				if (cl.hasOption(opt)) {
					// The option was present on the command line.
					// Make sure the types match, then set the value.
					validateOption(cl, opt, o);
					parsed.put(opt, cl.getOptionValue(opt));
				} else if (o.hasArg() && defaults.containsKey(opt)) {
					// If an argument was not parsed but does have a default value,
					// we set it here
					parsed.put(opt, defaults.get(opt));
				}
			}
		} catch (ParseException e) {
			printHelp(log, commandName, e.getMessage(), options);
			throw e;
		}
	}

	/**
	 * Helper method to ensure args conform to commons-cli format
	 */
	private static void formatArgs(String... args) {
		// Can't figure out how to configure CommandLineParser to allow arguments (but not flags) to
		// drop the "-" character... so, manually prepend it here.
		for (int i = 0; i < args.length; i++) {
			if (args[i].charAt(0) != '-') {
				StringBuilder sb = new StringBuilder().append('-').append(args[i]);
				args[i] = sb.toString();
			}
		}
	}

	/**
	 * Helper method to get the name of the given {@link Option} by checking both
	 * {@link Option#getOpt()} and {@link Option#getLongOpt()}.
	 */
	private static String getName(Option o) {
		String name = o.getOpt();
		if (name == null) {
			name = o.getLongOpt();
		}
		return name;
	}

	/**
	 * Helper method to print the help/usage dialog. If an optional {@code error} string is included,
	 * this message will be included after the standard usage text (helpful for reporting why a parse
	 * failed, for example).
	 */
	private static void printHelp(Logger log, String appName, String errorMessage, Options options) {
		HelpFormatter formatter = new HelpFormatter();
		formatter.setOptPrefix(FLAG_PREFIX);
		formatter.setLongOptPrefix(ARG_PREFIX);
		formatter.setLongOptSeparator(String.valueOf(ARG_SEPARATOR));
		formatter.setOptionComparator(null);

		StringBuilder sb = new StringBuilder().append(System.getProperty("line.separator"))
																					.append(errorMessage);
		String program = appName + " [-FLAG | ARG=value]...";
		PrintWriter pw = new PrintWriter(System.out);
		formatter.printHelp(pw, OUT_WIDTH, program, "", options, HelpFormatter.DEFAULT_LEFT_PAD,
												HelpFormatter.DEFAULT_DESC_PAD, sb.toString());
		pw.flush();

		String logfile = log.getFilename();
		if (logfile != null && !logfile.isEmpty()) {
			try {
				pw = new PrintWriter(new FileWriter(logfile));
				formatter.printHelp(pw, OUT_WIDTH, program, "", options, HelpFormatter.DEFAULT_LEFT_PAD,
														HelpFormatter.DEFAULT_DESC_PAD, sb.toString());
				pw.flush();
				pw.close();
			} catch (IOException e) {
				log.reportException(e);
			}
		}
	}

	/**
	 * Helper method to make sure an {@link Option} has a value in the given {@link CommandLine} that
	 * is compatible with its declared type.
	 *
	 * @throws ParseException if the value is not compatible with declared type
	 */
	private static void validateOption(CommandLine cl, String opt, Option o) throws ParseException {
		if (o.hasArg()) {
			try {
				cl.getParsedOptionValue(opt);
			} catch (ParseException e) {
				StringBuilder sb = new StringBuilder();
				sb.append("Argument \"").append(opt).append("\" with value \"")
					.append(cl.getOptionValue(opt)).append("\" does not conform to required type: ")
					.append(o.getType());
				throw new ParseException(sb.toString());
			}
		}
	}

	/**
	 * Helper method to forms a valid argument String
	 * 
	 * @param name Argument name as specified to
	 *        {@link #addArg(String, String, String, boolean, Class)}
	 * @param value Passed value for this argument
	 * @return Formatted command line argument
	 */
	public static String formCmdLineArg(String name, String value) {
		return ARG_PREFIX + name + ARG_SEPARATOR + value;
	}

	/**
	 * Helper method to form a valid flag String
	 * 
	 * @param name Flag name as specified to {@link #addFlag(String, String, String, boolean)}
	 * @return Formatted command line flag
	 */
	public static String formCmdLineFlag(String name) {
		return FLAG_PREFIX + name;
	}
}
