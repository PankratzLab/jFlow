package org.genvisis;

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
 * <li>Create an {@link Options} object with {@link #defaultOptions()}. This
 * ensures a standard set of help flags.</li>
 * <li>Build up this {@link Options} set by using {@link #addFlag} (e.g.
 * <code>-doThisOption</code>) and {@link #addArg} (e.g.,
 * <code>myOpt=someValue</code>)</li>
 * <li>(Optional) If your program includes mutually exclusive commands and you
 * want to ensure a user selects one (and only one), group these options with
 * {@link #addGroup}</li>
 * <li>Parse your {@link Options} with <b>either</b> {@link #parseWithExit}, if
 * you want your application to exit if the command-line usage is printed, or
 * {@link #parse} if you want to handle any exception cases yourself.</li>
 * <li>After parsing you can interrogate the returned {@link Map} of argument
 * names to values (or {@code null} for flags)</li>
 * </ol>
 */
public final class CLI {
	private static final String HELP_CHAR = "h";
	private static final String HELP_STRING = "help";
	private static final String ARG_VALUE = "value";
	private static final boolean REQUIRED = false;
	private static final Class<?> TYPE = PatternOptionBuilder.STRING_VALUE;
	private static final int OUT_WIDTH = 150;

	private CLI() {
		// Prevent construction of static utility class
	}

	/**
	 * @return A pre-configured {@link Options} that includes supports help flags
	 */
	public static Options defaultOptions() {
		Options options = new Options();
		addFlag(options, HELP_CHAR, "Print application usage", HELP_STRING);
		
		return options;
	}

	/**
	 * @see {@link #addFlag(Options, String, String, String, boolean)}
	 *
	 * @param options
	 *            {@link Options} collection, modified in-place by adding the
	 *            {@link Option} created by this call.
	 * @param name
	 *            Text to enable this flag on command line
	 * @param description
	 *            Text to use when printing command-line usage
	 */
	public static void addFlag(Options options, String name, String description) {
		addFlag(options, name, description, null);
	}

	/**
	 * @see {@link #addFlag(Options, String, String, String, boolean)}
	 *
	 * @param options
	 *            {@link Options} collection, modified in-place by adding the
	 *            {@link Option} created by this call.
	 * @param name
	 *            Text to enable this flag on command line
	 * @param description
	 *            Text to use when printing command-line usage
	 * @param required
	 *            If {@code true}, {@link #parse} will fail if this argument is
	 *            not set. <i>Default: {@code false}</i>
	 */
	public static void addFlag(Options options, String name, String description, boolean required) {
		addFlag(options, name, description, null, required);
	}

	/**
	 * @see {@link #addFlag(Options, String, String, String, boolean)}
	 *
	 * @param options
	 *            {@link Options} collection, modified in-place by adding the
	 *            {@link Option} created by this call.
	 * @param name
	 *            Text to enable this flag on command line
	 * @param description
	 *            Text to use when printing command-line usage
	 * @param longName
	 *            Alternate long-parse name. If present, this and {@code name} can be used interchangeably
	 */
	public static void addFlag(Options options, String name, String description, String longName) {
		addFlag(options, name, description, longName, REQUIRED);
	}


	/**
	 * Add a <b>flag</b> (e.g. "-flag") command-line option to the given
	 * {@link Options} object.
	 *
	 * @param options
	 *            {@link Options} collection, modified in-place by adding the
	 *            {@link Option} created by this call.
	 * @param name
	 *            Text to enable this flag on command line
	 * @param description
	 *            Text to use when printing command-line usage
	 * @param longName
	 *            Alternate long-parse name. If present, this and {@code name} can be used interchangeably
	 * @param required
	 *            If {@code true}, {@link #parse} will fail if this argument is
	 *            not set. <i>Default: {@code false}</i>
	 */
	public static void addFlag(Options options, String name, String description, String longName, boolean required) {
		Builder builder = Option.builder(name);

		add(options, builder, longName, description, required, null);
	}

	/**
	 * @see #addArg(Options, String, String, String, boolean, Class)
	 *
	 * @param options
	 *            {@link Options} collection, modified in-place by adding the
	 *            {@link Option} created by this call.
	 * @param name
	 *            Argument name
	 * @param description
	 *            Text to use when printing command-line usage
	 */
	public static void addArg(Options options, String name, String description) {
		addArg(options, name, description, null);
	}

	/**
	 * @see #addArg(Options, String, String, String, boolean, Class)
	 *
	 * @param options
	 *            {@link Options} collection, modified in-place by adding the
	 *            {@link Option} created by this call.
	 * @param name
	 *            Argument name
	 * @param description
	 *            Text to use when printing command-line usage
	 * @param required
	 *            If {@code true}, {@link #parse} will fail if this argument is
	 *            not explicitly set and no {@code argDefault} is provided.
	 *            <i>Default: {@code false}</i>
	 */
	public static void addArg(Options options, String name, String description, boolean required) {
		addArg(options, name, description, null, required);
	}

	/**
	 * @see #addArg(Options, String, String, String, boolean, Class)
	 *
	 * @param options
	 *            {@link Options} collection, modified in-place by adding the
	 *            {@link Option} created by this call.
	 * @param name
	 *            Argument name
	 * @param description
	 *            Text to use when printing command-line usage
	 * @param required
	 *            If {@code true}, {@link #parse} will fail if this argument is
	 *            not explicitly set and no {@code argDefault} is provided.
	 *            <i>Default: {@code false}</i>
	 * @param type
	 *            Any argument value must be convertible to this type (e.g.
	 *            through a string constructor, or {@code valueOf(String)}
	 *            method). <i>Default: {@link String}</i>
	 */
	public static void addArg(Options options, String name, String description, boolean required, Class<?> type) {
		addArg(options, name, description, null, required, type);
	}

	/**
	 * @see #addArg(Options, String, String, String, boolean, Class)
	 *
	 * @param options
	 *            {@link Options} collection, modified in-place by adding the
	 *            {@link Option} created by this call.
	 * @param name
	 *            Argument name
	 * @param description
	 *            Text to use when printing command-line usage
	 * @param argDefault
	 *            Default value of this argument. If specified, {@code required}
	 *            is set to {@code false}. If this argument is not present on
	 *            the command-line, the default value will be used.
	 *            <i>Default: {@code null}</i>
	 */
	public static void addArg(Options options, String name, String description, String argDefault) {
		addArg(options, name, description, argDefault, REQUIRED);
	}

	/**
	 * @see #addArg(Options, String, String, String, boolean, Class)
	 *
	 * @param options
	 *            {@link Options} collection, modified in-place by adding the
	 *            {@link Option} created by this call.
	 * @param name
	 *            Argument name
	 * @param description
	 *            Text to use when printing command-line usage
	 * @param argDefault
	 *            Default value of this argument. If specified, {@code required}
	 *            is set to {@code false}. If this argument is not present on
	 *            the command-line, the default value will be used.
	 *            <i>Default: {@code null}</i>
	 * @param required
	 *            If {@code true}, {@link #parse} will fail if this argument is
	 *            not explicitly set and no {@code argDefault} is provided.
	 *            <i>Default: {@code false}</i>
	 */
	public static void addArg(Options options, String name, String description, String argDefault, boolean required) {
		addArg(options, name, description, argDefault, required, TYPE);
	}

	/**
	 * Add an <b>argument</b> (e.g. "key=value") command-line option to the
	 * given {@link Options} object.
	 *
	 * @param options
	 *            {@link Options} collection, modified in-place by adding the
	 *            {@link Option} created by this call.
	 * @param name
	 *            Argument name
	 * @param description
	 *            Text to use when printing command-line usage
	 * @param argDefault
	 *            Default value of this argument. If specified, {@code required}
	 *            is set to {@code false}. If this argument is not present on
	 *            the command-line, the default value will be used. <i>Default:
	 *            {@code \<value\>}</i>
	 * @param required
	 *            If {@code true}, {@link #parse} will fail if this argument is
	 *            not explicitly set and no {@code argDefault} is provided.
	 *            <i>Default: {@code false}</i>
	 * @param type
	 *            Enforce typing of this argument. <b>Notice</b>: the specified
	 *            class must be one of the {@link PatternOptionBuilder} value
	 *            codes. <i>Default:
	 *            {@link PatternOptionBuilder#STRING_VALUE}</i>
	 */
	public static void addArg(Options options, String name, String description, String argDefault, boolean required, Class<?> type) {
		String argVal = ARG_VALUE;
		boolean reallyRequired = required;

		if (argDefault != null) {
			// If a default value is present, this arg is not really required.
			reallyRequired = false;
			argVal = argDefault;
		}

		Builder builder = Option.builder().hasArg().argName(argVal).valueSeparator('=');

		add(options, builder, name, description, reallyRequired, type);
	}

	/**
	 * Helper method to finish building an {@link Option}, whether from
	 * {@link #addFlag} or {@link #addArg}, and add it to the given
	 * {@link Options}
	 */
	private static void add(Options options, Builder builder, String longName, String desc, boolean required, Class<?> type) {
		if (longName != null) {
			builder = builder.longOpt(longName);
		}

		// Mark the description so required flags will be clearly indicated in usage text
		String description = (required ? "(required) " : "") + desc;
		builder = builder.type(type);
		builder = builder.desc(description).required(required);

		options.addOption(builder.build());
	}

	/**
	 * Create and add an {@link OptionGroup} to the given {@link Options}
	 * instance. An options group is a mutually-exclusive set of required
	 * options: <b>one and only one</b> of the options in this group will need
	 * to be set at runtime.
	 *
	 * @param options
	 *            {@link Options} collection containing the {@link Option}
	 *            instances to group.
	 * @param toGroup
	 *            The name of one or more {@link Option} instances to group
	 *            together.
	 */
	public static void addGroup(Options options, String... toGroup) {
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

			o.setDescription("(select one) " + o.getDescription());
			g.addOption(o);
		}
	
		g.addOption(options.getOption(HELP_STRING));
		
		options.addOptionGroup(g);
	}

	/**
	 * Note: unlike {@link #parse} methods, if the act of parsing does not
	 * return normally (e.g. due to a parsing exception or help flag passed)
	 * this method <b>will not return</b>. Use at your own risk.
	 * 
	 * @see {@link #parseWithExit(String, Options, Logger, String...)}
	 * 
	 * @param appClass
	 *            Class with {@code main} method being run (for usage message -
	 *            how a user should invoke from command line)
	 * @param options
	 *            {@link Options} instance to use in parsing
	 * @param args
	 *            Provided command-line arguments to parse
	 * @return A {@link Map} of {@link Option} names to parsed values.
	 */
	public static Map<String, String> parseWithExit(Class<?> appClass, Options options, String... args) {
		return parseWithExit(appClass, options, new Logger(), args);
	}

	/**
	 * Note: unlike {@link #parse} methods, if the act of parsing does not
	 * return normally (e.g. due to a parsing exception or help flag passed)
	 * this method <b>will not return</b>. Use at your own risk.
	 * 
	 * @see {@link #parseWithExit(String, Options, Logger, String...)}
	 * 
	 * @param appName
	 *            Name of the program being run (for usage message - how a user
	 *            should invoke from command line) * @param options
	 *            {@link Options} instance to use in parsing
	 * @param options
	 *            {@link Options} instance to use in parsing
	 * @param args
	 *            Provided command-line arguments to parse
	 * @return A {@link Map} of {@link Option} names to parsed values.
	 */
	public static Map<String, String> parseWithExit(String appName, Options options, String... args) {
		return parseWithExit(appName, options, new Logger(), args);
	}

	/**
	 * Note: unlike {@link #parse} methods, if the act of parsing does not
	 * return normally (e.g. due to a parsing exception or help flag passed)
	 * this method <b>will not return</b>. Use at your own risk.
	 * 
	 * @see {@link #parseWithExit(String, Options, Logger, String...)}
	 * 
	 * @param appClass
	 *            Class with {@code main} method being run (for usage message -
	 *            how a user should invoke from command line)
	 * @param options
	 *            {@link Options} instance to use in parsing
	 * @param log
	 *            {@link Logger} instance to report any information resulting
	 *            from parse
	 * @param args
	 *            Provided command-line arguments to parse
	 * @return A {@link Map} of {@link Option} names to parsed values.
	 */
	public static Map<String, String> parseWithExit(Class<?> appClass, Options options, Logger log, String... args) {
		return parseWithExit(appClass.getName(), options, log, args);
	}

	/**
	 * Performs the actual parsing of command-line options and return parsed
	 * values. This is the final step of command-line argument processing.
	 * <p>
	 * Note: unlike {@link #parse} methods, if the act of parsing does not
	 * return normally (e.g. due to a parsing exception or help flag passed)
	 * this method <b>will not return</b>. Use at your own risk.
	 * </p>
	 *
	 * @see {@link #parse(String, Options, Logger, String...)}
	 * 
	 * @param appName
	 *            Name of the program being run (for usage message - how a user
	 *            should invoke from command line) * @param options
	 *            {@link Options} instance to use in parsing
	 * @param options
	 *            {@link Options} instance to use in parsing
	 * @param log
	 *            {@link Logger} instance to report any information resulting
	 *            from parse
	 * @param args
	 *            Provided command-line arguments to parse
	 * @return A {@link Map} of {@link Option} names to parsed values.
	 */
	public static Map<String, String> parseWithExit(String appName, Options options, Logger log, String... args) {
		Map<String, String> parsed = null;
		try {
			parsed = parse(appName, options, log, args);
		}
		catch (ParseException e) {
			// If there's no error message then this exception was thrown because a help flag was passed.
			int errorCode = e.getMessage().isEmpty() ? 0 : 1;
			System.exit(errorCode);
		}

		return parsed;
	}

	/**
	 * @see {@link #parse(String, Options, Logger, String...)}
	 *
	 * @param appClass
	 *            Class with {@code main} method being run (for usage message -
	 *            how a user should invoke from command line)
	 * @param options
	 *            {@link Options} instance to use in parsing
	 * @param args
	 *            Provided command-line arguments to parse
	 * @return A {@link Map} of {@link Option} names to parsed values.
	 */
	public static Map<String, String> parse(Class<?> appClass, Options options, String... args) throws ParseException {
		return parse(appClass, options, new Logger(), args);
	}

	/**
	 * @see {@link #parse(String, Options, Logger, String...)}
	 *
	 * @param appName
	 *            Name of the program being run (for usage message - how a user
	 *            should invoke from command line)
	 * @param options
	 *            {@link Options} instance to use in parsing
	 * @param args
	 *            Provided command-line arguments to parse
	 * @return A {@link Map} of {@link Option} names to parsed values.
	 */
	public static Map<String, String> parse(String appName, Options options, String... args) throws ParseException {
		return parse(appName, options, new Logger(), args);
	}

	/**
	 * @see {@link #parse(String, Options, Logger, String...)}
	 *
	 * @param appClass
	 *            Class with {@code main} method being run (for usage message -
	 *            how a user should invoke from command line)
	 * @param options
	 *            {@link Options} instance to use in parsing
	 * @param log
	 *            {@link Logger} instance to report any information resulting
	 *            from parse
	 * @param args
	 *            Provided command-line arguments to parse
	 * @return A {@link Map} of {@link Option} names to parsed values.
	 */
	public static Map<String, String> parse(Class<?> appClass, Options options, Logger log, String... args) throws ParseException {
		return parse(appClass.getName(), options, log, args);
	}

	/**
	 * Performs the actual parsing of command-line options and return parsed
	 * values. This is the final step of command-line argument processing.
	 *
	 * @param appName
	 *            Name of the program being run (for usage message - how a user
	 *            should invoke from command line)
	 * @param options
	 *            {@link Options} instance to use in parsing
	 * @param log
	 *            {@link Logger} instance to report any information resulting
	 *            from parse. <i>Default: new Logger()</i>
	 * @param args
	 *            Provided command-line arguments to parse
	 * @return A {@link Map} of {@link Option} names to parsed values.
	 */
	public static Map<String, String> parse(String appName, Options options, Logger log, String... args) throws ParseException {
		// Can't figure out how to configure CommandLineParser to allow arguments (but not flags) to drop the "-"
		// character... so, manually prepend it here.
		for (int i=0; i<args.length; i++) {
			if (args[i].charAt(0) != '-') {
				StringBuilder sb = new StringBuilder().append('-').append(args[i]);
				args[i] = sb.toString();
			}
		}
		
		try {
			CommandLineParser parser = new DefaultParser();
			CommandLine cl = parser.parse(options, args);

			// First see if the help flags are present
			if (cl.hasOption(HELP_CHAR) || cl.hasOption(HELP_STRING)) {
				throw new ParseException("");
			}

			Map<String, String> parsed = new HashMap<String, String>();
			for (Option o : options.getOptions()) {
				String opt = getName(o);

				if (cl.hasOption(opt)) {
					// The option was present on the command line.
					// Make sure the types match, then set the value.
					validateOption(cl, opt, o);
					parsed.put(opt, cl.getOptionValue(opt));
				}
				else if (o.hasArg() && !ARG_VALUE.equals(o.getArgName())) {
					// If an argument was not parsed but does have a default value,
					// we set it here
					parsed.put(opt, o.getArgName());
				}
			}
			return parsed;
		} catch (ParseException e) {
			printHelp(log, appName, e.getMessage(), options);
			throw e;
		}
	}

	/**
	 * Helper method to get the name of the given {@link Option} by checking
	 * both {@link Option#getOpt()} and {@link Option#getLongOpt()}.
	 */
	private static String getName(Option o) {
		String name = o.getOpt();
		if (name == null) {
			name = o.getLongOpt();
		}
		return name;
	}

	/**
	 * Helper method to make sure an {@link Option} has a value in the given
	 * {@link CommandLine} that is compatible with its declared type.
	 *
	 * @throws ParseException if the value is not compatible with declared type
	 */
	private static void validateOption(CommandLine cl, String opt, Option o) throws ParseException {
		if (o.hasArg()) {
			try {
				cl.getParsedOptionValue(opt);
			}
			catch (ParseException e) {
				StringBuilder sb = new StringBuilder();
				sb.append("Argument \"").append(opt)
				.append("\" with value \"").append(cl.getOptionValue(opt))
				.append("\" does not conform to required type: ").append(o.getType());
				throw new ParseException(sb.toString());
			}
		}
	}

	/**
	 * Helper method to print the help/usage dialog. If an optional
	 * {@code error} string is included, this message will be included after the
	 * standard usage text (helpful for reporting why a parse failed, for
	 * example).
	 */
	private static void printHelp(Logger log, String appName, String errorMessage, Options options) {
		HelpFormatter formatter = new HelpFormatter();
		formatter.setOptPrefix("-");
		formatter.setLongOptPrefix("");
		formatter.setLongOptSeparator("=");
		formatter.setOptionComparator(null);

		StringBuilder sb = new StringBuilder().append(System.lineSeparator()).append(errorMessage);
		String program = appName + " [-FLAG | ARG=value]...";
		PrintWriter pw = new PrintWriter(System.out);
		formatter.printHelp(pw, OUT_WIDTH, program, "", options, HelpFormatter.DEFAULT_LEFT_PAD, HelpFormatter.DEFAULT_DESC_PAD, sb.toString());
		pw.flush();
		
		String logfile = log.getFilename();
		if (logfile != null && !logfile.isEmpty()) {
			try {
				pw = new PrintWriter(new FileWriter(logfile));
				formatter.printHelp(pw, OUT_WIDTH, program, "", options, HelpFormatter.DEFAULT_LEFT_PAD, HelpFormatter.DEFAULT_DESC_PAD, sb.toString());
				pw.flush();
				pw.close();
			} catch (IOException e) {
				log.reportException(e);
			}
		}
	}
}
